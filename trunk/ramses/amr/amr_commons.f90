module amr_commons
  use amr_parameters
  use hydro_parameters
  use hash
  
  logical::output_done=.false.                  ! Output just performed
  logical::init=.false.                         ! Set up or run
  logical::balance=.false.                      ! Load balance or run
  logical::shrink=.false.                       ! Shrink mesh or run
  integer::nstep=0                              ! Time step
  integer::nstep_coarse=0                       ! Coarse step
  integer::nstep_coarse_old=0                   ! Old coarse step
  integer::nflag,ncreate,nkill                  ! Refinements
  integer::ncoarse                              ! nx.ny.nz
  integer::ngrid_current                        ! Actual number of octs

  real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
  real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
  real(dp)::eint_tot=0.0D0                      ! Total internal energy
  real(dp)::epot_tot=0.0D0                      ! Total potential energy
  real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
  real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
  real(dp)::const=0.0D0                         ! Energy conservation
  real(dp)::aexp_old=1.0D0                      ! Old expansion factor
  real(dp)::rho_tot=0.0D0                       ! Mean density in the box
  real(dp)::t=0.0D0                             ! Time variable

  ! executable identification
  CHARACTER(LEN=80)::builddate,patchdir
  CHARACTER(LEN=80)::gitrepo,gitbranch,githash

  ! Save namelist filename
  CHARACTER(LEN=80)::namelist_file

  ! MPI variables
  integer::ncpu,ndomain,myid,overload=1

  ! Friedman model variables
  integer::n_frw
  real(dp),allocatable,dimension(:)::aexp_frw,hexp_frw,tau_frw,t_frw

  ! Initial conditions parameters from grafic
  integer                  ::nlevelmax_part
  real(dp)                 ::aexp_ini=10.
  real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
  real(dp),dimension(1:MAXLEVEL)::vfact
  real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
  integer ,dimension(1:MAXLEVEL)::n1,n2,n3

  ! Level related arrays
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
  real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level

  ! Oct structure
  type oct
     integer(kind=4)::lev
     integer(kind=4),dimension(1:ndim)::ckey
     integer(kind=8)::hkey
     integer(kind=4),dimension(1:twotondim)::flag1
     integer(kind=4),dimension(1:twotondim)::flag2
     logical,dimension(1:twotondim)::refined
     integer(kind=4)::superoct
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim)::rho
     real(kind=dp),dimension(1:twotondim)::phi
     real(kind=dp),dimension(1:twotondim)::phi_old
     real(kind=dp),dimension(1:twotondim,1:ndim)::f
#endif
#ifdef SOLVERhydro
     real(kind=dp),dimension(1:twotondim,1:nvar)::uold
     real(kind=dp),dimension(1:twotondim,1:nvar)::unew
#endif
#ifdef SOLVERmhd
     real(kind=dp),dimension(1:twotondim,1:nvar+3)::uold
     real(kind=dp),dimension(1:twotondim,1:nvar+3)::unew
#endif
#ifdef DUALENER
     real(kind=dp),dimension(1:twotondim)::divu
     real(kind=dp),dimension(1:twotondim)::enew
#endif
  end type oct

  ! Tile structure
  type tile
     ! Header information
     integer(kind=4)::lev    ! Level of the tile
     integer(kind=4)::cpu    ! CPU in which the tile sits
     integer(kind=4)::status ! Status of the tile
     ! AMR data
     ! Oct-based arrays
     integer(kind=4),dimension(1:ndim)::ckey ! Cartesian key of each oct
     integer(kind=8)::hkey                   ! Hilbert key of each oct
     ! Cell-based arrays
     integer,dimension(1:twotondim)::flag1   ! Flag1
     integer,dimension(1:twotondim)::flag2   ! Flag2
     logical,dimension(1:twotondim)::refined ! Refined flag
#ifdef GRAV
     ! Tile gravity data
     real(kind=dp),dimension(1:twotondim)::rho       ! Total mass density
     real(kind=dp),dimension(1:twotondim)::phi       ! Gravitationall potential
     real(kind=dp),dimension(1:twotondim)::phi_old   ! Previous potential
     real(kind=dp),dimension(1:twotondim,1:ndim)::f  ! Gravity acceleration
#endif
#ifdef SOLVERhydro
     ! Tile hydro data
     real(kind=dp),dimension(1:twotondim,1:nvar)::uold ! Old conservative variables
     real(kind=dp),dimension(1:twotondim,1:nvar)::unew ! New conservative variables
#endif
#ifdef SOLVERmhd
     ! Tile MHD data
     real(kind=dp),dimension(1:twotondim,1:nvar+3)::uold ! Old conservative variables
     real(kind=dp),dimension(1:twotondim,1:nvar+3)::unew ! New conservative variables
#endif
#ifdef DUALENER
     ! Tile dual energy forumlation data
     real(kind=dp),dimension(1:twotondim)::divu ! Velocity divergence
     real(kind=dp),dimension(1:twotondim)::enew ! New internal energy density
#endif
  end type tile

  type cache_cell
     integer(kind=4)::grid
     integer(kind=4)::cell
  end type cache_cell

  ! Persistent array for the AMR grid
  type(oct),dimension(:),allocatable::grid
  type(hash_table)::grid_dict   ! Oct hash table

  ! Starting index for each level 
  integer,allocatable,dimension(:)::head
  integer,allocatable,dimension(:)::tail
  integer,allocatable,dimension(:)::noct
  integer,allocatable,dimension(:)::ckey_max
  integer::noct_used

  ! Hilbert key
  integer(kind=8),allocatable,dimension(:,:)::bound_key_level

  ! Software cache array for the AMR grid
  type(tile),dimension(:),allocatable::cache
  type(hash_table)::cache_dict  ! Oct hash table
  integer::free_line
  integer::cache_type

  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro='hydro'
  CHARACTER(LEN=20)::type_accel='accel'
  CHARACTER(LEN=20)::type_flag='flag'

  ! Default units
  real(dp)::units_density=1.0 ! [g/cm^3]
  real(dp)::units_time=1.0    ! [seconds]
  real(dp)::units_length=1.0  ! [cm]

end module amr_commons


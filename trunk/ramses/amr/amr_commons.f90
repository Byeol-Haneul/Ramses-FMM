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
  integer::nlevelmax_part
  real(dp)::aexp_ini=10.
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
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar)::uold
     real(kind=dp),dimension(1:twotondim,1:nvar)::unew
#endif
#ifdef DUALENER
     real(kind=dp),dimension(1:twotondim)::divu
     real(kind=dp),dimension(1:twotondim)::enew
#endif
  end type oct

  ! Persistent array for the AMR grid
  type(oct),dimension(:),allocatable::grid
  type(hash_table)::grid_dict   ! Oct hash table

  ! Starting index for each level 
  integer,allocatable,dimension(:)::head
  integer,allocatable,dimension(:)::tail
  integer,allocatable,dimension(:)::noct
  integer,allocatable,dimension(:)::noct_min
  integer,allocatable,dimension(:)::noct_max
  integer,allocatable,dimension(:)::noct_tot
  integer,allocatable,dimension(:)::ckey_max
  integer::noct_used,noct_used_max,noct_used_tot

  ! Hilbert key
  integer(kind=8),allocatable,dimension(:,:)::bound_key_level

  ! Software cache parameters
  integer::cache_operation
  integer::operation_initflag=1,operation_upload=2,operation_godunov=3,operation_smooth=4
  integer::operation_hydro=5,operation_refine=6,operation_derefine=7,operation_loadbalance=8
  integer::operation_phi=9,operation_rho=10,operation_multipole=11,operation_conjgrad=12
  integer::cache_operation_type
  integer::operation_type_flag=1,operation_type_hydro=2,operation_type_poisson=3  
  integer::operation_type_refine=4

  ! Software cache array for the AMR grid
  logical,allocatable,dimension(:)::dirty
  logical,allocatable,dimension(:)::occupied
  logical,allocatable,dimension(:)::locked
  integer,allocatable,dimension(:)::parent_cpu
  integer::free_cache,ncache,ifree

  ! Software cache array for failed requests
  logical,allocatable,dimension(:)::occupied_null
  integer,allocatable,dimension(:)::lev_null
  integer,allocatable,dimension(:,:)::ckey_null
  integer::free_null,nnull

  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro='hydro'
  CHARACTER(LEN=20)::type_accel='accel'
  CHARACTER(LEN=20)::type_flag='flag'

  ! Default units
  real(dp)::units_density=1.0 ! [g/cm^3]
  real(dp)::units_time=1.0    ! [seconds]
  real(dp)::units_length=1.0  ! [cm]

  ! Communication-related objects
  integer::mail_counter=0
  integer::request_id,flush_id
  integer,allocatable,dimension(:)::reply_id
  logical::expect_flush=.false.

  ! New MPI derived types
  integer::new_mpi_int4_msg,new_mpi_realdp_msg,new_mpi_request
  integer::new_mpi_small_realdp_msg,new_mpi_large_realdp_msg
  integer::new_mpi_int4_flush,new_mpi_realdp_flush
  integer::new_mpi_small_realdp_flush,new_mpi_large_realdp_flush
  integer::flush_tag=1000,msg_tag=100,request_tag=10

  ! Request message buffer
  type request
     integer(kind=4)::lev
     integer(kind=4),dimension(1:ndim)::ckey
  end type request

  ! Response message buffer
  integer,parameter::ntilemax=16
  type int4_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
  end type int4_msg
  type realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar,1:ntilemax)::realdp
  end type realdp_msg
  type small_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     real(kind=dp),dimension(1:twotondim,1:ntilemax)::realdp
  end type small_realdp_msg
  type large_realdp_msg
     sequence
     integer(kind=4)::type
     integer(kind=4)::ntile
     integer(kind=4),dimension(1:ntilemax)::lev
     integer(kind=4),dimension(1:ndim,1:ntilemax)::ckey
     integer(kind=4),dimension(1:twotondim,1:ntilemax)::int4
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar,1:ntilemax)::realdp_hydro
#endif
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim,1:ndim+2,1:ntilemax)::realdp_poisson
#endif
  end type large_realdp_msg

  ! Fush message buffer
  integer,parameter::nflushmax=128
  type int4_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
  end type int4_flush
  type realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
     real(kind=dp),dimension(1:twotondim,1:nvar,1:nflushmax)::realdp
  end type realdp_flush
  type small_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     real(kind=dp),dimension(1:twotondim,1:nflushmax)::realdp
  end type small_realdp_flush
  type large_realdp_flush
     sequence
     integer(kind=4)::nflush
     integer(kind=4),dimension(1:nflushmax)::lev
     integer(kind=4),dimension(1:ndim,1:nflushmax)::ckey
     integer(kind=4),dimension(1:twotondim,1:nflushmax)::int4
#ifdef HYDRO
     real(kind=dp),dimension(1:twotondim,1:nvar,1:nflushmax)::realdp_hydro
#endif
#ifdef GRAV
     real(kind=dp),dimension(1:twotondim,1:ndim+2,1:nflushmax)::realdp_poisson
#endif
  end type large_realdp_flush

  ! Communication buffers
  type(request)::recv_request
  type(int4_msg),allocatable,dimension(:)::reply_flag
  type(realdp_msg),allocatable,dimension(:)::reply_hydro
  type(small_realdp_msg),allocatable,dimension(:)::reply_poisson
  type(large_realdp_msg),allocatable,dimension(:)::reply_refine
  type(int4_flush)::recv_flush_flag
  type(realdp_flush)::recv_flush_hydro
  type(small_realdp_flush)::recv_flush_poisson
  type(large_realdp_flush)::recv_flush_refine
  type(int4_flush),allocatable,dimension(:)::send_flush_flag
  type(realdp_flush),allocatable,dimension(:)::send_flush_hydro
  type(small_realdp_flush),allocatable,dimension(:)::send_flush_poisson
  type(large_realdp_flush),allocatable,dimension(:)::send_flush_refine

end module amr_commons


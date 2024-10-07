module amr_commons
  use amr_parameters
  use hydro_parameters
  use oct_commons
  use hydro_commons
  use hash
  use domain_m
  
  type multipole_t
    real(dp),dimension(1:ndim+1)::q
  end type multipole_t

  type run_t

     ! Run control
     logical::cosmo   =.false.   ! Cosmology activated
     logical::pic     =.false.   ! Particle In Cell activated
     logical::poisson =.false.   ! Poisson solver activated
     logical::hydro   =.false.   ! Hydro activated
     logical::star    =.false.   ! Stars and star formation activated
     logical::sink    =.false.   ! Sinks and sink formation activated
     logical::verbose =.false.   ! Write everything
     logical::debug   =.false.   ! Debug mode activated
     integer::nrestart=0         ! New run or backup file number
     integer::ncontrol=1         ! Write control variables
     integer::nstepmax=1000000   ! Maximum number of time steps
     integer,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level
     integer::nremap=10          ! Load balancing frequency (0: never)
     logical::static  =.false.   ! Static mode activated     
     integer::geom    =1         ! 1: cartesian, 2: cylindrical, 3: spherical
     integer::overload=1         ! MPI domain overloading
     integer::nsuperoct=0        ! Number of superoct levels
     
     ! Output parameters
     integer::noutput=1          ! Total number of outputs
     integer::foutput=1000000    ! Frequency of outputs
     real(dp),dimension(1:MAXOUT)::aout=1.1 ! Output expansion factors
     real(dp),dimension(1:MAXOUT)::tout=0.0 ! Output times
     integer::output_mode=0      ! Output mode (for hires runs)
     logical::gadget_output=.false. ! Output in gadget format
     real(kind=8)::bkp_time_hrs=2   ! Backup file frequency in hours
     real(kind=8)::run_time_hrs=0   ! Estimated run time in hrs
     real(kind=8)::bkp_last_min=10  ! Backup file before the end of run in min
     integer::bkp_modulo=0       ! Use modulo for backup file count
     integer::nfile=1            ! Number of file used per snapshot. Use -1 for nfile=ncpu

     ! Mesh parameters
     integer::levelmin=1         ! Full refinement up to levelmin
     integer::nlevelmax=1        ! Maximum number of level
     integer::ngridmax=0         ! Maximum number of grids
     integer::ncachemax=10000    ! Maximum number of cache lines
     integer::npartmax=0         ! Maximum number of DM particles
     integer::nstarmax=0         ! Maximum number of star particles
     integer::nsinkmax=0         ! Maximum number of sink particles
     integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
     real(dp)::boxlen=1.0D0      ! Cell size at level 0 (total box size)
     real(dp)::box_size=0.0D0    ! Box length of active domain along x direction
     integer::box_xmin,box_xmax  ! Min and max Cartesian keys at levelmin
     integer::box_ymin,box_ymax  ! Min and max Cartesian keys at levelmin
     integer::box_zmin,box_zmax  ! Min and max Cartesian keys at levelmin
          
     ! Poisson solver parameters
     real(dp)::epsilon=1.0D-4     ! Convergence criterion for Poisson solvers
     integer ::gravity_type=0     ! Type of force computation
     real(dp),dimension(1:10)::gravity_params=0.0 ! Gravity parameters
     integer :: cic_levelmax=0     ! Maximum level for CIC dark matter interpolation
     integer :: cg_levelmin=999    ! Min level for CG solver
     logical :: fast_solver = .false. ! Fast solver with MPI pre-fetch (memory intensive)

     ! Movie parameters
     integer::levelmax_frame=0
     integer::nw_frame=512 ! width of frame in pixels
     integer::nh_frame=512 ! height of frame in pixels
     integer::ivar_frame=1
     real(kind=8),dimension(1:20)::xcentre_frame=0d0
     real(kind=8),dimension(1:20)::ycentre_frame=0d0
     real(kind=8),dimension(1:20)::zcentre_frame=0d0
     real(kind=8),dimension(1:10)::deltax_frame=0d0
     real(kind=8),dimension(1:10)::deltay_frame=0d0
     real(kind=8),dimension(1:10)::deltaz_frame=0d0
     logical::movie=.false.
     logical::zoom_only=.false.
     integer::imovout=0    ! Increment for output times
     integer::imov=1       ! Initialize
     real(kind=8)::tendmov=0.
     real(kind=8)::aendmov=0.
     character(LEN=5)::proj_axis='z' ! x->x, y->y, projection along z
     integer,dimension(0:NVAR+2)::movie_vars=0
     character(len=5),dimension(0:NVAR+2)::movie_vars_txt=''
     
     ! Hydro solver parameters
     real(dp)::gamma=1.4d0
     real(dp)::courant_factor=0.5d0
     real(dp)::smallc=1.d-10
     real(dp)::smallr=1.d-10
     integer ::niter_riemann=10
     integer ::slope_type=1
     integer ::slope_mag_type=1
     real(dp)::difmag=0.0d0
     real(dp)::etamag=0.0d0
     real(dp),dimension(1:nener+1)::gamma_rad=1.33333333334d0
     logical ::induction=.false.
     logical ::entropy=.false.
     logical ::turb=.false.
     real(dp)::dual_energy=-1
     real(dp)::T2_fix=0d0
     character(LEN=10)::scheme='muscl'
     integer::riemann=0
     integer::riemann2d=0
     real(dp),dimension(1:3)::constant_gravity
     integer::inener,ientropy,imetal,iturb,ichem
     
     ! Physics parameters
     real(dp)::units_density=1.0 ! [g/cm^3]
     real(dp)::units_time=1.0    ! [seconds]
     real(dp)::units_length=1.0  ! [cm]

     ! Cosmological parameters (others are read from file)
     real(dp)::omega_b=0.0D0  ! Omega Baryon
     real(dp)::aexp_ini=10.   ! Starting expansion factor

     ! Refinement parameters for each level
     real(dp),dimension(1:MAXLEVEL)::m_refine = -1.0 ! Lagrangian threshold
     real(dp),dimension(1:MAXLEVEL)::r_refine = -1.0 ! Radius of refinement region
     real(dp),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
     real(dp),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
     real(dp),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
     real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
     real(dp),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
     real(dp),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
     real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0 ! Number of cells per Jeans length
     real(dp)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
     real(dp)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
     integer::ivar_refine=-1 ! Variable index for refinement
     logical::aexp_lock_refine=.false.
     logical::pic_lock_refine=.false.

     ! Refinement parameters for hydro
     integer ::interpol_var=0   ! Interpolated variables
     integer ::interpol_type=1  ! Interpolation scheme
     real(dp)::err_grad_d=-1.0  ! Density gradient
     real(dp)::err_grad_u=-1.0  ! Velocity gradient
     real(dp)::err_grad_p=-1.0  ! Pressure gradient
     real(dp)::floor_d=1.d-10   ! Density floor
     real(dp)::floor_u=1.d-10   ! Velocity floor
     real(dp)::floor_p=1.d-10   ! Pressure floor
     real(dp)::mass_sph=0.0D0   ! mass_sph
#ifdef MHD
     real(dp)::err_grad_b2=-1.0
     real(dp)::err_grad_A=-1.0
     real(dp)::err_grad_B=-1.0
     real(dp)::err_grad_C=-1.0
     real(dp)::floor_b2=1.d-10
     real(dp)::floor_A=1.d-10
     real(dp)::floor_B=1.d-10
     real(dp)::floor_C=1.d-10
#endif
#if NENER>0
     real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>5+NENER
     real(dp),dimension(1:NVAR-5-NENER)::err_grad_var=-1.0
#endif

     ! Initial condition regions parameters
     integer::nregion=0
     character(LEN=10),dimension(1:MAXREGION)::region_type='square'
     real(dp),dimension(1:MAXREGION)::x_center=0.
     real(dp),dimension(1:MAXREGION)::y_center=0.
     real(dp),dimension(1:MAXREGION)::z_center=0.
     real(dp),dimension(1:MAXREGION)::length_x=1.d10
     real(dp),dimension(1:MAXREGION)::length_y=1.d10
     real(dp),dimension(1:MAXREGION)::length_z=1.d10
     real(dp),dimension(1:MAXREGION)::exp_region=2.0

     ! Initial conditions hydro variables
     real(dp),dimension(1:MAXREGION)::d_region=0.
     real(dp),dimension(1:MAXREGION)::u_region=0.
     real(dp),dimension(1:MAXREGION)::v_region=0.
     real(dp),dimension(1:MAXREGION)::w_region=0.
     real(dp),dimension(1:MAXREGION)::p_region=0.
#if NENER>0
     real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>5+NENER
     real(dp),dimension(1:MAXREGION,1:NVAR-5-NENER)::var_region=0.0
#endif
#ifdef MHD
     real(dp),dimension(1:MAXREGION)::B_region=0.
     real(dp),dimension(1:MAXREGION)::C_region=0.
     real(dp)::A_ave=0.,B_ave=0.,C_ave=0.
#endif
     ! Initial condition files for each level
     character(LEN=20)::filetype='ascii'
     logical::multiple=.false.
     character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
     real(dp)::ic_scale_m=1.0d0

     ! Boundary conditions parameters
     logical,dimension(1:NDIM)::periodic=.true.
     integer::nbound=0
     logical::no_inflow=.false.
     integer,dimension(1:MAXBOUND)::bound_type=0
     integer,dimension(1:MAXBOUND)::bound_dir=0
     integer,dimension(1:MAXBOUND)::bound_shift=0
     integer,dimension(1:MAXBOUND)::bound_xmin=0
     integer,dimension(1:MAXBOUND)::bound_xmax=0
     integer,dimension(1:MAXBOUND)::bound_ymin=0
     integer,dimension(1:MAXBOUND)::bound_ymax=0
     integer,dimension(1:MAXBOUND)::bound_zmin=0
     integer,dimension(1:MAXBOUND)::bound_zmax=0
     real(dp),dimension(1:MAXBOUND)::d_bound=0
     real(dp),dimension(1:MAXBOUND)::u_bound=0
     real(dp),dimension(1:MAXBOUND)::v_bound=0
     real(dp),dimension(1:MAXBOUND)::w_bound=0
     real(dp),dimension(1:MAXBOUND)::p_bound=0
#if NENER>0
     real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0
#endif
#if NVAR>5+NENER
     real(dp),dimension(1:MAXBOUND,1:NVAR-5-NENER)::var_bound=0
#endif
#ifdef MHD
     real(dp),dimension(1:MAXBOUND)::B_bound=0
     real(dp),dimension(1:MAXBOUND)::C_bound=0
#endif
     ! Cooling parameters
     logical::cooling=.false.
     logical::cooling_ism=.false.
     logical::metal=.false.
     logical::isothermal=.false.
     logical::haardt_madau=.false.
     logical::self_shielding=.false.
     real(dp)::J21=0d0,a_spec=1d0,z_ave=0d0,z_reion=8.5d0
     integer::eos_type=1 ! 1=isothermal, 2=polytrope, 3=isothermal+polytrope
     real(dp)::eos_nH=1d50,eos_index=1d0,eos_T2=10d0
     real(dp)::T2max

     ! Star formation parameters
     real(dp)::T2_star=2e4
     real(dp)::n_star=0.1
     real(dp)::eps_star=0.01
     integer(kind=8),dimension(1:6)::seed
     real(dp)::m_star=1

     ! Supernovae feedback parameters
     real(dp)::M_SNII=10.
     real(dp)::E_SNII=1d51
     real(dp)::t_SNII=20.
     real(dp)::eta_SNII=0.1
     real(dp)::yield_SNII=0.1
     logical::thermal_feedback=.false.
     logical::mechanical_feedback=.false.
     
     ! Clump finder parameters
     logical::clump_finder=.false.
     logical::clump_info=.false.
     logical::output_clump=.false.
     logical::output_peak=.false.
     integer::rho_type_clump=1
     logical::output_peak_part=.false.
     logical::output_peak_star=.false.
     logical::output_peak_sink=.false.
     logical::output_halo_part=.false.
     logical::output_halo_star=.false.
     logical::output_halo_sink=.false.
     real(dp)::relevance_threshold=2
     real(dp)::density_threshold=-1
     real(dp)::saddle_threshold=-1
     real(dp)::mass_threshold=0

     ! Sink parameters
     integer::rho_type_sink=1
     logical::sink_descent=.false.
     real(dp)::fudge_descent=0.5d0

     ! Gadget initial conditions parameters
     character(len=flen)::ic_file, ic_format
     integer,dimension(1:6)::ic_skip_type=-1
     character(len=4)::ic_head_name = 'HEAD'
     character(len=4)::ic_pos_name = 'POS '
     character(len=4)::ic_vel_name = 'VEL '
     character(len=4)::ic_id_name = 'ID  '
     character(len=4)::ic_mass_name = 'MASS'
     character(len=4)::ic_u_name = 'U   '
     character(len=4)::ic_metal_name = 'Z   '
     character(len=4)::ic_age_name = 'AGE '
     real(dp)::gadget_scale_l = 3.085677581282D21 ! Gadget units in cgs
     real(dp)::gadget_scale_v = 1.0D5
     real(dp)::gadget_scale_m = 1.9891D43
     real(dp)::gadget_scale_t = 1.0D6*365*24*3600
     real(dp)::IG_rho = 1.0D-6
     real(dp)::IG_T2 = 1.0D7
     real(dp)::IG_metal = 0.01

  end type run_t
  
  type global_t

     ! MPI variables
     integer::ncpu, myid

     integer::iout=1             ! Increment for output times
     integer::ifout=1            ! Increment for output files
     integer::ifbkp=1            ! Increment for backup files

     logical::output_done=.false.                  ! Output just performed
     logical::init=.false.                         ! Set up or run
     integer::nstep=0                              ! Time step
     integer::nstep_coarse=0                       ! Coarse step
     integer::nstep_coarse_old=0                   ! Old coarse step
     integer::nflag,ncreate,nkill                  ! Refinements
     
     real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
     real(dp)::eint_tot=0.0D0                      ! Total internal energy
     real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
     real(dp)::epot_tot=0.0D0                      ! Total potential energy
     real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
     real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
     real(dp)::const=0.0D0                         ! Energy conservation
     real(dp)::aexp_old=1.0D0                      ! Old expansion factor
     real(dp)::rho_tot=0.0D0                       ! Mean density in the box
     real(dp)::t=0.0D0                             ! Time variable
     real(dp)::mass_tot=0.0D0                      ! Total gass mass
     real(dp)::mass_tot_0=0.0D0                    ! Initial total gas mass
     real(dp)::mass_star_tot=0.0D0                 ! Total mass in new stars
     real(dp)::mass_sink_tot=0.0D0                 ! Total mass in new sinks
     
     ! Level related arrays
     real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
     real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
     integer,dimension(1:MAXLEVEL)::isubcycle    ! Current subcycling step at each level

     ! Only one process can write at a time in an I/O group
     integer::IOGROUPSIZE=0           ! Main snapshot
     integer::IOGROUPSIZECONE=0       ! Lightcone
     integer::IOGROUPSIZEREP=0        ! Subfolder size
     logical::withoutmkdir=.false.    ! If true mkdir should be done before the run
     logical::print_when_io=.false.   ! If true print when IO
     logical::synchro_when_io=.false. ! If true synchronize when IO
     
     ! Lightcone parameters
     real(dp)::thetay_cone=12.5
     real(dp)::thetaz_cone=12.5
     real(dp)::zmax_cone=2.0
     
     ! Cosmology parameters
     real(dp)::boxlen_ini     ! Box size in h-1 Mpc
     real(dp)::omega_b=0.0D0  ! Omega Baryon
     real(dp)::omega_m=1.0D0  ! Omega Matter
     real(dp)::omega_l=0.0D0  ! Omega Lambda
     real(dp)::omega_k=0.0D0  ! Omega Curvature
     real(dp)::h0=1.0D0       ! Hubble constant in km/s/Mpc
     real(dp)::aexp=1.0D0     ! Current expansion factor
     real(dp)::hexp=0.0D0     ! Current Hubble parameter
     real(dp)::texp=0.0D0     ! Current proper time
     logical ::use_proper_time=.false.
     
     ! Executable identification
     CHARACTER(LEN=300)::builddate,buildcommand,patchdir
     CHARACTER(LEN=300)::gitrepo,gitbranch,githash
     
     ! Save namelist filename
     CHARACTER(LEN=300)::namelist_file
     
     ! Friedman model variables
     real(dp),dimension(0:n_frw)::aexp_frw,hexp_frw,tau_frw,t_frw
     
     ! Initial conditions parameters from grafic
     integer::nlevelmax_part
     real(dp)::aexp_ini=10.
     real(dp)::T2_start          ! Starting gas temperature     
     real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
     real(dp),dimension(1:MAXLEVEL)::vfact
     real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
     integer ,dimension(1:MAXLEVEL)::n1,n2,n3
          
     ! Minimum particle mass
     real(dp) :: mp_min=-1.0

     ! Minimum MG level
     integer :: levelmin_mg
     
     ! Multigrid safety switch
     logical, dimension(1:MAXLEVEL)::safe_mode=.false.
     
     ! Multipole coefficients
     !real(dp),dimension(1:ndim+1)::multipole
     type(multipole_t)::multipole
     
  end type global_t

  type mesh_t
     ! Level related arrays
     integer(kind=4),allocatable,dimension(:)::head      ! Starting index for each level
     integer(kind=4),allocatable,dimension(:)::tail      ! Final index for each level
     integer(kind=4),allocatable,dimension(:)::noct      ! Number of octs for each level
     integer(kind=4),allocatable,dimension(:)::noct_min  ! Min. number of octs across cpus
     integer(kind=4),allocatable,dimension(:)::noct_max  ! Max. number of octs across cpus
     integer(kind=8),allocatable,dimension(:)::noct_tot  ! Total number of octs across cpus
     integer(kind=4),allocatable,dimension(:)::ckey_max  ! Max. Cartesian key per level
     integer(kind=4),allocatable,dimension(:,:)::box_ckey_min  ! Min. Cartesian key per level for the domain
     integer(kind=4),allocatable,dimension(:,:)::box_ckey_max  ! Max. Cartesian key per level for the domain
     integer(kind=8),allocatable,dimension(:,:)::hkey_max ! Max. Hilbert key
     integer(kind=4),allocatable,dimension(:)::head_cache ! Starting index in the cache for each level
     integer(kind=4),allocatable,dimension(:)::tail_cache ! Final index in the cache for each level
     integer(kind=4)::noct_used,noct_used_max,noct_used_tot ! Total used octs in local memory
     integer(kind=4)::nx,ny,nz                   ! Size of mesh at levelmin
     real(kind=8),allocatable,dimension(:)::skip ! Coordinates of lower left corner of the box

     ! Persistent array for the AMR grid
!     type(oct),dimension(:),allocatable::grid
     type(oct),dimension(:),pointer::grid
     type(hash_table)::grid_dict   ! Oct hash table
     
     ! Arrays for the MG solver
     type(hash_table)::mg_dict     ! MG hash table
     integer(kind=4),allocatable,dimension(:)::head_mg ! Starting index for each level
     integer(kind=4),allocatable,dimension(:)::tail_mg ! Final index for each level
     integer(kind=4),allocatable,dimension(:)::noct_mg ! Number of octs for each level
     integer(kind=4)::ifree_mg ! Starting index in free memory

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
     
     ! Peano-Hilbert key boundaries for cpu domains
     type(domain_t),pointer,dimension(:)::domain,domain_mg
     type(domain_t),pointer,dimension(:)::domain_hilbert

     ! Hydro kernel workspace
     type(hydro_workspace_t)::hydro_w
     
  end type mesh_t

  ! Peano-Hilbert key boundaries for cpu domains
  type(domain_t), allocatable, target, dimension(:)::domain,domain_mg
  type(domain_t), pointer,             dimension(:)::domain_hilbert

contains

  subroutine print_run_parameters(run_p)
    class(run_t)::run_p
    type(run_t)::run_params
    namelist /run_parameters/ run_params
    run_params = run_p
    write(*,NML=run_parameters)
  end subroutine print_run_parameters

end module amr_commons


module params_module

contains

subroutine m_read_params(pst)
  use amr_parameters
  use hydro_parameters
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
  implicit none
  type(pst_t)::pst

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,levelmax
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  integer(kind=8)::nstartot=0
  integer(kind=8)::nsinktot=0
  integer(kind=8)::ntreetot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok

  !--------------------------------------------------
  ! Namelist variables
  !--------------------------------------------------

  ! Run control
  logical::cosmo   =.false.   ! Cosmology activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::hydro   =.false.   ! Hydro activated
  logical::star    =.false.   ! Stars and star formation activated
  logical::sink    =.false.   ! Sinks and sink formation activated
  logical::merger_tree=.false. ! Merger tree particles activated
  logical::orphan  =.false.   ! Orphan particles activated
  logical::verbose =.false.   ! Write everything
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::nremap=10          ! Particle load balancing frequency (0: never)

  ! Maximum number of allocatable particles
  integer::npartmax=0 
  integer::nstarmax=0
  integer::nsinkmax=0
  integer::ntreemax=0

  ! Number of superoct levels
  integer::nsuperoct=0
  
  ! MPI domain overloading
  integer::overload=1

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer::ncachemax=10000    ! Maximum number of cache lines
  real(dp)::boxlen=1.0D0      ! Cell sixe at level 0
  real(dp)::box_size=0.0D0    ! Box length in active domain along x direction
  integer::box_xmin=0         ! Min. Cartesian key for the box at levelmin in x direction
  integer::box_xmax=0         ! Max. Cartesian key for the box at levelmin in x direction
  integer::box_ymin=0         ! Min. Cartesian key for the box at levelmin in y direction
  integer::box_ymax=0         ! Max. Cartesian key for the box at levelmin in y direction
  integer::box_zmin=0         ! Min. Cartesian key for the box at levelmin in z direction
  integer::box_zmax=0         ! Max. Cartesian key for the box at levelmin in z direction

  ! Output parameters
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)
  logical::gadget_output=.false. ! Output in gadget format
  real(kind=8)::bkp_time_hrs=2   ! Backup file frequency in hours
  real(kind=8)::run_time_hrs=0   ! Estimated run time in hrs
  real(kind=8)::bkp_last_min=10  ! Backup file before the end of run in min
  integer::bkp_modulo=0       ! Use modulo for backup file count
  integer::nfile=1            ! Number of file per snapshot. Use -1 for nfile=ncpu

  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Movie
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tendmov=0.,aendmov=0.
  logical::movie=.false.
  logical::zoom_only=.false.
  integer::nw_frame=512 ! width of frame in pixels
  integer::nh_frame=512 ! height of frame in pixels
  integer::levelmax_frame=0
  integer::ivar_frame=1
  real(kind=8),dimension(1:20)::xcentre_frame=0d0
  real(kind=8),dimension(1:20)::ycentre_frame=0d0
  real(kind=8),dimension(1:20)::zcentre_frame=0d0
  real(kind=8),dimension(1:10)::deltax_frame=0d0
  real(kind=8),dimension(1:10)::deltay_frame=0d0
  real(kind=8),dimension(1:10)::deltaz_frame=0d0
  character(LEN=5)::proj_axis='z' ! x->x, y->y, projection along z
  integer,dimension(0:NVAR+2)::movie_vars=0
  character(len=5),dimension(0:NVAR+2)::movie_vars_txt=''

  ! Refinement parameters for each level
  integer ,dimension(1:MAXLEVEL)::nexpand = 1 ! Number of mesh expansion
  integer ,dimension(1:MAXLEVEL)::nsubcycle = 2 ! Subcycling at each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1.0 ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1.0 ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1 ! Variable index for refinement
  real(dp)::aexp_lock_refine=-1.0
  logical::pic_lock_refine=.false.

  ! Default units
  real(dp)::units_density=1.0 ! [g/cm^3]
  real(dp)::units_time=1.0    ! [seconds]
  real(dp)::units_length=1.0  ! [cm]

  ! Initial conditions parameters from grafic
  real(dp)::aexp_ini=10.
  real(dp)::omega_b=0.045

  ! Initial condition regions parameters
  integer::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)::x_center=0.
  real(dp),dimension(1:MAXREGION)::y_center=0.
  real(dp),dimension(1:MAXREGION)::z_center=0.
  real(dp),dimension(1:MAXREGION)::length_x=1.E10
  real(dp),dimension(1:MAXREGION)::length_y=1.E10
  real(dp),dimension(1:MAXREGION)::length_z=1.E10
  real(dp),dimension(1:MAXREGION)::exp_region=2.0

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=20)::filetype='ascii'
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  real(dp)::ic_scale_m=1.0d0

  ! Refinement parameters for hydro
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
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

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
  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  integer ::slope_mag_type=-1
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::etamag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
  character(LEN=10)::riemann2d='none'
  logical ::induction=.false.
  logical ::entropy=.false.
  logical ::turb=.false.
  real(dp)::dual_energy=-1
  real(dp)::T2_fix=0d0
  real(dp),dimension(1:3)::constant_gravity=0.0d0

  ! Non-thernal energies and passive scalars index
  integer ::inener,ientropy,imetal,iturb,ichem

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Poisson solver parameters
  real(dp)::epsilon=1.0D-4 ! Convergence criterion
  real(dp),dimension(1:10)::gravity_params=0.0 ! Gravity parameters
  integer :: gravity_type=0 ! Type of gravity calculations (see user guide)
  integer :: cic_levelmax=0 ! Maximum level for CIC dark matter interpolation
  integer :: cg_levelmin=999   ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  logical :: fast_solver = .false.   ! Fast solver with MPI pre-fetch (memory intensive)
  integer :: part_mass_deposition_scheme=1     ! part mass deposition schemes (CIC 1, TSC 2, PCS 3)
  integer :: part_force_interpolation_scheme=1 ! part force interpolation schemes (CIC 1, TSC 2, PCS 3)
  integer :: star_mass_deposition_scheme=1     ! star mass deposition schemes
  integer :: star_force_interpolation_scheme=1 ! star force interpolation schemes
  integer :: sink_mass_deposition_scheme=1     ! sink mass deposition schemes
  integer :: sink_force_interpolation_scheme=1 ! sink force interpolation schemes
  integer :: tree_mass_deposition_scheme=1     ! tree mass deposition schemes
  integer :: tree_force_interpolation_scheme=1 ! tree force interpolation schemes

  ! Boundary conditions parameters
  integer::nbound=0
  logical::no_inflow=.false.
  logical,dimension(1:NDIM)::periodic=.true.
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
  real(dp),dimension(1:MAXBOUND)::p_bound=0
  real(dp),dimension(1:MAXBOUND)::u_bound=0
  real(dp),dimension(1:MAXBOUND)::v_bound=0
  real(dp),dimension(1:MAXBOUND)::w_bound=0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0
#endif
#if NVAR>5+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-5-NENER)::var_bound=0
#endif

  ! Cooling parameters
  logical::cooling=.false.
  logical::cooling_ism=.false.
  logical::metal=.false.
  logical::isothermal=.false. ! Force temperature to eos value
  logical::haardt_madau=.false.
  logical::self_shielding=.false.
  real(dp)::J21=0d0,a_spec=1d0,z_ave=0d0,z_reion=8.5d0
  integer::eos_type=1 ! 1=isothermal, 2=polytrope, 3=isothermal+polytrope
  real(dp)::eos_nH=1d50,eos_index=1d0,eos_T2=10d0
  real(dp)::T2max=1d50

  ! Star formation parameters
  integer::sf_model=1
  real(dp)::T2_star=2e4
  real(dp)::n_star=0.1
  real(dp)::eps_star=0.01
  integer(kind=8),dimension(1:6)::seed=(/123,456,789,1,1,1/)
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
  logical::output_peak_grid=.false.
  logical::output_peak_part=.false.
  logical::output_peak_star=.false.
  logical::output_peak_sink=.false.
  logical::output_peak_tree=.false.
  integer::rho_type_clump=1 ! 1: DM, 2: stars, 3: sinks, 4: gas
  real(dp)::relevance_threshold=2
  real(dp)::density_threshold=-1
  real(dp)::saddle_threshold=-1
  real(dp)::mass_threshold=0
  real(dp)::purity_threshold=-1
  real(dp)::fraction_threshold=0.1d0

  ! Sink parameters
  integer::rho_type_sink=1
  logical::sink_descent=.false.
  real(dp)::fudge_descent=0.5d0
  real(dp)::sink_relevance_threshold=2
  real(dp)::sink_density_threshold=-1
  real(dp)::sink_saddle_threshold=-1
  real(dp)::sink_mass_threshold=0
  real(dp)::sink_purity_threshold=-1
  real(dp)::sink_fraction_threshold=2d0
  logical::sink_form=.false.
  logical::sink_refine=.false.
  logical::sink_dump=.false.

  ! Black hole parameters
  integer::accretion_type = 0 ! 0: None, 1: Bondi
  character(len=10)::accretion_method = 'mass' ! Whether to mass-weigh the accretion 
  real(dp)::acc_sink_boost = 1.0d0 ! Boost for bondi accretion
  logical::bondi_use_vrel = .true. ! Whether to use the relative sink velocity for BHL accretion
  real(dp)::eddington_cap = -1 ! Factor of Eddington rate to cap accretion at
  integer::sink_b_spline_order = 4 ! Order of B-spline interpolation used for sink accretion and dynamics
  logical::verbose_sink = .false. ! Whether to print verbose statements for sink particles
  logical::bondi_use_gas_mass = .true. ! Whether to include the local gas mass in the Bondi calculation
  logical::use_local_bondi_rate = .false. ! Switch to average after (true) or before (false) computing the Bondi rate

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

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  ! Global run parameter
  namelist/run_params/cosmo,pic,poisson,hydro,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap &
       & ,static,geom,overload,nsuperoct
  ! Output parameters
  namelist/output_params/foutput,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout,gadget_output &
       & ,run_time_hrs,bkp_time_hrs,bkp_last_min,bkp_modulo,nfile
  ! AMR grid basic parameters
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen,box_size &
       & ,box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax
  ! Poisson solver parameters
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax,fast_solver,part_mass_deposition_scheme &
       & ,part_force_interpolation_scheme,star_mass_deposition_scheme &
       & ,star_force_interpolation_scheme,sink_mass_deposition_scheme &
       & ,sink_force_interpolation_scheme,tree_mass_deposition_scheme &
       & ,tree_force_interpolation_scheme 
  ! Movies parameters
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only &
       & ,imovout,imov,tendmov,aendmov,proj_axis,movie_vars,movie_vars_txt
  ! Initial conditions parameters
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini,omega_b &
       & ,length_x,length_y,length_z,exp_region &
       & ,ic_scale_m &
#if NENER>0
       & ,prad_region &
#endif
#if NVAR>5+NENER
       & ,var_region &
#endif
#ifdef MHD
       & ,B_region,C_region,A_ave,B_ave,C_ave &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  ! Hydro solver parameters
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,slope_mag_type,difmag,etamag,gamma_rad &
       & ,dual_energy,T2_fix,induction,entropy,turb,scheme,riemann,riemann2d,constant_gravity
  ! Grid refinement parameters
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
#ifdef MHD
       & ,err_grad_b2,err_grad_A,err_grad_B,err_grad_C &
       & ,floor_b2,floor_A,floor_B,floor_C &
#endif
#if NENER>0
       & ,err_grad_prad &
#endif
#if NVAR>5+NENER
       & ,err_grad_var &
#endif
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type &
       & ,aexp_lock_refine,pic_lock_refine,sink_refine
  ! Units parameters
  namelist/units_params/units_density,units_time,units_length
  ! Boundary conditions parameters
  namelist/boundary_params/periodic,nbound,bound_type,bound_dir,bound_shift &
       & ,bound_xmin,bound_xmax,bound_ymin,bound_ymax,bound_zmin,bound_zmax &
#if NENER>0
       & ,prad_bound &
#endif
#if NVAR>5+NENER
       & ,var_bound &
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound
  ! Cooling / basic chemistry parameters
  namelist/cooling_params/cooling,metal,isothermal,haardt_madau,J21 &
       & ,eos_type,eos_nH,eos_index,eos_T2 &
       & ,a_spec,self_shielding,z_ave,z_reion,T2max,cooling_ism
  ! Star particles and star formation recipe
  namelist/star_params/star,nstarmax,nstartot,T2_star,n_star,eps_star,seed,m_star,sf_model
  ! Sink particles and black hole parameters
  namelist/sink_params/sink,nsinkmax,nsinktot,rho_type_sink,sink_descent,fudge_descent &
       & ,sink_relevance_threshold,sink_density_threshold,sink_saddle_threshold &
       & ,sink_mass_threshold,sink_purity_threshold,sink_fraction_threshold &
       & ,accretion_type,acc_sink_boost,bondi_use_vrel,accretion_method &
       & ,eddington_cap,sink_form,sink_b_spline_order,verbose_sink,bondi_use_gas_mass &
       & ,use_local_bondi_rate,sink_dump
  ! Supernovae feedback parameters
  namelist/feedback_params/M_SNII,E_SNII,t_SNII,eta_SNII,yield_SNII,thermal_feedback,mechanical_feedback
  ! Clump finder parameters
  namelist/clump_params/clump_finder,clump_info &
       & ,output_clump,output_peak_grid,output_peak_part,output_peak_star,output_peak_sink,output_peak_tree &
       & ,relevance_threshold,density_threshold,saddle_threshold &
       & ,mass_threshold,purity_threshold,fraction_threshold &
       & ,merger_tree,orphan,ntreemax,ntreetot,rho_type_clump
  ! Gadget initial conditions parameters
  namelist/gadget_params/ic_file,ic_format,IG_rho,IG_T2,IG_metal &
       & ,ic_head_name,ic_pos_name,ic_vel_name,ic_id_name,ic_mass_name &
       & ,ic_u_name,ic_metal_name,ic_age_name &
       & ,gadget_scale_l, gadget_scale_v, gadget_scale_m ,gadget_scale_t &
       & ,ic_skip_type

  associate(s=>pst%s)

  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (Princeton University)       '
  write(*,*)'        (c) CEA 1999-2007, UZH 2008-2021, PU 2022-2023         '
  write(*,*)' '

  ! Check nvar is not too small
  write(*,'(" Using solver = hydro with nvar = ",I2," and ndim = ",I1)')nvar,ndim
  if(nvar<5)then
     write(*,*)'You should have: nvar>=5'
     write(*,'(" Please recompile with -DNVAR=5")')
     call mdl_abort(s%mdl)
  endif
  
  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = command_argument_count()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call mdl_abort(s%mdl)
  END IF
  CALL getarg(1,infile)

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'File '//TRIM(infile)//' does not exist'
     call mdl_abort(s%mdl)
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (narg==2) then
    CALL getarg(2,cmdarg)
    read(cmdarg,*) nrestart
  endif

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=max(count(aout>0.and.aout<=1),count(tout>0))
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  if(nfile==-1)nfile=s%g%ncpu
  nfile=max(nfile,1)
  nfile=min(nfile,s%g%ncpu)

  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nsuperoct=MIN(nsuperoct,5)
  nml_ok=.true.
  if(levelmin<1)then
     write(*,*)'Error in the namelist:'
     write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     write(*,*)'Error in the namelist:'
     write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        write(*,*)'Error in the namelist:'
        write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=int(ngridtot/int(s%g%ncpu,kind=8),kind=4)
     endif
  end if
  if(npartmax==0)then
     npartmax=int(nparttot/int(s%g%ncpu,kind=8),kind=4)
  endif
  if(nstarmax==0)then
     nstarmax=int(nstartot/int(s%g%ncpu,kind=8),kind=4)
  endif
  if(nsinkmax==0)then
     nsinkmax=int(nsinktot/int(s%g%ncpu,kind=8),kind=4)
  endif
  if(ntreemax==0)then
     ntreemax=int(ntreetot/int(s%g%ncpu,kind=8),kind=4)
     if(ntreemax==0)then
        ntreemax=int(npartmax/100)
     endif
  endif
#ifdef HYDRO
  if(.not. hydro)then
     write(*,*)'You are not using the hydro solver but'
     write(*,*)'the code was compiled with HYDRO=1'
     write(*,*)'This might not be optimal but I am still running.'
  endif
#else
  if(hydro)then
     write(*,*)'You are using the hydro solver but'
     write(*,*)'the code was compiled with HYDRO=0'
     write(*,*)'Please recompile with HYDRO=1'
     call mdl_abort(s%mdl)
  endif  
#endif
#ifdef GRAV
  if(.not. poisson)then
     write(*,*)'You are not using the poisson solver but'
     write(*,*)'the code was compiled with GRAV=1'
     write(*,*)'This might not be optimal but I am still running.'
  endif
#else
  if(poisson)then
     write(*,*)'You are using the poisson solver but'
     write(*,*)'the code was compiled with GRAV=0'
     write(*,*)'Please recompile with GRAV=1'
     call mdl_abort(s%mdl)
  endif
#endif

  !----------------------------
  ! Read hydro parameters 
  !----------------------------
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call mdl_abort(s%mdl)
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=units_params,END=105)
105 continue
  rewind(1)
  read(1,NML=boundary_params,END=106)
106 continue
  rewind(1)
  read(1,NML=cooling_params,END=107)
107 continue
  rewind(1)
  read(1,NML=star_params,END=108)
108 continue
  rewind(1)
  read(1,NML=feedback_params,END=109)
109 continue
  rewind(1)
  read(1,NML=clump_params,END=110)
110 continue
  rewind(1)
  read(1,NML=gadget_params,END=111)
111 continue
  rewind(1)
  read(1,NML=sink_params,END=112)
112 continue
  close(1)

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call mdl_abort(s%mdl)
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call mdl_abort(s%mdl)
  end if

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
     jeans_refine(i)=-1.0
  end do

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
  if(nvar<5+nener)then
     write(*,*)'Error: non-thermal energy need nvar >= 5+nener'
     write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check for metals
  !--------------------------------------------------
  if(metal.and.nvar<6)then
     write(*,*)'Error: metal=.true. need nvar >= 6'
     write(*,*)'Modify NVAR and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for entropy
  !--------------------------------------------------
  if(entropy.and.nvar<6)then
     write(*,*)'Error: entropy=.true. need nvar >= 6'
     write(*,*)'Modify NVAR and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Compute indices for passive scalars
  ! and non-thermal energies
  !--------------------------------------------------
  inener=6
  ientropy=inener+nener
  imetal=ientropy
  if(entropy)imetal=ientropy+1
  iturb=imetal
  if(metal)then
     iturb=imetal+1
  endif
  ichem=iturb
  if(turb)then
     ichem=iturb+1
  endif
  if(hydro.and.(nvar>5)) then
     write(*,'(A50)')"__________________________________________________"
     write(*,*) 'Hydro var extra indices:'
#if NENER>0
                 write(*,*) '   inener   = ',inener
#endif
     if(entropy) write(*,*) '   ientropy = ',ientropy
     if(metal)   write(*,*) '   imetal   = ',imetal
     if(turb)    write(*,*) '   iturb    = ',iturb
     if(ichem.LE.nvar)then
                 write(*,*) '   ichem    = ',ichem
     endif
     write(*,'(A50)')"__________________________________________________"
  endif

  if(.not. nml_ok)then
     write(*,*)'Too many errors in the namelist'
     write(*,*)'Aborting...'
     call mdl_abort(s%mdl)
  end if

  ! Fill in all run parameters in corresponding structure

  s%r%cosmo=cosmo
  s%r%pic=pic
  s%r%poisson=poisson
  s%r%hydro=hydro
  s%r%star=star
  s%r%sink=sink
  s%r%tree=merger_tree
  s%r%orphan=orphan
  s%r%verbose=verbose
  s%r%debug=debug
  s%r%nrestart=nrestart
  s%r%ncontrol=ncontrol
  s%r%nstepmax=nstepmax
  s%r%nsubcycle=nsubcycle
  s%r%nremap=nremap
  s%r%static=static
  s%r%geom=geom
  s%r%overload=overload
  s%r%nsuperoct=nsuperoct

  s%r%noutput=noutput
  s%r%foutput=foutput
  s%r%aout=aout
  s%r%tout=tout
  s%r%output_mode=output_mode
  s%r%gadget_output=gadget_output
  s%r%run_time_hrs=run_time_hrs
  s%r%bkp_time_hrs=bkp_time_hrs
  s%r%bkp_last_min=bkp_last_min
  s%r%bkp_modulo=bkp_modulo
  s%r%nfile=nfile

  s%r%levelmin=levelmin
  s%r%nlevelmax=nlevelmax
  s%r%ngridmax=ngridmax
  s%r%ncachemax=ncachemax
  s%r%npartmax=npartmax
  s%r%nstarmax=nstarmax
  s%r%nsinkmax=nsinkmax
  s%r%ntreemax=ntreemax
  s%r%nexpand=nexpand
  s%r%boxlen=boxlen
  s%r%box_size=box_size
  s%r%box_xmin=box_xmin
  s%r%box_xmax=box_xmax
  s%r%box_ymin=box_ymin
  s%r%box_ymax=box_ymax
  s%r%box_zmin=box_zmin
  s%r%box_zmax=box_zmax

  s%r%epsilon=epsilon
  s%r%gravity_type=gravity_type
  s%r%gravity_params=gravity_params
  s%r%cic_levelmax=cic_levelmax
  s%r%cg_levelmin=cg_levelmin
  s%r%fast_solver=fast_solver
  s%r%part_mass_deposition_scheme=part_mass_deposition_scheme
  s%r%part_force_interpolation_scheme=part_force_interpolation_scheme
  s%r%star_mass_deposition_scheme=star_mass_deposition_scheme
  s%r%star_force_interpolation_scheme=star_force_interpolation_scheme
  s%r%sink_mass_deposition_scheme=sink_mass_deposition_scheme
  s%r%sink_force_interpolation_scheme=sink_force_interpolation_scheme
  s%r%tree_mass_deposition_scheme=tree_mass_deposition_scheme
  s%r%tree_force_interpolation_scheme=tree_force_interpolation_scheme

  s%r%nw_frame=nw_frame
  s%r%nh_frame=nh_frame
  s%r%levelmax_frame=levelmax_frame
  s%r%ivar_frame=ivar_frame
  s%r%xcentre_frame=xcentre_frame
  s%r%ycentre_frame=ycentre_frame
  s%r%zcentre_frame=zcentre_frame
  s%r%deltax_frame=deltax_frame
  s%r%deltay_frame=deltay_frame
  s%r%deltaz_frame=deltaz_frame
  s%r%movie=movie
  s%r%zoom_only=zoom_only
  s%r%imovout=imovout
  s%r%imov=imov
  s%r%tendmov=tendmov
  s%r%aendmov=aendmov
  s%r%proj_axis=proj_axis
  s%r%movie_vars_txt=movie_vars_txt
  if(s%r%movie)call set_movie_vars(s%r)

  s%r%gamma=gamma
  s%r%courant_factor=courant_factor
  s%r%smallc=smallc
  s%r%smallr=smallr
  s%r%niter_riemann=niter_riemann
  s%r%slope_type=slope_type
  if(slope_mag_type<0)then
     s%r%slope_mag_type=slope_type
  else
     s%r%slope_mag_type=slope_mag_type
  endif
  s%r%difmag=difmag
  s%r%etamag=etamag
  s%r%gamma_rad=gamma_rad(1:nener+1)
  s%r%dual_energy=dual_energy
  s%r%T2_fix=T2_fix
  s%r%induction=induction
  s%r%entropy=entropy
  s%r%turb=turb
  s%r%inener=inener
  s%r%ientropy=ientropy
  s%r%imetal=imetal
  s%r%ichem=ichem
  s%r%iturb=iturb
  s%r%scheme=scheme
  s%r%constant_gravity=constant_gravity
#ifndef MHD
  if(riemann=='llf')s%r%riemann=solver_llf
  if(riemann=='hll')s%r%riemann=solver_hll
  if(riemann=='hllc')s%r%riemann=solver_hllc
#endif
#ifdef MHD
  if(riemann=='llf')s%r%riemann=solver_llf
  if(riemann=='hll')s%r%riemann=solver_hll
  if(riemann=='hlld')s%r%riemann=solver_hlld
  if(riemann=='roe')s%r%riemann=solver_roe
  if(riemann=='upwind')s%r%riemann=solver_upwind
#endif
#ifdef MHD
  if(riemann2d=='none')riemann2d=riemann
  if(riemann2d=='llf')s%r%riemann2d=solver2d_llf
  if(riemann2d=='hll')s%r%riemann2d=solver2d_hllf
  if(riemann2d=='hllf')s%r%riemann2d=solver2d_hllf
  if(riemann2d=='hlla')s%r%riemann2d=solver2d_hlla
  if(riemann2d=='hlld')s%r%riemann2d=solver2d_hlld
  if(riemann2d=='roe')s%r%riemann2d=solver2d_roe
  if(riemann2d=='upwind')s%r%riemann2d=solver2d_upwind
#endif

  s%r%units_density=units_density
  s%r%units_time=units_time
  s%r%units_length=units_length

  s%r%m_refine=m_refine
  s%r%r_refine=r_refine
  s%r%x_refine=x_refine
  s%r%y_refine=y_refine
  s%r%z_refine=z_refine
  s%r%exp_refine=exp_refine
  s%r%a_refine=a_refine
  s%r%b_refine=b_refine
  s%r%jeans_refine=jeans_refine
  s%r%var_cut_refine=var_cut_refine
  s%r%mass_cut_refine=mass_cut_refine
  s%r%ivar_refine=ivar_refine
  s%r%aexp_lock_refine=aexp_lock_refine
  s%r%pic_lock_refine=pic_lock_refine

  s%r%interpol_var=interpol_var
  s%r%interpol_type=interpol_type
  s%r%err_grad_d=err_grad_d
  s%r%err_grad_u=err_grad_u
  s%r%err_grad_p=err_grad_p
  s%r%floor_d=floor_d
  s%r%floor_u=floor_u
  s%r%floor_p=floor_p
  s%r%mass_sph=mass_sph
#ifdef MHD
  s%r%err_grad_b2=err_grad_b2
  s%r%err_grad_A=err_grad_A
  s%r%err_grad_B=err_grad_B
  s%r%err_grad_C=err_grad_C
  s%r%floor_b2=floor_b2
  s%r%floor_A=floor_A
  s%r%floor_B=floor_B
  s%r%floor_C=floor_C
#endif
#if NENER>0
  s%r%err_grad_prad=err_grad_prad
#endif
#if NVAR>5+NENER
  s%r%err_grad_var=err_grad_var
#endif

  if(nrestart>0)filetype='restart'
  s%r%filetype=filetype
  s%r%initfile=initfile
  s%r%multiple=multiple
  s%r%aexp_ini=aexp_ini
  s%r%omega_b=omega_b
  s%r%ic_scale_m=ic_scale_m

  s%r%nregion=nregion
  s%r%region_type=region_type
  s%r%x_center=x_center
  s%r%y_center=y_center
  s%r%z_center=z_center
  s%r%length_x=length_x
  s%r%length_y=length_y
  s%r%length_z=length_z
  s%r%exp_region=exp_region
  s%r%d_region=d_region
  s%r%u_region=u_region
  s%r%v_region=v_region
  s%r%w_region=w_region
  s%r%p_region=p_region
#if NENER>0
  s%r%prad_region=prad_region
#endif
#if NVAR>5+NENER
  s%r%var_region=var_region
#endif
#ifdef MHD
  s%r%B_region=B_region
  s%r%C_region=C_region
  s%r%A_ave=A_ave
  s%r%B_ave=B_ave
  s%r%C_ave=C_ave
#endif

  s%r%periodic=periodic
  s%r%nbound=nbound
  s%r%no_inflow=no_inflow
  s%r%bound_dir=bound_dir
  s%r%bound_type=bound_type
  s%r%bound_shift=bound_shift
  s%r%bound_xmin=bound_xmin
  s%r%bound_xmax=bound_xmax
  s%r%bound_ymin=bound_ymin
  s%r%bound_ymax=bound_ymax
  s%r%bound_zmin=bound_zmin
  s%r%bound_zmax=bound_zmax
  s%r%d_bound=d_bound
  s%r%u_bound=u_bound
  s%r%v_bound=v_bound
  s%r%w_bound=w_bound
  s%r%p_bound=p_bound
#if NENER>0
  s%r%prad_bound=prad_bound
#endif
#if NVAR>5+NENER
  s%r%var_bound=var_bound
#endif

  s%r%cooling=cooling
  s%r%cooling_ism=cooling_ism
  s%r%metal=metal
  s%r%isothermal=isothermal
  s%r%haardt_madau=haardt_madau
  s%r%self_shielding=self_shielding
  s%r%J21=J21
  s%r%a_spec=a_spec
  s%r%z_ave=z_ave
  s%r%z_reion=z_reion
  s%r%eos_type=eos_type
  s%r%eos_nH=eos_nH
  s%r%eos_index=eos_index
  s%r%eos_T2=eos_T2
  s%r%T2max=T2max

  s%r%T2_star=T2_star
  s%r%n_star=n_star
  s%r%eps_star=eps_star
  s%r%seed=seed
  s%r%m_star=m_star
  s%r%sf_model=sf_model

  s%r%M_SNII=M_SNII
  s%r%E_SNII=E_SNII
  s%r%t_SNII=t_SNII
  s%r%eta_SNII=eta_SNII
  s%r%yield_SNII=yield_SNII
  s%r%thermal_feedback=thermal_feedback
  s%r%mechanical_feedback=mechanical_feedback

  s%r%clump_finder=clump_finder
  s%r%clump_info=clump_info
  s%r%output_clump=output_clump
  s%r%output_peak_grid=output_peak_grid
  s%r%output_peak_part=output_peak_part
  s%r%output_peak_star=output_peak_star
  s%r%output_peak_sink=output_peak_sink
  s%r%output_peak_tree=output_peak_tree
  s%r%relevance_threshold=relevance_threshold
  s%r%density_threshold=density_threshold
  s%r%saddle_threshold=saddle_threshold
  s%r%mass_threshold=mass_threshold
  s%r%purity_threshold=purity_threshold
  s%r%fraction_threshold=fraction_threshold
  s%r%rho_type_clump=rho_type_clump

  s%r%rho_type_sink=rho_type_sink
  s%r%sink_descent=sink_descent
  s%r%fudge_descent=fudge_descent
  s%r%sink_relevance_threshold=sink_relevance_threshold
  s%r%sink_density_threshold=sink_density_threshold
  s%r%sink_saddle_threshold=sink_saddle_threshold
  s%r%sink_mass_threshold=sink_mass_threshold
  s%r%sink_purity_threshold=sink_purity_threshold
  s%r%sink_fraction_threshold=sink_fraction_threshold

  s%r%accretion_type = accretion_type
  s%r%accretion_method = accretion_method
  s%r%acc_sink_boost = acc_sink_boost
  s%r%bondi_use_vrel = bondi_use_vrel
  s%r%eddington_cap = eddington_cap
  s%r%sink_form = sink_form
  s%r%sink_refine = sink_refine
  s%r%sink_dump = sink_dump
  s%r%sink_b_spline_order = sink_b_spline_order
  s%r%verbose_sink = verbose_sink
  s%r%bondi_use_gas_mass = bondi_use_gas_mass
  s%r%use_local_bondi_rate = use_local_bondi_rate

  s%r%ic_file=ic_file
  s%r%ic_format=ic_format
  s%r%IG_rho=IG_rho
  s%r%IG_T2=IG_T2
  s%r%IG_metal=IG_metal
  s%r%ic_head_name=ic_head_name
  s%r%ic_pos_name=ic_pos_name
  s%r%ic_vel_name=ic_vel_name
  s%r%ic_id_name=ic_id_name
  s%r%ic_mass_name=ic_mass_name
  s%r%ic_u_name=ic_u_name
  s%r%ic_metal_name=ic_metal_name
  s%r%ic_age_name=ic_age_name
  s%r%gadget_scale_l=gadget_scale_l
  s%r%gadget_scale_v=gadget_scale_v
  s%r%gadget_scale_m=gadget_scale_m
  s%r%gadget_scale_t=gadget_scale_t
  s%r%ic_skip_type=ic_skip_type

  ! Broadcast parameters to all CPUs.
  call m_broadcast_params(pst)

  end associate

end subroutine m_read_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_params(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------

  ! Broadcast parameters to all CPUs.
  call r_broadcast_params(pst,pst%s%r,storage_size(pst%s%r)/32)
end subroutine m_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_params(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use amr_commons, only: run_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(run_t)::input

  integer::rID

  if(pst%nLower>0)then
    rID = mdl_send_request(pst%s%mdl,MDL_BCAST_PARAMS,pst%iUpper+1,input_size,0,input)
    call r_broadcast_params(pst%pLower,input,input_size)
    call mdl_get_reply(pst%s%mdl,rID,0)
  else
    pst%s%r=input
  endif

end subroutine r_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_global(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------
  integer::input_size
  integer,dimension(1:1)::dummy
  integer,dimension(:),allocatable::input_array

  ! Broadcast parameters to all CPUs.
  input_size=storage_size(pst%s%g)/32
  allocate(input_array(1:input_size))
  input_array=transfer(pst%s%g,input_array)
  call r_broadcast_global(pst,input_array,input_size,dummy,0)
  deallocate(input_array)

end subroutine m_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_global(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BCAST_GLOBAL,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_global(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     pst%s%g=transfer(input_array,pst%s%g)
     pst%s%g%myid=mdl_self(pst%s%mdl)
  endif

end subroutine r_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module params_module

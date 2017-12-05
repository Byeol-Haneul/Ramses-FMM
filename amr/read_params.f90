subroutine m_read_params(r,g,m,p,mdl)
  use amr_parameters
  use hydro_parameters
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,levelmax
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok

  !--------------------------------------------------
  ! Namelist variables
  !--------------------------------------------------
  ! Maximum number of allocatable particles
  integer::npartmax=0 

  ! Number of superoct levels
  integer::nsuperoct=0
  
  ! MPI domain overloading
  integer::overload=1

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer::ncachemax=10000    ! Maximum number of cache lines
  real(dp)::boxlen=1.0D0      ! Box length along x direction

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::nremap=0           ! Load balancing frequency (0: never)

  ! Output parameters
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)
  logical::gadget_output=.false. ! Output in gadget format

  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Physics parameters
  logical ::pressure_fix=.false.

  ! Movie
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tendmov=0.,aendmov=0.
  real(kind=8),dimension(1:10000)::amovout
  real(kind=8),dimension(1:10000)::tmovout
  logical::movie=.false.
  logical::zoom_only=.false.
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
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
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level
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

  ! Default units
  real(dp)::units_density=1.0 ! [g/cm^3]
  real(dp)::units_time=1.0    ! [seconds]
  real(dp)::units_length=1.0  ! [cm]

  ! Initial conditions parameters from grafic
  real(dp)::aexp_ini=10.

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

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:NVAR-NDIM-2)::err_grad_var=-1.0
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
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-NDIM-2-NENER)::var_region=0.0
#endif

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Other hydro solver parameters
  real(dp)::T2_star=10.
  real(dp)::g_star=1.0
  real(dp)::n_star=1d100
  logical::isothermal
  logical::cooling
  
  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0.0

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Fast solver with MPI pre-fetch (memory intensive)
  logical :: fast_solver = .false.

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/cosmo,pic,poisson,hydro,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap &
       & ,static,geom,overload,nsuperoct
  namelist/output_params/noutput,foutput,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout,gadget_output
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax,fast_solver
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only &
       & ,imovout,imov,tendmov,aendmov,proj_axis,movie_vars,movie_vars_txt
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
#if NENER>0
       & ,prad_region &
#endif
#if NVAR>NDIM+2+NENER
       & ,var_region &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag,gamma_rad &
       & ,pressure_fix,scheme,riemann
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
#if NENER>0
       & ,err_grad_prad &
#endif
#if NVAR>NDIM+2+NENER
       & ,err_grad_var &
#endif
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type
  namelist/physics_params/cooling,units_density,units_time,units_length &
       & ,T2_star,g_star,n_star,isothermal
  
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
  write(*,*)'       written by Romain Teyssier (University of Zurich)       '
  write(*,*)'               (c) CEA 1999-2007, UZH 2008-2014                '
  write(*,*)' '

  ! Check nvar is not too small
  write(*,'(" Using solver = hydro with nvar = ",I2," and ndim = ",I1)')nvar,ndim
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call mdl_abort
  endif

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = command_argument_count()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call mdl_abort
  END IF
  CALL getarg(1,infile)

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'File '//TRIM(infile)//' does not exist'
     call mdl_abort
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
  noutput=MIN(noutput,MAXOUT)
  tmovout=1d100
  amovout=1d100
  if(imovout>0) then
     if(tendmov>0)then
        do i=1,imovout
           tmovout(i)=tendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(aendmov>0)then
        do i=1,imovout
           amovout(i)=aendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  
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
        ngridmax=int(ngridtot/int(g%ncpu,kind=8),kind=4)
     endif
  end if
  if(npartmax==0)then
     npartmax=int(nparttot/int(g%ncpu,kind=8),kind=4)
  endif

  !----------------------------
  ! Read hydro parameters 
  !----------------------------
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call mdl_abort
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue
  close(1)
  
  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call mdl_abort
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call mdl_abort
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
  if(nvar<(ndim+2+nener))then
     write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif
  
  if(.not. nml_ok)then
     write(*,*)'Too many errors in the namelist'
     write(*,*)'Aborting...'
     call mdl_abort
  end if

  ! Fill in all run parameters in corresponding structure

  r%cosmo=cosmo
  r%pic=pic
  r%poisson=poisson
  r%hydro=hydro
  r%verbose=verbose
  r%debug=debug
  r%nrestart=nrestart
  r%ncontrol=ncontrol
  r%nstepmax=nstepmax
  r%nsubcycle=nsubcycle
  r%nremap=nremap
  r%static=static
  r%geom=geom
  r%overload=overload
  r%nsuperoct=nsuperoct

  r%noutput=noutput
  r%foutput=foutput
  r%aout=aout
  r%tout=tout
  r%output_mode=output_mode
  r%gadget_output=gadget_output

  r%levelmin=levelmin
  r%nlevelmax=nlevelmax
  r%ngridmax=ngridmax
  r%ncachemax=ncachemax
  r%npartmax=npartmax
  r%nexpand=nexpand
  r%boxlen=boxlen

  r%epsilon=epsilon
  r%gravity_type=gravity_type
  r%gravity_params=gravity_params
  r%cic_levelmax=cic_levelmax
  r%cg_levelmin=cg_levelmin
  r%fast_solver=fast_solver

  r%nw_frame=nw_frame
  r%nh_frame=nh_frame
  r%levelmax_frame=levelmax_frame
  r%ivar_frame=ivar_frame
  r%xcentre_frame=xcentre_frame
  r%ycentre_frame=ycentre_frame
  r%zcentre_frame=zcentre_frame
  r%deltax_frame=deltax_frame
  r%deltay_frame=deltay_frame
  r%deltaz_frame=deltaz_frame
  r%movie=movie
  r%zoom_only=zoom_only
  r%imovout=imovout
  r%imov=imov
  r%tendmov=tendmov
  r%aendmov=aendmov
  r%amovout=amovout
  r%tmovout=tmovout
  r%proj_axis=proj_axis
  r%movie_vars_txt=movie_vars_txt
  if(r%movie)call set_movie_vars(r)

  r%gamma=gamma
  r%courant_factor=courant_factor
  r%smallc=smallc
  r%smallr=smallr
  r%niter_riemann=niter_riemann
  r%slope_type=slope_type
  r%difmag=difmag
  r%gamma_rad=gamma_rad(1:nener)
  r%pressure_fix=pressure_fix
  r%scheme=scheme
  if(riemann=='llf')r%riemann=solver_llf
  if(riemann=='hll')r%riemann=solver_hll
  if(riemann=='hllc')r%riemann=solver_hllc

  r%cooling=cooling
  r%units_density=units_density
  r%units_time=units_time
  r%units_length=units_length
  r%T2_star=T2_star
  r%g_star=g_star
  r%n_star=n_star
  r%isothermal=isothermal

  r%m_refine=m_refine
  r%r_refine=r_refine
  r%x_refine=x_refine
  r%y_refine=y_refine
  r%z_refine=z_refine
  r%exp_refine=exp_refine
  r%a_refine=a_refine
  r%b_refine=b_refine
  r%jeans_refine=jeans_refine
  r%var_cut_refine=var_cut_refine
  r%mass_cut_refine=mass_cut_refine
  r%ivar_refine=ivar_refine

  r%interpol_var=interpol_var
  r%interpol_type=interpol_type
  r%err_grad_d=err_grad_d
  r%err_grad_u=err_grad_u
  r%err_grad_p=err_grad_p
  r%floor_d=floor_d
  r%floor_u=floor_u
  r%floor_p=floor_p
  r%mass_sph=mass_sph
#if NENER>0
  r%err_grad_prad=err_grad_prad
#endif
#if NVAR>NDIM+2+NENER
  r%err_grad_var=err_grad_var
#endif

  if(nrestart>0)filetype='restart'
  r%filetype=filetype
  r%initfile=initfile
  r%multiple=multiple
  r%nregion=nregion
  r%region_type=region_type
  r%x_center=x_center
  r%y_center=y_center
  r%z_center=z_center
  r%length_x=length_x
  r%length_y=length_y
  r%length_z=length_z
  r%exp_region=exp_region
  r%d_region=d_region
  r%d_region=u_region
  r%d_region=v_region
  r%d_region=w_region
  r%d_region=p_region
#if NENER>0
  r%prad_region=prad_region
#endif
#if NVAR>NDIM+2+NENER
  r%var_region=var_region
#endif

  ! Broadcast parameters to all CPUs.
  call m_broadcast_params(r,g,m,p,mdl)
  
end subroutine m_read_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_params(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------
  integer,dimension(:),allocatable::input_array
  integer,dimension(:),allocatable::output_array

  ! Broadcast parameters to all CPUs.
  allocate(input_array(1:storage_size(r)/32))
  input_array=transfer(r,input_array)
  allocate(output_array(1:1))
  output_array=0
  write(*,*)"BROADCASTING parameters"
  call r_broadcast_params(r,g,m,p,mdl,g%ncpu,storage_size(r)/32,1,input_array,output_array)
  write(*,*)"READ_PARAMS completed with output ",output_array(1)
  deallocate(input_array,output_array)

end subroutine m_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_params(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array,output_array)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::next_range,next_cpu
  integer,dimension(1:output_size)::next_output

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_BCAST_PARAMS,next_cpu,next_range,input_size,output_size,input_array)
     call r_broadcast_params(r,g,m,p,mdl,next_range,input_size,output_size,input_array,output_array)
     call mdl_get_reply(mdl,next_cpu,output_size,next_output)
     output_array(1)=output_array(1)+next_output(1)
  else
     r=transfer(input_array,r)
     output_array(1)=1
  endif

end subroutine r_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_global(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------
  integer,dimension(:),allocatable::input_array

  ! Broadcast parameters to all CPUs.
  allocate(input_array(1:storage_size(g)/32))
  input_array=transfer(g,input_array)
  call r_broadcast_global(r,g,m,p,mdl,g%ncpu,storage_size(g)/32,0,input_array)
  deallocate(input_array)

end subroutine m_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_global(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_BCAST_GLOBAL,next_cpu,next_range,input_size,output_size,input_array)
     call r_broadcast_global(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     g=transfer(input_array,g)
  endif

end subroutine r_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

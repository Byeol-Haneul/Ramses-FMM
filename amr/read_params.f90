subroutine read_params(run_p,global_v)
  use amr_commons
  use pm_commons
  use poisson_parameters
  use hydro_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  type(run_t)::run_p
  type(global_t)::global_v

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok
  real(kind=4)::real_mem,real_mem_tot

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

  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
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
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

#ifndef WITHOUTMPI
  call getmem(real_mem)
  call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(myid==1)then
     write(*,*)'Diagnostic right at start-up'
     call writemem(real_mem_tot)
  endif
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
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
  if (myid==1 .and. narg == 2) then
    CALL getarg(2,cmdarg)
    read(cmdarg,*) nrestart
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

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
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  if(npartmax==0)then
     npartmax=nparttot/int(ncpu,kind=8)
  endif
  if(myid>1)verbose=.false.
  call read_hydro_params(nml_ok)

  close(1)

  if (movie)call set_movie_vars

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
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
  end do
     
  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

  ! Fill in all run parameters in corresponding structure

  run_p%cosmo=cosmo
  run_p%pic=pic
  run_p%poisson=poisson
  run_p%hydro=hydro
  run_p%verbose=verbose
  run_p%debug=debug
  run_p%nrestart=nrestart
  run_p%ncontrol=ncontrol
  run_p%nstepmax=nstepmax
  run_p%nsubcycle=nsubcycle
  run_p%nremap=nremap
  run_p%static=static
  run_p%geom=geom
  run_p%overload=overload
  run_p%nsuperoct=nsuperoct

  run_p%noutput=noutput
  run_p%foutput=foutput
  run_p%aout=aout
  run_p%tout=tout
  run_p%output_mode=output_mode
  run_p%gadget_output=gadget_output

  run_p%levelmin=levelmin
  run_p%nlevelmax=nlevelmax
  run_p%ngridmax=ngridmax
  run_p%ncachemax=ncachemax
  run_p%npartmax=npartmax
  run_p%nexpand=nexpand
  run_p%boxlen=boxlen

  run_p%epsilon=epsilon
  run_p%gravity_type=gravity_type
  run_p%gravity_params=gravity_params
  run_p%cic_levelmax=cic_levelmax
  run_p%cg_levelmin=cg_levelmin
  run_p%fast_solver=fast_solver

  run_p%nw_frame=nw_frame
  run_p%nh_frame=nh_frame
  run_p%levelmax_frame=levelmax_frame
  run_p%ivar_frame=ivar_frame
  run_p%xcentre_frame=xcentre_frame
  run_p%ycentre_frame=ycentre_frame
  run_p%zcentre_frame=zcentre_frame
  run_p%deltax_frame=deltax_frame
  run_p%deltay_frame=deltay_frame
  run_p%deltaz_frame=deltaz_frame
  run_p%movie=movie
  run_p%zoom_only=zoom_only
  run_p%imovout=imovout
  run_p%imov=imov
  run_p%tendmov=tendmov
  run_p%aendmov=aendmov
  run_p%amovout=amovout
  run_p%tmovout=tmovout
  run_p%proj_axis=proj_axis
  run_p%movie_vars=movie_vars
  run_p%movie_vars_txt=movie_vars_txt

  run_p%gamma=gamma
  run_p%courant_factor=courant_factor
  run_p%smallc=smallc
  run_p%smallr=smallr
  run_p%niter_riemann=niter_riemann
  run_p%slope_type=slope_type
  run_p%difmag=difmag
  run_p%gamma_rad=gamma_rad
  run_p%pressure_fix=pressure_fix
  run_p%scheme=scheme
  run_p%riemann=riemann

  run_p%cooling=cooling
  run_p%units_density=units_density
  run_p%units_time=units_time
  run_p%units_length=units_length
  run_p%T2_star=T2_star
  run_p%g_star=g_star
  run_p%n_star=n_star
  run_p%isothermal=isothermal

  run_p%m_refine=m_refine
  run_p%r_refine=r_refine
  run_p%x_refine=x_refine
  run_p%y_refine=y_refine
  run_p%z_refine=z_refine
  run_p%exp_refine=exp_refine
  run_p%a_refine=a_refine
  run_p%b_refine=b_refine
  run_p%jeans_refine=jeans_refine
  run_p%var_cut_refine=var_cut_refine
  run_p%mass_cut_refine=mass_cut_refine
  run_p%ivar_refine=ivar_refine

  run_p%interpol_var=interpol_var
  run_p%interpol_type=interpol_type
  run_p%err_grad_d=err_grad_d
  run_p%err_grad_u=err_grad_u
  run_p%err_grad_p=err_grad_p
  run_p%floor_d=floor_d
  run_p%floor_u=floor_u
  run_p%floor_p=floor_p
  run_p%mass_sph=mass_sph
#if NENER>0
  run_p%err_grad_prad=err_grad_prad
#endif
#if NVAR>NDIM+2+NENER
  run_p%err_grad_var=err_grad_var
#endif

  run_p%filetype=filetype
  run_p%initfile=initfile
  run_p%multiple=multiple
  run_p%nregion=nregion
  run_p%region_type=region_type
  run_p%x_center=x_center
  run_p%y_center=y_center
  run_p%z_center=z_center
  run_p%length_x=length_x
  run_p%length_y=length_y
  run_p%length_z=length_z
  run_p%exp_region=exp_region
  run_p%d_region=d_region
  run_p%d_region=u_region
  run_p%d_region=v_region
  run_p%d_region=w_region
  run_p%d_region=p_region
#if NENER>0
  run_p%prad_region=prad_region
#endif
#if NVAR>NDIM+2+NENER
  run_p%var_region=var_region
#endif

  global_v%ncpu=ncpu
  global_v%myid=myid
  

end subroutine read_params


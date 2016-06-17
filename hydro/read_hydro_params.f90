subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  integer::i,idim
  real(dp)::scale

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
#if NENER>0
       & ,prad_region &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
#if NENER>0
       & ,gamma_rad &
#endif
       & ,pressure_fix,scheme,riemann
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type
  namelist/physics_params/cooling,units_density,units_time,units_length &
       & ,T2_star,g_star,n_star,isothermal

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
  if(nvar<(ndim+2+nener))then
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do

end subroutine read_hydro_params


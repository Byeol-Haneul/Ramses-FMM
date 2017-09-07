module hydro_parameters
  use amr_parameters

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=ndim+2+nener
#else
  integer,parameter::nvar=NVAR
#endif

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

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6

end module hydro_parameters

module const
  use amr_parameters

  ! Some useful constant
  real(dp),parameter ::bigreal = 1.0d+30
  real(dp),parameter ::zero = 0.0
  real(dp),parameter ::one = 1.0
  real(dp),parameter ::two = 2.0
  real(dp),parameter ::three = 3.0
  real(dp),parameter ::four = 4.0
  real(dp),parameter ::two3rd = 0.6666666666666667
  real(dp),parameter ::half = 0.5
  real(dp),parameter ::third = 0.33333333333333333
  real(dp),parameter ::forth = 0.25
  real(dp),parameter ::sixth = 0.16666666666666667

end module const

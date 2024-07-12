module hydro_parameters
  use amr_parameters, only: ndim

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=5+nener
#else
  integer,parameter::nvar=NVAR
#endif

#ifdef MHD
  integer,parameter::nprim=NVAR+3
  integer,parameter::ie=8
#else
  integer,parameter::nprim=NVAR
  integer,parameter::ie=5
#endif

! Ionizattion fractions
#ifdef NIONS
  integer,parameter::nIons=NIONS
#else
  integer,parameter::nIons=0
#endif

  integer,parameter::solver_llf=1
  integer,parameter::solver_hll=2
  integer,parameter::solver_hllc=3
  integer,parameter::solver_hlld=4
  integer,parameter::solver_roe=5
  integer,parameter::solver_upwind=6

  integer,parameter::solver2d_llf=1
  integer,parameter::solver2d_hllf=2
  integer,parameter::solver2d_hlla=3
  integer,parameter::solver2d_hlld=4
  integer,parameter::solver2d_roe=5
  integer,parameter::solver2d_upwind=6

end module hydro_parameters

module const
  use amr_parameters, only: dp

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

#ifdef RT

module rt_parameters
  use amr_parameters, only: ndim, clight, dp

#ifdef NRTGROUPS
  integer,parameter::nrtgroups=NRTGROUPS          ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgroups=1
#endif
  integer,parameter::nrtvar=nrtgroups*(1+ndim) ! # of RT variables (photon density and flux)
  real(dp)::rt_c=1.                  !   Reduced lightspeed in code units
  real(dp)::rt_c_cgs=clight          !   Reduced lightspeed [cm s-1]

end module rt_parameters

module rt_const
  use amr_parameters, only: dp

  ! Some useful constant
  !real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; SI

end module rt_const

#endif
subroutine units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_parameters, only: dp,ndim
  use constants, only: kB,mH,X_H,rhoc
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  real(kind=8)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

#define COSMO 1
#define COEUR 2
#define MERGER 3

#if UNITS==COSMO

  ! Units for cosmological simulations are hard-coded
  scale_d = g%omega_m * rhoc *(g%h0/100.)**2 / g%aexp**3
  scale_t = g%aexp**2 / (g%h0*1d5/3.0856776d24)
  scale_l = g%aexp * g%boxlen_ini * 3.0856776d24 / (g%h0/100)
  scale_v = scale_l / scale_t    ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = X_H/mH * scale_d    ! scale_nH converts rho in user units into nH in H/cc

#elif UNITS==COEUR

  ! Units for molecular core collapse
  scale_d = 0.443201646421875d-17  ! scale_d converts mass density from user units into g/cc
  scale_t = 0.183887441975074d+13  ! 1/(G*d0)**0.5
  scale_l = 25*0.478712000000000d+18  ! 25*32000 AU
  scale_v = scale_l / scale_t     ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2   ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = 1.0/mH * scale_d     ! scale_nH converts rho in user units into nH in H/cc

#elif UNITS==MERGER

  ! Units are in kpc, Gyr, km/s and time is 0.98 Gyr
  scale_l = 3.086568025d21
  scale_d = 1.573391528d-26
  scale_t = 1.0/sqrt(6.67d-8*scale_d)
  scale_v = scale_l / scale_t
  scale_T2 = mH/kB * scale_v**2
  scale_nH = X_H/mH * scale_d

#else

  ! Default units from namelist
  scale_d = r%units_density      ! scale_d converts mass density from user units into g/cc
  scale_t = r%units_time         ! scale_t converts time from user units into seconds
  scale_l = r%units_length       ! scale_l converts distance from user units into cm
  scale_v = scale_l / scale_t    ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = X_H/mH * scale_d    ! scale_nH converts rho in user units into nH in H/cc

#endif

end subroutine units

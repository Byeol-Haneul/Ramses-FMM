subroutine units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_parameters, only: dp,ndim,kB,mH,X_H,rhoc
  use amr_commons, only: run_t,global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! Default units from namelist
#if UNITS==default
  scale_d = r%units_density      ! scale_d converts mass density from user units into g/cc
  scale_t = r%units_time         ! scale_t converts time from user units into seconds
  scale_l = r%units_length       ! scale_l converts distance from user units into cm
  scale_v = scale_l / scale_t    ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = X_H/mH * scale_d    ! scale_nH converts rho in user units into nH in H/cc
#endif

  ! Units for cosmological simulations are hard-coded
#if UNITS==cosmo
  scale_d = g%omega_m * rhoc *(g%h0/100.)**2 / g%aexp**3
  scale_t = g%aexp**2 / (g%h0*1d5/3.08d24)
  scale_l = g%aexp * g%boxlen_ini * 3.08d24 / (g%h0/100)
  scale_v = scale_l / scale_t    ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = X_H/mH * scale_d    ! scale_nH converts rho in user units into nH in H/cc
#endif

#if UNITS==coeur
  scale_d = 0.443201646421875d-17  ! scale_d converts mass density from user units into g/cc
  scale_t = 0.183887441975074d+13  !1/(G*d0)**0.5
  scale_l = 25*0.478712000000000d+18  !25*32000 AU
  scale_v = scale_l / scale_t    ! scale_v converts velocity in user units into cm/s
  scale_T2 = mH/kB * scale_v**2  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_nH = 1.0/mH * scale_d    ! scale_nH converts rho in user units into nH in H/cc
#endif

end subroutine units

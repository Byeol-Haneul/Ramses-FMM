subroutine rt_units(r,g,scale_Np,scale_Fp)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  real(kind=8)::scale_Np, scale_Fp
  real(kind=8)::scale_nH, scale_T2, scale_t, scale_v, scale_d, scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units for RT variables.
  ! scale_Np <= photon number density scale [# cm-3]
  ! scale_Fp <= photon flux density scale [# cm-2 s-1]
  !-----------------------------------------------------------------------

#define COSMO 1

  call units(r, g, scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

#if UNITS==COSMO
  ! Units for cosmological simulations are hard-coded
  scale_Np = 1./g%aexp**3
  scale_Fp = scale_Np * scale_v
#else
  ! Default units from namelist
  scale_Np = r%units_Np          ! scale_Np converts photon density from user units into #/cc
  scale_Fp = scale_Np * scale_v  ! scale_Fp converts photon flux from user units into #/cm2/s
#endif

end subroutine rt_units

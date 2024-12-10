#ifdef RT

module rt_parameters
  use amr_parameters, only: ndim, clight, dp

#ifdef NRTGROUPS
  integer,parameter::nrtgroups=NRTGROUPS          ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgroups=1
#endif
  integer,parameter::nrtvar=nrtgroups*(1+ndim) ! # of RT variables (photon density and flux)
  real(dp)::rt_c=1.,rt_c2                      !   Reduced lightspeed in code units
  real(dp)::rt_c_cgs=clight          !   Reduced lightspeed [cm s-1]

end module rt_parameters

module rt_const
  use amr_parameters, only: dp

  ! Some useful constant
  !real(dp),parameter ::c_cgs        = 2.9979246d+10 ! Speed of light [cm s-1]; SI

end module rt_const

#endif
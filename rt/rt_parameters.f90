module rt_parameters
  use amr_parameters, only: ndim
#ifdef NRTGROUPS
  integer,parameter::nrtgroups=NRTGROUPS       ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgroups=1
#endif
  integer,parameter::nrtvar=nrtgroups*(1+ndim) ! # of RT variables (photon density and flux)
end module rt_parameters

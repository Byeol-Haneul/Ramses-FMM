module rt_parameters
  use amr_parameters, only: ndim, dp

#ifdef NRTGRP
  integer,parameter::nrtgrp=NRTGRP             ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgrp=1
#endif
  integer,parameter::nrtvar=nrtgrp*(1+ndim) ! # of RT variables (photon density and flux)

#ifndef NPRE
  real(dp),parameter::smallNp=1d-30            !               Minimum photon density
#else
#if NPRE==4
  real(dp),parameter::smallNp=1d-30            !               Minimum photon density
#else
  real(dp),parameter::smallNp=1d-50            !               Minimum photon density
#endif
#endif

end module rt_parameters

module rt_parameters
  use amr_parameters, only: ndim

#ifdef NRTGRP
  integer,parameter::nrtgrp=NRTGRP             ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgrp=1
#endif
  integer,parameter::nrtvar=nrtgrp*(1+ndim)    ! # of RT variables (photon density and flux)

  real(kind=8),parameter::smallNp=1d-50        !  Minimum photon density

end module rt_parameters

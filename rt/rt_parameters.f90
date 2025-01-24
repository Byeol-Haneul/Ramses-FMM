module rt_parameters
  use amr_parameters, only: ndim, dp

#ifdef NRTGRP
  integer,parameter::nrtgrp=NRTGRP             ! # of photon groups (set in Makefile)
#else
  integer,parameter::nrtgrp=1
#endif
  integer,parameter::nrtvar=nrtgrp*(1+ndim) ! # of RT variables (photon density and flux)

#ifdef NION
  integer,parameter::nion=NION                 !   # of ion species (set in Makefile)
#else
  integer,parameter::nion=3                    !   HII and optionally HI, HeII, HeIII
#endif

#ifndef NPRE
  real(dp),parameter::smallNp=1d-30            !               Minimum photon density
#else
#if NPRE==4
  real(dp),parameter::smallNp=1d-30            !               Minimum photon density
#else
  real(dp),parameter::smallNp=1d-50            !               Minimum photon density
#endif
#endif

  real(dp),parameter::ionEv_HI    = 11.20d0
  real(dp),parameter::ionEv_HII   = 13.60d0
  real(dp),parameter::ionEv_HeII  = 24.59d0
  real(dp),parameter::ionEv_HeIII = 54.42d0
  
end module rt_parameters

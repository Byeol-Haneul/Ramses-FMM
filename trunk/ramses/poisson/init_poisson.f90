subroutine init_poisson
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(verbose)write(*,*)'Entering init_poisson'

end subroutine init_poisson





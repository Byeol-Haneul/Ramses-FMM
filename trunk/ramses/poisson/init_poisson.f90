subroutine init_poisson
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(verbose)write(*,*)'Entering init_poisson'

  ! Allocate multigrid parameters
  allocate(safe_mode(1:nlevelmax))
  safe_mode = .false.

end subroutine init_poisson





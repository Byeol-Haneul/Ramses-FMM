subroutine init_poisson
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(verbose)write(*,*)'Entering init_poisson'

  ! Allocate multigrid parameters
  if(fast_solver)then
     allocate(buffer_mg(1:nlevelmax))
  endif

end subroutine init_poisson





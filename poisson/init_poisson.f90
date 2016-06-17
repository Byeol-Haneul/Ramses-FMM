subroutine init_poisson
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(verbose)write(*,*)'Entering init_poisson'

  ! Allocate Poisson oct-based data
  allocate(grav(1:ngridmax+ncachemax))

  ! Allocate multigrid parameters
  allocate(buffer_mg(1:nlevelmax))
  allocate(safe_mode(1:nlevelmax))
  safe_mode = .false.

end subroutine init_poisson





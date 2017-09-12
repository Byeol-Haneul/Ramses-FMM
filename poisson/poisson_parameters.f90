module poisson_parameters

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters

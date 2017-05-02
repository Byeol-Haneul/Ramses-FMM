module poisson_parameters
  use amr_parameters
  use hash

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

  ! Poisson variables

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0.0

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Fast solver with MPI pre-fetch (memory intensive)
  logical :: fast_solver = .false.

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, dimension(1:MAXLEVEL)::safe_mode=.false.

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

  ! MG hash table
  type(hash_table)::mg_dict

  ! MG grid
  integer,allocatable,dimension(:)::head_mg,tail_mg,noct_mg,noct_tot_mg
  integer::ifree_mg

end module poisson_parameters

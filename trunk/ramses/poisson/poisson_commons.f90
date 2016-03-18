module poisson_commons 
  use amr_commons
  use hash
  use poisson_parameters

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Minimum MG level
  integer :: levelmin_mg

  ! MG hash table
  type(hash_table)::mg_dict

  ! MG grid
  integer,allocatable,dimension(:)::head_mg,tail_mg,noct_mg,noct_tot_mg
  integer::ifree_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons

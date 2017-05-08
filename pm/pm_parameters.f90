module pm_parameters
  use amr_parameters, ONLY: dp

  real(dp) :: mp_min=-1.0 ! Minimum particle mass

  logical :: part_memory=.true. ! Optimize particle memory distribution

  integer :: action_kick_only = 1
  integer :: action_kick_drift = 2  

end module pm_parameters

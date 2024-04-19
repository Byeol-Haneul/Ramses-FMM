module pm_parameters

  logical,parameter :: part_memory=.true. ! Optimize particle memory distribution

  integer,parameter :: action_kick_only = 1
  integer,parameter :: action_kick_drift = 2  

  integer,parameter :: DM_TYPE = 0
  integer,parameter :: STAR_TYPE = 1

end module pm_parameters

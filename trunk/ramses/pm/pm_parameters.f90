module pm_parameters
  use amr_parameters, ONLY: dp
  integer  :: npartmax=0                  ! Maximum number of allocatable particles
  integer  :: npart=0                     ! Actual number of particles in processor
  integer  :: npart_tot=0                 ! Total number of particles in all processors
  integer  :: npart_max=0                 ! Maximum number of particles in all processors

end module pm_parameters

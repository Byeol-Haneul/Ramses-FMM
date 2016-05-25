module pm_parameters
  use amr_parameters, ONLY: dp
  integer  :: npartmax=0                  ! Maximum number of particles
  integer  :: npart=0                     ! Actual number of particles in processor
  integer  :: npart_tot=0                 ! Total number of particles in all processors
  real(dp) :: n_dump_parts_direct = 100000000   ! If a leaf cell contains more than n_dump_parts_direct
                                          ! the histogrammed quantities are used for to compute rho

end module pm_parameters

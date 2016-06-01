module pm_commons
  use amr_parameters
  use pm_parameters

  real(dp),allocatable,dimension(:,:)       ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)       ::vp       ! Velocities
  real(dp),allocatable,dimension(:)         ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)         ::phip     ! Potential
#endif
  integer ,allocatable,dimension(:)         ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)     ::idp      ! Particle unique identifier
  integer ,allocatable,dimension(:)         ::sortp    ! Sorted indices
  integer ,allocatable,dimension(:)         ::workp    ! Work space

  integer ,allocatable,dimension(:)         ::headp    ! Particle levels head
  integer ,allocatable,dimension(:)         ::tailp    ! Particle levels tail

  real(dp)                                  ::mp_min   ! Minimum particle mass

  integer :: action_kick_only = 1
  integer :: action_kick_drift = 2  

end module pm_commons

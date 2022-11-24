module pm_commons
  use amr_parameters, only: dp, i8b

  type part_t

     integer :: npart=0     ! Actual number of particles in processor
     integer(kind=8):: npart_tot=0 ! Total number of particles in all processors
     integer :: npart_max=0 ! Maximum number of particles in all processors
     
     ! Particle dependent arrays
     real(dp),allocatable,dimension(:,:)   ::xp       ! Positions
     real(dp),allocatable,dimension(:,:)   ::vp       ! Velocities
     real(dp),allocatable,dimension(:)     ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
     real(dp),allocatable,dimension(:)     ::phip     ! Potential
#endif
     integer ,allocatable,dimension(:)     ::levelp   ! Current level of particle
     integer(i8b),allocatable,dimension(:) ::idp      ! Particle unique identifier
     integer ,allocatable,dimension(:)     ::sortp    ! Sorted indices
     integer ,allocatable,dimension(:)     ::workp    ! Work space
     
     ! Level dependent arrays
     integer ,allocatable,dimension(:)::headp    ! Particle levels head
     integer ,allocatable,dimension(:)::tailp    ! Particle levels tail
     
  end type part_t

end module pm_commons

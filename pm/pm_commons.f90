module pm_commons

  use amr_parameters, only: dp, i8b

  type part_t

     integer :: type ! Particle type
     integer :: npart=0 ! Actual number of particles in processor
     integer(kind=8):: npart_tot=0 ! Total number of particles in all processors
     integer :: npart_max=0 ! Maximum number of particles in all processors
     integer :: nvaralloc ! Number of allocated variables
     integer :: norphan_peak ! Number of orphan particles outside of peak patch
     integer :: norphan_halo ! Number of orphan particles outside of halo patch

     ! Particle dependent arrays
     real(dp),allocatable,dimension(:,:)   ::xp       ! Position
     real(dp),allocatable,dimension(:,:)   ::vp       ! Velocity
     real(dp),allocatable,dimension(:,:)   ::fp       ! Acceleration
     real(dp),allocatable,dimension(:)     ::mp       ! Mass
     real(dp),allocatable,dimension(:)     ::zp       ! Metallicity
     real(dp),allocatable,dimension(:)     ::tp       ! Formation time
     real(dp),allocatable,dimension(:)     ::tm       ! Merging time
     real(dp),allocatable,dimension(:)     ::up       ! Specific energy
     real(dp),allocatable,dimension(:)     ::phip     ! Potential
     integer ,allocatable,dimension(:)     ::levelp   ! Current level of particle
     integer(i8b),allocatable,dimension(:) ::idp      ! Particle unique identifier
     integer(i8b),allocatable,dimension(:) ::idm      ! Merging particle id
     integer ,allocatable,dimension(:)     ::sortp    ! Sorted index
     integer ,allocatable,dimension(:)     ::workp    ! Work space
     integer ,allocatable,dimension(:)     ::pid      ! Peak ID
     integer ,allocatable,dimension(:)     ::hid      ! Halo ID
     
     ! Level dependent arrays
     integer ,allocatable,dimension(:)::headp    ! First particle in level
     integer ,allocatable,dimension(:)::tailp    ! Last particle in level
     
  end type part_t

end module pm_commons

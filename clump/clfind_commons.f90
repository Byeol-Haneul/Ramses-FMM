module clfind_commons
    use amr_parameters, only: dp

    integer::ntest,itest
    real(dp),allocatable,dimension(:)::denp ! Density of the cells
    integer,allocatable,dimension(:)::testp_sort ! Sort indices

    type peak_t

        integer :: npart=0     ! Actual number of particles in processor
        integer(kind=8):: npart_tot=0 ! Total number of particles in all processors
        !integer :: npart_max=0 ! Maximum number of particles in all processors
        
        ! Particle dependent arrays
        real(dp),allocatable,dimension(:,:)   ::xp       ! Positions
        real(dp),allocatable,dimension(:)     ::denp       ! Density

        integer ,allocatable,dimension(:)     ::levelp   ! Current level of particle
        integer(kind=8),allocatable,dimension(:) ::idp      ! Particle unique identifier
        integer ,allocatable,dimension(:)     ::sortp    ! Sorted indices
        logical ,allocatable,dimension(:)    ::peak        ! whether is a peak
        real(dp),allocatable,dimension(:,:)   ::maxxp       ! Positions of neighbor-peaks
        
        !integer ,allocatable,dimension(:)     ::workp    ! Work space
        
        ! Level dependent arrays
        !integer ,allocatable,dimension(:)::headp    ! Particle levels head
        !integer ,allocatable,dimension(:)::tailp    ! Particle levels tail
        
    end type peak_t
end module clfind_commons
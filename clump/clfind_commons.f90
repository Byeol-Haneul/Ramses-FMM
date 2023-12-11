module clfind_commons
    use amr_parameters, only: dp
    use sparse_matrix
    !integer::ntest,itest
    real(dp),allocatable,dimension(:)::denp ! Density of the cells
    integer,allocatable,dimension(:)::testp_sort ! Sort indices
    integer,allocatable,dimension(:)::npeak_cum !number of peak in all processors
    !integer(kind=8),dimension(0:g%ncpu)::npeak_cum !number of peak in all processors
    logical::clinfo=.false.
    integer::npeaks_max
    ! Spare matrix for saddle points densities
    type(sparse_mat)::sparse_saddle_dens

    ! Hash table variables
    integer::nhash,hfree,hcollision
    integer,dimension(:),allocatable::gkey,nkey,hkey


    real(dp)::relevance_threshold=2
    real(dp)::density_threshold=-1
    real(dp)::saddle_threshold=-1
    real(dp)::mass_threshold=0

    type peak_t

        integer :: npart=0     ! Actual number of particles in processor
        integer(kind=8):: npart_tot=0 ! Total number of particles in all processors
        integer(kind=8):: npeak_tot=0 ! Total number of peaks in all processors
        !integer :: npart_max=0 ! Maximum number of particles in all processors
        
        ! Particle dependent arrays
        real(dp),allocatable,dimension(:,:)   ::xp       ! Positions
        real(dp),allocatable,dimension(:)     ::denp       ! Density
        real(dp),allocatable,dimension(:)     ::denpm       ! Density of the peak

        integer ,allocatable,dimension(:)     ::levelp   ! Current level of particle
        integer ,allocatable,dimension(:)     ::levelpm   ! Current level of peak
        integer(kind=8),allocatable,dimension(:) ::idp      ! Particle unique identifier
        integer,allocatable,dimension(:) ::pid     ! the peak id in all processors
        integer ,allocatable,dimension(:)     ::sortp    ! Sorted indices
        logical ,allocatable,dimension(:)    ::peak        ! whether is a peak
        real(dp),allocatable,dimension(:,:)   ::maxxp       ! Positions of neighbor-peaks
        
        !integer ,allocatable,dimension(:)     ::workp    ! Work space
        
        ! Level dependent arrays
        integer ,allocatable,dimension(:)::headp    ! Particle levels head
        integer ,allocatable,dimension(:)::tailp    ! Particle levels tail
        
    end type peak_t


    type clump_t

        integer :: npart=0     ! Actual number of particles in processor
        integer(kind=8):: npeak_tot=0 ! Total number of peaks in all processors
        !integer :: npart_max=0 ! Maximum number of particles in all processors
        
        ! Particle dependent arrays
        real(dp),allocatable,dimension(:,:)   ::xp       ! Positions
        real(dp),allocatable,dimension(:)     ::denp       ! Density

        integer ,allocatable,dimension(:)     ::levelp   ! Current level of particle
        integer(kind=8),allocatable,dimension(:) ::idp      ! the peak id in all processors
        integer(kind=8),allocatable,dimension(:) ::idc      ! the clump id in all processors
        integer,allocatable,dimension(:) ::pid     ! the reference id of the particle
        integer ,allocatable,dimension(:)     ::sortp    ! Sorted indices
        integer,allocatable,dimension(:) ::alive     ! whether this peak is alive
        integer,allocatable,dimension(:) ::lev_peak   ! peak levels
        integer,allocatable,dimension(:) ::ind_halo   ! ind of halo
        integer,allocatable,dimension(:) ::n_cell_halo   ! ind of halo
        integer,allocatable,dimension(:) ::n_cells   ! ind of halo
        real(dp),allocatable,dimension(:)     ::halo_mass       ! halo mass
        real(dp),allocatable,dimension(:)     ::clump_mass       ! halo mass

        real(dp),allocatable,dimension(:) ::relevance     ! temp relevance list
        
        !integer ,allocatable,dimension(:)     ::workp    ! Work space
        
        ! Level dependent arrays
        integer ,allocatable,dimension(:)::headp    ! Particle levels head
        integer ,allocatable,dimension(:)::tailp    ! Particle levels tail
        
    end type clump_t
end module clfind_commons
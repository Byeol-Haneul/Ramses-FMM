module clfind_commons
    use amr_parameters, only: dp
    use hash

    type clump_t

       real(dp)::relevance_threshold=2
       real(dp)::density_threshold=-1
       real(dp)::saddle_threshold=-1
       real(dp)::mass_threshold=0
       real(dp)::purity_threshold=-1
       real(dp)::bound_threshold=2
       real(dp)::fraction_threshold=0.1d0

       integer :: ntest=0 ! Actual number of test particles in current processor
       integer(kind=8) :: ntest_tot=0 ! Total number of test particles across all processors
       integer,allocatable,dimension(:) :: cell ! Cell index of test particle
       integer,allocatable,dimension(:) :: grid ! Grid index of test particle
       integer,allocatable,dimension(:) :: level ! Level of test particle
       integer,allocatable,dimension(:,:) :: hash ! Hash key of densest neighbor

       integer :: npeak=0 ! Actual number of density peaks in current processor
       integer :: npeak_max ! Maximum number of peaks per processor including ghost peaks
       integer(kind=8) :: npeak_tot=0 ! Total number of density peaks across all processors
       integer(kind=8),allocatable,dimension(:) :: npeak_cum ! Cumulative number of peak per processor
       integer :: merge_levelmax ! Max level of sub-halo merging tree

       integer,allocatable,dimension(:) :: peak_cell ! Cell index of peak
       integer,allocatable,dimension(:) :: peak_grid ! Grid index of peak
       integer,allocatable,dimension(:) :: peak_level ! Level of peak
       integer,allocatable,dimension(:) :: saddle_nbor ! Neighboring densest saddle peak global id
       integer,allocatable,dimension(:) :: n_cells ! Number of AMR cells per peak patch
       integer,allocatable,dimension(:) :: lev_peak ! AMR level of the peak
       integer,allocatable,dimension(:) :: new_peak ! Peak ID in which current peak has merged
       real(dp),allocatable,dimension(:) :: max_dens ! Density value at the peak
       real(dp),allocatable,dimension(:) :: clump_mass ! Mass inside peak patch
       real(dp),allocatable,dimension(:) :: relevance ! Relevance (prominence) of the peak
       real(dp),allocatable,dimension(:) :: tidal_dens ! Density at the tidal radius
       real(dp),allocatable,dimension(:) :: min_dens ! min density of the clump
       real(dp),allocatable,dimension(:) :: clump_vol ! volume of the clump
       real(dp),allocatable,dimension(:) :: particle_mass ! clump mass using directly dark matter particles
       real(dp),allocatable,dimension(:,:) :: peak_pos ! position of the peak
       real(dp),allocatable,dimension(:,:) :: peak_vel ! velocity of the peak
       real(dp),allocatable,dimension(:,:) :: peak_acc ! acceleration of the peak
       real(dp),allocatable,dimension(:,:) :: peak_com ! center of mass of the peak

       integer,allocatable,dimension(:) :: ind_halo ! Peak index of halo patch
       integer,allocatable,dimension(:) :: npart ! number of particles per clump
       integer,allocatable,dimension(:) :: pid ! peak id of the parent (bound to) clump
       real(dp),allocatable,dimension(:) :: halo_mass ! Mass inside halo patch
       real(dp),allocatable,dimension(:) :: saddle_dens ! Density of the densest saddle point
       real(dp),allocatable,dimension(:,:) :: mass_bin ! cumulative mass profile of halo
       real(dp),allocatable,dimension(:,:) :: phi ! gravitational self-potential

       integer,allocatable,dimension(:) :: nsink ! number of sinks per clump
       integer,allocatable,dimension(:) :: form_sink ! does peak form a new sink particle?
       integer,allocatable,dimension(:) :: min_sink_id ! minimum sink id in clump

       ! Software cache array for peaks
       integer :: ncachemax=1000
       integer :: free_cache,ncache,ifree
       logical,allocatable,dimension(:) :: dirty
       logical,allocatable,dimension(:) :: occupied
       logical,allocatable,dimension(:) :: locked
       integer,allocatable,dimension(:) :: parent_cpu
       integer(kind=8),allocatable,dimension(:) :: gid

       type(hash_simple)::peak_dict ! Peak hash table

    end type clump_t

end module clfind_commons

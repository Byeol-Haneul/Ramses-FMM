module clfind_commons
    use amr_parameters, only: dp
    use sparse_matrix

    type clump_t

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

       integer,allocatable,dimension(:) :: peak_cell ! Cell index of peak
       integer,allocatable,dimension(:) :: peak_grid ! Grid index of peak
       integer,allocatable,dimension(:) :: peak_level ! Level of peak
       integer,allocatable,dimension(:) :: n_cells ! Number of AMR cells per peak patch
       integer,allocatable,dimension(:) :: n_cells_halo ! Number of AMR cells per halo
       integer,allocatable,dimension(:) :: lev_peak ! AMR level of the peak
       integer,allocatable,dimension(:) :: new_peak ! Peak ID in which current peak has merged
       integer,allocatable,dimension(:) :: ind_halo ! Peak ID of the halo densest peak
       real(dp),allocatable,dimension(:) :: max_dens ! Density value at the peak
       real(dp),allocatable,dimension(:) :: halo_mass ! Total halo mass
       real(dp),allocatable,dimension(:) :: clump_mass ! Mass inside peak patch
       real(dp),allocatable,dimension(:) :: relevance ! Relevance (prominence) of the peak
       real(dp),allocatable,dimension(:,:) :: clump_size ! Size of clump
       real(dp),allocatable,dimension(:,:) :: peak_pos ! position of the peak
       real(dp),allocatable,dimension(:,:) :: center_of_mass ! position of the center of mass
       real(dp),allocatable,dimension(:) :: min_dens ! min density of the clump
       real(dp),allocatable,dimension(:) :: av_dens ! average density of the clump
       real(dp),allocatable,dimension(:) :: clump_vol ! volume of the clump

       integer::peak_recv_tot,peak_send_tot ! Peak communicator arrays
       integer,allocatable,dimension(:)::peak_send_cnt,peak_send_oft
       integer,allocatable,dimension(:)::peak_recv_cnt,peak_recv_oft
       integer,allocatable,dimension(:)::peak_send_buf,peak_recv_buf

       type(sparse_mat) :: sparse_saddle_dens ! Spare matrix for saddle points densities
       
       integer :: nhash,hfree,hcollision  ! Hash table variables
       integer,dimension(:),allocatable :: gkey,nkey,hkey
        
       ! Prime numbers for hash table
       integer,dimension(0:30)::prime=(/2,3,7,13,23,53,97,193,389,769,1543,&
            & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
            & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
            & 402653189,805306457,1610612741/)
       
    end type clump_t
    
end module clfind_commons

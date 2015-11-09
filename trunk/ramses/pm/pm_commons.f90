module pm_commons
  use amr_parameters
  use pm_parameters

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)       ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)       ::vp       ! Velocities
  real(dp),allocatable,dimension(:)         ::mp       ! Masses
  real(dp),allocatable,dimension(:,:)       ::ap       ! Accelerations (more convenient communication if allocated fully
                                                                   ! for all particles
  integer(kind=8),allocatable,dimension(:)  ::part_ref_mask        ! mask for refinement
  integer(kind=8),allocatable,dimension(:,:)::part_hkey
  integer(kind=4),allocatable,dimension(:)  ::current_state
  integer(kind=4),allocatable,dimension(:)  ::sorted_particle_index
  integer(kind=4),allocatable,dimension(:)  ::sort_index
  integer(kind=4),allocatable,dimension(:)  ::part_ind_permutation, part_ind_permutation2
  ! Particle histogram related variables and arrays
  integer(kind=8),allocatable,dimension(:,:)::bin_keys,particle_histogram_keys
  real(dp),allocatable,dimension(:)         ::bin_mass,particle_histogram_mass, bin_count
  integer                                   ::nbins

#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  integer ,allocatable,dimension(:)  :: part_level_offset    
  integer ,allocatable,dimension(:)  :: bin_start_offset


  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid

  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0

  real(dp),allocatable,dimension(:)::dp_part_send_buf,dp_part_recv_buf

end module pm_commons

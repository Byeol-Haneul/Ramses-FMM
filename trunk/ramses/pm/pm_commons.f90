module pm_commons
  use amr_parameters
  use pm_parameters

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)       ::xp, xp_andreas       ! Positions
  real(dp),allocatable,dimension(:,:)       ::vp, vp_andreas       ! Velocities
  real(dp),allocatable,dimension(:)         ::mp, mp_andreas       ! Masses
  integer(kind=8),allocatable,dimension(:)  ::part_ref_mask        ! mask for refinement
  integer(kind=8),allocatable,dimension(:,:)::part_hkey
  integer(kind=4),allocatable,dimension(:)  ::current_state
  integer(kind=4),allocatable,dimension(:)  ::sorted_particle_index
  integer(kind=4),allocatable,dimension(:)  ::sort_index
  integer(kind=4),allocatable,dimension(:)  ::part_ind_permutation, part_ind_permutation2
  ! Particle histogram related variables and arrays
  integer(kind=8),allocatable,dimension(:,:)::bin_keys,particle_histogram_keys
  real(dp),allocatable,dimension(:)         ::bin_mass,particle_histogram_mass
  integer,allocatable,dimension(:)         ::bin_count
  integer                                   ::nbins

#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp,levelp_andreas   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp,idp_andreas    ! Identity of particle

  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid

  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0

  ! Particle communicator arrays
  integer,allocatable,dimension(:)::part_send_cnt,part_send_oft
  integer,allocatable,dimension(:)::part_recv_cnt,part_recv_oft
  integer::part_recv_tot,part_send_tot
  integer(kind=8),allocatable,dimension(:,:)::receive_keys, send_keys
  real(dp),allocatable,dimension(:)::dp_part_send_buf,dp_part_recv_buf
  ! Particle histogram communicator arrays
  integer,allocatable,dimension(:)::bin_send_cnt,bin_send_oft
  integer,allocatable,dimension(:)::bin_recv_cnt,bin_recv_oft
  integer(kind=8),allocatable,dimension(:,:)::recv_bin_keys, send_bin_keys
  real(dp),allocatable,dimension(:)::recv_bin_mass, send_bin_mass
  integer::bin_recv_tot,bin_send_tot
  integer :: mybins, mybins_offset



  
end module pm_commons

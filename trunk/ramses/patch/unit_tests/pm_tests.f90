subroutine pm_tests(all_ok)
  implicit none
  logical::all_ok
  logical::pm_ok=.true.

  call sort_particle_tests(pm_ok)
!  call other_test(pm_ok)
!  call other_test2(pm_ok)

  if(.not. pm_ok)then
     write(*,*)'PM_TESTS FAILED'
     all_ok=.false.
  end if
end subroutine pm_tests


! =====================================================================================
! =====================================================================================
! ADD UNIT TESTS HERE
! =====================================================================================
! =====================================================================================
subroutine sort_particle_tests(all_ok)
  use amr_commons
  use amr_parameters
  use pm_commons
  use hilbert,       only: hilbert3d_for_particle
  use sort,          only: lsd_radix_sort_particles, gt_3keys
  implicit none

  ! A simple test which transforms integer coordinates into a 3 integer hilbert key
  ! and the hilbert key back into integer coordinate. Results must be equal to input.

  integer, parameter :: size=1000
  integer, parameter :: offs=317
  integer::ilevel,i
  logical::all_ok
  real(dp),dimension(1:size, 1:3)::xfloat


  ! a little bit of init_part
  allocate(xp_andreas    (size,3))
  allocate(vp_andreas    (size,3))
  vp_andreas=0.
  allocate(mp_andreas    (size))
  mp_andreas=1
  allocate(levelp_andreas(size))
  levelp_andreas=1
  allocate(idp_andreas   (size))
  allocate(part_hkey(size,0:2))
  allocate(current_state(size))
  allocate(part_ind_permutation(size))
  allocate(part_ind_permutation2(size))
  allocate(bin_mass(1:2))
  allocate(bin_count(1:2))
  allocate(bin_keys(1:2,0:2))
  boxlen=1

  do ilevel=1,63

     call random_number(xfloat)
     xp_andreas(1:size,1:3)=xfloat(1:size,1:3)

     call hilbert3d_for_particle(0, size,0,ilevel)
     call lsd_radix_sort_particles(0, size, ilevel, ilevel)

     do i=1,size-1
        if (gt_3keys(part_hkey(i,0:2),part_hkey(i+1,0:2)))then
           write(*,*)'particle sort test FAILED for level', ilevel, ilevel
           all_ok=.false.
        end if
     end do
  end do
  
  do ilevel=1,63

     call random_number(xfloat)
     xp_andreas(1:size,1:3)=xfloat(1:size,1:3)

     call hilbert3d_for_particle(offs, size-offs,0,ilevel)
     call lsd_radix_sort_particles(offs, size-offs, ilevel, ilevel)

     do i=offs+1,size-1
        if (gt_3keys(part_hkey(i,0:2),part_hkey(i+1,0:2)))then
           write(*,*)'particle sort test FAILED for level', ilevel, ilevel
           write(*,*)part_hkey(i,0:2), i
           write(*,*)part_hkey(i+1,0:2), i
           all_ok=.false.
        end if
     end do
  end do
  
end subroutine sort_particle_tests

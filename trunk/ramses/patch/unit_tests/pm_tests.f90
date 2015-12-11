subroutine pm_tests(all_ok)
  use pm_commons
  use pm_parameters
  use amr_parameters
  implicit none
  logical::all_ok
  logical::pm_ok=.true.
  integer, parameter :: size=16



  ! a little bit of init_part
  allocate(xp    (size,ndim))
  allocate(vp    (size,ndim))
  vp=0.
  allocate(mp    (size))
  mp=1
  allocate(levelp(size))
  levelp=1
  allocate(idp   (size))
  allocate(part_hkey(size,0:ndim-1))
  allocate(current_state(size))
  allocate(part_ind_permutation(size))
  allocate(part_ind_permutation2(size))
  allocate(bin_mass(1:2))
  allocate(bin_count(1:2))
  allocate(bin_keys(1:2,0:ndim-1))
  npart=size
  boxlen=3





  call sort_particle_tests(pm_ok)
!  call other_test(pm_ok)

  if(.not. pm_ok)then
     write(*,*)'PM_TESTS FAILED'
     all_ok=.false.
  end if

contains

! =====================================================================================
! =====================================================================================
! ADD UNIT TESTS HERE
! =====================================================================================
! =====================================================================================
  subroutine sort_particle_tests(all_ok)
    use amr_commons
    use amr_parameters
    use pm_commons
    use hilbert,       only: hilbert_for_particle
    use sort,          only: lsd_radix_sort_particles, gt_3keys, gt_2keys, apply_particle_permutation
    implicit none
    
    ! A simple test which transforms integer coordinates into a 3 integer hilbert key
    ! and the hilbert key back into integer coordinate. Results must be equal to input.

    integer, parameter :: offs=317
    integer::ilevel,i
    logical::all_ok
    real(dp),dimension(1:size, 1:ndim)::xfloat


#if NDIM > 1
#if NDIM == 2
    do ilevel=1,62
#endif
#if NDIM == 3
    do ilevel=1,63
#endif 

       call random_number(xfloat)
       xp(1:size,1:ndim)=xfloat(1:size,1:ndim)
    
       call hilbert_for_particle(0, size,0,ilevel)
       call lsd_radix_sort_particles(0, size, ilevel, ilevel, .true.)
       call apply_particle_permutation(0,size,ilevel)
       
       do i=1,size-1
#if NDIM == 2
          if (gt_2keys(part_hkey(i,0:1),part_hkey(i+1,0:1)))then
#endif
#if NDIM == 3
          if (gt_3keys(part_hkey(i,0:2),part_hkey(i+1,0:2)))then
#endif
             write(*,*)'particle sort test FAILED for level', ilevel, ilevel
             all_ok=.false.
          end if
       end do
    end do

#if NDIM == 2
    do ilevel=1,62
#endif
#if NDIM == 3
    do ilevel=1,63
#endif 

       call random_number(xfloat)
       xp(1:size,1:ndim)=xfloat(1:size,1:ndim)
       
       call hilbert_for_particle(offs, size-offs,0,ilevel)
       call lsd_radix_sort_particles(offs, size-offs, ilevel, ilevel, .true.)
       call apply_particle_permutation(offs,size-offs,ilevel)
       
       do i=offs+1,size-1
#if NDIM == 2
          if (gt_2keys(part_hkey(i,0:1),part_hkey(i+1,0:1)))then
#endif
#if NDIM == 3
          if (gt_3keys(part_hkey(i,0:2),part_hkey(i+1,0:2)))then
#endif
             write(*,*)'particle sort test FAILED for level', ilevel, ilevel
             all_ok=.false.
          end if
       end do
    end do
    
#endif
    
  end subroutine sort_particle_tests
end subroutine pm_tests
            

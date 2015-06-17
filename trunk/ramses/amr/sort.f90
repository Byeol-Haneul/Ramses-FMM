module sort
  
  ! module containing various routines related to sorting

  ! quick sort routines
  ! subroutine quick_sort(list, order, n)  
  ! subroutine quick_sort_dp(list, order, n)

  ! not used, might be useful at some point...
  ! subroutine swap_2parts_in_mem(i,j)

  ! composite hilbert key comparison functions
  ! function ge_2keys(key_a, key_b)
  ! function ge_3keys(key_a, key_b)
  ! function gt_2keys(key_a, key_b)
  ! function gt_3keys(key_a, key_b)

  ! radix sorts:
  ! recursive msd radix sort
  ! subroutine lsd_radix_sort_particles(offset, np, initial_level, final_level, key_level)

  ! sort paricles by lsd radix sort
  ! subroutine lsd_radix_sort_particles(offset, np, final_level, key_level)
  ! uses
  ! subroutine lsd_counting_sort_3digits(np, ilevel, key_level, sigma1, sigma2, hkey2, hkey1, hkey0)      

  ! reshuffle particles in memory according
  ! subroutine apply_particle_permutation(offset, np, key_level)

contains

  SUBROUTINE quick_sort(list, order, n)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    use amr_parameters, ONLY: qdp
    IMPLICIT NONE
    INTEGER :: n
    REAL(qdp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

    ! Local variable
    INTEGER :: i

    DO i = 1, n
       order(i) = i
    END DO

    CALL quick_sort_1(1, n)

  CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

      use amr_parameters, ONLY: qdp
      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER             :: i, j, itemp
      REAL(qdp)              :: reference, temp
      INTEGER, PARAMETER  :: max_simple_sort_size = 6

      IF (right_end < left_end + max_simple_sort_size) THEN
         ! Use interchange sort for small lists
         CALL interchange_sort(left_end, right_end)

      ELSE
         ! Use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1; j = right_end + 1

         DO
            ! Scan list from left end until element >= reference is found
            DO
               i = i + 1
               IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
               j = j - 1
               IF (list(j) <= reference) EXIT
            END DO


            IF (i < j) THEN
               ! Swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            ELSE IF (i == j) THEN
               i = i + 1
               EXIT
            ELSE
               EXIT
            END IF
         END DO

         IF (left_end < j) CALL quick_sort_1(left_end, j)
         IF (i < right_end) CALL quick_sort_1(i, right_end)
      END IF

    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)

      use amr_parameters, ONLY: qdp
      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER             :: i, j, itemp
      REAL(qdp)           :: temp

      DO i = left_end, right_end - 1
         DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            END IF
         END DO
      END DO

    END SUBROUTINE interchange_sort

  END SUBROUTINE quick_sort
  !########################################################################
  !########################################################################
  !########################################################################
  !########################################################################
  !########################################################################
  SUBROUTINE quick_sort_dp(list, order, n)
    use amr_parameters, ONLY: dp
    IMPLICIT NONE
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.


    INTEGER :: n
    REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

    ! Local variable
    INTEGER :: i

    DO i = 1, n
       order(i) = i
    END DO

    CALL quick_sort_1_dp(1, n)

  CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)

      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER             :: i, j, itemp
      REAL(kind=8)              :: reference, temp
      INTEGER, PARAMETER  :: max_simple_sort_size = 6

      IF (right_end < left_end + max_simple_sort_size) THEN
         ! Use interchange sort for small lists
         CALL interchange_sort_dp(left_end, right_end)

      ELSE
         ! Use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1; j = right_end + 1

         DO
            ! Scan list from left end until element >= reference is found
            DO
               i = i + 1
               IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
               j = j - 1
               IF (list(j) <= reference) EXIT
            END DO


            IF (i < j) THEN
               ! Swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            ELSE IF (i == j) THEN
               i = i + 1
               EXIT
            ELSE
               EXIT
            END IF
         END DO
         IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
         IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
      END IF

    END SUBROUTINE quick_sort_1_dp


    SUBROUTINE interchange_sort_dp(left_end, right_end)

      INTEGER, INTENT(IN) :: left_end, right_end

      !     Local variables                                                                                                                                                                           
      INTEGER             :: i, j, itemp
      REAL(kind=8)           :: temp

      DO i = left_end, right_end - 1
         DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            END IF
         END DO
      END DO

    END SUBROUTINE interchange_sort_dp

  END SUBROUTINE quick_sort_dp

!   subroutine swap_2parts_in_mem(i,j)
!     use amr_parameters, only:i8b
!     use pm_commons

!     integer, intent(in) :: i, j
    
!     ! temporary storage for one element
!     integer(kind=8), save :: i8_temp
!     integer(kind=4), save :: i4_temp
!     integer(i8b),    save :: i8b_temp
!     real(dp),        save :: dp_temp
    
!     ! hilbert keys
!     i8_temp = part_hkey(i,0); part_hkey(i,0) = part_hkey(j,0); part_hkey(j,0) = i8_temp
!     i8_temp = part_hkey(i,1); part_hkey(i,1) = part_hkey(j,1); part_hkey(j,1) = i8_temp
!     i8_temp = part_hkey(i,2); part_hkey(i,2) = part_hkey(j,2); part_hkey(j,2) = i8_temp
!     ! current_state of hilbert state diagram
!     i4_temp = current_state(i); current_state(i) = current_state(j); current_state(j) = i4_temp
!     ! position
!     dp_temp = xp_andreas(i,1); xp_andreas(i,1) = xp_andreas(j,1); xp_andreas(j,1) = dp_temp
!     dp_temp = xp_andreas(i,2); xp_andreas(i,2) = xp_andreas(j,2); xp_andreas(j,2) = dp_temp
!     dp_temp = xp_andreas(i,3); xp_andreas(i,3) = xp_andreas(j,3); xp_andreas(j,3) = dp_temp
!     ! velocity
!     dp_temp = vp_andreas(i,1); vp_andreas(i,1) = vp_andreas(j,1); vp_andreas(j,1) = dp_temp
!     dp_temp = vp_andreas(i,2); vp_andreas(i,2) = vp_andreas(j,2); vp_andreas(j,2) = dp_temp
!     dp_temp = vp_andreas(i,3); vp_andreas(i,3) = vp_andreas(j,3); vp_andreas(j,3) = dp_temp
!     ! mass
!     dp_temp = mp_andreas(i); mp_andreas(i) = mp_andreas(j); mp_andreas(j) = dp_temp
!     ! level
!     i4_temp = levelp_andreas(i); levelp_andreas(i) = levelp_andreas(j); levelp_andreas(j) = i4_temp
!     ! particle index
!     i8b_temp = idp_andreas(i); idp_andreas(i) = idp_andreas(j); idp_andreas(j) = i8b_temp
    
! #ifdef OUTPUT_PARTICLE_POTENTIAL
!     print*, 'add particle potential here'
!     stop
! #endif
    
!   end SUBROUTINE swap_2parts_in_mem


  !########################################################################
  !########################################################################
  !########################################################################
  function ge_2keys(key_a, key_b)
    implicit none
    integer(kind=8), intent(in), dimension (0:1):: key_a, key_b
    logical::ge_2keys
    ! Function to test wether a >= b for a and b two-integer hilbert keys 
    if     (key_a(1) > key_b(1)) then
       ge_2keys  =  .true.
    elseif (key_a(1) < key_b(1)) then
       ge_2keys = .false.
    elseif (key_a(0) > key_b(0)) then
       ge_2keys = .true.
    elseif (key_a(0) < key_b(0)) then
       ge_2keys = .false.
    else
       ge_2keys = .true.
    end if
  end function ge_2keys

  function ge_3keys(key_a, key_b)
    implicit none
    integer(kind=8), intent(in), dimension (0:2):: key_a, key_b
    logical::ge_3keys
    ! Function to test wether a >= b for a and b three-integer hilbert keys 
    if     (key_a(2) > key_b(2)) then
       ge_3keys = .true.
    elseif (key_a(2) < key_b(2)) then
       ge_3keys = .false.
    elseif (key_a(1) > key_b(1)) then
       ge_3keys = .true.
    elseif (key_a(1) < key_b(1)) then
       ge_3keys = .false.
    elseif (key_a(0) > key_b(0)) then
       ge_3keys = .true.
    elseif (key_a(0) < key_b(0)) then
       ge_3keys = .false.
    else
       ge_3keys = .true.
    end if
  end function ge_3keys

  function gt_2keys(key_a, key_b)
    implicit none
    integer(kind=8), intent(in), dimension (0:1):: key_a, key_b
    logical::gt_2keys
    ! Function to test wether a >= b for a and b two-integer hilbert keys 
    if     (key_a(1) > key_b(1)) then
       gt_2keys  =  .true.
    elseif (key_a(1) < key_b(1)) then
       gt_2keys = .false.
    elseif (key_a(0) > key_b(0)) then
       gt_2keys = .true.
    elseif (key_a(0) < key_b(0)) then
       gt_2keys = .false.
    else
       gt_2keys = .false.
    end if
  end function gt_2keys

  function gt_3keys(key_a, key_b)
    implicit none
    integer(kind=8), intent(in), dimension (0:2):: key_a, key_b
    logical::gt_3keys
    ! Function to test wether a > b for a and b three-integer hilbert keys 
    if     (key_a(2) > key_b(2)) then
       gt_3keys = .true.
    elseif (key_a(2) < key_b(2)) then
       gt_3keys = .false.
    elseif (key_a(1) > key_b(1)) then
       gt_3keys = .true.
    elseif (key_a(1) < key_b(1)) then
       gt_3keys = .false.
    elseif (key_a(0) > key_b(0)) then
       gt_3keys = .true.
    elseif (key_a(0) < key_b(0)) then
       gt_3keys = .false.
    else
       gt_3keys = .false.
    end if
  end function gt_3keys
  
  !########################################################################
  !########################################################################
  !########################################################################
  
  subroutine msd_radix_sort_particles(offset, np, initial_level, final_level, key_level)
    implicit none
    integer, intent(in) :: offset, np, initial_level, final_level, key_level

    ! -> Probalby not necessary to implement msd_radix sort at all...
    ! for sorting down to levelmin, lsd_radix is more efficient. After that, only one 
    ! additional level will be sorted at the time. In this case, simply call lsd_radix sort
    ! on shorter sequences (of equal coarse level hilbert keys) for level 1 ...

    ! use a level-dependent bucket_offset because allocation at each recursion is slow
    integer, dimension(0:7, 1:final_level) :: bucket_offset, bucket_count
    integer, dimension(     1:final_level) :: ibucket_level
    

    if (final_level > key_level) then
       write(*,*)'you are trying to sort the hilbert keys to a too high level'
       call clean_stop
    end if

    call msd_counting_sort_3digits(offset, np, initial_level + 1, key_level)
  contains    

    recursive subroutine  msd_counting_sort_3digits(offset, np, ilevel, key_level)
      use pm_commons,  only: part_hkey, part_ind_permutation
      integer, intent(in) :: offset, np, ilevel, key_level

      ! Recursive subroutine that sorts three bits and calls itself 8 times to
      ! sort the shorter particle sequences by the next three (less significant) bits.

      integer,         save :: ind, ibit1, ikey
      integer,         save :: ipart, inext, ip
      integer,         save :: dummy_arg1, dummy_arg2
      integer(kind=8), save :: ibucket

      ! get bit and key to read from 
      ibit1 = (key_level-ilevel)*3       
      ikey = ibit1/63
      ibit1 = mod(ibit1,63)

      ! Count particles per bucket
      bucket_count(0:7, ilevel)=0
      do ipart = offset+1, offset+np
         ibucket=ibits(part_hkey(ipart, ikey),ibit1,3)
         bucket_count(ibucket, ilevel) = bucket_count(ibucket, ilevel) + 1
      end do

      ! "Prefix sum"
      bucket_offset(0, ilevel) = offset
      do ibucket = 1, 7
         bucket_offset(ibucket, ilevel) = bucket_offset(ibucket-1, ilevel) &
                                + bucket_count(ibucket-1, ilevel)
      end do

      ! Update permutation
      bucket_count( 0:7, ilevel) = 0
      do ipart = offset+1, offset+np
         ibucket=ibits(part_hkey(ipart, ikey),ibit1,3)
         bucket_count(ibucket, ilevel) = bucket_count(ibucket, ilevel) + 1
         ind = bucket_offset(ibucket, ilevel) + bucket_count(ibucket, ilevel)
         part_ind_permutation(ind) = ipart
      end do
      
      call apply_particle_permutation(offset, np, key_level)
      
      if (ilevel < final_level) then
         ibucket = 0
         do while (ibucket < 8)
            ! Store level dependent loop variable
            ibucket_level(ilevel) = ibucket
            dummy_arg1 = bucket_offset(ibucket, ilevel)
            dummy_arg2 = bucket_count(ibucket, ilevel)
            if (bucket_count(ibucket, ilevel) > 0) then
               call msd_counting_sort_3digits(dummy_arg1, dummy_arg2, ilevel + 1, key_level)
            end if
            ibucket = ibucket_level(ilevel)
            ibucket = ibucket + 1
         end do
      end if      
    end subroutine msd_counting_sort_3digits
  end subroutine msd_radix_sort_particles

  !########################################################################
  !########################################################################
  !########################################################################

  subroutine lsd_radix_sort_particles(offset, np, final_level, key_level)
    use pm_commons,     only: part_ind_permutation, part_ind_permutation2, part_hkey
    use amr_parameters, only: ndim
    implicit none
    integer, intent(in) :: offset, np, final_level, key_level

    ! This routine sorts a contiguous sequence of the pariticle 
    ! array (specified by offset and np) by ascending hilbert keys

    integer, save       :: ilevel, ip
    
    if (ndim ==1)then
       print*, 'no 1D particle sorting case yet'
       stop
    end if

    if (final_level > key_level) then
       write(*,*)'you are trying to sort the hilbert keys to a too high level'
       call clean_stop
    end if
    
    do ip = 1, np
       part_ind_permutation(offset+ip) = ip
    end do

    do ilevel=final_level,1,-1
#if NDIM == 3
       call lsd_counting_sort_onelevel(np, ilevel, key_level, &
            part_ind_permutation(offset+1:offset+np), part_ind_permutation2(offset+1:offset+np), &
            part_hkey(offset+1:offset+np, 2), part_hkey(offset+1:offset+np, 1), part_hkey(offset+1:offset+np, 0))
#endif 
#if NDIM == 2
       call lsd_counting_sort_onelevel(np, ilevel, key_level, &
            part_ind_permutation(offset+1:offset+np), part_ind_permutation2(offset+1:offset+np), &
            part_hkey(offset+1:offset+np, 1), part_hkey(offset+1:offset+np, 0))
#endif
    end do

    do ip = 1, np
       part_ind_permutation(offset+ip) = part_ind_permutation(offset+ip) + offset
    end do

    call apply_particle_permutation(offset,np, key_level)

  end subroutine lsd_radix_sort_particles

  !########################################################################
  !########################################################################
  !########################################################################
#if NDIM == 3
  subroutine lsd_counting_sort_onelevel(np, ilevel, key_level, sigma1, sigma2, hkey2, hkey1, hkey0)      
    implicit none
    integer(kind=8), dimension(1:np), intent(in), target :: hkey2, hkey1, hkey0
#endif

#if NDIM == 2
  subroutine lsd_counting_sort_onelevel(np, ilevel, key_level, sigma1, sigma2, hkey1, hkey0)      
    implicit none
    integer(kind=8), dimension(1:np), intent(in), target :: hkey1, hkey0
#endif

#if NDIM == 1
  subroutine lsd_counting_sort_onelevel(np, ilevel, key_level, sigma1, sigma2, hkey0)
    implicit none
    integer(kind=8), dimension(1:np), intent(in), target :: hkey0
#endif

    integer,                          intent(in)         :: np, ilevel, key_level
    integer,         dimension(1:np), intent(inout)      :: sigma1, sigma2

    ! Update the given permuations sigma1, sigma2 such that sigma1  will 
    ! sort the 3 bits belonging to level ilevel of the input hilbert keys.
    ! For the 2D version, only two bits (corresponding to one level of refinement
    ! are read at the time.

    integer(kind=8),      pointer :: use_key(:)
    integer,                 save :: ipart, ip, ibit1, ikey
    integer(kind=8),         save :: ibucket

#if NDIM==3
    integer, dimension(0:7), save :: bucket_offset, bucket_count
    integer, parameter :: nbucket = 7
    integer, parameter :: nbits_read = 3

    ! get bit and key to read from 
    ibit1 = (key_level-ilevel)*3       
    ikey = ibit1/63
    ibit1 = mod(ibit1,63)

    ! a bit ugly, but move this line here as 
    ! hkey2 does not exist in 2d version of the routine
    if (ikey==2) use_key => hkey2

#endif

#if NDIM==2
    integer, dimension(0:3), save :: bucket_offset, bucket_count
    integer, parameter :: nbucket = 3
    integer, parameter :: nbits_read = 2
 
    ! get bit and key to read from 
    ibit1 = (key_level-ilevel)*2       
    ikey = ibit1/62
    ibit1 = mod(ibit1,62)
#endif

    ! use a pointer here to define which of the three integer keys                          
    ! must be accessed.                                                                     
    if (ikey==0) use_key => hkey0
    if (ikey==1) use_key => hkey1


    ! Count particles per bucket
    bucket_count=0
    do ipart = 1, np
       ibucket=ibits(use_key(ipart),ibit1,nbits_read)
       bucket_count(ibucket) = bucket_count(ibucket) + 1
    end do

    ! "Prefix sum"
    bucket_offset(0) = 0
    do ibucket = 1, nbucket
       bucket_offset(ibucket) = bucket_offset(ibucket-1) &
            + bucket_count(ibucket-1)
    end do

    ! Build up index permutation that will sort array
    do ip = 1, np
       ipart = sigma1(ip)
       ibucket=ibits(use_key(ipart),ibit1,nbits_read)
       bucket_offset(ibucket) = bucket_offset(ibucket) + 1
       sigma2(bucket_offset(ibucket))=ipart
    end do
    sigma1(1:np) = sigma2(1:np)
    
  end subroutine  lsd_counting_sort_onelevel
  !########################################################################
  !########################################################################
  !########################################################################
  
  subroutine apply_particle_permutation(offset, np, key_level)
    use amr_commons,    only: dp, ndim
    use pm_commons
    implicit none
    integer, intent(in)                          :: offset, np, key_level
    integer, save                                :: ipart, ikey, nkey, idim

    real(dp),        allocatable, dimension(:)   :: extra_storage_dp
    integer(kind=8), allocatable, dimension(:)   :: extra_storage_i8


    ! The permutation sigma built during the radix sort sorts the index array I such that
    ! key(sigma(i)) is sorted. We thus have to invert sigma and apply
    ! this to the key array to get (sigma^-1(key))(i) sorted in memory

    do ipart = offset + 1, offset + np
       part_ind_permutation2( part_ind_permutation(ipart) ) = ipart
    end do

    ! Apply  particle_permutation2

    ! Rearrange float arrays
    allocate(extra_storage_dp(offset + 1 : offset + np))

    do idim = 1, ndim       
       do ipart = offset + 1, offset + np
          extra_storage_dp(part_ind_permutation2(ipart)) = xp_andreas(ipart, idim) 
       end do
       xp_andreas(offset + 1 : offset + np, idim) = extra_storage_dp(offset + 1 : offset + np)
    end do

    do idim = 1, ndim       
       do ipart = offset + 1, offset + np
          extra_storage_dp(part_ind_permutation2(ipart)) = vp_andreas(ipart, idim) 
       end do
       vp_andreas(offset + 1 : offset + np, idim) = extra_storage_dp(offset + 1 : offset + np)
    end do

    do ipart = offset + 1, offset + np
       extra_storage_dp(part_ind_permutation2(ipart)) = mp_andreas(ipart) 
    end do
    mp_andreas(offset + 1 : offset + np) = extra_storage_dp(offset + 1 : offset + np)

    deallocate(extra_storage_dp)

    ! Rearrange long integer arrays
    allocate(extra_storage_i8(offset+1 : offset + np))

    nkey = 0
#if NDIM==3
    if (key_level > 21) nkey = nkey + 1
    if (key_level > 42) nkey = nkey + 1
#endif
#if NDIM==2
    if (key_level > 31) nkey = nkey + 1
#endif

    do ikey = 0, nkey  
       do ipart = offset + 1, offset + np
          extra_storage_i8(part_ind_permutation2(ipart)) = part_hkey(ipart, ikey) 
       end do
       part_hkey(offset + 1 : offset + np, ikey) = extra_storage_i8(offset + 1 : offset + np)
    end do
    
#ifdef LONGINT
    do ipart = offset + 1, offset + np
       extra_storage_i8(part_ind_permutation2(ipart)) = idp_andreas(ipart)
    end do
    idp_andreas(offset + 1 : offset + np) = extra_storage_i8(offset + 1 : offset + np)
#endif

    
    deallocate(extra_storage_i8)

    ! Rearrange short integer arrays
    ! use permanent permutation array as tmp storage space    
#ifndef LONGINT
    do ipart = offset + 1, offset + np
       part_ind_permutation(part_ind_permutation2(ipart)) = idp_andreas(ipart)
    end do
    idp_andreas(offset + 1 : offset + np) = part_ind_permutation(offset + 1 : offset + np)
#endif

    do ipart = offset + 1, offset + np
       part_ind_permutation(part_ind_permutation2(ipart)) = levelp_andreas(ipart)
    end do
    levelp_andreas(offset + 1 : offset + np) = part_ind_permutation(offset + 1 : offset + np)
    
    do ipart = offset + 1, offset + np
       part_ind_permutation(part_ind_permutation2(ipart)) = current_state(ipart)
    end do
    current_state(offset + 1 : offset + np) = part_ind_permutation(offset + 1 : offset + np)

  end subroutine apply_particle_permutation

  !########################################################################
  !########################################################################
  !########################################################################

end module sort

module hilbert

  ! Some parameters used in many of the modules routines

#if NDIM == 3
  integer, parameter :: big_shift = -60, left_shift = 4, right_shift = -1
#endif

#if NDIM == 2
  integer, parameter :: big_shift = -60, left_shift = 4, right_shift = -2
#endif
  
#if NDIM == 1
  integer, parameter :: big_shift = -62, left_shift = 2, right_shift = -1
#endif

  ! How many levels of refinement can be used per 64bit integer Hilbert key
  integer(kind=4), dimension(1:3), parameter :: levels_per_key = (/63, 31, 21/)
  integer(kind=4), dimension(1:3), parameter :: bits_per_int = (/63, 62, 63/)

  ! State diagrams taken from:
  ! J K Lawder, Using State Diagrams for Hilbert Curve Mappings
  ! http://www.dcs.bbk.ac.uk/TriStarp/pubs/JL2_00.pdf
  !
  ! Usage of state diagrams: 
  !
  ! - cartesian index: position of a cell inside 
  !   its oct in cartesian order (0 to 7)
  !
  ! - current state: integer encoding the "orientation" 
  !   of the hilbert curve inside the oct
  !
  ! - three_digit_diagram(cartesion index + current state * 8)  
  !   hilbert index of the cell given by cartesion index
  !
  ! - next_state_diagram(cartesian index + current state * 8)  
  !   orientation of curve inside the cell given by cartesian index

  !================================================================
  !================================================================
  !================================================================
  !================================================================

#if NDIM==3
  integer(kind=8),parameter,dimension(0:95)::next_digits_diagram=(/&
       &   0, 1, 3, 2, 7, 6, 4, 5,&
       &   0, 7, 1, 6, 3, 4, 2, 5,&
       &   0, 3, 7, 4, 1, 2, 6, 5,&
       &   2, 3, 1, 0, 5, 4, 6, 7,&
       &   4, 3, 5, 2, 7, 0, 6, 1,&
       &   6, 5, 1, 2, 7, 4, 0, 3,&
       &   4, 7, 3, 0, 5, 6, 2, 1,&
       &   6, 7, 5, 4, 1, 0, 2, 3,&
       &   2, 5, 3, 4, 1, 6, 0, 7,&
       &   2, 1, 5, 6, 3, 0, 4, 7,&
       &   4, 5, 7, 6, 3, 2, 0, 1,&
       &   6, 1, 7, 0, 5, 2, 4, 3 /)

  integer(kind=4),parameter,dimension(0:95)::next_state_diagram=(/&
       &   1, 2, 3, 2, 4, 5, 3, 5,&
       &   2, 6, 0, 7, 8, 8, 0, 7,&
       &   0, 9,10, 9, 1, 1,11,11,&
       &   6, 0, 6,11, 9, 0, 9, 8,&
       &  11,11, 0, 7, 5, 9, 0, 7,&
       &   4, 4, 8, 8, 0, 6,10, 6,&
       &   5, 7, 5, 3, 1, 1,11,11,&
       &   6, 1, 6,10, 9, 4, 9,10,&
       &  10, 3, 1, 1,10, 3, 5, 9,&
       &   4, 4, 8, 8, 2, 7, 2, 3,&
       &   7, 2,11, 2, 7, 5, 8, 5,&
       &  10, 3, 2, 6,10, 3, 4, 4 /)
    
  integer(kind=4),parameter,dimension(0:95, 1:3) :: one_digit_diagram=reshape((/&
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  

       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 1,  1,  1,  1,  0,  0,  0,  0,&  
       & 0,  0,  1,  1,  1,  1,  0,  0,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 1,  0,  0,  1,  1,  0,  0,  1,&  
       & 0,  0,  0,  0,  1,  1,  1,  1,&  
       & 1,  1,  0,  0,  0,  0,  1,  1,&  
       & 0,  1,  1,  0,  0,  1,  1,  0,&  
       & 1,  1,  1,  1,  0,  0,  0,  0/), (/96, 3/))  
  
  ! Next state diagram for reverse (hilbert to cartesian key) conversion
  integer(kind=4),parameter,dimension(0:95)::next_state_diagram_reverse=(/&
       & 1,   2,   2,   3,   3,   5,   5,   4,&  
       & 2,   0,   0,   8,   8,   7,   7,   6,&  
       & 0,   1,   1,   9,   9,  11,  11,  10,&  
       &11,   6,   6,   0,   0,   9,   9,   8,&  
       & 9,   7,   7,  11,  11,   0,   0,   5,&  
       &10,   8,   8,   6,   6,   4,   4,   0,&  
       & 3,  11,  11,   5,   5,   1,   1,   7,&  
       & 4,   9,   9,  10,  10,   6,   6,   1,&  
       & 5,  10,  10,   1,   1,   3,   3,   9,&  
       & 7,   4,   4,   2,   2,   8,   8,   3,&  
       & 8,   5,   5,   7,   7,   2,   2,  11,&  
       & 6,   3,   3,   4,   4,  10,  10,   2/)  

#endif

  !================================================================
  !================================================================
  !================================================================
  !================================================================

#if NDIM==2
  ! State diagrams for 2D case
  integer(kind=4),parameter,dimension(0:15)::next_state_diagram=(/&
       & 1, 0, 2, 0, &
       & 0, 3, 1, 1, &
       & 2, 2, 0, 3, &
       & 3, 1, 3, 2/)
  
  integer(kind=8),parameter,dimension(0:15)::next_digits_diagram=(/&
       & 0, 1, 3, 2, &
       & 0, 3, 1, 2, &
       & 2, 1, 3, 0, &
       & 2, 3, 1, 0/)

  integer(kind=4),parameter,dimension(0:15)::next_state_diagram_reverse=(/&
       & 1, 0, 0, 2, &
       & 0, 1, 1, 3, &
       & 3, 2, 2, 0, &
       & 2, 3, 3, 1/)

  integer(kind=4),parameter,dimension(0:15, 1:2)::one_digit_diagram=reshape((/&
       & 0, 0, 1, 1, &
       & 0, 1 ,1, 0, &
       & 1, 0, 0, 1, &
       & 1, 1, 0, 0, &
       & 0, 1, 1, 0, &
       & 0, 0 ,1, 1, &
       & 1, 1, 0, 0, &
       & 1, 0, 0, 1/), (/16,2/))
#endif

  !================================================================
  !================================================================
  !================================================================
  !================================================================

#if NDIM==1
  ! State diagrams for 1D case - a bit silly, but...
  integer(kind=4),parameter,dimension(0:1)::next_state_diagram=(/0, 0/)
  
  integer(kind=8),parameter,dimension(0:1)::next_digits_diagram=(/0, 1/)

  integer(kind=4),parameter,dimension(0:1)::next_state_diagram_reverse=(/0, 0/)

  integer(kind=4),parameter,dimension(0:1, 1:1)::one_digit_diagram=reshape((/0, 1/), (/2,1/))
#endif

contains
  
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function hilbert_key(ix, level) result(hkey)
    use amr_parameters, only: nhilbert, ndim, twotondim
    implicit none
    integer,                            intent(in) :: level
    integer(kind=8), dimension(1:ndim), intent(in) :: ix
    integer(kind=8), dimension(1:nhilbert) :: hkey

    ! Compute 3-integer hilbert keys from the cartesian keys ix

    integer(kind=4) :: cstate, sdigit, ind
    integer         :: ibit, add_digit, idim, ikey, nkey_local

    cstate = 0
    hkey = 0
    nkey_local = ceiling(1.d0 * level / levels_per_key(ndim))

    do ibit = level - 1, 0, -1
       if (.true.) then ! This prevents annoying bogus out of bounds message
       do ikey = nkey_local, 2, -1
          hkey(ikey) = ISHFT(hkey(ikey), left_shift)
          hkey(ikey) = ISHFT(hkey(ikey), right_shift)
          hkey(ikey) = hkey(ikey) + ISHFT(hkey(ikey - 1), big_shift)
       end do
       endif
       hkey(1) = ISHFT(hkey(1), left_shift)
       hkey(1) = ISHFT(hkey(1), right_shift)

       sdigit=0
       do idim = 1, ndim
          add_digit = 2 ** (ndim - idim)
          if(btest(ix(idim),ibit)) sdigit = sdigit + add_digit
       end do

       ind = cstate * twotondim + sdigit
       cstate = next_state_diagram(ind)
       hkey(1) = hkey(1) + next_digits_diagram(ind)
    enddo
    
  end function hilbert_key

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function hilbert_reverse(hkey, key_level) result(ix)
    use amr_parameters, only: nhilbert, ndim, twotondim
    implicit none
    integer,                                intent(in) :: key_level
    integer(kind=8), dimension(1:nhilbert), intent(in) :: hkey
    integer(kind=8), dimension(1:ndim) :: ix

    ! Compute cartesian keys from the corresponding 3-integer hilbert keys.

    integer         :: ibit, ikey, ilevel, idim
    integer(kind=4) :: cstate, nstate, sdigit, ind

    cstate = 0
    ix = 0

    do ilevel = 1, key_level
       ibit = (key_level-ilevel) * ndim
       ikey = ibit / bits_per_int(ndim) + 1
       ibit = mod(ibit, bits_per_int(ndim))
       do idim = 1, ndim
          ix(idim) = ISHFT(ix(idim),1)
       end do
       sdigit = int(ibits(hkey(ikey), ibit, ndim),kind=4)
       ind = cstate * twotondim + sdigit
       nstate = next_state_diagram_reverse(ind)
       do idim = 1, ndim
          ix(idim) = ix(idim) + one_digit_diagram(ind,idim)
       end do
       cstate = nstate
    enddo

  end function hilbert_reverse

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine hilbert_key_vec(ix, hkey, cstate, initial_level, final_level, npoint)
    use amr_parameters, only: nvector, nhilbert, ndim, twotondim
    implicit none
    integer, intent(in) :: initial_level, final_level, npoint
    integer(kind=8), dimension(:,:), intent(in)    :: ix
    integer(kind=4), dimension(:  ), intent(inout) :: cstate
    integer(kind=8), dimension(:,:), intent(inout) :: hkey
    
    ! Compute nvector 3-integer hilbert keys from the cartesian keys ix
        
    ! Local vars
    integer :: ibit, ip, add_digit, idim, ikey, nkey_local
    integer(kind=4), dimension(1:nvector) :: nstate, sdigit, ind

    ! if no keys present yet
    if (initial_level == 0) then
       cstate(1:npoint) = 0
       hkey(1:npoint, 1:nhilbert) = 0
    end if

    nkey_local = ceiling(1.d0 * final_level / levels_per_key(ndim))
    
    do ibit = final_level - initial_level - 1, 0, -1
       
       do ikey = nkey_local, 1, -1
          do ip = 1, npoint
             hkey(ip, ikey) = ISHFT(hkey(ip, ikey), left_shift)
             hkey(ip, ikey) = ISHFT(hkey(ip, ikey), right_shift)
          end do
          if (ikey > 1)then
             do ip=1,npoint
                hkey(ip, ikey) = hkey(ip, ikey) + ISHFT(hkey(ip, ikey - 1), big_shift)
             end do
          end if
       end do

       sdigit=0
       do idim = 1, ndim
          add_digit = 2 ** (ndim - idim)
          do ip = 1, npoint
             if(btest(ix(ip, idim),ibit)) sdigit(ip) = sdigit(ip) + add_digit
          end do
       end do

       do ip=1,npoint
          ind(ip) = cstate(ip) * twotondim + sdigit(ip)
       end do

       do ip=1,npoint
          nstate(ip) = next_state_diagram(ind(ip))
       end do

       do ip=1,npoint
          hkey(ip, 1) = hkey(ip, 1) + next_digits_diagram(ind(ip))
       end do

       do ip=1,npoint
          cstate(ip) = nstate(ip)
       end do
    enddo
  end subroutine hilbert_key_vec

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine hilbert_reverse_vec(ix, hkey, key_level, npoint)
    use amr_parameters, only: nvector, ndim, twotondim
    implicit none

    ! Inpu/Output variables:
    integer        , intent(in)                  :: key_level, npoint
    integer(kind=8), intent(out), dimension(:,:) :: ix
    integer(kind=8), intent(in),  dimension(:,:) :: hkey

    ! Descripton:
    ! Compute nvector cartesian keys from the corresponding 3-integer hilbert keys.

    ! pointer to one of the three hkey arrays
    integer                               :: ip, ibit1, ikey, ilevel, idim
    integer(kind=4), dimension(1:nvector) :: cstate, nstate, ind
    integer(kind=4), dimension(1:nvector) :: sdigit

    ! Build the cartesian key using the state diagrams
    cstate = 0; ix = 0

    do ilevel = 1, key_level
       ibit1 = (key_level-ilevel) * ndim
       ikey = ibit1 / bits_per_int(ndim) + 1
       ibit1 = mod(ibit1, bits_per_int(ndim))

       ! leftshift the cartesian keys by one position
       do idim = 1, ndim
          do ip = 1, npoint
             ix(ip, idim) = ISHFT(ix(ip, idim),1)
          end do
       end do

       ! read the next ndim bits from the hilbert key
       do ip = 1, npoint
          sdigit(ip) = int(ibits(hkey(ip, ikey), ibit1, ndim),kind=4)
       end do

       ! Compute lookup index in state diagrams
       do ip=1,npoint
          ind(ip) = cstate(ip) * twotondim + sdigit(ip)
       end do

       ! save next state
       do ip=1,npoint
          nstate(ip) = next_state_diagram_reverse(ind(ip))
       end do

       ! add one integer key digit each
       do idim = 1, ndim
          do ip=1,npoint
             ix(ip, idim) = ix(ip, idim) + one_digit_diagram(ind(ip), idim)
          end do
       end do

       do ip=1,npoint
          cstate(ip) = nstate(ip)
       end do
    enddo

  end subroutine hilbert_reverse_vec

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function ge_keys(key_a, key_b)
    implicit none
    logical::ge_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b

#if NHILBERT <= 1
    ge_keys = (key_a(1) >= key_b(1))
#endif
  
#if NHILBERT == 2
    if     (key_a(2) > key_b(2)) then
       ge_keys  =  .true.
    elseif (key_a(2) < key_b(2)) then
       ge_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       ge_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       ge_keys = .false.
    else
       ge_keys = .true.
    end if
#endif
    
#if NHILBERT == 3  
    if     (key_a(3) > key_b(3)) then
       ge_keys = .true.
    elseif (key_a(3) < key_b(3)) then
       ge_keys = .false.
    elseif (key_a(2) > key_b(2)) then
       ge_keys = .true.
    elseif (key_a(2) < key_b(2)) then
       ge_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       ge_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       ge_keys = .false.
    else
       ge_keys = .true.
    end if
#endif

  end function ge_keys
 
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function average_keys(key_a, key_b)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension(1:nhilbert) :: average_keys
    integer(kind=8), dimension(1:nhilbert) :: sum_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b
    ! This function assumes that key_a and key_b contain only positive integers

#if NHILBERT <= 1
    sum_keys(1) = key_a(1) + key_b(1)
    average_keys(1) = ishft(sum_keys(1),-1) ! Avoid overflow
#endif

#if NHILBERT == 2
    ! Multi-precision arithmetics
    sum_keys(1) = key_a(1) + key_b(1)
    sum_keys(2) = key_a(2) + key_b(2)
    if(sum_keys(2)==0)then
       average_keys(2) = 0
       average_keys(1) = ishft(sum_keys(1),-1) ! Avoid overflow
    else
       sum_keys(2) = sum_keys(2) + ishft(sum_keys(1),-63) ! Add the carry to the MS sum
       sum_keys(1) = ibits(sum_keys(1),0,63)              ! Compute the LS sum
       average_keys(2) = ishft(sum_keys(2),-1)            ! Avoid overflow
       sum_keys(1) = sum_keys(1) + ishft(sum_keys(2), 63) ! Add back the shift
       average_keys(1) = ishft(sum_keys(1),-1)            ! Avoid overflow
    endif
#endif

  end function average_keys

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function difference_keys(key_a, key_b)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension(1:nhilbert) :: difference_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b
    ! This function assumes that key_a is always greater or equal than key_b
    ! and that key_a and key_b contain only positive integers

#if NHILBERT <= 1
    difference_keys(1) = key_a(1) - key_b(1)
#endif

#if NHILBERT == 2
    ! Multi-precision arithmetics
    difference_keys(1) = key_a(1) - key_b(1)
    if(difference_keys(1).GE.0)then
       difference_keys(2) = key_a(2) - key_b(2)
    else
       difference_keys(1) = difference_keys(1) - ishft(1_8, 63)
       difference_keys(2) = key_a(2) - key_b(2) - 1_8
    endif
#endif

  end function difference_keys

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function gt_keys(key_a, key_b)
    implicit none
    logical::gt_keys
    integer(kind=8), intent(in), dimension (:):: key_a, key_b

#if NHILBERT <= 1
    gt_keys = (key_a(1) > key_b(1))
#endif

#if NHILBERT == 2
    if     (key_a(2) > key_b(2)) then
       gt_keys  =  .true.
    elseif (key_a(2) < key_b(2)) then
       gt_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       gt_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       gt_keys = .false.
    else
       gt_keys = .false.
    end if
#endif

#if NHILBERT == 3
    if     (key_a(3) > key_b(3)) then
       gt_keys = .true.
    elseif (key_a(3) < key_b(3)) then
       gt_keys = .false.
    elseif (key_a(2) > key_b(2)) then
       gt_keys = .true.
    elseif (key_a(2) < key_b(2)) then
       gt_keys = .false.
    elseif (key_a(1) > key_b(1)) then
       gt_keys = .true.
    elseif (key_a(1) < key_b(1)) then
       gt_keys = .false.
    else
       gt_keys = .false.
    end if
#endif
  end function gt_keys

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function eq_keys(key_a, key_b)
    implicit none
    logical::eq_keys
    integer(kind=8), intent(in), dimension (:):: key_a, key_b

#if NHILBERT <= 1
    eq_keys = (key_a(1) == key_b(1))
#endif

#if NHILBERT == 2
    eq_keys = (key_a(1) == key_b(1)) .AND. (key_a(2) == key_b(2))
#endif

#if NHILBERT == 3
    eq_keys = (key_a(1) == key_b(1)) .AND. (key_a(2) == key_b(2)) .AND. (key_a(3) == key_b(3))
#endif

  end function eq_keys
  
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function refine_key(key_in,key_level)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension (1:nhilbert) :: refine_key
    integer(kind=8), intent(in) , dimension (1:nhilbert) :: key_in
    integer(kind=4), intent(in) :: key_level

    integer :: ikey, nkey_local

    refine_key(1:nhilbert) = key_in(1:nhilbert) 
    nkey_local = ceiling(1.d0 * key_level / levels_per_key(ndim))
    do ikey = nkey_local, 1, -1
       refine_key(ikey) = ISHFT(refine_key(ikey), left_shift)
       refine_key(ikey) = ISHFT(refine_key(ikey), right_shift)
       if(ikey>1)then
          refine_key(ikey) = refine_key(ikey) + ISHFT(refine_key(ikey-1),big_shift)
       endif
    end do
    
  end function refine_key

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  pure function coarsen_key(key_in,key_level)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension (1:nhilbert) :: coarsen_key
    integer(kind=8), intent(in) , dimension (1:nhilbert) :: key_in
    integer(kind=4), intent(in) :: key_level

    integer :: ikey, nkey_local

    coarsen_key(1:nhilbert) = key_in(1:nhilbert) 
    nkey_local = ceiling(1.d0 * key_level / levels_per_key(ndim))
    do ikey = 1, nkey_local
       coarsen_key(ikey) = ISHFT(coarsen_key(ikey),-ndim)
       if(ikey<nkey_local)then
          coarsen_key(ikey) = coarsen_key(ikey) + ISHFT(ISHFT(coarsen_key(ikey+1),64-ndim),bits_per_int(ndim)-64)
       endif
    end do
    
  end function coarsen_key

  !================================================================
  !================================================================
  !================================================================
  !================================================================

end module hilbert

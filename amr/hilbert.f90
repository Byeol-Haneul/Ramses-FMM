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

  ! Usage of state diagrams: 
  
  ! - cartesian index: position of a cell inside 
  !   its oct in cartesian order (0 to 7)
  
  ! - current state: integer encoding the "orientation" 
  !   of the hilbert curve inside the oct
  
  ! - three_digit_diagram(cartesion index + current state * 8)  
  !   hilbert index of the cell given by cartesion index
  
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

  subroutine hilbert_key(ix, hkey, cstate, initial_level, final_level, npoint)
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
  end subroutine hilbert_key

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  subroutine hilbert_reverse(ix, hkey, key_level, npoint)
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
    integer(kind=8), dimension(1:nvector) :: sdigit

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
          sdigit(ip) = ibits(hkey(ip, ikey), ibit1, ndim)
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

  end subroutine hilbert_reverse

  !================================================================
  !================================================================
  !================================================================
  !================================================================

!!$  subroutine hilbert_for_particle(offset, nparts, initial_level, final_level)
!!$    use amr_parameters, only: nvector, boxlen, dp, ndim, nhilbert
!!$    use pm_commons,     only: part_hkey, current_state, xp
!!$    implicit none
!!$
!!$    integer, intent(in) :: initial_level, final_level
!!$    integer, intent(in) :: offset, nparts
!!$
!!$    ! Description:
!!$    ! This subroutine computes 3D hilbert keys for particles 
!!$    ! It assumes that the particles are stored as contiguous 
!!$    ! arrays in memory and that positions, 3-integer hilbert keys
!!$    ! and next_state are allocated as particle-based quantities.
!!$
!!$    ! Iputs: 
!!$    ! - Starting offset in particle arrays and number of particles
!!$    !   to process (np)
!!$    ! - Level of already computed hilbert key 
!!$    ! - Desired level of hilbert key on exit
!!$
!!$    ! Example: 
!!$    ! call hilbert_for_particle(0, npart_levelmin, nlevelmax-1, nlevelmax) 
!!$    ! will compute the last 3 bits of the hilbert key for the
!!$    ! levelmin particles (resulting in a total of 3 * nlevelmax bits)
!!$
!!$    ! Local variables
!!$    integer :: ibit, ip, ind_part, idim, np, ioft
!!$    integer(kind=8), dimension(1:nvector, 1:ndim) :: ix
!!$    real(dp) :: ckey_factor
!!$    
!!$    ! Compute particle position to cartesian key factor
!!$    ckey_factor = 2.0**final_level / dble(boxlen)
!!$
!!$    do ioft = offset, offset + nparts - 1, nvector
!!$       np = min(nvector, offset + nparts - ioft)
!!$
!!$       ! compute cartesian keys
!!$       do idim = 1, ndim
!!$          do ip = 1, np
!!$             ix(ip, idim) = floor(xp(ioft + ip, idim) * ckey_factor, kind=8)
!!$          end do
!!$       end do
!!$
!!$       ! Passing in array slices is ok (no copying) if the dummy argument has assumed shape and the interface is explicit!
!!$       call hilbert_key(ix, part_hkey(ioft + 1: ioft + np, 1:nhilbert), current_state(ioft + 1: ioft + np), initial_level, final_level, np)
!!$    end do
!!$
!!$  end subroutine hilbert_for_particle

  !================================================================
  !================================================================
  !================================================================
  !================================================================

  recursive subroutine sort_hilbert(head_part, tail_part, ix_coarse, cstate_coarse, ilevel, final_level)
    use amr_commons, only: boxlen, dp, ndim, twotondim, myid
    use pm_commons,  only: workp, sortp, xp
    implicit none

    integer, intent(in) :: ilevel, final_level
    integer, intent(in) :: head_part, tail_part
    integer, dimension(1:ndim), intent(in) :: ix_coarse
    integer, intent(in) :: cstate_coarse

    ! Description:
    ! This subroutine sort particles along the Hilbert key at the resolution
    ! set by final_level. It should be called first with ilevel=1
    ! arrays in memory and that positions, 3-integer hilbert keys
    ! and next_state are allocated as particle-based quantities.

    ! Iputs: 
    ! - Head_part and tail_part are head and tail of particle distribution to work on.
    ! - Array sortp must be initialized with sortp(i)=i between head_part and tail_part.
    ! - Cartesian key of coarse cell in which these particles are contained.
    ! - State of the coarse cell for Hilbert ordering
    ! - Current and final level

    ! Example: 
    ! ix=(/0,0,0/)
    ! call sort_hilbert(1, npart, ix, 0, 1, nlevelmax) 
    ! will sort all particles according to their Hilbert key at levelmax.
    ! On output, array sortp is modified.

    ! Local variables
    integer :: ibit, ip, ind_part, idim, np, ioft, ipart, new_ipart
    integer :: ckey_max, cstate_fine, ind_cart_part, head_fine, tail_fine
    real(dp) :: ckey_factor
    integer, dimension(1:ndim) :: ix_fine, ix_ref, ix_part
    integer, dimension(0:twotondim-1,1:ndim) :: ix, ix_child
    integer, dimension(0:twotondim-1) :: nstate, sdigit, ind, ind_cart, ind_hilbert
    integer, dimension(0:twotondim-1) :: numb_part, offset

    ! Compute particle position to cartesian key factor
    ckey_max = 2**ilevel
    ckey_factor = 2.0**ilevel / dble(boxlen)

    ! Initial Cartesian offset for fine cells
    do idim = 1, ndim
       ix_ref(idim) = ISHFT(ix_coarse(idim),1)
    end do
    
    ! Compute the Hilbert index for fine cells
    do ip = 0, twotondim-1
       sdigit(ip) = ip
    end do

    ! Compute lookup index in state diagrams
    do ip = 0, twotondim-1
       ind(ip) = cstate_coarse * twotondim + sdigit(ip)
    end do

    ! Save next state
    do ip = 0, twotondim-1
       nstate(ip) = next_state_diagram_reverse(ind(ip))
    end do

    ! Add one integer key digit each
    do idim = 1, ndim
       do ip = 0, twotondim-1
          ix(ip, idim) = one_digit_diagram(ind(ip), idim)
       end do
    end do

    ! Compute Cartesian index for children cells
    ind_cart = 0
    do idim = 1, ndim
       do ip = 0, twotondim-1
          ix_child(ip, idim) = ix_ref(idim) + ix(ip, idim)
          ind_cart(ip) = ind_cart(ip) + ix(ip, idim) * 2**(idim-1)
       end do
    end do

    ! Compute mapping from Cartesian to Hilbert order
    ind_hilbert = 0
    do ip = 0, twotondim-1
       ind_hilbert(ind_cart(ip))=ip
    end do

    ! Count particles per children cell
    numb_part = 0
    do ipart = head_part, tail_part
       ind_part = sortp(ipart)
       ind_cart_part = 0
       do idim = 1,ndim
          ix_part(idim) = int(xp(ind_part,idim)*ckey_factor) - ix_ref(idim)
          ind_cart_part = ind_cart_part + ix_part(idim) * 2**(idim-1)
       end do
       ip = ind_hilbert(ind_cart_part)
       numb_part(ip) = numb_part(ip) + 1
    end do
    
    offset = head_part-1
    do ip = 1, twotondim-1
       offset(ip) = offset(ip-1) + numb_part(ip-1)
    end do

    ! Compute new sortp array
    numb_part = 0
    do ipart = head_part, tail_part
       ind_part = sortp(ipart)
       ind_cart_part = 0
       do idim = 1,ndim
          ix_part(idim) = int(xp(ind_part,idim)*ckey_factor) - ix_ref(idim)
          ind_cart_part = ind_cart_part + ix_part(idim) * 2**(idim-1)
       end do
       ip = ind_hilbert(ind_cart_part)
       numb_part(ip) = numb_part(ip) + 1
       new_ipart = offset(ip) + numb_part(ip)
       workp(new_ipart) = ind_part
    end do
    do ipart = head_part,tail_part
       sortp(ipart) = workp(ipart)
    end do

    ! Recursive call
    if(ilevel < final_level)then
       do ip = 0, twotondim-1
          if(numb_part(ip) > 0)then
             head_fine = offset(ip) + 1
             tail_fine = offset(ip) + numb_part(ip)
             ix_fine(1:ndim) = ix_child(ip,1:ndim)
             cstate_fine = nstate(ip)
             call sort_hilbert(head_fine,tail_fine,ix_fine,cstate_fine,ilevel+1,final_level)
          endif
       end do
    endif

  end subroutine sort_hilbert
  
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function ge_keys(key_a, key_b)
    implicit none
    logical::ge_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b

#if NHILBERT == 1
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

  function average_keys(key_a, key_b)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension(1:nhilbert) :: average_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b

#if NHILBERT == 1
    average_keys(1) = (key_a(1) + key_b(1)) / 2
#endif
  
  end function average_keys
 
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function difference_keys(key_a, key_b)
    use amr_parameters, only: ndim, nhilbert
    implicit none
    integer(kind=8), dimension(1:nhilbert) :: difference_keys
    integer(kind=8), intent(in), dimension(:) :: key_a, key_b

#if NHILBERT == 1
    difference_keys(1) = (key_a(1) - key_b(1))
#endif
  
  end function difference_keys
 
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  function gt_keys(key_a, key_b)
    implicit none
    logical::gt_keys
    integer(kind=8), intent(in), dimension (:):: key_a, key_b

#if NHILBERT == 1
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

  function eq_keys(key_a, key_b)
    implicit none
    logical::eq_keys
    integer(kind=8), intent(in), dimension (:):: key_a, key_b

#if NHILBERT == 1
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

  function refine_key(key_in,key_level)
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

  function coarsen_key(key_in,key_level)
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

end module hilbert

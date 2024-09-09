! Hash table module for the use inside RAMSES.

! - KEY: A tuple (ilevel, ix, iy, iz) acts as hash key.

! - VALUE: Integer (typically grid indices) are stored in the hash table.

! - HASH FUNCTION: Simple hash function based on multiplication with constants.
! - interface to murmur3 hash exists but is not used anymore

! - COLLISIONS: A linked list is used to deal with collisions.

module hash
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
  use amr_parameters, only: ndim, nvector
  implicit none

  ! General module parameters
  integer, dimension(0:3), parameter :: constants = (/5, -1640531527, 97, 1003313/)

  ! Define a bucket as a derived type (sequence statement!) for better
  ! cache efficiency.
  type,bind(c) :: bucket
     !sequence
     integer(kind=8), dimension(0:ndim) :: key
     type(C_PTR) :: valuep
     integer :: next_ibucket
  end type bucket

  type bucket_simple
     integer(kind=8) :: key
     integer :: valuep
     integer :: next_ibucket
  end type bucket_simple

  ! The actual hash table is an array of buckets
  type hash_table
     type(bucket), allocatable, dimension(:)  :: data
     integer         :: total_size, head_free, nfree_chain, nfree
     integer(kind=8) :: size
     integer(kind=8) :: bitmask
     integer, allocatable, dimension(:) :: next_free
  end type hash_table  

  type hash_simple
     type(bucket_simple), allocatable, dimension(:)  :: data
     integer :: total_size, head_free, nfree_chain, nfree
     integer :: size
     integer, allocatable, dimension(:) :: next_free
  end type hash_simple

  ! Prime numbers for hash table
  integer,dimension(0:30) :: prime=(/2,3,7,13,23,53,97,193,389,769,1543, &
       & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
       & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
       & 402653189,805306457,1610612741/)

contains

  ! ============================================================================= 
  pure function hash_func(key)
    integer(kind=8), dimension(0:ndim), intent(in) :: key
    integer(kind=8)                                :: hash_func
    hash_func = dot_product(key(0:ndim), constants(0:ndim))
  end function hash_func
  ! =============================================================================

  ! =============================================================================
  subroutine init_empty_hash(htable, req_size, hash_type)
    implicit none
    type(hash_table), intent(inout) :: htable
    integer         , intent(in)    :: req_size
    character(6)    , intent(in)    :: hash_type

    ! Allocate all hash table arrays and variables.
    ! Chose size (excluding the chaining space) as the smallest
    ! power of two >= the required_size.

    if (hash_type .ne. 'simple') then
       print*, 'only simple hash is currently supported'
       stop
    end if

    htable%size = 2
    do while (htable%size < req_size)
       htable%size = htable%size * 2
    end do

    call reset_entire_hash(htable, .false.)

  end subroutine init_empty_hash

  subroutine init_empty_hash_simple(htable, req_size)
    implicit none
    type(hash_simple), intent(inout) :: htable
    integer          , intent(in)    :: req_size

    ! Allocate all hash table arrays and variables.

    integer :: bit_length, ncode

    ncode = req_size
    do bit_length=1,32
       ncode=ncode/2
       if(ncode<=1) exit
    end do
    htable%size=prime(bit_length+1)

    call reset_entire_hash_simple(htable)

  end subroutine init_empty_hash_simple
  ! =============================================================================

  ! =============================================================================
  subroutine reset_entire_hash(htable, resize)
    implicit none
    logical, intent(in)             :: resize
    type(hash_table), intent(inout) :: htable
    
    ! Subroutine to reset the entire hash table
    ! IMPORTANT: The new size of the hash table is adapted based on the
    ! load factor before resetting the hash table.

    integer :: i
    real :: load_factor

    if (resize) then
       load_factor = (htable%size - htable%nfree) * 1.0 / htable%size    
       if (load_factor > 0.6) then
          htable%size = htable%size * 2
          deallocate(htable%data, htable%next_free)
       else if (load_factor < 0.2 .and. htable%size > 2)then
          htable%size = htable%size / 2
          deallocate(htable%data, htable%next_free)
       end if
    end if
    
    ! Compute sizes and allocate arrays
    htable%total_size = int(htable%size/4,kind=4) + int(htable%size,kind=4)
    htable%nfree = int(htable%size,kind=4)
    htable%nfree_chain = htable%total_size - int(htable%size,kind=4)
    htable%head_free = int(htable%size,kind=4) + 1
    htable%bitmask = htable%size - 1

    if (.not. allocated(htable%data))then
       allocate(htable%data(1:htable%total_size))
       allocate(htable%next_free(htable%size+1:htable%total_size))
    end if

    ! Initialize data
    do i = 1, int(htable%total_size,kind=4)
       call reset_bucket(htable%data(i))
    end do
    do i = int(htable%size,kind=4) + 1, htable%total_size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%total_size) = 0

  end subroutine reset_entire_hash

  subroutine reset_entire_hash_simple(htable)
    implicit none
    type(hash_simple), intent(inout) :: htable

    ! Subroutine to reset the entire hash table

    integer :: i

    ! Compute sizes and allocate arrays
    htable%total_size = int(htable%size/4,kind=4) + int(htable%size,kind=4)
    htable%nfree = int(htable%size,kind=4)
    htable%nfree_chain = htable%total_size - int(htable%size,kind=4)
    htable%head_free = int(htable%size,kind=4) + 1

    if (.not. allocated(htable%data))then
       allocate(htable%data(1:htable%total_size))
       allocate(htable%next_free(htable%size+1:htable%total_size))
    end if

    ! Initialize data
    do i = 1, htable%total_size
       call reset_bucket_simple(htable%data(i))
    end do
    do i = htable%size + 1, htable%total_size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%total_size) = 0

  end subroutine reset_entire_hash_simple
  ! =============================================================================

  ! =============================================================================
  subroutine reset_bucket(buck)
    implicit none
    type(bucket), intent(inout) :: buck

    ! Reset the content of a bucket

    buck%next_ibucket = -1
    buck%key = 0
  end subroutine reset_bucket

  subroutine reset_bucket_simple(buck)
    implicit none
    type(bucket_simple), intent(inout) :: buck

    ! Reset the content of a bucket

    buck%next_ibucket = -1
    buck%key = 0
  end subroutine reset_bucket_simple
  ! =============================================================================

  ! =============================================================================
  subroutine hash_setp(htable, key, val)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_NULL_PTR
    implicit none
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key
    type(*),optional,target,             intent(in)    :: val
    
    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer(kind=8) :: full_hash
    integer(kind=8) :: ibucket

    ! Compute ibucket
    full_hash = hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

    if (htable%data(ibucket)%next_ibucket < 0) then          

       ! Bucket is empty, simply insert value       
       htable%data(ibucket)%next_ibucket = 0
       if (present(val)) then
         htable%data(ibucket)%valuep       = C_LOC(val)
       else
         htable%data(ibucket)%valuep       = C_NULL_PTR
       endif
       htable%data(ibucket)%key(0:ndim) = key(0:ndim)
       htable%nfree = htable%nfree - 1

    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%data(ibucket)%next_ibucket .ne. 0)
          ! Check if key already exists - abort if so
          if (same_keys(htable%data(ibucket)%key(0:ndim),key(0:ndim)))then
             write(*,*) "trying to insert already existing key: ",key !htable%data(ibucket)%valuep
             write(*,*) "existing key: ", htable%data(ibucket)%key(0:ndim)
             stop
          end if
          ibucket = htable%data(ibucket)%next_ibucket
       end do

       ! Check again (at the end of linked list)
       if (same_keys(htable%data(ibucket)%key(0:ndim),key(0:ndim)))then
          write(*,*) "trying to insert already existing key: ",key
          CALL ABORT
          stop
       end if
       
       ! Have reached end of chain, val not present yet -> add
       htable%data(ibucket)%next_ibucket = htable%head_free
       ibucket = htable%head_free
       htable%data(ibucket)%next_ibucket = 0
       if (present(val)) then
         htable%data(ibucket)%valuep       = C_LOC(val)
       else
         htable%data(ibucket)%valuep       = C_NULL_PTR
       endif
       htable%data(ibucket)%key(0:ndim) = key(0:ndim)

       ! remove bucket from head of free linked list
       htable%head_free   = htable%next_free(htable%head_free)
       htable%nfree_chain = htable%nfree_chain - 1

    else
       write(*,*)"hash chaining space full "
       stop
    end if
  end subroutine hash_setp

  subroutine hash_setp_simple(htable, key, val)
    implicit none
    type(hash_simple), intent(inout) :: htable
    integer(kind=8),   intent(in)    :: key
    integer,           intent(in)    :: val

    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer :: ibucket

    ! Compute ibucket
    ibucket = MOD(key,int(htable%size,kind=8)) + 1

    if (htable%data(ibucket)%next_ibucket < 0) then

       ! Bucket is empty, simply insert value
       htable%data(ibucket)%valuep = val
       htable%data(ibucket)%next_ibucket = 0
       htable%data(ibucket)%key = key
       htable%nfree = htable%nfree - 1

    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%data(ibucket)%next_ibucket .ne. 0)
          ! Check if key already exists - abort if so
          if (htable%data(ibucket)%key == key)then
             write(*,*) "trying to insert already existing key: ",key
             write(*,*) "existing key: ", htable%data(ibucket)%key
             stop
          end if
          ibucket = htable%data(ibucket)%next_ibucket
       end do

       ! Check again (at the end of linked list)
       if (htable%data(ibucket)%key == key)then
          write(*,*) "trying to insert already existing key: ",key
          CALL ABORT
          stop
       end if

       ! Have reached end of chain, val not present yet -> add
       htable%data(ibucket)%next_ibucket = htable%head_free
       ibucket = htable%head_free
       htable%data(ibucket)%valuep = val
       htable%data(ibucket)%next_ibucket = 0
       htable%data(ibucket)%key = key

       ! remove bucket from head of free linked list
       htable%head_free   = htable%next_free(htable%head_free)
       htable%nfree_chain = htable%nfree_chain - 1

    else
       write(*,*)"hash chaining space full "
       stop
    end if
  end subroutine hash_setp_simple
  ! =============================================================================

  ! =============================================================================
  function hash_getp(htable, key, absent)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR, C_ASSOCIATED
    USE oct_commons
    implicit none
    type(hash_table),                    intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in) :: key
    type(C_PTR)                                     :: hash_getp
    logical,optional,                    intent(out):: absent
    ! Function which retrieves the hash table value for a given key. If no entry exists, return 0.
    integer(kind=8) :: ibucket, full_hash

    if(present(absent)) absent = .false.
    
    full_hash = hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

    if (same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))then
       hash_getp = htable%data(ibucket)%valuep
       if(present(absent) .and. .not.C_ASSOCIATED(hash_getp)) absent = .true.
       return
    end if
    
    ! Walk linked list until key is found or to the end is reached
    do while( htable%data(ibucket)%next_ibucket > 0)
       ibucket = htable%data(ibucket)%next_ibucket
       if (same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))then
          hash_getp = htable%data(ibucket)%valuep
          if(present(absent) .and. .not.C_ASSOCIATED(hash_getp)) absent = .true.
          return
       end if
    end do

    ! Nothing found...
    hash_getp = C_NULL_PTR

  end function hash_getp

  function hash_getp_simple(htable, key)
    USE oct_commons
    implicit none
    type(hash_simple), intent(in) :: htable
    integer(kind=8),   intent(in) :: key
    integer :: hash_getp_simple

    ! Function which retrieves the hash table value for a given key. If no entry exists, return 0.

    integer :: ibucket

    ibucket = MOD(key,int(htable%size,kind=8)) + 1

    if (htable%data(ibucket)%key == key)then
       hash_getp_simple = htable%data(ibucket)%valuep
       return
    end if

    ! Walk linked list until key is found or to the end is reached
    do while( htable%data(ibucket)%next_ibucket > 0)
       ibucket = htable%data(ibucket)%next_ibucket
       if (htable%data(ibucket)%key == key)then
          hash_getp_simple = htable%data(ibucket)%valuep
          return
       end if
    end do

    ! Nothing found...
    hash_getp_simple = 0

  end function hash_getp_simple
  ! =============================================================================

  ! =============================================================================`
  function hash_getp_vec(htable, keys, absent, n)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR, C_ASSOCIATED
    use oct_commons
    implicit none
    type(hash_table),                                 intent(in) :: htable
    integer(kind = 8) , dimension(1:nvector, 0:ndim), intent(in) :: keys
    type(C_PTR),dimension(1:nvector)                             :: hash_getp_vec
    logical, optional,dimension(1:nvector)       , intent(inout) :: absent
    integer                                                      :: n
    ! Function to obtain up to nvector values from the hash key at once.
    ! This subroutine is only valid if the simple hash is used.
    integer(kind = 8), dimension(1:nvector)         :: ibucket, full_hash
    integer(kind = 8), dimension(1:nvector, 0:ndim) :: bucket_keys
    logical,           dimension(1:nvector)         :: ok
    integer :: i, idim, n_coll

    if(present(absent)) absent = .false.

    full_hash = 0
    do idim = 0, ndim
       do i = 1, n
          full_hash(i) = full_hash(i) + keys(i, idim) * constants(idim)
       end do
    end do

    do i = 1, n
       ibucket(i) = IAND(full_hash(i), htable%bitmask) + 1
    end do
    
    do idim = 0, ndim
       do i = 1, n
          bucket_keys(i, idim) = htable%data(ibucket(i))%key(idim)
       end do
    end do

    ok = .true.
    do idim = 0, ndim
       do i = 1, n
          ok(i) = ok(i) .and. (bucket_keys(i, idim) == keys(i, idim))
       end do
    end do

    n_coll = 0
    do i = 1, n
       if (ok(i)) then
          hash_getp_vec(i) = htable%data(ibucket(i))%valuep
       else
          n_coll = n_coll + 1
       endif
    end do

    if(present(absent))then
       do i = 1, n
          if(ok(i) .and. .not. C_ASSOCIATED(hash_getp_vec(i))) absent(i) = .true.
       end do
    end if

    if(n_coll == 0)return

    do i = 1, n
       if (.not. ok(i)) then
          ! Walk linked list until key is found or to the end is reached
          do while( htable%data(ibucket(i))%next_ibucket > 0)
             ibucket(i) = htable%data(ibucket(i))%next_ibucket
             if (same_keys(htable%data(ibucket(i))%key(0:ndim), keys(i,0:ndim)))then
                hash_getp_vec(i) = htable%data(ibucket(i))%valuep
                ok(i) = .true.
                if(present(absent) .and. .not. C_ASSOCIATED(hash_getp_vec(i))) absent(i)=.true.
             end if
          end do
          ! Nothing found...
          if (.not. ok(i))then
             hash_getp_vec(i) = C_NULL_PTR
          end if
       end if
    end do

  end function hash_getp_vec
  ! =============================================================================  
  
  ! =============================================================================  
  subroutine hash_free(htable, key)
    implicit none
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key

    ! Remove the hash table entry for a given key 

    integer(kind=8) :: ibucket, previous_ibucket=0, full_hash

    full_hash = hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

    ! No collision case
    if (htable%data(ibucket)%next_ibucket == 0) then     
       htable%data(ibucket)%next_ibucket = -1
       htable%data(ibucket)%key(0:ndim) = 0
       htable%nfree = htable%nfree + 1
    else
       ! Collision case
       do while (.not. same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))
          previous_ibucket=ibucket
          ibucket=htable%data(ibucket)%next_ibucket
       end do
       if (ibucket <= htable%size) then           
          ! It's the first element we need to erase: Move first element from chaning 
          ! space into bucket and do as if the value to remove had been in the chaning space
          htable%data(ibucket)%valuep = htable%data(htable%data(ibucket)%next_ibucket)%valuep
          htable%data(ibucket)%key = htable%data(htable%data(ibucket)%next_ibucket)%key
          previous_ibucket = ibucket
          ibucket = htable%data(ibucket)%next_ibucket
       end if
       ! fill the hole and reconnect linked list
       htable%data(previous_ibucket)%next_ibucket = htable%data(ibucket)%next_ibucket
       htable%next_free(ibucket) = htable%head_free
       htable%head_free = int(ibucket,kind=4)
       htable%nfree_chain = htable%nfree_chain + 1
    end if
  end subroutine hash_free

  subroutine hash_free_simple(htable, key)
    implicit none
    type(hash_simple), intent(inout) :: htable
    integer(kind=8), intent(in)      :: key

    ! Remove the hash table entry for a given key

    integer :: ibucket, previous_ibucket=0

    ibucket = MOD(key,int(htable%size,kind=8)) + 1

    ! No collision case
    if (htable%data(ibucket)%next_ibucket == 0) then
       htable%data(ibucket)%next_ibucket = -1
       htable%data(ibucket)%key = 0
       htable%nfree = htable%nfree + 1
    else
       ! Collision case
       do while (.not. (htable%data(ibucket)%key == key))
          previous_ibucket=ibucket
          ibucket=htable%data(ibucket)%next_ibucket
       end do
       if (ibucket <= htable%size)then
          ! It's the first element we need to erase: Move first element from chaining
          ! space into bucket and do as if the value to remove had been in the chaining space
          htable%data(ibucket)%valuep = htable%data(htable%data(ibucket)%next_ibucket)%valuep
          htable%data(ibucket)%key = htable%data(htable%data(ibucket)%next_ibucket)%key
          previous_ibucket = ibucket
          ibucket = htable%data(ibucket)%next_ibucket
       end if
       ! fill the hole and reconnect linked list
       htable%data(previous_ibucket)%next_ibucket = htable%data(ibucket)%next_ibucket
       htable%next_free(ibucket) = htable%head_free
       htable%head_free = ibucket
       htable%nfree_chain = htable%nfree_chain + 1
    end if
  end subroutine hash_free_simple
! =============================================================================
  
  ! =============================================================================
  pure function same_keys(key1, key2)
    logical :: same_keys
    integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2
    logical, dimension(0:ndim)                     :: ok
    integer                                        :: i
    do i = 0, ndim
       ok(i) = (key1(i)==key2(i))
    end do
    same_keys = ALL(ok)
  end function same_keys
  ! =============================================================================
  
  ! =============================================================================
  subroutine hash_stats(htable)
    implicit none
    type(hash_table)::htable

    write(*,*)"Total values stored in hash table: "&
         ,htable%total_size - htable%nfree - htable%nfree_chain
    write(*,*)"Size of hash table (without chaning space): "&
         ,htable%size
    write(*,*)"Load factor: "&
         ,(htable%size - htable%nfree) * 1.D0 / (htable%size + tiny(0.D0))
    write(*,*)"Total collisions in hash table: "&
         ,htable%total_size - htable%size - htable%nfree_chain
    write(*,*)"Collision fraction: "&
         ,(htable%total_size - htable%size - htable%nfree_chain)&
         *1./(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Perfect collision fraction (assuming perfect randomness): "&
         ,(htable%total_size - htable%nfree - htable%nfree_chain - &
         htable%size * (1.d0 - ((htable%size - 1.d0)/(htable%size)) &
         **(htable%total_size - htable%nfree - htable%nfree_chain))) & 
         *1./(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Fraction of collision space used: "&
         ,(htable%total_size - htable%size - htable%nfree_chain)&
         * 1.D0 / (htable%total_size - htable%size + tiny(0.D0))
  end subroutine hash_stats

  subroutine hash_stats_simple(htable)
    implicit none
    type(hash_simple)::htable

    write(*,*)"Total values stored in hash table: "&
         ,htable%total_size - htable%nfree - htable%nfree_chain
    write(*,*)"Size of hash table (without chaning space): "&
         ,htable%size
    write(*,*)"Load factor: "&
         ,(htable%size - htable%nfree) * 1.D0 / (htable%size + tiny(0.D0))
    write(*,*)"Total collisions in hash table: "&
         ,htable%total_size - htable%size - htable%nfree_chain
    write(*,*)"Collision fraction: "&
         ,(htable%total_size - htable%size - htable%nfree_chain)&
         *1./(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Perfect collision fraction (assuming perfect randomness): "&
         ,(htable%total_size - htable%nfree - htable%nfree_chain - &
         htable%size * (1.d0 - ((htable%size - 1.d0)/(htable%size)) &
         **(htable%total_size - htable%nfree - htable%nfree_chain))) &
         *1./(htable%total_size - htable%nfree - htable%nfree_chain + tiny(0.D0))
    write(*,*)"Fraction of collision space used: "&
         ,(htable%total_size - htable%size - htable%nfree_chain)&
         * 1.D0 / (htable%total_size - htable%size + tiny(0.D0))
  end subroutine hash_stats_simple
  ! =============================================================================
end module hash

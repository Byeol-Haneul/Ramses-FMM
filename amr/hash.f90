! Hash table module for the use inside RAMSES.

! - KEY: A tuple (ilevel, ix, iy, iz) acts as hash key.

! - VALUE: Integer (typically grid indices) are stored in the hash table.

! - HASH FUNCTION: Either the murmur3 (A. Appleby,
!   https://sites.google.com/site/murmurhash/, currently only in 3d)
!   or a simpler hash function based on simple multiplication with constants.

! - COLLISIONS: A linked list is used to deal with collisions.

! - UNSAFE_HASH: Avoid comparing 4-integer keys (slow) and compare the full
!   64 (or maybe 128 bit) hash of the keys instead. EXPERIMENTAL!!!

module hash
  use amr_parameters, only: ndim
  implicit none
  
  ! General module parameters
  integer, parameter :: key_length = (ndim + 1) * 8
  integer, dimension(0:3), parameter :: constants = (/5, -1640531527, 97, 1003313/)
  
  ! Define a bucket as a derived type (sequence statement!) for better
  ! cache efficiency.
  type bucket
     sequence
#ifndef UNSAFE_HASH
     integer(kind=8), dimension(0:ndim) :: key
#else
     integer(kind=8) :: full_hash
#endif 
     integer :: value
     integer :: next_ibucket
  end type bucket     

  ! The actual hash table. Contains a procedure pointer to the desired
  ! hash function which is set at initialization of the hash table.
  type hash_table
     type(bucket), allocatable, dimension(:)  :: data
     integer         :: total_size, head_free, nfree_chain, nfree
     integer(kind=8) :: size
     integer(kind=8) :: bitmask
     integer, allocatable, dimension(:) :: next_free
     procedure(hfunc), nopass, pointer :: hash_func => null()
  end type hash_table  

  ! Interface to a general hash function - necessary for the procedure
  ! pointer used above.
  abstract interface
     pure function hfunc(key)
       use amr_parameters, only: ndim
       integer(kind=8), dimension(0:ndim), intent(in) :: key
       integer(kind=8)                                 :: hfunc
     end function hfunc
  end interface

contains

  ! Note: All these functions and subroutines could actually be contained in the
  ! hash_table type to optain a more object oriented fashion (Fortran 2003 standard).
  
  ! ============================================================================= 
  pure function simple_hash_func(key)
    integer(kind=8), dimension(0:ndim), intent(in) :: key
    integer(kind=8)                                 :: simple_hash_func
    
    simple_hash_func = dot_product(key(0:ndim), constants(0:ndim))
  end function simple_hash_func
  ! =============================================================================
  
  ! ============================================================================= 
  pure function murmur3_hash_func(key)
    integer(kind=8), dimension(0:ndim), intent(in) :: key
    integer(kind=8)                                 :: murmur3_hash_func
    integer(kind=4), parameter :: seed=42

    ! Murmur3 hash adapted for a key which is a multiple of 32 bits
    
    ! Explicit interface for the c subroutine (needed because the interface must
    ! be known at compilation time of the module)
    interface
       pure subroutine murmurhash3_x64_128(key, key_length, seed, hash_func)
         use amr_parameters, only: ndim
         integer(kind=8) , dimension(0:ndim), intent(in) :: key
         integer, intent(in)                              :: seed, key_length
         integer(kind=8), intent(inout)                   :: hash_func
       end subroutine murmurhash3_x64_128
    end interface
    
    call murmurhash3_x64_128(key, key_length, seed, murmur3_hash_func)    
    
  end function murmur3_hash_func
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
    
    if (hash_type == 'simple') then
       htable%hash_func => simple_hash_func
    else if (hash_type == 'murmur') then
       if (ndim .ne. 3)then
          print*, 'murmur3 hash currently only in 3d'
          stop
       end if
       htable%hash_func => murmur3_hash_func
    else
       htable%hash_func => simple_hash_func
    end if
    
    htable%size = 2
    do while (htable%size < req_size)
       htable%size = htable%size * 2
    end do

    call reset_entire_hash(htable, .false.)

  end subroutine init_empty_hash
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
    htable%total_size = htable%size / 4 + htable%size
    htable%nfree = htable%size
    htable%nfree_chain = htable%total_size - htable%size
    htable%head_free = htable%size + 1
    htable%bitmask = htable%size - 1

    if (.not. allocated(htable%data))then
       allocate(htable%data(1: htable%total_size))
       allocate(htable%next_free (htable%size + 1: htable%total_size))
    end if

    ! Initialize data
    do i = 1, htable%total_size
       call reset_bucket(htable%data(i))
    end do
    do i = htable%size + 1, htable%total_size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%total_size) = 0

  end subroutine reset_entire_hash
  ! =============================================================================
  
  ! =============================================================================
  subroutine reset_bucket(buck)
    implicit none
    type(bucket), intent(inout) :: buck
    
    ! Reset the content of a bucket
    buck%next_ibucket = -1
#ifndef UNSAFE_HASH
    buck%key = 0
#else
    buck%full_hash = 0
#endif
  end subroutine reset_bucket
  ! =============================================================================

  ! =============================================================================
  subroutine hash_set(htable, key, val)
    implicit none
    type(hash_table),                     intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key
    integer,                              intent(in)    :: val    
    
    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer(kind=8) :: ibucket, full_hash    

    if (val == 0)then
       write(*,*) "trying to insert 0 (0 is used to indicate absence of a value) "
       stop
    end if

    ! Compute ibucket
    full_hash = htable%hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

    if (htable%data(ibucket)%next_ibucket < 0) then          

       ! Bucket is empty, simply insert value       
       htable%data(ibucket)%next_ibucket = 0
       htable%data(ibucket)%value       = val
#ifndef UNSAFE_HASH
       htable%data(ibucket)%key(0:ndim) = key(0:ndim)
#else
       htable%data(ibucket)%full_hash   = full_hash
#endif       
       htable%nfree = htable%nfree - 1
       
    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%data(ibucket)%next_ibucket .ne. 0)
          ! Check if key already exists - abort if so
#ifndef UNSAFE_HASH
          if (same_keys(htable%data(ibucket)%key(0:ndim),key(0:ndim)))then
#else
          if (htable%data(ibucket)%full_hash == full_hash)then
#endif
             write(*,*) "trying to insert already existing key: ",key
             write(*,*) "existing key: ", htable%data(ibucket)%key(0:ndim)
             stop
          end if
          ibucket = htable%data(ibucket)%next_ibucket
       end do

       ! Check again (at the end of linked list)
#ifndef UNSAFE_HASH
       if (same_keys(htable%data(ibucket)%key(0:ndim),key(0:ndim)))then
#else
       if (htable%data(ibucket)%full_hash == full_hash)then
#endif
          write(*,*) "trying to insert already existing key: ",key
          stop
       end if
       
       ! Have reached end of chain, val not present yet -> add
       htable%data(ibucket)%next_ibucket = htable%head_free
       ibucket = htable%head_free
       htable%data(ibucket)%next_ibucket = 0
       htable%data(ibucket)%value = val
#ifndef UNSAFE_HASH
       htable%data(ibucket)%key(0:ndim) = key(0:ndim)
#else
       htable%data(ibucket)%full_hash = full_hash
#endif
       ! remove bucket from head of free linked list
       htable%head_free   = htable%next_free(htable%head_free)
       htable%nfree_chain = htable%nfree_chain - 1

    else
       write(*,*)"hash chaining space full "
       stop
    end if
  end subroutine hash_set
  ! =============================================================================

  ! =============================================================================
  pure function hash_get(htable, key)
    implicit none
    type(hash_table),                     intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in) :: key
    integer                                          :: hash_get
    
    ! Function (not subroutine, could also be changed...? ) which retrieves the 
    ! hash table value for a given key. If no entry exists, return 0
    integer(kind=8) :: ibucket, full_hash
    
    full_hash = htable%hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

#ifndef UNSAFE_HASH
    if (same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))then
#else
    if (htable%data(ibucket)%full_hash == full_hash)then
#endif
       hash_get = htable%data(ibucket)%value
       return
    end if
    
    ! Walk linked list until key is found or to the end is reached
    do while( htable%data(ibucket)%next_ibucket > 0)
       ibucket = htable%data(ibucket)%next_ibucket
#ifndef UNSAFE_HASH
       if (same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))then
#else
       if (htable%data(ibucket)%full_hash == full_hash)then
#endif
          hash_get = htable%data(ibucket)%value
          return
       end if
    end do

    ! Nothing found...
    hash_get = 0

  end function hash_get
  ! =============================================================================

  ! =============================================================================  
  subroutine hash_free(htable, key)
    implicit none
    type(hash_table),                     intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key

    ! Remove the hash table entry for a given key 

    integer(kind=8) :: ibucket, previous_ibucket, full_hash

    full_hash = htable%hash_func(key)
    ibucket = IAND(full_hash, htable%bitmask) + 1

    ! No collision case
    if (htable%data(ibucket)%next_ibucket == 0) then     
       htable%data(ibucket)%next_ibucket = -1
#ifndef UNSAFE_HASH
       htable%data(ibucket)%key(0:ndim) = 0
#else
       htable%data(ibucket)%full_hash = 0
#endif       
       htable%nfree = htable%nfree + 1
    else
       ! Collision case
#ifndef UNSAFE_HASH
       do while (.not. same_keys(htable%data(ibucket)%key(0:ndim), key(0:ndim)))
#else
       do while (htable%data(ibucket)%full_hash .ne. full_hash)
#endif
          previous_ibucket=ibucket
          ibucket=htable%data(ibucket)%next_ibucket
       end do
       if (ibucket <= htable%size) then           
          ! It's the first element we need to erase: Move first element from chaning 
          ! space into bucket and do as if the value to remove had been in the chaning space
          htable%data(ibucket)%value = htable%data(htable%data(ibucket)%next_ibucket)%value
#ifndef UNSAFE_HASH
          htable%data(ibucket)%key = htable%data(htable%data(ibucket)%next_ibucket)%key
#else
          htable%data(ibucket)%full_hash = htable%data(htable%data(ibucket)%next_ibucket)%full_hash
#endif
          previous_ibucket = ibucket
          ibucket = htable%data(ibucket)%next_ibucket
       end if
       ! fill the hole and reconnect linked list
       htable%data(previous_ibucket)%next_ibucket = htable%data(ibucket)%next_ibucket
       htable%next_free(ibucket) = htable%head_free
       htable%head_free = ibucket
       htable%nfree_chain = htable%nfree_chain + 1
    end if
  end subroutine hash_free
  ! =============================================================================

  ! =============================================================================
  pure function same_keys(key1, key2)
    logical :: same_keys
    integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2       
    
    ! Function to test the equality of two provided keys
    ! using the c standar library function memcmp
    interface
       pure function memcmp(key1, key2, key_length)
         use amr_parameters, only: ndim
         integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2     
         integer, intent(in) :: key_length
         integer :: memcmp
       end function memcmp
    end interface
    same_keys =  memcmp(key1, key2, key_length) == 0_4
  end function same_keys

  ! ALTERNATIVE VERSION - CAN BE USED INSTEAD OF THE C CODE.
  ! function same_keys(key1, key2)
  !   logical :: same_keys
  !   integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2
  !   logical, dimension(0:ndim), save :: ok
  !   do i = 0, ndmin
  !      ok(i) = (key1(i)==key2(i))
  !   end do
  !   same_keys = ALL(ok)
  ! end function same_keys
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
  ! =============================================================================
end module hash

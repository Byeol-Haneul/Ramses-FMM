! Hash table for the use inside ramses. Prime murmur3 hash, double-linked-list 
! for chaining, assumes a NDIM+1 -integer hilbert-key as key.
! TODO: test this 

module hash
  use amr_parameters, only: ndim, nlevelmax
  type hash_table
     integer        , allocatable, dimension(:)   :: value
     integer        , allocatable, dimension(:)   :: next_bucket
     integer        , allocatable, dimension(:)   :: next_free
     integer(kind=8), allocatable, dimension(:,:) :: key
     integer         :: size, head_free, nfree_chain, nfree
     integer         :: c1, c2, c3
     integer(kind=8) :: prime     
  end type hash_table
contains

  ! ============================================================================= 
  function hash_func(htable, key)
    type(hash_table),                    intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in) :: key
    integer(kind=4)                                 :: hash_func
    integer(kind=4), dimension(1:2),save            :: hash
    integer(kind=4), parameter :: seed=42, len = 32
    integer(kind=4), save :: tablesize

    ! compute the "bucket" as a function of the nkey-integer key.
#if NDIM==1
    hash_func = key(0)+htable%c1*key(1)
#endif    
#if NDIM==2
    hash_func = key(0)+htable%c1*key(1)+htable%c2*key(2)
#endif    
#if NDIM==3
    hash_func = key(0)+htable%c1*key(1)+htable%c2*key(2)+htable%c3*key(3) 
#endif    
    hash_func = MODULO(hash_func,htable%prime) + 1

!    tablesize = htable%prime-1
    !call murmurhash3_x64_128(key, len, tablesize, seed, hash)
    ! TODO: maybe remove this by allocating the buckets starting from )...
    !hash_func = hash(1) + 1

  end function hash_func
  ! =============================================================================

  ! =============================================================================
  subroutine init_empty_hash(htable, req_size)
    use amr_parameters, only:ndim
    implicit none
    type(hash_table), intent(inout) :: htable
    integer         , intent(in)    :: req_size

    ! Allocate all hash table arrays and variables, choose appropriate prime
    ! based on the required size of the hash table.

    integer                  :: ncode, bit_length, i
    integer, dimension(0:30) :: prime=(/2,3,7,13,23,53,97,193,389,769,1543,&
         & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
         & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
         & 402653189,805306457,1610612741/)

    ! TODO: rename prime since it's not a prime anymore...
    ! Compute prime number
    ncode=req_size
    do bit_length=1,32
       ncode=ncode/2
       if(ncode<=1) exit
    end do

    !htable%prime = 2
    !do while (htable%prime < req_size)
    !   htable%prime = htable%prime * 2
    !end do

    ! Allocate and initialize arrays
    htable%prime = prime(bit_length)
    htable%size = htable%prime / 4 + htable%prime
    htable%nfree = htable%prime
    allocate(htable%value      (1:htable%size))
    allocate(htable%key (0:ndim,1:htable%size))
    htable%key = 0
    allocate(htable%next_bucket(1:htable%size))
    htable%next_bucket = -1

    ! Build linked list of free slots in the chaning part of the array
    allocate(htable%next_free (htable%prime + 1 : htable%size))
    do i = htable%prime + 1, htable%size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%size) = 0
    htable%head_free = htable%prime + 1
    htable%nfree_chain = htable%size - htable%prime

    ! build constants
    htable%c1 = nlevelmax
    htable%c2 = nlevelmax*2**nlevelmax
    htable%c3 = nlevelmax*4**nlevelmax

  end subroutine init_empty_hash
  ! =============================================================================

  ! =============================================================================
  subroutine reset_entire_hash(htable)
    implicit none
    type(hash_table), intent(inout) :: htable

    ! Subroutine to reset the entire hash table

    integer :: i
    
    ! Reinitialize arrays
    htable%nfree = htable%prime
    htable%next_bucket = -1
    htable%key = 0
    do i = htable%prime + 1, htable%size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%size) = 0
    htable%head_free = htable%prime + 1
    htable%nfree_chain = htable%size - htable%prime
  end subroutine reset_entire_hash
  ! =============================================================================

  ! =============================================================================
  subroutine hash_set(htable, key, val)
    implicit none
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key
    integer,                             intent(in)    :: val    
    
    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer :: bucket

    if (val == 0)then
       write(*,*) "trying to insert 0 (0 is used to indicate absence of a value) "
       stop
    end if
    
    ! Compute bucket
    bucket = hash_func(htable, key)
    
    if (htable%next_bucket(bucket) < 0) then          

       ! Bucket is empty, simply insert value       
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key (0:ndim,bucket) = key(0:ndim)
       htable%nfree = htable%nfree - 1

    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%next_bucket(bucket) .ne. 0)

          ! Check if key already exists
          if (same_keys(htable%key(0:ndim,bucket),key(0:ndim)))then
             write(*,*) "trying to insert already existing key: ",key
             stop
          end if
          bucket = htable%next_bucket(bucket)
       end do

       ! Check if key is already there
       if (same_keys(htable%key(0:ndim,bucket),key(0:ndim)))then
          write(*,*) "trying to insert already existing key: ",key
          stop
       end if
       
       ! Have reached end of chain, val not present yet -> add
       htable%next_bucket(bucket) = htable%head_free
       bucket = htable%head_free
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key  (0:ndim,bucket) = key(0:ndim)

       ! remove bucket from head of free linked list
       htable%head_free   = htable%next_free(htable%head_free)
       htable%nfree_chain = htable%nfree_chain - 1

    else
       write(*,*)"hash chaining space full on process "
       stop
    end if
  end subroutine hash_set
  ! =============================================================================

  ! =============================================================================
  function hash_get(htable, key)
    implicit none
    type(hash_table),                    intent(in) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in) :: key
    integer                                         :: hash_get

    ! Function (not subroutine, could also be changed...? ) which retrieves the 
    ! hash table value for a given key. If no entry exists, return 0
    integer :: bucket

    bucket = hash_func(htable, key)

    if (same_keys(htable%key(0:ndim,bucket), key(0:ndim)))then
       hash_get = htable%value(bucket)
       return
    end if
    
    ! Walk linked list until key is found or to the end is reached
    do while( htable%next_bucket(bucket) > 0)
       bucket = htable%next_bucket(bucket)
       if (same_keys(htable%key(0:ndim,bucket), key(0:ndim)))then
          hash_get = htable%value(bucket)
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
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:ndim), intent(in)    :: key
    
    ! Remove the hash table entry for a given key 

    integer :: bucket, previous_bucket
    
    bucket = hash_func(htable, key)

    if (htable%next_bucket(bucket) == 0) then     ! No collision case
       htable%next_bucket(bucket)  = -1
       htable%nfree = htable%nfree + 1
       htable%key(0:ndim, bucket) = 0
    else                                          ! Collision case
       do while (.not. same_keys(htable%key(0:ndim,bucket), key(0:ndim)))
          previous_bucket=bucket
          bucket=htable%next_bucket(bucket)
       end do
       if (bucket <= htable%prime) then           
          ! It's the first element we need to erase: Move first element from chaning 
          ! space into bucket and do as if the value to remove had been in the chaning space
          htable%value(       bucket) = htable%value(       htable%next_bucket(bucket))
          htable%key  (0:ndim,bucket) = htable%key  (0:ndim,htable%next_bucket(bucket))
          previous_bucket = bucket
          bucket = htable%next_bucket(bucket)
       end if
       ! fill the hole and reconnect linked list
       htable%next_bucket(previous_bucket) = htable%next_bucket(bucket)
       htable%next_free(bucket) = htable%head_free
       htable%head_free = bucket
       htable%nfree_chain = htable%nfree_chain+1
    end if
  end subroutine hash_free
  ! =============================================================================

  ! =============================================================================
  subroutine hash_stats(htable)
    implicit none
    type(hash_table)::htable

    write(*,*)"Total values stored in hash table: "&
         ,htable%size-htable%nfree-htable%nfree_chain
    write(*,*)"Total collisions in hash table: "&
         ,htable%size-htable%prime-htable%nfree_chain
    write(*,*)"Collision fraction: "&
         ,(htable%size-htable%prime-htable%nfree_chain)&
         *1./(htable%size-htable%nfree-htable%nfree_chain+tiny(0.D0))
    write(*,*)"Fraction of collision space used: "&
         ,(htable%size-htable%prime-htable%nfree_chain)&
         *1./ (htable%size-htable%prime+tiny(0.D0))
    write(*,*)"Fraction of proper space used: "&
         ,(htable%prime-htable%nfree)*1./(htable%prime+tiny(0.D0))
  end subroutine hash_stats
  ! =============================================================================

  ! function same_keys(key1, key2)
  !   logical :: same_keys
  !   integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2     
  !   same_keys =  ( IOR(IEOR(key1(3), key2(3)), &
  !        IOR(IEOR(key1(2), key2(2)), &
  !        IOR(IEOR(key1(1), key2(1)), &
  !        IEOR(key1(0), key2(0)))))) == 0_8
  ! end function same_keys
!!$  function same_keys(key1, key2)
!!$    logical :: same_keys
!!$    integer, parameter :: thirtytwo=32
!!$    integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2     
!!$    same_keys =  memcmp(key1, key2, thirtytwo) == 0_4
!!$  end function same_keys
  function same_keys(key1, key2)
    logical :: same_keys
    integer(kind=8), dimension(0:ndim), intent(in) :: key1, key2     
#if NDIM==1
    same_keys =  (key1(0)==key2(0) .and. key1(1)==key2(1))
#endif
#if NDIM==2
    same_keys =  (key1(0)==key2(0) .and. key1(1)==key2(1) .and. key1(2)==key2(2))
#endif
#if NDIM==3
    same_keys =  (key1(0)==key2(0) .and. key1(1)==key2(1) .and. key1(2)==key2(2) .and. key1(3)==key2(3))
#endif
    return
  end function same_keys
end module hash

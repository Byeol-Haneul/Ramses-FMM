! Basic hash table for the use inside ramses. Prime number hash, double-linked-list 
! for chaining, assumes a NDIM-integer hilber-key  as key!!!
! The hash statistics optained in some simple tests suggest that the hash function
! might have to be changed to something more elaborate at some point...

module hash
  use amr_parameters, only: ndim, nlevelmax
  type hash_table
     integer        , allocatable, dimension(:)   :: value
     integer        , allocatable, dimension(:)   :: next_bucket
     integer        , allocatable, dimension(:)   :: next_free
     integer(kind=8), allocatable, dimension(:,:) :: key
     integer(kind=4), allocatable, dimension(:)   :: key_level
     integer         :: size, head_free, nfree_chain, nfree
     integer         :: c1, c2, c3
     integer(kind=8) :: prime     
  end type hash_table
  integer, parameter :: nkey = ndim - 1                   
contains

  ! ============================================================================= 
  function hash_func(htable, key, key_level)
    type(hash_table),                    intent(in) :: htable
    integer(kind=8) , dimension(0:nkey), intent(in) :: key
    integer,                             intent(in) :: key_level
    integer                                         :: hash_func
    ! compute the "bucket" as a function of the nkey-integer key.
    hash_func = MOD(MOD(key(0),htable%prime) + htable%c1 * MOD(key(1),htable%prime)&
         + htable%c2 * MOD(key(2),htable%prime) + htable%c3 * key_level, htable%prime) + 1
    return
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
    
    ! Compute prime number
    ncode=req_size
    do bit_length=1,32
       ncode=ncode/2
       if(ncode<=1) exit
    end do

    ! Allocate and initialize arrays
    htable%prime = prime(bit_length + 1)
    htable%size = htable%prime / 4 + htable%prime
    htable%nfree = htable%prime
    allocate(htable%value      (1:htable%size))
    allocate(htable%key (0:nkey,1:htable%size))
    htable%key = 0
    allocate(htable%key_level (1:htable%size))
    htable%key_level = 0
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
    htable%c1 = 1
    do i = 1, nlevelmax
       htable%c1 = mod(2 * htable%c1, htable%prime)
    end do

    htable%c2 = htable%c1
    do i = 1, nlevelmax
       htable%c2 = mod(2 * htable%c2, htable%prime)
    end do

    htable%c3 = htable%c2
    do i = 1, nlevelmax
       htable%c3 = mod(2 * htable%c3, htable%prime)
    end do
    

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
    htable%key_level = 0    
    do i = htable%prime + 1, htable%size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%size) = 0
    htable%head_free = htable%prime + 1
    htable%nfree_chain = htable%size - htable%prime
  end subroutine reset_entire_hash
  ! =============================================================================

  ! =============================================================================
  subroutine hash_set(htable, key, key_level, val)
    implicit none
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:nkey), intent(in)    :: key
    integer,                             intent(in)    :: key_level
    integer,                             intent(in)    :: val    
    
    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer :: bucket

    if (val == 0)then
       write(*,*) "trying to insert 0 (0 is used to indicate absence of a value) "
       stop
    end if
    
    ! Compute bucket
    bucket = hash_func(htable, key, key_level)
    
    if (htable%next_bucket(bucket) < 0) then          

       ! Bucket is empty, simply insert value       
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key (0:nkey,bucket) = key(0:nkey)
       htable%key_level(bucket)   = key_level
       htable%nfree = htable%nfree - 1

    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%next_bucket(bucket) .ne. 0)

          ! a bit ugly: check if htable%key == key
          if (htable%key(0,bucket) == key(0) .and. htable%key_level(bucket) == key_level &
#if NDIM>1
               .and. htable%key(1,bucket) == key(1) &
#if NDIM>2
               .and. htable%key(2,bucket) == key(2) &
#endif
#endif
               ) then
             write(*,*) "trying to insert already existing key: ",key, key_level
!             print*, 1/mod(2,2)
             stop
          end if
          bucket = htable%next_bucket(bucket)
       end do

       ! a bit ugly: check if htable%key == key
       if (htable%key(0,bucket) == key(0) .and. htable%key_level(bucket) == key_level &
#if NDIM>1
            .and. htable%key(1,bucket) == key(1) &
#if NDIM>2
            .and. htable%key(2,bucket) == key(2) &
#endif
#endif
            ) then
          !          print*, 1/mod(2,2)
          write(*,*) "trying to insert already existing key: ",key, key_level
          stop
       end if

       ! Have reached end of chain, val not present yet -> add
       htable%next_bucket(bucket) = htable%head_free
       bucket = htable%head_free
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key  (0:nkey,bucket) = key(0:nkey)
       htable%key_level   (bucket) = key_level

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
  function hash_get(htable, key, key_level)
    implicit none
    type(hash_table),                    intent(in) :: htable
    integer(kind=8) , dimension(0:nkey), intent(in) :: key
    integer,                             intent(in) :: key_level
    integer                                         :: hash_get

    ! Function (not subroutine, could also be changed...? ) which retrieves the 
    ! hash table value for a given key. If no entry exists, return 0
    integer :: bucket

    bucket = hash_func(htable, key, key_level)

    ! Walk linked list until key is found or to the end

    ! a bit ugly: do while htable%key .ne. key .and. next_bucket > 0
    do while(((htable%key(0,bucket) .ne. key(0)) .or. (htable%key_level(bucket) .ne. key_level) &
#if NDIM>1
         .or. (htable%key(1,bucket) .ne. key(1)) &
#if NDIM>2
         .or. (htable%key(2,bucket) .ne. key(2)) &
#endif
#endif
         ) .and. htable%next_bucket(bucket) > 0)       
       bucket = htable%next_bucket(bucket)
    end do
    ! a bit ugly: check if htable%key == key
    if (htable%key(0,bucket) == key(0) .and. htable%key_level(bucket) == key_level &
#if NDIM>1
         .and. htable%key(1,bucket) == key(1) &
#if NDIM>2
         .and. htable%key(2,bucket) == key(2) &
#endif
#endif
         ) then
       hash_get = htable%value(bucket)
       return
    else
       hash_get = 0
       return 
    end if
  end function hash_get
  ! =============================================================================

  ! =============================================================================  
  subroutine hash_free(htable, key, key_level)
    implicit none
    type(hash_table),                    intent(inout) :: htable
    integer(kind=8) , dimension(0:nkey), intent(in)    :: key
    integer,                             intent(in)    :: key_level
    
    ! Remove the hash table entry for a given key 

    integer :: bucket, previous_bucket
    
    bucket = hash_func(htable, key, key_level)

    if (htable%next_bucket(bucket) == 0) then     ! No collision case
       htable%next_bucket(bucket)  = -1
       htable%nfree = htable%nfree + 1
    else                                          ! Collision case
       do while ((htable%key(0,bucket) .ne. key(0)) .or. (htable%key_level(bucket) .ne. key_level) &            
#if NDIM>1
            .or. htable%key(1,bucket) .ne. key(1) &
#if NDIM>2
            .or. htable%key(2,bucket) .ne. key(2) &
#endif
#endif            
            )
          previous_bucket=bucket
          bucket=htable%next_bucket(bucket)
       end do
       if (bucket <= htable%prime) then           
          ! It's the first element we need to erase: Move first element from chaning 
          ! space into bucket and do as if the value to remove had been in the chaning space
          htable%value(       bucket) = htable%value(       htable%next_bucket(bucket))
          htable%key  (0:nkey,bucket) = htable%key  (0:nkey,htable%next_bucket(bucket))
          htable%key_level   (bucket) = htable%key_level   (htable%next_bucket(bucket))
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

end module hash

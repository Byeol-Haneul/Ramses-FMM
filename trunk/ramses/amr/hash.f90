! Basic hash table for the use inside ramses. Prime number hash, double-linked-list 
! for chaining, assumes integers as keys.
! The hash statistics optained in some simple tests suggest that the hash function
! might have to be changed to something more elaborate at some point...

module hash
  type hash_table
     integer        , allocatable, dimension(:) :: value
     integer        , allocatable, dimension(:) :: next_bucket
     integer        , allocatable, dimension(:) :: next_free
     integer(kind=8), allocatable, dimension(:) :: key     
     integer         :: size, head_free, nfree_chain, nfree
     integer(kind=8) :: prime     
  end type hash_table
contains

  ! =============================================================================
  subroutine init_empty_hash(htable, req_size)
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
    htable%prime         = prime(bit_length + 1)
    htable%size          = htable%prime / 4 + htable%prime
    htable%nfree         = htable%prime
    allocate(htable%value      (1:htable%size))
    allocate(htable%key        (1:htable%size))
    allocate(htable%next_bucket(1:htable%size))
    htable%next_bucket   = -1

    ! Build linked list of free slots in the chaning part of the array
    allocate(htable%next_free (htable%prime + 1 : htable%size))
    do i = htable%prime + 1, htable%size - 1
       htable%next_free(i) = i + 1
    end do
    htable%next_free(htable%size) = 0
    htable%head_free = htable%prime + 1
    htable%nfree_chain = htable%size - htable%prime
  end subroutine init_empty_hash
  ! =============================================================================

  ! =============================================================================
  subroutine hash_set(htable, key, val)
    use amr_commons, only: myid
    implicit none
    type(hash_table), intent(inout) :: htable
    integer(kind=8) , intent(in)    :: key
    integer         , intent(in)    :: val    
    
    ! Add a key/value pair to the hash table. If there is already a key/value
    ! pair stored for this key, return an error message.

    integer :: bucket
    
    ! Compute bucket
    bucket = MOD(key,htable%prime) + 1
    
    if (htable%next_bucket(bucket) < 0) then          

       ! Bucket is empty, simply insert value       
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key        (bucket) = key
       htable%nfree = htable%nfree - 1

    else if (htable%nfree_chain>0)then

       ! Bucket is not empty, walk through linked list
       do while (htable%next_bucket(bucket) .ne. 0)
          if (htable%key(bucket) == key) then
             write(*,*) "trying to insert already existing key: ",key
             stop
          end if
          bucket = htable%next_bucket(bucket)
       end do
       if (htable%key(bucket) == key) then
          write(*,*) "trying to insert already existing key: ",key
          stop
       end if

       ! Have reached end of chain, val not present yet -> add
       htable%next_bucket(bucket) = htable%head_free
       bucket = htable%head_free
       htable%next_bucket(bucket) = 0
       htable%value      (bucket) = val
       htable%key        (bucket) = key

       ! remove bucket from head of free linked list
       htable%head_free = htable%next_free(htable%head_free)
       htable%nfree_chain = htable%nfree_chain - 1

    else
       write(*,*)"hash chaining space full on process ", myid
       stop
    end if
  end subroutine hash_set
  ! =============================================================================

  ! =============================================================================
  function hash_get(htable, key)
    implicit none
    type(hash_table), intent(in) :: htable
    integer(kind=8) , intent(in) :: key
    integer                      :: hash_get

    ! Function (not subroutine, could also be changed...? ) which retrieves the 
    ! hash table value for a given key. If no entry exists, print an error and abort.
    ! This could be changed to return 0 if no entry exists. This return value could
    ! trigger insertion a new entry for this key.

    integer :: bucket

    bucket = MOD(key, htable%prime) + 1

    ! Walk linked list until key is found or to the end
    do while ((htable%key(bucket) .ne. key) .and. htable%next_bucket(bucket) > 0)       
       bucket = htable%next_bucket(bucket)
    end do
    if (htable%key(bucket) == key)then
       hash_get = htable%value(bucket)
       return
    else
       write(*,*)"no entry found for key: ", key
       stop
    end if
  end function hash_get
  ! =============================================================================

  ! =============================================================================  
  subroutine hash_free(htable, key)
    implicit none
    type(hash_table), intent(inout) :: htable
    integer(kind=8) , intent(in)    :: key
    
    ! Remove the hash table entry for a given key 

    integer :: bucket, previous_bucket
    
    bucket = MOD(key,htable%prime) + 1

    if (htable%next_bucket(bucket) == 0) then     ! No collision case
       htable%next_bucket(bucket)  = -1
       htable%nfree = htable%nfree + 1
    else                                          ! Collision case
       do while (htable%key(bucket) .ne. key)
          previous_bucket=bucket
          bucket=htable%next_bucket(bucket)
       end do
       if (bucket <= htable%prime) then           
          ! It's the first element we need to erase: Move first element from chaning 
          ! space into bucket and do as if the value to remove had been in the chaning space
          htable%value(bucket) = htable%value(htable%next_bucket(bucket))
          htable%key  (bucket) = htable%key  (htable%next_bucket(bucket))
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

! basic hash table for the use inside ramses. Prime number hash, double-linked-list for chaining.
! assumes integers as keys.

module hash

  type hash_table
     integer,allocatable,dimension(:)::value
     integer,allocatable,dimension(:)::next_bucket
     integer,allocatable,dimension(:)::next_free
     integer(kind=8),allocatable,dimension(:)::key
     
     integer(kind=8)::prime
     integer::size,head_free,nfree_chain,nfree
     
  end type hash_table
  
contains
  
  subroutine init_empty_hash(htable,req_size)
    implicit none
    type(hash_table)::htable
    integer::req_size
    integer,dimension(0:30)::prime=(/2,3,7,13,23,53,97,193,389,769,1543,&
         & 3079,6151,12289,24593,49157,98317,196613,393241,786433,1572869, &
         & 3145739,6291469,12582917,25165843,50331653,100663319,201326611, &
         & 402653189,805306457,1610612741/)
    integer::ncode,bit_length,i

    ! Compute prime of appropriate size
    ncode=req_size
    do bit_length=1,32
       ncode=ncode/2
       if(ncode<=1) exit
    end do

    ! Allocate and initialize arrays
    htable%prime=prime(bit_length+1)
    htable%size=htable%prime/4+htable%prime
    allocate(htable%value(1:htable%size))
    htable%nfree=htable%prime
    allocate(htable%key(1:htable%size))
    allocate(htable%next_bucket(1:htable%size))
    htable%next_bucket=-1
    allocate(htable%next_free(htable%prime+1:htable%size))
    do i=htable%prime+1,htable%size-1
       htable%next_free(i)=i+1
    end do
    htable%next_free(htable%size)=0
    htable%head_free=htable%prime+1
    htable%nfree_chain=htable%size-htable%prime

  end subroutine init_empty_hash


  subroutine hash_set(htable,key,val)
    use amr_commons, only:myid
    implicit none
    type(hash_table)::htable
    integer(kind=8)::key
    integer::val
    
    integer::bucket
    
    bucket=MOD(key,htable%prime)+1
    
    if (htable%next_bucket(bucket) < 0)then
       htable%next_bucket(bucket)=0
       htable%value(bucket)=val
       htable%key(bucket)=key
       htable%nfree=htable%nfree-1
    else if (htable%nfree_chain>0)then
       do while(htable%next_bucket(bucket) .ne. 0)
          if (htable%key(bucket)==key)then
             print*, 'trying to insert already existing key'
             stop
          end if
          bucket=htable%next_bucket(bucket)
       end do
       if (htable%key(bucket)==key)then
          print*, 'trying to insert already existing key'
          stop
       end if       
       htable%next_bucket(bucket)=htable%head_free
       bucket=htable%head_free
       htable%next_bucket(bucket)=0
       htable%value(bucket)=val
       htable%key(bucket)=key
       ! reset head of free
       htable%head_free=htable%next_free(htable%head_free)
       htable%nfree_chain=htable%nfree_chain-1
    else
       print*,'hash chaining space full', myid
       stop
    end if

  end subroutine hash_set


  function hash_get(htable,key)
    implicit none
    type(hash_table)::htable
    integer(kind=8)::key
    integer::hash_get
    
    integer::bucket

    bucket=MOD(key,htable%prime)+1
!    print*,key,htable%prime,MOD(key,htable%prime),bucket
    do while ((htable%key(bucket) .ne. key) .and. htable%next_bucket(bucket)>0)       
       bucket=htable%next_bucket(bucket)
    end do
    if (htable%key(bucket) == key)then
       hash_get=htable%value(bucket)
       return
    else
       print*,'no entry found for key: ',key
       stop
    end if
  end function hash_get

  
  subroutine hash_free(htable,key)
    implicit none
    type(hash_table)::htable
    integer(kind=8)::key
        
    integer::bucket,previous_bucket
    
    bucket=MOD(key,htable%prime)+1
    if (htable%next_bucket(bucket)==0)then     ! No collision case
       htable%next_bucket(bucket)=-1
       htable%nfree=htable%nfree+1
    else                                      ! Collision case
       do while (htable%key(bucket) .ne. key)
          previous_bucket=bucket
          bucket=htable%next_bucket(bucket)
       end do
       if (bucket <= htable%prime)then        ! It's the first element you need to erase
          htable%value(bucket)=htable%value(htable%next_bucket(bucket))
          htable%key(bucket)=htable%key(htable%next_bucket(bucket))
!          htable%next_bucket(bucket)=htable%next_bucket(nbucket)
          previous_bucket=bucket
          bucket=htable%next_bucket(bucket)
       end if
       ! fill the hole
       htable%next_bucket(previous_bucket)=htable%next_bucket(bucket)
       htable%next_free(bucket)=htable%head_free
       htable%head_free=bucket
       htable%nfree_chain=htable%nfree_chain+1
    end if

  end subroutine hash_free

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
         ,(htable%size-htable%prime-htable%nfree_chain)*1./(htable%size-htable%prime+tiny(0.D0))

    write(*,*)"Fraction of proper space used: "&
         ,(htable%prime-htable%nfree)*1./(htable%prime+tiny(0.D0))
  end subroutine hash_stats

end module hash

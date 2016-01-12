subroutine init_amr
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  integer(kind=4)::ilevel,icpu
  integer(kind=8)::maxkey
  logical::ok

  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0

  ! Allocate main oct array
  allocate(grid(1:ngridmax))

  ! Allocate oct hash table
  if(verbose.and.myid==1)write(*,*)'Initialize empty hash'
  call init_empty_hash(grid_dict,2*ngridmax)

  ! Set initial cpu boundaries
  ! Set maximum Cartesian key per level
  if(verbose.and.myid==1)write(*,*)'Initialize level cpu boundaries'
  allocate(bound_key_level(0:ncpu,levelmin:nlevelmax))
  allocate(ckey_max(levelmin:nlevelmax))
  do ilevel=levelmin,nlevelmax
     ckey_max(ilevel)=2**(ilevel-1)
     maxkey=2**((ilevel-1)*ndim)
     do icpu=1,ncpu-1
        bound_key_level(icpu,ilevel) = (icpu*maxkey)/ncpu
     end do
     bound_key_level(0,ilevel) = 0
     bound_key_level(ncpu,ilevel) = maxkey
  end do

  ! Allocate head, tail and numbers for each level
  if(verbose.and.myid==1)write(*,*)'Initialize oct decomposition'
  allocate(head(levelmin:nlevelmax))
  allocate(tail(levelmin:nlevelmax))
  allocate(noct(levelmin:nlevelmax))
  head=0       ! Head oct in the level
  tail=0       ! Tail oct in the level
  noct=0       ! Number of oct in the level
  noct_used=0  ! Total number of oct used

end subroutine init_amr


subroutine hash_tests(all_ok)
  use hash
  implicit none

  type(hash_table)::htable
  integer::i,nfree_store,nfree_chain_store,ipos
  logical::ok,all_ok
  real,dimension(1:3000)::val_float
  real,dimension(0:ndim,1:3000)::key_float
  integer,dimension(1:3000)::val
  integer(kind=8),dimension(0:ndim,1:3000)::key

  ok=.true.

  write(*,*)'Entering hash_tests'

  call random_number(key_float)
  call random_number(val_float)

  do i=1,3000
     key(0,i)=int(key_float(0,i)*nlevelmax,kind=8)
     key(1:ndim,i)=int(key_float(1:ndim,i)*2.0**nlevelmax,kind=8)
     val(i)=int(val_float(i)*2000,kind=4) + 1
  end do

  call init_empty_hash(htable,3000)
  write(*,*)htable%prime

  call hash_stats(htable)

  nfree_store=htable%nfree
  nfree_chain_store=htable%nfree_chain

  do i=1,2000
     call hash_set(htable,key(0:ndim,i),val(i))
  end do

  call hash_stats(htable)

  do i=2000,1,-1
     ipos=hash_get(htable,key(0:ndim,i))
     ok=ok .and. (val(i)==ipos)
  end do

  do i=1001,2000
     call hash_free(htable,key(0:ndim ,i) )
  end do

  call hash_stats(htable)

  do i=2001,3000
     call hash_set(htable,key(0:ndim ,i) ,val(i))
  end do

  call hash_stats(htable)

  do i=1,1000
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
  end do

  do i=2001,3000
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
     call hash_free(htable,key(0:ndim ,i) )
  end do

  call hash_stats(htable)

  do i=1000,1,-1
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
     call hash_free(htable,key(0:ndim ,i) )
  end do

  call hash_stats(htable)

  ok=ok .and. (nfree_store==htable%nfree)
  ok=ok .and. (nfree_chain_store==htable%nfree_chain)

  print*,ok

  if (.not. ok)then
     write(*,*)'hash test FAILED '
     all_ok=.false.
  end if

end subroutine hash_tests



subroutine amr_tests(all_ok)
  implicit none
  logical::all_ok
  logical::amr_ok=.true.

  call hilbert_tests(amr_ok)
  call hash_tests(amr_ok)
!  call other_test2(amr_ok)

  if(.not. amr_ok)then
     write(*,*)'AMR_TESTS FAILED'
     all_ok=.false.
  end if
end subroutine amr_tests


! =====================================================================================
! =====================================================================================
! ADD UNIT TESTS HERE
! =====================================================================================
! =====================================================================================

subroutine hilbert_tests(all_ok)
  use amr_commons
  use amr_parameters
  use hilbert,       only: hilbert3d_reverse, hilbert3d
  implicit none

  ! A simple test which transforms integer coordinates into a 3 integer hilbert key
  ! and the hilbert key back into integer coordinate. Results must be equal to input.

  integer::bit_length,i
  logical::ok,all_ok
  real,dimension(1:nvector)::xfloat,yfloat,zfloat
  integer(kind=8),dimension(1:nvector)::xint,yint,zint
  integer(kind=8),dimension(1:nvector)::hkey2,hkey1,hkey0
  integer(kind=4),dimension(1:nvector)::cstate
  integer(kind=8),dimension(1:nvector)::xint_store,yint_store,zint_store

  do bit_length=1,63
     ok=.true.

     call random_number(xfloat)
     call random_number(yfloat)
     call random_number(zfloat)
     
     do i=1,nvector
        xint(i)=int(xfloat(i)*2.0**bit_length,kind=8)
        yint(i)=int(yfloat(i)*2.0**bit_length,kind=8)
        zint(i)=int(zfloat(i)*2.0**bit_length,kind=8)
     end do

     xint_store(1:nvector)=xint(1:nvector)
     yint_store(1:nvector)=yint(1:nvector)
     zint_store(1:nvector)=zint(1:nvector)

     call hilbert3d(xint,yint,zint,hkey2,hkey1,hkey0,cstate,0,bit_length,nvector)
     call hilbert3d_reverse(xint,yint,zint,hkey2,hkey1,hkey0,bit_length,nvector)

     do i=1,nvector
        if( xint_store(i) .ne. xint(i) )ok=.false.
        if( yint_store(i) .ne. yint(i) )ok=.false.
        if( zint_store(i) .ne. zint(i) )ok=.false.
     end do

     if (.not. ok)then
        write(*,*)'hilbert test FAILED for bit_length ',bit_length
        all_ok=.false.
     end if

  end do

end subroutine hilbert_tests
! =====================================================================================
! =====================================================================================
subroutine hash_tests(all_ok)
  use hash
  implicit none


  type(hash_table)::htable
  integer::i,nfree_store,nfree_chain_store
  logical::ok,all_ok
  real,dimension(1:3000)::val_float,key_float
  integer,dimension(1:3000)::val
  integer(kind=8),dimension(1:3000)::key

  ok=.true.
 
  call random_number(key_float)
  call random_number(val_float)
  
  do i=1,3000
     key(i)=int(key_float(i)*2.0**62,kind=8)
     val(i)=int(val_float(i)*2000,kind=4)
  end do

  call init_empty_hash(htable,3000)

!  call hash_stats(htable)

  nfree_store=htable%nfree
  nfree_chain_store=htable%nfree_chain
  
  do i=1,2000
     call hash_set(htable,key(i),val(i))     
  end do

!  call hash_stats(htable)


  do i=2000,1,-1
     ok=ok .and. (val(i)==hash_get(htable,key(i)))
  end do


  do i=1001,2000
     call hash_free(htable,key(i))
  end do

!  call hash_stats(htable)

  do i=2001,3000
     call hash_set(htable,key(i),val(i))
  end do
  
  call hash_stats(htable)

  do i=1,1000
     ok=ok .and. (val(i)==hash_get(htable,key(i)))
  end do


  do i=2001,3000
     ok=ok .and. (val(i)==hash_get(htable,key(i)))
     call hash_free(htable,key(i))
  end do

!  call hash_stats(htable)

  do i=1000,1,-1
     ok=ok .and. (val(i)==hash_get(htable,key(i)))
     call hash_free(htable,key(i))
  end do

!  call hash_stats(htable)
  
  ok=ok .and. (nfree_store==htable%nfree)
  ok=ok .and. (nfree_chain_store==htable%nfree_chain)

  
  if (.not. ok)then
     write(*,*)'hash test FAILED '
     all_ok=.false.
  end if
  
end subroutine hash_tests

module cache
contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache(s,hash_dict)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  type(hash_table)::hash_dict
#ifndef MDL2
  !
  ! This routine closes all cache operations.
  ! It purges all remaining flush messages.
  !
  integer::info,icache,igrid,icpu,iskip
  integer::send_flush_id,nflush
  integer::dummy_int,close_tag=7,close_id
  integer(kind=8),dimension(0:ndim)::hash_child
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status
#endif
  
#ifndef WITHOUTMPI

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
#ifdef MDL20
    call mdl_cache_close(mdl%mdl2,0)
#else
  ! EMPTY AND CLEAN THE CACHE
  do icache=1,m%ncache
     igrid=r%ngridmax+icache
     m%locked(icache)=.false.
     if(m%occupied(icache))call destage(s,igrid,hash_dict)
     m%occupied(icache)=.false.
     m%dirty(icache)=.false.
  end do
  m%free_cache=1
  m%ncache=0

  do icache=1,m%nnull
     if(m%occupied_null(icache))then
        hash_child(0)=m%lev_null(icache)
        hash_child(1:ndim)=m%ckey_null(1:ndim,icache)
        call hash_free(hash_dict,hash_child)
     endif
     m%occupied_null(icache)=.false.
  end do
  m%free_null=1
  m%nnull=0

  ! COMPLETE THE LAST FLUSH
  do icpu=1,g%ncpu
     iskip=mdl%size_flush_array*(icpu-1)+1
     nflush=mdl%send_flush_array(iskip)
     if(nflush>0)then
        ! Post send
        call MPI_ISSEND(mdl%send_flush_array(iskip),mdl%size_flush_array,MPI_INTEGER,icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
        ! While waiting for completion, check on incoming messages and perform actions
        call check_mail(s,send_flush_id,hash_dict)
        mdl%send_flush_array(iskip)=0
     endif
  end do
  
  ! CHECK-IN CHECK-OUT
  if(g%myid.NE.1)then
     call MPI_ISEND(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(s,close_id,hash_dict)
     call MPI_IRECV(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(s,close_id,hash_dict)
  else
     do icpu=2,g%ncpu
        call MPI_IRECV(dummy_int,1,MPI_INTEGER,MPI_ANY_SOURCE,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(s,close_id,hash_dict)
     end do
     do icpu=2,g%ncpu
        call MPI_ISEND(dummy_int,1,MPI_INTEGER,icpu-1,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(s,close_id,hash_dict)
     end do
  endif

  ! Barrier to get the last flush message
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  call check_mail(s,MPI_REQUEST_NULL,hash_dict)

  ! Finally CANCEL THE 2 RECV
  call MPI_CANCEL(mdl%request_id,info)
  call MPI_CANCEL(mdl%flush_id,info)

  ! Test to free memory in corresponding MPI buffer
  call MPI_WAIT(mdl%request_id,request_status,info)
  call MPI_WAIT(mdl%flush_id,flush_status,info)
  do icpu=1,g%ncpu
     call MPI_WAIT(mdl%reply_id(icpu),reply_status,info)
  end do
  
  ! Barrier to prevent interference with the next cache
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  end associate
  
#endif
#endif
end subroutine close_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module cache

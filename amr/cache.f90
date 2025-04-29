module cache
  type cache_key_ptr
     integer(kind=8),dimension(:),pointer::p
  end type cache_key_ptr
contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_tile(s,child,nkey,keys,grid,ntile)
  use ramses_commons, only: ramses_t
  use amr_commons, only: nbor,oct
  use amr_parameters, only: ndim
  use cache_commons, only: ntilemax
  implicit none
  type(ramses_t)::s
  type(oct),pointer,intent(in)::child
  integer,value::nkey
  type(cache_key_ptr),dimension(1:nkey),intent(inout)::keys
  type(nbor),dimension(1:nkey),intent(out)::grid
  integer,intent(out)::ntile
  integer::igrid,itile,ilevel,ipos,i

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
    ilevel = child%lev
    igrid=(loc(child)-loc(m%grid(1)))/(loc(m%grid(2))-loc(m%grid(1)))+1
    if(igrid.GE.m%head_cache(ilevel).and.igrid.LE.m%tail_cache(ilevel))then
       itile=(igrid-m%head_cache(ilevel))/ntilemax
       ntile=MIN(m%tail_cache(ilevel)-itile*ntilemax-m%head_cache(ilevel)+1,ntilemax)
       do i=1,ntile
          ipos=m%head_cache(ilevel)+itile*ntilemax+i-1
          keys(i)%p(0) = m%grid(ipos)%lev
          keys(i)%p(1:ndim) = m%grid(ipos)%ckey(1:ndim)
          grid(i)%p => m%grid(ipos)
       end do
    else
       ntile=1
       keys(1)%p(0) = child%lev
       keys(1)%p(1:ndim) = child%ckey(1:ndim)
       grid(1)%p => child
    endif
  end associate

end subroutine get_tile
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache(s,hash_dict)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
  use mdl_module
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(hash_table)::hash_dict
  !
  ! This routine closes all cache operations.
  ! It purges all remaining flush messages.
  !
  integer::info,icache,igrid,icpu,ibuf,iskip,ipeak
  integer::send_flush_id,send_flush_id_clump,nflush
  integer::dummy_int,close_tag=7,close_id
  integer(kind=8),dimension(0:ndim)::hash_child
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status
  integer,dimension(MPI_STATUS_SIZE)::reply_status_clump,request_status_clump,flush_status_clump
#endif
  
  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  ! EMPTY AND CLEAN THE GRID CACHE
  if(mdl%cache_opened)then
     do icache=1,m%ncache
        igrid=r%ngridmax+icache
        m%locked(icache)=.false.
        if(m%occupied(icache))call destage(s,igrid,hash_dict)
        m%occupied(icache)=.false.
        m%dirty(icache)=.false.
        m%ghost_parent_grid(icache)=0
        m%ghost_parent_cell(icache)=0
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
  endif

  ! EMPTY AND CLEAN THE CLUMP CACHE
  if(mdl%cache_opened_clump)then
     do icache=1,c%ncache
        ipeak=c%npeak+icache
        c%locked(icache)=.false.
        if(c%occupied(icache))call destage_clump(s,ipeak,hash_dict)
        c%occupied(icache)=.false.
        c%dirty(icache)=.false.
     end do
     c%free_cache=1
     c%ncache=0
  endif

#ifndef WITHOUTMPI

  if(mdl%cache_opened)then
     ! COMPLETE THE LAST GRID FLUSH
     do icpu=1,g%ncpu
        ibuf=mdl%cpu2buf_flush(icpu)
        if(ibuf>0)then
           iskip=1
           nflush=mdl%send_flush(ibuf)%array(iskip)
           if(nflush>0)then
              ! Post send
              call MPI_ISSEND(mdl%send_flush(ibuf)%array(iskip),mdl%size_flush_array,MPI_INTEGER,icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)
              ! While waiting for completion, check on incoming messages and perform actions
              call check_mail(s,send_flush_id,hash_dict)
              mdl%send_flush(ibuf)%array(iskip)=0
           endif
        endif
     end do
  endif

  if(mdl%cache_opened_clump)then
     ! COMPLETE THE LAST CLUMP FLUSH
     do icpu=1,g%ncpu
        ibuf=mdl%cpu2buf_flush_clump(icpu)
        if(ibuf>0)then
           iskip=1
           nflush=mdl%send_flush_clump(ibuf)%array(iskip)
           if(nflush>0)then
              ! Post send
              call MPI_ISSEND(mdl%send_flush_clump(ibuf)%array(iskip),mdl%size_flush_array_clump,MPI_INTEGER,icpu-1,flush_tag_clump,MPI_COMM_WORLD,send_flush_id_clump,info)
              ! While waiting for completion, check on incoming messages and perform actions
              call check_mail(s,send_flush_id_clump,hash_dict)
              mdl%send_flush_clump(ibuf)%array(iskip)=0
           endif
        endif
     end do
  endif

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

  ! Finally CANCEL THE 2 GRID RECV
  if(mdl%cache_opened)then

     call MPI_CANCEL(mdl%request_id,info)
     call MPI_CANCEL(mdl%flush_id,info)

     ! Test to free memory in corresponding MPI buffer
     call MPI_WAIT(mdl%request_id,request_status,info)
     call MPI_WAIT(mdl%flush_id,flush_status,info)
     do icpu=1,g%ncpu
        call MPI_WAIT(mdl%reply_id(icpu),reply_status,info)
     end do

     ! Reset cpu mapping to flush and fetch buffers
     do icpu=1,g%ncpu
        mdl%cpu2buf_fetch(icpu)=0
        mdl%cpu2buf_flush(icpu)=0
     end do
     mdl%ibuffer_fetch=0
     mdl%ibuffer_flush=0
  endif

  ! Finally CANCEL THE 2 CLUMP RECV
  if(mdl%cache_opened_clump)then

     call MPI_CANCEL(mdl%request_id_clump,info)
     call MPI_CANCEL(mdl%flush_id_clump,info)

     ! Test to free memory in corresponding MPI buffer
     call MPI_WAIT(mdl%request_id_clump,request_status_clump,info)
     call MPI_WAIT(mdl%flush_id_clump,flush_status_clump,info)
     do icpu=1,g%ncpu
        call MPI_WAIT(mdl%reply_id_clump(icpu),reply_status_clump,info)
     end do

     ! Reset cpu mapping to flush and fetch buffers
     do icpu=1,g%ncpu
        mdl%cpu2buf_fetch_clump(icpu)=0
        mdl%cpu2buf_flush_clump(icpu)=0
     end do
     mdl%ibuffer_fetch_clump=0
     mdl%ibuffer_flush_clump=0
  endif

  ! Barrier to prevent interference with the next cache
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  mdl%cache_opened=.false.
  mdl%cache_opened_clump=.false.

  end associate
  
end subroutine close_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
integer function get_thread_id(s,uhash,size,hash_key) result(grid_cpu)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use ramses_commons, only: ramses_t
  use hilbert
  use mdl_module
  implicit none
  type(ramses_t)::s
  integer,value::uhash, size
  integer(kind=8),dimension(0:ndim),intent(in)::hash_key
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer::ilevel
  logical::in_rank

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
    ilevel=hash_key(0)
    ix(1:ndim)=hash_key(1:ndim)
    hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

    grid_cpu = mdl_self(mdl)
    in_rank = ge_keys(hk,m%domain_hilbert(ilevel)%b(1:nhilbert,mdl_self(mdl)-1)).and. &
       &    gt_keys(m%domain_hilbert(ilevel)%b(1:nhilbert,grid_cpu),hk)

    ! Determine parent processor
    if (.not.in_rank) grid_cpu = m%domain_hilbert(ilevel)%get_rank(hk)

    grid_cpu = grid_cpu - 1 ! Fortran is 1 based
  end associate

end function get_thread_id
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache(s,table,data_size,hilbert,pack_size,pack,unpack,init,flush,combine,bound)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use ramses_commons, only: ramses_t
  use cache_commons
  use call_back, only: cache_function, cache_f
  use domain_m, only: domain_t
  use hash
  use mdl_module
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(hash_table)::table
  type(domain_t),pointer,dimension(:)::hilbert
  integer::data_size,pack_size
  procedure(cache_function_unpack)::unpack
  procedure(cache_function)::pack
  procedure(cache_function),optional::flush
  procedure(cache_function_init),optional::init
  procedure(cache_function_unpack),optional::combine
  procedure(cache_function_bound),optional::bound
  integer::info,icpu,iskip

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    m%domain_hilbert => hilbert

    if (loc(m%domain_hilbert) .eq. loc(m%domain)) then
       m%head_cache(r%levelmin:r%nlevelmax)=m%head
       m%tail_cache(r%levelmin:r%nlevelmax)=m%tail
    else if (loc(m%domain_hilbert) .eq. loc(m%domain_mg)) then
       m%head_cache(1:r%nlevelmax)=m%head_mg
       m%tail_cache(1:r%nlevelmax)=m%tail_mg
    else
       write(*,*) 'Unknown domain decomposition scheme'
       stop
    end if

    mdl%size_msg_array = pack_size

    pack_fetch%proc => pack
    unpack_fetch%proc => unpack
    init_flush%proc => null()
    pack_flush%proc => null()
    init_bound%proc => null()
    unpack_flush%proc => null()
    if (present(init))      init_flush%proc => init
    if (present(flush))     pack_flush%proc => flush
    if (present(bound))     init_bound%proc => bound
    if (present(combine)) unpack_flush%proc => combine
    mdl%cache_opened=.true.

#ifndef WITHOUTMPI
    do icpu=1,g%ncpu
      mdl%reply_id(icpu)=MPI_REQUEST_NULL
    end do

    mdl%mail_counter=0

    if (.not.present(init) .and. present(flush)) then
      mdl%combiner_rule = COMBINER_CREATE
    else
      mdl%combiner_rule = COMBINER_EXIST
    end if
    mdl%size_request_array=1+ndim
    mdl%size_flush_array=1+(1+ndim+mdl%size_msg_array)*nflushmax
    mdl%size_fetch_array=2+(1+ndim+mdl%size_msg_array)*ntilemax

    ! Set communication counters to zero
    do icpu=1,mdl%nbuffer_flush
      mdl%send_flush(icpu)%array(1)=0
    end do
  
    ! Post the first RECV for request
    call MPI_IRECV(mdl%recv_request_array,mdl%size_request_array,MPI_INTEGER,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,mdl%request_id,info)
  
    ! Post the first RECV for flush
    call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)
#endif

  end associate

end subroutine open_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache_clump(s,pack_size,pack,unpack,init,flush,combine)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use ramses_commons, only: ramses_t
  use cache_commons
  use call_back, only: cache_function, cache_f
  use domain_m, only: domain_t
  use hash
  use mdl_module
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer::pack_size
  procedure(cache_function_clump),optional::pack
  procedure(cache_function_unpack_clump),optional::unpack
  procedure(cache_function_clump),optional::flush
  procedure(cache_function_init_clump),optional::init
  procedure(cache_function_unpack_clump),optional::combine
  integer::info,icpu,iskip

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  mdl%size_msg_array_clump = pack_size

  pack_fetch_clump%proc => null()
  unpack_fetch_clump%proc => null()
  init_flush_clump%proc => null()
  pack_flush_clump%proc => null()
  unpack_flush_clump%proc => null()
  if (present(pack))      pack_fetch_clump%proc => pack
  if (present(unpack))  unpack_fetch_clump%proc => unpack
  if (present(init))      init_flush_clump%proc => init
  if (present(flush))     pack_flush_clump%proc => flush
  if (present(combine)) unpack_flush_clump%proc => combine

  mdl%cache_opened_clump = .true.

#ifndef WITHOUTMPI
  do icpu=1,g%ncpu
     mdl%reply_id_clump(icpu)=MPI_REQUEST_NULL
  end do

  mdl%mail_counter=0

  mdl%size_request_array_clump=2
  mdl%size_flush_array_clump=1+(2+mdl%size_msg_array_clump)*nflushmax
  mdl%size_fetch_array_clump=2+(2+mdl%size_msg_array_clump)*ntilemax

  ! Set communication counters to zero
  do icpu=1,mdl%nbuffer_flush_clump
     mdl%send_flush_clump(icpu)%array(1)=0
  end do

  ! Post the first RECV for request
  call MPI_IRECV(mdl%recv_request_array_clump,mdl%size_request_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,request_tag_clump,MPI_COMM_WORLD,mdl%request_id_clump,info)\

  ! Post the first RECV for flush
  call MPI_IRECV(mdl%recv_flush_array_clump,mdl%size_flush_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag_clump,MPI_COMM_WORLD,mdl%flush_id_clump,info)
#endif

  end associate

end subroutine open_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module cache

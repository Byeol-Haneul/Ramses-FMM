!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine check_mail(s,comm_id,hash_dict)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_F_POINTER, C_ASSOCIATED
  use mdl_module
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use amr_commons, only: oct,nbor
  use hilbert
  use hash
  use cache, only:cache_key_ptr,get_tile
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(hash_table)::hash_dict
  integer::comm_id
  !
  ! This routine checks for incoming messages.
  ! It can be a cache request, and the routine
  ! assembles a fetch message and sends it back.
  ! It can be a flush message, and the routine
  ! unpacks it and combine it in the local memory.
  !
  integer::i,ind,ivar,idim,info,ipos,iskip,igrid,ichild,grid_cpu,ilevel,itile,ntile_reply,nflush
  integer::buffer_size_fetch_array,buffer_size_msg_array,ibuf
  logical::comm_completed,request_received,flush_received=.false.
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status,comm_status
#endif
  integer(kind=8), dimension(1:nhilbert)::hk
  integer(kind=8), dimension(1:ndim)::ix
  integer(kind=8), dimension(0:ndim)::hash_key,hash_child
  type(oct),pointer::child
  integer(kind=8), dimension(0:ndim,1:ntilemax),target::raw_keys
  type(cache_key_ptr),dimension(1:ntilemax)::keys
  type(nbor),dimension(1:ntilemax)::grid
  type(msg_large_realdp)::dummy_large_realdp

#ifndef WITHOUTMPI

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  comm_completed=.false.
  do while (.not. comm_completed)

     ! USE A NESTED DO WHILE LOOP FOR THE 2 ADD. TESTS

     !===========================
     ! Check for incoming request
     !===========================
     request_received=.true.
     do while(request_received)
        call MPI_TEST(mdl%request_id,request_received,request_status,info)
        if(request_received)then
           
           ! Assemble a reply and send it back
           ilevel=mdl%recv_request_array(1)
           hash_key(0)=ilevel
           hash_key(1:ndim)=mdl%recv_request_array(2:ndim+1)
           call c_f_pointer(hash_getp(hash_dict,hash_key),child)
           grid_cpu=request_status(MPI_SOURCE)+1

           ! Identify the corresponding communication buffer
           ibuf=mdl%cpu2buf_fetch(grid_cpu)

           ! If none is available, get the next free one
           if(ibuf==0)then
              mdl%ibuffer_fetch=mdl%ibuffer_fetch+1
              ! If there is no communication buffer left, allocate a new one
              if(mdl%ibuffer_fetch>mdl%nbuffer_fetch)then
                 mdl%nbuffer_fetch=mdl%nbuffer_fetch+1
                 buffer_size_msg_array=storage_size(dummy_large_realdp)/32
                 buffer_size_fetch_array=2+(1+ndim+buffer_size_msg_array)*ntilemax
                 allocate(mdl%send_fetch(mdl%nbuffer_fetch)%array(1:buffer_size_fetch_array))
              endif
              mdl%cpu2buf_fetch(grid_cpu)=mdl%ibuffer_fetch
              ibuf=mdl%ibuffer_fetch
           endif

           ! If grid does not exist, send a null reply
           if(.not.ASSOCIATED(child))then
              
              ! Store type corresponding to a null reply
              iskip=1
              mdl%send_fetch(ibuf)%array(iskip)=-1

           ! Otherwise, assemble a proper reply with a complete tile
           else
              iskip=1
              ! Record the location of where we want the keys to be stored, then ask for a tile
              do i=1,ntilemax
!                 ipos = iskip + 2 + (i-1)*(1+ndim+mdl%size_msg_array)
!                 keys(i)%p(0:ndim) => mdl%send_fetch(ibuf)%array(ipos:ipos+ndim)
                 keys(i)%p(0:ndim) => raw_keys(0:ndim,i)
              end do
              call get_tile(s,child,ntilemax,keys,grid,ntile_reply)

              ! Store type of reply and number of entries
              mdl%send_fetch(ibuf)%array(iskip)=1
              iskip=iskip+1
              mdl%send_fetch(ibuf)%array(iskip)=ntile_reply
              iskip=iskip+1

              ! Store data, depending on reply type
              do i=1,ntile_reply
                 mdl%send_fetch(ibuf)%array(iskip:iskip+ndim)=raw_keys(0:ndim,i)
                 iskip=iskip+1+ndim ! Skip over the key already present
                 call pack_fetch%proc(grid(i)%p,mdl%size_msg_array,mdl%send_fetch(ibuf)%array(iskip:iskip+mdl%size_msg_array-1))
                 iskip=iskip+mdl%size_msg_array
              end do
           endif
           
           ! Wait for the old SEND to free memory in corresponding MPI buffer
           call MPI_WAIT(mdl%reply_id(grid_cpu),reply_status,info)
           
           ! Send back the reply
           iskip=1
           call MPI_ISEND(mdl%send_fetch(ibuf)%array(iskip),mdl%size_fetch_array,MPI_INTEGER,grid_cpu-1,msg_tag,MPI_COMM_WORLD,mdl%reply_id(grid_cpu),info)
           
           !=================================
           ! Post a new RECV for request
           !=================================
           call MPI_IRECV(mdl%recv_request_array,mdl%size_request_array,MPI_INTEGER,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,mdl%request_id,info)

        endif
     end do

     !=========================
     ! Check for incoming flush
     !=========================
     flush_received=.true.
     do while(flush_received)
        call MPI_TEST(mdl%flush_id,flush_received,flush_status,info)
        if(flush_received)then
           
           ! Combine received data to local memory only if grid exists
           if(mdl%combiner_rule.eq.COMBINER_EXIST)then
              iskip=1
              nflush=mdl%recv_flush_array(iskip)
              iskip=iskip+1

              do i=1,nflush
                 ilevel=mdl%recv_flush_array(iskip)
                 hash_child(0)=ilevel
                 hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
                 iskip=iskip+ndim+1
                 ! Get grid from hash table
                 call c_f_pointer(hash_getp(hash_dict,hash_child),child)
                 if(ASSOCIATED(child))then
                    call unpack_flush%proc(child,mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)
                 endif

                 iskip=iskip+mdl%size_msg_array

              end do
              
           endif
           
           ! Combine received data to local memory only if grid does not exist
           if(mdl%combiner_rule.eq.COMBINER_CREATE)then
              iskip=1
              nflush=mdl%recv_flush_array(iskip)
              iskip=iskip+1

              do i=1,nflush
                 ilevel=mdl%recv_flush_array(iskip)
                 hash_child(0)=ilevel
                 hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
                 iskip=iskip+ndim+1

                 ! Create new grid if grid does not exist
                 if(.not.C_ASSOCIATED(hash_getp(hash_dict,hash_child)))then
                    
                    ! Compute Hilbert keys of new octs
                    ix(1:ndim)=hash_child(1:ndim)
                    hk(1:nhilbert)=hilbert_key(ix,ilevel-1)
                    
                    ! Set grid index to a virtual grid in local main memory
                    ichild=m%ifree
                    
                    ! Go to next main memory free line
                    m%ifree=m%ifree+1
                    if(m%ifree.GT.r%ngridmax)then
                       write(*,*)'No more free memory'
                       write(*,*)'Increase ngridmax'
                       call mdl_abort(mdl)
                    endif
                    ! TODO: THIS IS NOW DONE IN UNPACK / INIT and seems to work (but should be checked)
                    !m%grid(ichild)%lev=hash_child(0)
                    !m%grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
                    m%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
                    m%grid(ichild)%superoct=1
                    m%grid(ichild)%flag1(1:twotondim)=0
                    m%grid(ichild)%flag2(1:twotondim)=0
                    
                    ! Insert new grid in hash table
                    call hash_setp(hash_dict,hash_child,m%grid(ichild))

                    ! Unpack message content
                    call unpack_flush%proc(m%grid(ichild),mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)
                    
                 endif
                 
                 iskip=iskip+mdl%size_msg_array
                 
              end do

           endif

           !=================================
           ! Post a new RECV for flush
           !=================================
           call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)
           
        endif
     end do
     
     !=================================
     ! Check for input comm. completion
     !=================================
     if(comm_id==MPI_REQUEST_NULL)then
        comm_completed=.true.
     else
        call MPI_TEST(comm_id,comm_completed,comm_status,info)
     endif
  end do

  end associate
  
#endif
#endif
end subroutine check_mail
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine destage(s,igrid,hash_dict)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_F_POINTER, C_ASSOCIATED
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  type(hash_table)::hash_dict
  integer::igrid
  !
  ! This routine frees the cache memory.
  ! It assembles flush messages, and when the message
  ! buffer is full, it sends it to the target CPU.
  !
  integer::ind,ivar,idim,info,icache,iflush,grid_cpu,ibuf
  integer::send_flush_id,iskip,nflush
  integer::buffer_size_flush_array,buffer_size_msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_large_realdp)::dummy_large_realdp

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  hash_key(0)=m%grid(igrid)%lev
  hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)

  if(.not.C_ASSOCIATED(hash_getp(hash_dict,hash_key)))then
     write(*,*)'PE ',g%myid,' trying to free non existing grid'
     stop
  endif

  call hash_free(hash_dict,hash_key)

#ifndef WITHOUTMPI

  ! Check if the destage requires a flush
  icache=igrid-r%ngridmax

  if(m%dirty(icache))then
     grid_cpu=m%parent_cpu(icache)
     m%dirty(icache)=.false.
  
     ! Identify the corresponding communication buffer
     ibuf=mdl%cpu2buf_flush(grid_cpu)
 
     ! If none is available, get the next free one
     if(ibuf==0)then
        mdl%ibuffer_flush=mdl%ibuffer_flush+1
        ! If there is no communication buffer left, allocate a new one
        if(mdl%ibuffer_flush>mdl%nbuffer_flush)then
           mdl%nbuffer_flush=mdl%nbuffer_flush+1
           buffer_size_msg_array=storage_size(dummy_large_realdp)/32
           buffer_size_flush_array=1+(1+ndim+buffer_size_msg_array)*nflushmax
           allocate(mdl%send_flush(mdl%nbuffer_flush)%array(1:buffer_size_flush_array))
        endif
        mdl%cpu2buf_flush(grid_cpu)=mdl%ibuffer_flush
        ibuf=mdl%ibuffer_flush
     endif

     ! Filling the flush buffer
     iskip=1
     nflush=mdl%send_flush(ibuf)%array(iskip)

     ! If buffer full, send it to remote CPU.
     if(nflush==nflushmax)then
        ! Post send
        call MPI_ISSEND(mdl%send_flush(ibuf)%array(iskip),mdl%size_flush_array,MPI_INTEGER,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)
        ! While waiting for completion, check on incoming messages and perform actions
        call check_mail(s,send_flush_id,hash_dict)
        ! Reset counter
        mdl%send_flush(ibuf)%array(iskip)=0
     endif
     
     ! Increment counter
     mdl%send_flush(ibuf)%array(iskip)=mdl%send_flush(ibuf)%array(iskip)+1

     ! Skip to the last available position
     iflush=mdl%send_flush(ibuf)%array(iskip)
     iskip=iskip+(1+ndim+mdl%size_msg_array)*(iflush-1)+1
     
     ! Pack message header
     mdl%send_flush(ibuf)%array(iskip)=m%grid(igrid)%lev
     mdl%send_flush(ibuf)%array(iskip+1:iskip+ndim)=m%grid(igrid)%ckey(1:ndim)
     iskip=iskip+ndim+1

     ! Pack message content
     call pack_flush%proc(m%grid(igrid),mdl%size_msg_array,mdl%send_flush(ibuf)%array(iskip:iskip+mdl%size_msg_array-1))

  endif
#endif

  end associate

#endif
end subroutine destage
!##############################################################
!##############################################################
!##############################################################
!##############################################################

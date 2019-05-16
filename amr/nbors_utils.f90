!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_threetondim_nbor_parent_cell(s,hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,twotondim,threetondim
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(1:threetondim)::igrid_nbor,ind_nbor
  !
  ! This routine computes and acquire the 3**ndim neighboring father cells 
  ! for the input hash_key. The output arrays are the father cells
  ! parent oct indices and their associated cell indices within the oct.
  ! The corresponding data can be accessed using: grid(igrid)%data(ind).
  ! If the grid index is zero, it means that this oct does not exist.
  ! Note that the 2**ndim grids are all locked if remote.
  !
  integer,dimension(1:twotondim)::igrid_twotondim_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor,hash_ref,hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer,dimension(1:3,1:8)::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer::i1,j1,k1
  integer::i1min=-1
  integer::i1max=+1
  integer::j1min=0*(1-ndim/2)-1*(ndim/2)
  integer::j1max=0*(1-ndim/2)+1*(ndim/2)
  integer::k1min=0*(1-ndim/3)-1*(ndim/3)
  integer::k1max=0*(1-ndim/3)+1*(ndim/3)
  integer::ind,ipos,idim,ilevel,inbor
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ilevel=hash_key(0)

  hash_father(0)=hash_key(0)-1

  ! Gather twotondim neighboring father grids
  do inbor=1,twotondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift_oct(1:ndim,inbor)
     ! Periodic boundary conditons
     do idim=1,ndim
        if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
        if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
     enddo
     hash_father(1:ndim)=hash_nbor(1:ndim)/2
     ! Store lower left neighbor coordinates 
     if(inbor==1)hash_ref(1:ndim)=hash_father(1:ndim)
     ! Get grid into memory and lock it if remote 
     call get_grid(s,hash_father,hash_dict,ipos,flush_cache,fetch_cache)
     call lock_cache(s,ipos)
     igrid_twotondim_nbor(inbor)=ipos
  end do
     
  ! Deal with neighboring father cells
  inbor=0
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max           
           inbor=inbor+1
#if NDIM>0
           hash_nbor(1)=hash_key(1)+i1
#endif
#if NDIM>1
           hash_nbor(2)=hash_key(2)+j1
#endif
#if NDIM>2
           hash_nbor(3)=hash_key(3)+k1
#endif
           ! Periodic boundary conditons
           do idim=1,ndim
              if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
              if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
           enddo
           ! Compute neighboring cell index
           hash_father(1:ndim)=hash_nbor(1:ndim)/2
           ii(1:ndim)=hash_nbor(1:ndim)-2*hash_father(1:ndim)
           ind=1
           do idim=1,ndim
              ind=ind+2**(idim-1)*ii(idim)
           end do
           ind_nbor(inbor)=ind
           ! Compute neighboring grid index
           ii(1:ndim)=hash_father(1:ndim)-hash_ref(1:ndim)
           ! Periodic boundary conditons
           do idim=1,ndim
              if(ii(idim)<0)ii(idim)=ii(idim)+m%ckey_max(ilevel-1)
           enddo
           ind=1
           do idim=1,ndim
              ind=ind+2**(idim-1)*ii(idim)
           end do
           igrid_nbor(inbor)=igrid_twotondim_nbor(ind)
        end do
     end do
  end do

  end associate
  
end subroutine get_threetondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_twondim_nbor_parent_cell(s,hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,twotondim,twondim
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(0:twondim)::igrid_nbor,ind_nbor
  !
  ! This routine computes and acquires the 2xndim neighboring father cells 
  ! for the input hash_key. The output arrays are the father cells
  ! parent oct indices and their associated cell indices within the oct.
  ! The corresponding data can be accessed using: grid(igrid)%data(ind).
  ! The first element (0) stands for the central father cell.
  ! If the grid index is zero, it means that this oct does not exist.
  ! Note that the parent grids are all locked if remote.
  !
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer(kind=8),dimension(0:ndim)::hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer,dimension(1:3,1:6),save::shift=reshape((/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::ind,ipos,idim,ilevel,inbor
  
  associate(r=>s%r,g=>s%g,m=>s%m)

  ilevel=hash_key(0)

  ! Deal with central parent cell first
  hash_father(0)=hash_key(0)-1
  hash_father(1:ndim)=hash_key(1:ndim)/2
  ii(1:ndim)=hash_key(1:ndim)-2*hash_father(1:ndim)
  ind=1
  do idim=1,ndim
     ind=ind+2**(idim-1)*ii(idim)
  end do

  ! Get grid into memory and lock it if remote 
  call get_grid(s,hash_father,hash_dict,ipos,flush_cache,fetch_cache)
  call lock_cache(s,ipos)
  igrid_nbor(0)=ipos
  ind_nbor(0)=ind
  
  ! Deal with neighboring father cells
  do inbor=1,twondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift(1:ndim,inbor)

     ! Periodic boundary conditons
     do idim=1,ndim
        if(hash_nbor(idim)<0)hash_nbor(idim)=m%ckey_max(ilevel)-1
        if(hash_nbor(idim)==m%ckey_max(ilevel))hash_nbor(idim)=0
     enddo
     hash_father(1:ndim)=hash_nbor(1:ndim)/2
     ii(1:ndim)=hash_nbor(1:ndim)-2*hash_father(1:ndim)
     ind=1
     do idim=1,ndim
        ind=ind+2**(idim-1)*ii(idim)
     end do

     ! Get grid into memory and lock it if remote 
     call get_grid(s,hash_father,hash_dict,ipos,flush_cache,fetch_cache)
     call lock_cache(s,ipos)
     igrid_nbor(inbor)=ipos
     ind_nbor(inbor)=ind
  end do

  end associate
  
end subroutine get_twondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_parent_cell(s,hash_key,hash_dict,igrid,ind,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,twotondim
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer::igrid,ind
  !
  ! This routine acquires the parent cell of the grid 
  ! corresponding to the input hash key.
  !
  integer(kind=8),dimension(0:ndim)::hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer::idim
  hash_father(0)=hash_key(0)-1
  hash_father(1:ndim)=hash_key(1:ndim)/2
  ii(1:ndim)=hash_key(1:ndim)-2*hash_father(1:ndim)
  ind=1
  do idim=1,ndim
     ind=ind+2**(idim-1)*ii(idim)
  end do
  call get_grid(s,hash_father,hash_dict,igrid,flush_cache,fetch_cache)
end subroutine get_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache(s,child_grid)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::child_grid
  !
  ! This routine locks a cache line because
  ! it will be updated later.
  !
  integer::icache
  if(child_grid>s%r%ngridmax)then
     icache=child_grid-s%r%ngridmax
     s%m%locked(icache)=.true.
  endif
end subroutine lock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine unlock_cache(s,child_grid)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer::child_grid
  !
  ! This routine unlocks a cache line because
  ! it has been updated and can be flushed.
  !
  integer::icache
  if(child_grid>s%r%ngridmax)then
     icache=child_grid-s%r%ngridmax
     s%m%locked(icache)=.false.
  endif
end subroutine unlock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_grid(s,hash_key,hash_dict,child_grid,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  integer::child_grid
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
#ifndef MDL2
  !
  ! This routine acquires the grid 
  ! corresponding to the input hash key.
  ! Two logicals set the type of acquire.
  ! "fetch" means read-only, "flush" means write-only.
  ! and both used together means read-write.
  !
#ifndef WITHOUTMPI
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_child
  integer::i,ind,idim,ivar,iskip,ichild,ilevel,info,grid_cpu,ntile_response,icounter
  integer::send_request_id,response_id  
  logical::failed_request,in_rank
  integer,dimension(MPI_STATUS_SIZE)::send_request_status
#endif
 
  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(s,MPI_REQUEST_NULL,hash_dict)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Access hash table
  child_grid=hash_get(hash_dict,hash_key)

#ifndef WITHOUTMPI

  ! If grid index is positive, then return
  if(child_grid>0)then
     return
  endif

  ! If grid index is -1, then set it to 0 and return
  ! This means we already know the remote grid does not exist
  if(child_grid.EQ.-1)then
     child_grid=0
     return
  endif

  ! Now we know child_grid=0

  ! Compute the Hilbert key
  ilevel=hash_key(0)
  ix(1:ndim)=hash_key(1:ndim)
  hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

  ! Check if grid sits inside processor boundaries
!  if (m%domain_hilbert(ilevel)%in_rank(hk)) return
  in_rank = ge_keys(hk,m%domain_hilbert(ilevel)%b(1:nhilbert,mdl%myid-1)).and. &
       &    gt_keys(m%domain_hilbert(ilevel)%b(1:nhilbert,mdl%myid),hk)
  if (in_rank) return

  ! Determine parent processor
  grid_cpu = m%domain_hilbert(ilevel)%get_rank(hk)

  !============================================
  ! We have a fetch and possibly a flush cache
  ! Fetch alone means read-only cache operations.
  ! Both together means read-write cache.
  !============================================
  if(fetch_cache)then

     ! Send a request to the relevant cpu
     mdl%send_request_array(1)=ilevel
     mdl%send_request_array(2:ndim+1)=hash_key(1:ndim)
     
     ! Post RECV for the expected response
     call MPI_IRECV(mdl%recv_fetch_array,mdl%size_fetch_array,MPI_INTEGER,grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  

     ! Post SEND for the request
     call MPI_ISEND(mdl%send_request_array,mdl%size_request_array,MPI_INTEGER,grid_cpu-1,request_tag,MPI_COMM_WORLD,send_request_id,info)

     ! While waiting for reply, check on incoming messages and perform actions
     call check_mail(s,response_id,hash_dict)

     ! Wait for ISEND completion to free memory in corresponding MPI buffer
     call MPI_WAIT(send_request_id,send_request_status,info)

     ! Check header for type of response
     iskip=1
     failed_request=(mdl%recv_fetch_array(iskip)==-1)
     iskip=iskip+1

     ! If grid does not exist, store -1 in the cache
     ! The output grid index is still zero
     if(failed_request)then

        ! Delete old null grid if occupied
        if(m%occupied_null(m%free_null))then
           hash_child(0)=m%lev_null(m%free_null)
           hash_child(1:ndim)=m%ckey_null(1:ndim,m%free_null)
           call hash_free(hash_dict,hash_child)
        endif
        call hash_set(hash_dict,hash_key,-1)
        m%occupied_null(m%free_null)=.true.
        m%lev_null(m%free_null)=ilevel
        m%ckey_null(1:ndim,m%free_null)=hash_key(1:ndim)

        ! Go to next free cache line
        m%free_null=m%free_null+1
        m%nnull=m%nnull+1
        if(m%free_null.GT.r%ncachemax)then
           m%free_null=1
        endif
        if(m%nnull.GT.r%ncachemax)m%nnull=r%ncachemax

     ! If grid exists, store incoming tile in the cache
     else

        ! Number of tiles in the response buffer
        ntile_response=mdl%recv_fetch_array(iskip)
        iskip=iskip+1

        ! Loop over tiles
        do i=1,ntile_response

           ! If next cache line is occupied, free it.
           if(m%locked(m%free_cache))then
              icounter=0
              do while(m%locked(m%free_cache))
                 m%free_cache=m%free_cache+1
                 icounter=icounter+1
                 if(m%free_cache>r%ncachemax)m%free_cache=1
                 if(icounter>r%ncachemax)then
                    write(*,*)'PE ',g%myid,'cache entirely locked'
                    stop
                 endif
              end do
           end if

           ! Next available grid in memory
           ichild=r%ngridmax+m%free_cache

           ! Get grid coordinates from message header
           hash_child(0)=mdl%recv_fetch_array(iskip)
           hash_child(1:ndim)=mdl%recv_fetch_array(iskip+1:iskip+ndim)
           iskip=iskip+ndim+1

           ! If grid does not already exist, create it in local memory
           if(hash_get(hash_dict,hash_child).EQ.0)then

              if(m%occupied(m%free_cache))call destage(s,r%ngridmax+m%free_cache,hash_dict)

              call hash_set(hash_dict,hash_child,ichild)
              
              m%occupied(m%free_cache)=.true.
              m%parent_cpu(m%free_cache)=grid_cpu
              m%dirty(m%free_cache)=.false.
              
              ! Set the grid index of the requested grid
              if(same_keys(hash_key,hash_child))then
                 child_grid=ichild
              endif
              
              ! Store the grid coordinates for the entire tile
              m%grid(ichild)%lev=hash_child(0)
              m%grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
             
              ! Unpack response to fetch request
              call unpack_fetch%proc(m%grid(ichild),mdl%size_msg_array,mdl%recv_fetch_array(iskip:iskip+mdl%size_msg_array-1))
              
              ! If we also have also a flush cache...
              ! This is for combined read-write cache operations
              if(flush_cache)then
                 m%dirty(m%free_cache)=.true.
                 
                 ! Set initialisation rule for combiner operations
                 call init_flush%proc(m%grid(ichild),0)
                 
              endif

              ! Go to next free cache line
              m%free_cache=m%free_cache+1
              m%ncache=m%ncache+1
              if(m%free_cache.GT.r%ncachemax)then
                 m%free_cache=1
              endif
              if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax

           endif

           ! Go to next tile
           iskip=iskip+mdl%size_msg_array

        end do
        ! End loop over tiles

     endif

     !=================================================
     ! If we have only a flush cache (write-only cache) 
     !=================================================
  else if(flush_cache)then   

     ! If next cache line is occupied, free it.
     if(m%locked(m%free_cache))then
        do while(m%locked(m%free_cache))
           m%free_cache=m%free_cache+1
           if(m%free_cache>r%ncachemax)m%free_cache=1
        end do
     end if
     if(m%occupied(m%free_cache))call destage(s,r%ngridmax+m%free_cache,hash_dict)

     ! Set grid index to a virtual grid in local memory
     child_grid=r%ngridmax+m%free_cache
     call hash_set(hash_dict,hash_key,child_grid)

     ! Store the grid coordinates
     m%grid(child_grid)%lev=hash_key(0)
     m%grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
     m%occupied(m%free_cache)=.true.
     m%parent_cpu(m%free_cache)=grid_cpu
     m%dirty(m%free_cache)=.true.

     ! Set initialisation rule for combiner operation
     call init_flush%proc(m%grid(child_grid),0)

     ! Go to next free cache line
     m%free_cache=m%free_cache+1
     m%ncache=m%ncache+1
     if(m%free_cache.GT.r%ncachemax)then
        m%free_cache=1
     endif
     if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax

  endif

#endif

  end associate
#endif
end subroutine get_grid
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine check_mail(s,comm_id,hash_dict)
  use mdl_module
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  logical::comm_completed,request_received,flush_received=.false.
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status,comm_status
#endif
  integer(kind=8), dimension(1:nhilbert)::hk
  integer(kind=8), dimension(1:ndim)::ix
  integer(kind=8), dimension(0:ndim)::hash_key,hash_child

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
           igrid=hash_get(hash_dict,hash_key)
           grid_cpu=request_status(MPI_SOURCE)+1
           
           ! If grid does not exist, send a null reply
           if(igrid.EQ.0)then
              
              ! Store type corresponding to a null reply
              iskip=mdl%size_fetch_array*(grid_cpu-1)+1
              mdl%send_fetch_array(iskip)=-1

           ! Otherwise, assemble a proper reply with a complete tile
           else
              
              itile=(igrid-m%head_cache(ilevel))/ntilemax
              ntile_reply=MIN(m%tail_cache(ilevel)-itile*ntilemax-m%head_cache(ilevel)+1,ntilemax)              

              ! Store type of reply and number of entries
              iskip=mdl%size_fetch_array*(grid_cpu-1)+1
              mdl%send_fetch_array(iskip)=1
              iskip=iskip+1
              mdl%send_fetch_array(iskip)=ntile_reply
              iskip=iskip+1
              
              ! Store data, depending on reply type
              do i=1,ntile_reply
                 ipos=m%head_cache(ilevel)+itile*ntilemax+i-1
                 
                 ! Store message header
                 mdl%send_fetch_array(iskip)=m%grid(ipos)%lev
                 mdl%send_fetch_array(iskip+1:iskip+ndim)=m%grid(ipos)%ckey(1:ndim)
                 iskip=iskip+1+ndim

                 ! Store message content
                 call pack_fetch%proc(m%grid(ipos),mdl%size_msg_array,mdl%send_fetch_array(iskip:iskip+mdl%size_msg_array-1))

                 iskip=iskip+mdl%size_msg_array
              end do
           endif
           
           ! Wait for the old SEND to free memory in corresponding MPI buffer
           call MPI_WAIT(mdl%reply_id(grid_cpu),reply_status,info)
           
           ! Send back the reply
           iskip=mdl%size_fetch_array*(grid_cpu-1)+1
           call MPI_ISEND(mdl%send_fetch_array(iskip),mdl%size_fetch_array,MPI_INTEGER,grid_cpu-1,msg_tag,MPI_COMM_WORLD,mdl%reply_id(grid_cpu),info)
           
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

                 ! Get grid index from hash table
                 ichild=hash_get(hash_dict,hash_child)

                 ! Unpack message content only if grid exists
                 if(ichild>0)then
                    call unpack_flush%proc(m%grid(ichild),mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1))
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
                 if(hash_get(hash_dict,hash_child).EQ.0)then
                    
                    ! Compute Hilbert keys of new octs
                    ix(1:ndim)=hash_child(1:ndim)
                    hk(1:nhilbert)=hilbert_key(ix,ilevel-1)
                    
                    ! Set grid index to a virtual grid in local main memory
                    ichild=m%ifree
                    
                    ! Go to next main memory free line
                    m%ifree=m%ifree+1
                    if(m%ifree.GT.r%ngridmax)then
                       write(*,*)'No more free memory'
                       write(*,*)'while refining...'
                       write(*,*)'Increase ngridmax'
                       call mdl_abort(mdl)
                    endif
                    
                    m%grid(ichild)%lev=hash_child(0)
                    m%grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
                    m%grid(ichild)%hkey(1:nhilbert)=hk(1:nhilbert)
                    m%grid(ichild)%superoct=1
                    m%grid(ichild)%flag1(1:twotondim)=0
                    m%grid(ichild)%flag2(1:twotondim)=0
                    
                    ! Insert new grid in hash table
                    call hash_set(hash_dict,hash_child,ichild)

                    ! Unpack message content
                    call unpack_flush%proc(m%grid(ichild),mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1))
                    
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
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  type(hash_table)::hash_dict
  integer::igrid
  !
  ! This routine frees the cache memory.
  ! It assembles flush messages, and when the message
  ! buffer is full, it sends it to the target CPU.
  !
  integer::ind,ivar,idim,ipos,info,icache,iflush,grid_cpu
  integer::send_flush_id,iskip,nflush
  integer(kind=8),dimension(0:ndim)::hash_key

#ifndef WITHOUTMPI

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  hash_key(0)=m%grid(igrid)%lev
  hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)
  ipos=hash_get(hash_dict,hash_key)

  if(hash_get(hash_dict,hash_key).EQ.0)then
     write(*,*)'PE ',g%myid,' trying to free non existing grid'
     stop
  endif

  call hash_free(hash_dict,hash_key)

  ! Check if the destage requires a flush
  icache=igrid-r%ngridmax

  if(m%dirty(icache))then
     grid_cpu=m%parent_cpu(icache)
     m%dirty(icache)=.false.
  
     ! Filling the flush buffer
     iskip=mdl%size_flush_array*(grid_cpu-1)+1
     nflush=mdl%send_flush_array(iskip)

     ! If buffer full, send it to remote CPU.
     if(nflush==nflushmax)then
        ! Post send
        call MPI_ISSEND(mdl%send_flush_array(iskip),mdl%size_flush_array,MPI_INTEGER,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
        ! While waiting for completion, check on incoming messages and perform actions
        call check_mail(s,send_flush_id,hash_dict)
        ! Reset counter
        mdl%send_flush_array(iskip)=0
     endif
     
     ! Increment counter
     mdl%send_flush_array(iskip)=mdl%send_flush_array(iskip)+1

     ! Skip to the last available position
     iflush=mdl%send_flush_array(iskip)
     iskip=iskip+(1+ndim+mdl%size_msg_array)*(iflush-1)+1
     
     ! Pack message header
     mdl%send_flush_array(iskip)=m%grid(igrid)%lev
     mdl%send_flush_array(iskip+1:iskip+ndim)=m%grid(igrid)%ckey(1:ndim)
     iskip=iskip+ndim+1

     ! Pack message content
     call pack_flush%proc(m%grid(igrid),mdl%size_msg_array,mdl%send_flush_array(iskip:iskip+mdl%size_msg_array-1))

  endif

  end associate

#endif
#endif
end subroutine destage
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache(s,hash_dict)
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  type(hash_table)::hash_dict
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

  end associate
  
#endif
#endif
end subroutine close_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache(s,cache_operation,domain_decompos)
#ifndef MDL2
  use amr_parameters, only: ndim,nhilbert,twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use cache_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(ramses_t)::s
  integer::cache_operation
  integer::domain_decompos
  !
  ! This routine opens a new cache operation.
  ! Depending on the cache operation, it initializes
  ! the size of the communication buffers and the cache functions.
  ! The type of combiner operation is also set.
  !
  integer::info,icpu,iskip
  type(msg_int4)::dummy_int4
  type(msg_realdp)::dummy_realdp
  type(msg_small_realdp)::dummy_small_realdp
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_three_realdp)::dummy_three_realdp
  type(msg_twin_realdp)::dummy_twin_realdp
  
#ifndef WITHOUTMPI

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  do icpu=1,g%ncpu
     mdl%reply_id(icpu)=MPI_REQUEST_NULL
  end do

  mdl%mail_counter=0

  ! Domain decomposition to use
  if(domain_decompos==domain_decompos_amr)then
     m%domain_hilbert => m%domain
     m%head_cache(r%levelmin:r%nlevelmax)=m%head
     m%tail_cache(r%levelmin:r%nlevelmax)=m%tail
  endif
  if(domain_decompos==domain_decompos_mg )then
     m%domain_hilbert => m%domain_mg
     m%head_cache(1:r%nlevelmax)=m%head_mg
     m%tail_cache(1:r%nlevelmax)=m%tail_mg
  endif

  ! Default combiner rule
  mdl%combiner_rule=COMBINER_EXIST
  
  ! Operations of type "flag"
  if(cache_operation.EQ.operation_initflag)then
     mdl%size_msg_array = storage_size(dummy_int4)/32
     init_flush%proc => init_flush_initflag
     pack_fetch%proc => pack_fetch_flag
     unpack_fetch%proc => unpack_fetch_flag
     pack_flush%proc => pack_flush_initflag
     unpack_flush%proc => unpack_flush_initflag
  endif
  if(cache_operation.EQ.operation_smooth)then
     mdl%size_msg_array = storage_size(dummy_int4)/32
     init_flush%proc => null()
     pack_fetch%proc => pack_fetch_flag
     unpack_fetch%proc => unpack_fetch_flag
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif
  if(cache_operation.EQ.operation_derefine)then
     mdl%size_msg_array = storage_size(dummy_int4)/32
     init_flush%proc => init_flush_derefine
     pack_fetch%proc => pack_fetch_flag
     unpack_fetch%proc => unpack_fetch_flag
     pack_flush%proc => pack_flush_derefine
     unpack_flush%proc => unpack_flush_derefine
  endif
  if(cache_operation.EQ.operation_split)then
     mdl%size_msg_array = storage_size(dummy_int4)/32
     init_flush%proc => null()
     pack_fetch%proc => pack_fetch_split
     unpack_fetch%proc => unpack_fetch_split
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif

  ! Operations of type "hydro"
  if(cache_operation.EQ.operation_hydro)then
     mdl%size_msg_array = storage_size(dummy_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_hydro
     unpack_fetch%proc => unpack_fetch_hydro
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif
  if(cache_operation.EQ.operation_upload)then
     mdl%size_msg_array = storage_size(dummy_realdp)/32
     init_flush%proc => init_flush_upload
     pack_fetch%proc => pack_fetch_hydro
     unpack_fetch%proc => unpack_fetch_hydro
     pack_flush%proc => pack_flush_upload
     unpack_flush%proc => unpack_flush_upload
  endif
  if(cache_operation.EQ.operation_multipole)then
     mdl%size_msg_array = storage_size(dummy_realdp)/32
     init_flush%proc => init_flush_multipole
     pack_fetch%proc => pack_fetch_hydro
     unpack_fetch%proc => unpack_fetch_hydro
     pack_flush%proc => pack_flush_multipole
     unpack_flush%proc => unpack_flush_multipole
  endif

  ! Operations of type "poisson"
  if(cache_operation.EQ.operation_cg)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_cg
     unpack_fetch%proc => unpack_fetch_cg
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif
  if(cache_operation.EQ.operation_phi)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_phi
     unpack_fetch%proc => unpack_fetch_phi
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif
  if(cache_operation.EQ.operation_rho)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => init_flush_rho
     pack_fetch%proc => pack_fetch_phi
     unpack_fetch%proc => unpack_fetch_phi
     pack_flush%proc => pack_flush_rho
     unpack_flush%proc => unpack_flush_rho
  endif
  if(cache_operation.EQ.operation_build_mg)then
     mdl%combiner_rule = COMBINER_CREATE
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_phi
     unpack_fetch%proc => unpack_fetch_phi
     pack_flush%proc => pack_flush_build_mg
     unpack_flush%proc => unpack_flush_build_mg
  endif
  if(cache_operation.EQ.operation_restrict_mask)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => init_flush_restrict_mask
     pack_fetch%proc => pack_fetch_phi
     unpack_fetch%proc => unpack_fetch_phi
     pack_flush%proc => pack_flush_restrict_mask
     unpack_flush%proc => unpack_flush_restrict_mask
  endif
  if(cache_operation.EQ.operation_restrict_res)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => init_flush_restrict_res
     pack_fetch%proc => pack_fetch_restrict_res
     unpack_fetch%proc => unpack_fetch_restrict_res
     pack_flush%proc => pack_flush_restrict_res
     unpack_flush%proc => unpack_flush_restrict_res
  endif
  if(cache_operation.EQ.operation_scan)then
     mdl%size_msg_array = storage_size(dummy_small_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_scan
     unpack_fetch%proc => unpack_fetch_scan
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif

  ! Operations of type "refine"
  if(cache_operation.EQ.operation_godunov)then
     mdl%size_msg_array = storage_size(dummy_large_realdp)/32
     init_flush%proc => init_flush_godunov
     pack_fetch%proc => pack_fetch_refine
     unpack_fetch%proc => unpack_fetch_refine
     pack_flush%proc => pack_flush_godunov
     unpack_flush%proc => unpack_flush_godunov
  endif
  if(cache_operation.EQ.operation_refine)then
     mdl%combiner_rule = COMBINER_CREATE
     mdl%size_msg_array = storage_size(dummy_large_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_refine
     unpack_fetch%proc => unpack_fetch_refine
     pack_flush%proc => pack_flush_refine
     unpack_flush%proc => unpack_flush_refine
  endif
  if(cache_operation.EQ.operation_loadbalance)then
     mdl%combiner_rule = COMBINER_CREATE
     mdl%size_msg_array = storage_size(dummy_large_realdp)/32
     init_flush%proc => null()
     pack_fetch%proc => pack_fetch_refine
     unpack_fetch%proc => unpack_fetch_refine
     pack_flush%proc => pack_flush_loadbalance
     unpack_flush%proc => unpack_flush_loadbalance
  endif

  ! Operations of type "mg"
  if(cache_operation.EQ.operation_mg)then
     mdl%size_msg_array = storage_size(dummy_twin_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_mg
     unpack_fetch%proc => unpack_fetch_mg
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif

  ! Operations of type "interpol"
  if(cache_operation.EQ.operation_interpol)then
     mdl%size_msg_array = storage_size(dummy_three_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_interpol
     unpack_fetch%proc => unpack_fetch_interpol
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif
  if(cache_operation.EQ.operation_kick)then
     mdl%size_msg_array = storage_size(dummy_three_realdp)/32
     init_flush%proc => null()     
     pack_fetch%proc => pack_fetch_kick
     unpack_fetch%proc => unpack_fetch_kick
     pack_flush%proc => null()
     unpack_flush%proc => null()
  endif

  ! Compute useful size of communication buffers
  mdl%size_request_array=1+ndim
  mdl%size_flush_array=1+(1+ndim+mdl%size_msg_array)*nflushmax
  mdl%size_fetch_array=2+(1+ndim+mdl%size_msg_array)*ntilemax

  ! Set communication counters to zero
  do icpu=1,g%ncpu
     iskip=mdl%size_flush_array*(icpu-1)+1
     mdl%send_flush_array(iskip)=0
  end do
  
  ! Post the first RECV for request
  call MPI_IRECV(mdl%recv_request_array,mdl%size_request_array,MPI_INTEGER,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,mdl%request_id,info)
  
  ! Post the first RECV for flush
  call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)

  end associate
  
#endif
#endif
end subroutine open_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################

!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache(child_grid)
  use amr_commons
  implicit none
  integer::child_grid
  integer::icache
  if(child_grid>ngridmax)then
     icache=child_grid-ngridmax
     locked(icache)=.true.
  endif
end subroutine lock_cache
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine unlock_cache(child_grid)
  use amr_commons
  implicit none
  integer::child_grid
  integer::icache
  if(child_grid>ngridmax)then
     icache=child_grid-ngridmax
     locked(icache)=.false.
  endif
end subroutine unlock_cache
!###############################################################
!###############################################################
!###############################################################
!###############################################################
integer function get_parent_cell(hash_key,flush_cache,fetch_cache) result(parent_cell)
  use amr_commons
  use hash
  implicit none
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  !
  integer(kind=8),dimension(0:ndim)::hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer::ind,ipos,idim,get_grid

  hash_father(0)=hash_key(0)-1
  hash_father(1:ndim)=hash_key(1:ndim)/2
  ii(1:ndim)=hash_key(1:ndim)-2*hash_father(1:ndim)
  ind=1
  do idim=1,ndim
     ind=ind+2**(idim-1)*ii(idim)
  end do
  ipos=get_grid(hash_father,flush_cache,fetch_cache)
  parent_cell=0
  if(ipos>0)parent_cell=(ipos-1)*twotondim+ind
end function get_parent_cell
!##############################################################
!##############################################################
!##############################################################
!##############################################################
integer function get_grid(hash_key,flush_cache,fetch_cache) result(child_grid)
  use amr_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  !
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector),save::hk0,hk1,hk2
  integer(kind=8),dimension(1:nvector),save::ix,iy,iz
  integer(kind=8),dimension(0:ndim)::hash_child
  integer::i,ind,ichild,ilevel,info,icpu,grid_cpu,ntile_response,icounter
  integer::send_request_id
  type(request),save::send_request
  integer::response_id
  type(int4_msg),save::response_flag
  type(realdp_msg),save::response_hydro
  logical::failed_request,send_request_completed
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status
#endif
  !
#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mail_counter==32)then
     call check_mail(MPI_REQUEST_NULL)
     mail_counter=0
  endif
  mail_counter=mail_counter+1
#endif

  ! Access hash table
  child_grid=hash_get(grid_dict,hash_key)

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
#if NDIM==1
  ix(1)=hash_key(1)
  call hilbert1d(ix,hk0,1)
#endif
#if NDIM==2
  ix(1)=hash_key(1)
  iy(1)=hash_key(2)
  call hilbert2d(ix,iy,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
#if NDIM==3
  ix(1)=hash_key(1)
  iy(1)=hash_key(2)
  iz(1)=hash_key(3)
  call hilbert3d(ix,iy,iz,hk2,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif

  ! Check if grid sits inside processor boundaries
  if(    hk0(1).ge.bound_key_level(myid-1,ilevel).AND. &
       & hk0(1).lt.bound_key_level(myid  ,ilevel))then
     return
  endif

  ! Determine parent processor
  do icpu=1,ncpu
     if(    hk0(1).ge.bound_key_level(icpu-1,ilevel).AND. &
          & hk0(1).lt.bound_key_level(icpu  ,ilevel))then
        grid_cpu=icpu
        exit
     end if
  end do

  !============================================
  ! We have a fetch and possibly a flush cache
  !============================================
  if(fetch_cache)then

     ! Send a request to the relevant cpu
     send_request%lev=ilevel
     send_request%ckey(1:ndim)=hash_key(1:ndim)
     
     ! Post RECV for the expected response
     if(cache_operation_type.EQ.operation_type_flag)then
        call MPI_IRECV(response_flag,1,new_mpi_int4_msg,grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_hydro)then
        call MPI_IRECV(response_hydro,1,new_mpi_realdp_msg,grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     call MPI_ISEND(send_request,1,new_mpi_request,grid_cpu-1,request_tag,MPI_COMM_WORLD,send_request_id,info)

     ! While waiting for reply, check on incoming messages and perform actions
     call check_mail(response_id)

     ! Test for ISEND completion to free memory in corresponding MPI buffer
!     call MPI_TEST(send_request_id,send_request_completed,send_request_status,info)
     call MPI_WAIT(send_request_id,send_request_status,info)
     
     if(cache_operation_type.EQ.operation_type_flag)then
        failed_request=response_flag%type==-1
        ntile_response=response_flag%ntile
     endif
     if(cache_operation_type.EQ.operation_type_hydro)then
        failed_request=response_hydro%type==-1
        ntile_response=response_hydro%ntile
     endif

     ! If grid does not exist, store -1 in the cache
     ! The output grid index is still zero
     if(failed_request)then

        ! Delete old null grid if occupied
        if(occupied_null(free_null))then
           hash_child(0)=lev_null(free_null)
           hash_child(1:ndim)=ckey_null(1:ndim,free_null)
           call hash_free(grid_dict,hash_child)
        endif
        call hash_set(grid_dict,hash_key,-1)
        occupied_null(free_null)=.true.
        lev_null(free_null)=ilevel
        ckey_null(1:ndim,free_null)=hash_key(1:ndim)

        ! Go to next free cache line
        free_null=free_null+1
        nnull=nnull+1
        if(free_null.GT.ncachemax)then
           free_null=1
!           write(*,*)'cache null full'
        endif
        if(nnull.GT.ncachemax)nnull=ncachemax

     ! If grid exists, store incoming tile in the cache
     else    
        do i=1,ntile_response

           ! If next cache line is occupied, free it.
           if(locked(free_cache))then
              icounter=0
              do while(locked(free_cache))
                 free_cache=free_cache+1
                 icounter=icounter+1
                 if(free_cache>ncachemax)free_cache=1
                 if(icounter>ncachemax)then
                    write(*,*)'PE ',myid,'cache entirely locked'
                    stop
                 endif
              end do
           end if

           ! Create the grid in local memory
           if(cache_operation_type.EQ.operation_type_flag)then
              hash_child(0)=response_flag%lev(i)
              hash_child(1:ndim)=response_flag%ckey(1:ndim,i)
           endif
           if(cache_operation_type.EQ.operation_type_hydro)then
              hash_child(0)=response_hydro%lev(i)
              hash_child(1:ndim)=response_hydro%ckey(1:ndim,i)
           endif
           ichild=ngridmax+free_cache

           if(hash_get(grid_dict,hash_child).EQ.0)then

              if(occupied(free_cache))call destage(ngridmax+free_cache)

              call hash_set(grid_dict,hash_child,ichild)
              
              occupied(free_cache)=.true.
              parent_cpu(free_cache)=grid_cpu
              dirty(free_cache)=.false.
              
              ! Set the grid index of the requested grid
              if(same_keys(hash_key,hash_child))then
                 child_grid=ichild
              endif
              
              ! Store the grid coordinates for the entire tile
              grid(ichild)%lev=hash_child(0)
              grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
              
              ! Depends on the type of cache operations
              if(cache_operation_type.EQ.operation_type_flag)then
                 grid(ichild)%flag1(1:twotondim)=response_flag%int4(1:twotondim,i)
              endif
              if(cache_operation_type.EQ.operation_type_hydro)then
                 do ind=1,twotondim
                    if(response_hydro%int4(ind,i)==1)then
                       grid(ichild)%refined(ind)=.true.
                    else
                       grid(ichild)%refined(ind)=.false.
                    end if
                 end do
                 grid(ichild)%uold(1:twotondim,1:nvar)=response_hydro%realdp(1:twotondim,1:nvar,i)
              endif
              
              ! If we also have a flush cache...
              if(flush_cache)then
                 dirty(free_cache)=.true.
                 
                 !================================================
                 ! Set initialisation rule for combiner operations
                 !================================================
                 
                 ! Initialisation rule for Godunov update
                 if(cache_operation.EQ.operation_godunov)then           
                    grid(ichild)%unew(1:twotondim,1:nvar)=0.0
                 endif
                 
                 ! Initialisation rule for derefine
                 if(cache_operation.EQ.operation_derefine)then           
                    grid(ichild)%refined(1:twotondim)=.true.
                 endif
                 
              endif

              ! Go to next free cache line
              free_cache=free_cache+1
              ncache=ncache+1
              if(free_cache.GT.ncachemax)then
                 free_cache=1
              endif
              if(ncache.GT.ncachemax)ncache=ncachemax
           endif
        end do
     endif

     !=================================
     ! If we have only a flush cache
     !=================================
  else if(flush_cache)then   

     ! If next cache line is occupied, free it.
     if(locked(free_cache))then
        do while(locked(free_cache))
           free_cache=free_cache+1
           if(free_cache>ncachemax)free_cache=1
        end do
     end if
     if(occupied(free_cache))call destage(ngridmax+free_cache)

     ! Set grid index to a virtual grid in local memory
     child_grid=ngridmax+free_cache
     call hash_set(grid_dict,hash_key,child_grid)

     ! Store the grid coordinates
     grid(child_grid)%lev=hash_key(0)
     grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
     occupied(free_cache)=.true.
     parent_cpu(free_cache)=grid_cpu
     dirty(free_cache)=.true.

     !===============================================
     ! Set initialisation rule for combiner operation
     !===============================================

     ! Initialisation rule for initflag
     if(cache_operation.EQ.operation_initflag)then
        grid(child_grid)%flag1(1:twotondim)=0
     endif

     ! Initialisation rule for hydro upload
     if(cache_operation.EQ.operation_upload)then           
        grid(child_grid)%uold(1:twotondim,1:nvar)=0.0
     endif

     ! Go to next free cache line
     free_cache=free_cache+1
     ncache=ncache+1
     if(free_cache.GT.ncachemax)then
!        write(*,*)'cache 2 full'
        free_cache=1
     endif
     if(ncache.GT.ncachemax)ncache=ncachemax
  endif

#endif
end function get_grid
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine check_mail(comm_id)
  use amr_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::comm_id
  !
  integer::i,ind,ivar,info,ipos,igrid,ichild,grid_cpu,ilevel,itile,ntile_reply
  logical::comm_completed,request_received,flush_received=.false.,reply_sent
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status,comm_status
#endif
  integer(kind=4), dimension(1:nvector),save::dummy_state
  integer(kind=8), dimension(1:nvector),save::hk0,hk1,hk2
  integer(kind=8), dimension(1:nvector),save::ix,iy,iz
  integer(kind=8),dimension(0:ndim)::hash_key,hash_child
  !
#ifndef WITHOUTMPI
  comm_completed=.false.
  do while (.not. comm_completed)

     ! USE A NESTED DO WHILE LOOP FOR THE 2 ADD. TESTS

     !===========================
     ! Check for incoming request
     !===========================
     call MPI_Test(request_id,request_received,request_status,info)
     if(request_received)then
        
        ! Assemble a reply and send it back
        ilevel=recv_request%lev
        hash_key(0)=ilevel
        hash_key(1:ndim)=recv_request%ckey(1:ndim)
        igrid=hash_get(grid_dict,hash_key)
        grid_cpu=request_status(MPI_SOURCE)+1

        ! If grid does not exist, send a null reply
        if(igrid.EQ.0)then
           if(cache_operation_type.EQ.operation_type_flag)then
              reply_flag(grid_cpu)%type=-1
           endif
           if(cache_operation_type.EQ.operation_type_hydro)then
              reply_hydro(grid_cpu)%type=-1
           endif

        ! Otherwise, assemble a proper reply with a complete tile
        else
           itile=(igrid-head(ilevel))/ntilemax
           ntile_reply=MIN(tail(ilevel)-itile*ntilemax-head(ilevel)+1,ntilemax)

           ! Store type of reply and number of entries
           if(cache_operation_type.EQ.operation_type_flag)then
              reply_flag(grid_cpu)%type=1
              reply_flag(grid_cpu)%ntile=ntile_reply
           endif
           if(cache_operation_type.EQ.operation_type_hydro)then
              reply_hydro(grid_cpu)%type=1
              reply_hydro(grid_cpu)%ntile=ntile_reply
           endif

           ! Store data, depending on reply type
           do i=1,ntile_reply
              ipos=head(ilevel)+itile*ntilemax+i-1

              ! Reply of type flag
              if(cache_operation_type.EQ.operation_type_flag)then
                 reply_flag(grid_cpu)%lev(i)=grid(ipos)%lev
                 reply_flag(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
                 reply_flag(grid_cpu)%int4(1:twotondim,i)=grid(ipos)%flag1(1:twotondim)
              endif

              ! Reply of type hydro
              if(cache_operation_type.EQ.operation_type_hydro)then
                 reply_hydro(grid_cpu)%lev(i)=grid(ipos)%lev
                 reply_hydro(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
                 do ind=1,twotondim
                    if(grid(ipos)%refined(ind))then
                       reply_hydro(grid_cpu)%int4(ind,i)=1
                    else
                       reply_hydro(grid_cpu)%int4(ind,i)=0
                    endif
                 end do
                 reply_hydro(grid_cpu)%realdp(1:twotondim,1:nvar,i)=grid(ipos)%uold(1:twotondim,1:nvar)
              endif

           end do
        endif

        ! Test the old SEND to free memory in corresponding MPI buffer
!        call MPI_Test(reply_id(grid_cpu),reply_sent,reply_status,info)
        call MPI_WAIT(reply_id(grid_cpu),reply_status,info)

        ! Send back the reply
        if(cache_operation_type.EQ.operation_type_flag)then
           call MPI_ISEND(reply_flag(grid_cpu),1,new_mpi_int4_msg,grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
        endif
        if(cache_operation_type.EQ.operation_type_hydro)then
           call MPI_ISEND(reply_hydro(grid_cpu),1,new_mpi_realdp_msg,grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
        endif

        ! Post a new RECV for request
        call MPI_IRECV(recv_request,1,new_mpi_request,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,request_id,info)
     endif

     !=========================
     ! Check for incoming flush
     !=========================
     call MPI_Test(flush_id,flush_received,flush_status,info)
     if(flush_received)then

        !===========================================================
        ! Combine received data to local memory using combiner rules
        !===========================================================

        ! Combiner rules for initflag
        if(cache_operation.EQ.operation_initflag)then
           do i=1,recv_flush_flag%nflush
              hash_child(0)=recv_flush_flag%lev(i)
              hash_child(1:ndim)=recv_flush_flag%ckey(1:ndim,i)
              ichild=hash_get(grid_dict,hash_child)
              do ind=1,twotondim
                 grid(ichild)%flag1(ind)=MAX(grid(ichild)%flag1(ind),recv_flush_flag%int4(ind,i))
              end do
           end do
        endif

        ! Combiner rules for derefine
        if(cache_operation.EQ.operation_derefine)then
           do i=1,recv_flush_flag%nflush
              hash_child(0)=recv_flush_flag%lev(i)
              hash_child(1:ndim)=recv_flush_flag%ckey(1:ndim,i)
              ichild=hash_get(grid_dict,hash_child)
              if(ichild>0)then ! Since we are in the process of derefining,
                 do ind=1,twotondim ! we need to check if the grid is still here.
                    if(grid(ichild)%refined(ind))then
                       if(recv_flush_flag%int4(ind,i).EQ.0)then
                          grid(ichild)%refined(ind)=.false.
                       endif
                    endif
                 end do
              endif
           end do
        endif

        ! Combiner rules for hydro upload
        if(cache_operation.EQ.operation_upload)then
           do i=1,recv_flush_hydro%nflush
              hash_child(0)=recv_flush_hydro%lev(i)
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)
              ichild=hash_get(grid_dict,hash_child)
              do ivar=1,nvar
                 do ind=1,twotondim
                    if(grid(ichild)%refined(ind))then
                       grid(ichild)%uold(ind,ivar)=grid(ichild)%uold(ind,ivar)&
                            & +recv_flush_hydro%realdp(ind,ivar,i)
                    endif
                 end do
              end do
           end do
        endif

        ! Combiner rules for godunov update
        if(cache_operation.EQ.operation_godunov)then
           do i=1,recv_flush_hydro%nflush
              hash_child(0)=recv_flush_hydro%lev(i)
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)
              ichild=hash_get(grid_dict,hash_child)
              do ivar=1,nvar
                 do ind=1,twotondim
                    grid(ichild)%unew(ind,ivar)=grid(ichild)%unew(ind,ivar)+recv_flush_hydro%realdp(ind,ivar,i)
                 end do
              end do
           end do
        endif

        ! Combiner rules for refinements
        if(cache_operation.EQ.operation_refine)then

           do i=1,recv_flush_hydro%nflush
              ilevel=recv_flush_hydro%lev(i)
              hash_child(0)=ilevel
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)

              ! Compute Hilbert keys of new octs
#if NDIM==1
              ix(1)=hash_child(1)
              call hilbert1d(ix,hk0,1)
#endif
#if NDIM==2
              ix(1)=hash_child(1)
              iy(1)=hash_child(2)
              call hilbert2d(ix,iy,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
#if NDIM==3
              ix(1)=hash_child(1)
              iy(1)=hash_child(2)
              iz(1)=hash_child(3)
              call hilbert3d(ix,iy,iz,hk2,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
              ! Set grid index to a virtual grid in local main memory
              ichild=ifree

              ! Go to next main memory free line
              ifree=ifree+1
              if(ifree.GT.ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'Increase ngridmax'
                 call clean_abort
              endif

              grid(ichild)%lev=hash_child(0)
              grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
              grid(ichild)%hkey=hk0(1)
              grid(ichild)%refined(1:twotondim)=.false.
              grid(ichild)%flag1(1:twotondim)=0
              grid(ichild)%flag2(1:twotondim)=0
              grid(ichild)%superoct=1
              
              ! Insert new grid in hash table
              call hash_set(grid_dict,hash_child,ichild)

              ! Flush hydro variables
              do ivar=1,nvar
                 do ind=1,twotondim
                    grid(ichild)%uold(ind,ivar)=recv_flush_hydro%realdp(ind,ivar,i)
                 end do
              end do
           end do
        endif

        ! Combiner rules for refinements
        if(cache_operation.EQ.operation_loadbalance)then

           do i=1,recv_flush_hydro%nflush
              ilevel=recv_flush_hydro%lev(i)
              hash_child(0)=ilevel
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)

              ! Compute Hilbert keys of new octs
#if NDIM==1
              ix(1)=hash_child(1)
              call hilbert1d(ix,hk0,1)
#endif
#if NDIM==2
              ix(1)=hash_child(1)
              iy(1)=hash_child(2)
              call hilbert2d(ix,iy,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
#if NDIM==3
              ix(1)=hash_child(1)
              iy(1)=hash_child(2)
              iz(1)=hash_child(3)
              call hilbert3d(ix,iy,iz,hk2,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
              ! Set grid index to a virtual grid in local main memory
              ichild=ifree

              ! Go to next main memory free line
              ifree=ifree+1
              if(ifree.GT.ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'Increase ngridmax'
                 call clean_abort
              endif

              grid(ichild)%lev=hash_child(0)
              grid(ichild)%ckey(1:ndim)=hash_child(1:ndim)
              grid(ichild)%hkey=hk0(1)
              do ind=1,twotondim
                 if(recv_flush_hydro%int4(ind,i)==1)then
                    grid(ichild)%refined(ind)=.true.
                 else
                    grid(ichild)%refined(ind)=.false.
                 endif
              end do
              grid(ichild)%flag1(1:twotondim)=0
              grid(ichild)%flag2(1:twotondim)=0
              grid(ichild)%superoct=1
              
              ! Insert new grid in hash table
              call hash_set(grid_dict,hash_child,ichild)

              ! Flush hydro variables
              do ivar=1,nvar
                 do ind=1,twotondim
                    grid(ichild)%uold(ind,ivar)=recv_flush_hydro%realdp(ind,ivar,i)
                 end do
              end do
           end do
        endif

        ! Post a new RECV for flush
        if(cache_operation_type.EQ.operation_type_flag)then
           call MPI_IRECV(recv_flush_flag,1,new_mpi_int4_flush,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_hydro)then
           call MPI_IRECV(recv_flush_hydro,1,new_mpi_realdp_flush,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
     endif

     !=================================
     ! Check for input comm. completion
     !=================================
     if(comm_id==MPI_REQUEST_NULL)then
        comm_completed=.true.
     else
        call MPI_Test(comm_id,comm_completed,comm_status,info)
     endif
  end do
#endif
end subroutine check_mail
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine destage(igrid)
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid
  !
  integer::ind,ipos,info,icache,iflush,grid_cpu
  integer::send_flush_id
  integer(kind=8),dimension(0:ndim)::hash_key
  !
#ifndef WITHOUTMPI
  hash_key(0)=grid(igrid)%lev
  hash_key(1:ndim)=grid(igrid)%ckey(1:ndim)
  ipos=hash_get(grid_dict,hash_key)

  if(hash_get(grid_dict,hash_key).EQ.0)then
     write(*,*)'PE ',myid,' trying to free non existing grid'
     stop
  endif

  call hash_free(grid_dict,hash_key)

  ! Check if the destage requires a flush
  icache=igrid-ngridmax

  if(dirty(icache))then
     grid_cpu=parent_cpu(icache)
     dirty(icache)=.false.
  
     if(cache_operation.EQ.operation_initflag)then
        send_flush_flag(grid_cpu)%nflush=send_flush_flag(grid_cpu)%nflush+1
        if(send_flush_flag(grid_cpu)%nflush>nflushmax)then
           send_flush_flag(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_flag(grid_cpu),1,new_mpi_int4_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_flag(grid_cpu)%nflush=1
        endif
        iflush=send_flush_flag(grid_cpu)%nflush
        send_flush_flag(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_flag(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_flag(grid_cpu)%int4(1:twotondim,iflush)=grid(igrid)%flag1(1:twotondim)
     endif
     
     if(cache_operation.EQ.operation_derefine)then
        send_flush_flag(grid_cpu)%nflush=send_flush_flag(grid_cpu)%nflush+1
        if(send_flush_flag(grid_cpu)%nflush>nflushmax)then
           send_flush_flag(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_flag(grid_cpu),1,new_mpi_int4_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_flag(grid_cpu)%nflush=1
        endif
        iflush=send_flush_flag(grid_cpu)%nflush
        send_flush_flag(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_flag(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        do ind=1,twotondim
           if(grid(igrid)%refined(ind))then
              send_flush_flag(grid_cpu)%int4(ind,iflush)=1.0
           else
              send_flush_flag(grid_cpu)%int4(ind,iflush)=0.0
           endif
        end do
     endif
     
     if(cache_operation.EQ.operation_upload)then
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        if(send_flush_hydro(grid_cpu)%nflush>nflushmax)then
           send_flush_hydro(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_hydro(grid_cpu)%nflush=1
        endif
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_hydro(grid_cpu)%realdp(1:twotondim,1:nvar,iflush)=grid(igrid)%uold(1:twotondim,1:nvar)
     endif

     if(cache_operation.EQ.operation_godunov)then
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        if(send_flush_hydro(grid_cpu)%nflush>nflushmax)then
           send_flush_hydro(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_hydro(grid_cpu)%nflush=1
        endif
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_hydro(grid_cpu)%realdp(1:twotondim,1:nvar,iflush)=grid(igrid)%unew(1:twotondim,1:nvar)
     endif

     if(cache_operation.EQ.operation_refine)then
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        if(send_flush_hydro(grid_cpu)%nflush>nflushmax)then
           send_flush_hydro(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_hydro(grid_cpu)%nflush=1
        endif
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_hydro(grid_cpu)%realdp(1:twotondim,1:nvar,iflush)=grid(igrid)%uold(1:twotondim,1:nvar)
     endif

     if(cache_operation.EQ.operation_loadbalance)then
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        if(send_flush_hydro(grid_cpu)%nflush>nflushmax)then
           send_flush_hydro(grid_cpu)%nflush=nflushmax
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_hydro(grid_cpu)%nflush=1
        endif
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_hydro(grid_cpu)%realdp(1:twotondim,1:nvar,iflush)=grid(igrid)%uold(1:twotondim,1:nvar)
        do ind=1,twotondim
           if(grid(igrid)%refined(ind))then
              send_flush_hydro(grid_cpu)%int4(ind,iflush)=1.0
           else
              send_flush_hydro(grid_cpu)%int4(ind,iflush)=0.0
           endif
        end do
     endif

  endif
#endif
end subroutine destage
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !
  integer::info,icache,igrid,icpu
  integer::send_flush_id,ndebug
  integer::dummy_int,close_tag=7,close_id
  integer(kind=8),dimension(0:ndim)::hash_child
  logical::request_received,flush_received
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status
#endif
  !
#ifndef WITHOUTMPI

  ! EMPTY AND CLEAN THE CACHE
  do icache=1,ncache
     igrid=ngridmax+icache
     locked(icache)=.false.
     call destage(igrid)
     occupied(icache)=.false.
     dirty(icache)=.false.
  end do
  free_cache=1
  ncache=0

  do icache=1,nnull
     if(occupied_null(icache))then
        hash_child(0)=lev_null(icache)
        hash_child(1:ndim)=ckey_null(1:ndim,icache)
        call hash_free(grid_dict,hash_child)
     endif
     occupied_null(icache)=.false.
  end do
  free_null=1
  nnull=0

  ! COMPLETE THE LAST FLUSH
  do icpu=1,ncpu
     if(cache_operation_type.EQ.operation_type_flag)then
        if(send_flush_flag(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_flag(icpu),1,new_mpi_int4_flush,icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_flag(icpu)%nflush=0
        endif
     endif     
     if(cache_operation_type.EQ.operation_type_hydro)then

        if(send_flush_hydro(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_hydro(icpu),1,new_mpi_realdp_flush,icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id)
           send_flush_hydro(icpu)%nflush=0
        endif
     endif
  end do
  
    ! CHECK-IN CHECK-OUT
  if(myid.NE.1)then
     call MPI_ISEND(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(close_id)
     call MPI_IRECV(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(close_id)
  else
     do icpu=2,ncpu
        call MPI_IRECV(dummy_int,1,MPI_INTEGER,MPI_ANY_SOURCE,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(close_id)
     end do
     do icpu=2,ncpu
        call MPI_ISEND(dummy_int,1,MPI_INTEGER,icpu-1,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(close_id)
     end do
  endif

  ! Barrier to get the last flush message
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  call check_mail(MPI_REQUEST_NULL)

  ! Finally CANCEL THE 2 RECV
  call MPI_CANCEL(request_id,info)
  call MPI_CANCEL(flush_id,info)

  ! Test to free memory in corresponding MPI buffer
  call MPI_Wait(request_id,request_status,info)
  call MPI_Wait(flush_id,flush_status,info)
  do icpu=1,ncpu
     call MPI_WAIT(reply_id(icpu),reply_status,info)
  end do
  
  ! Barrier to prevent interference with the next cache
  call MPI_BARRIER(MPI_COMM_WORLD,info)

#endif
end subroutine close_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine open_cache
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !
  integer::info,icpu
  !
#ifndef WITHOUTMPI

  do icpu=1,ncpu
     reply_id(icpu)=MPI_REQUEST_NULL
  end do
  mail_counter=0

  do icpu=1,ncpu
     send_flush_flag(icpu)%nflush=0
     send_flush_hydro(icpu)%nflush=0
  end do

  cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_initflag)cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_smooth  )cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_derefine)cache_operation_type=operation_type_flag

  if(cache_operation.EQ.operation_hydro   )cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_upload  )cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_godunov )cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_refine  )cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_loadbalance)cache_operation_type=operation_type_hydro

  ! Post the first RECV for request
  call MPI_IRECV(recv_request,1,new_mpi_request,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,request_id,info)
  
  ! Post the first RECV for flush
  if(cache_operation_type.EQ.operation_type_flag)then
     call MPI_IRECV(recv_flush_flag,1,new_mpi_int4_flush,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_hydro)then
     call MPI_IRECV(recv_flush_hydro,1,new_mpi_realdp_flush,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif

#endif

end subroutine open_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################

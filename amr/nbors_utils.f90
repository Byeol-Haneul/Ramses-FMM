module nbors_utils
contains
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_threetondim_nbor_parent_cell(s,hash_key,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(1:threetondim)::igrid_nbor
  integer,dimension(1:threetondim)::ind_nbor
  logical::flush_cache,fetch_cache
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
  integer,dimension(1:3,1:8)::start_oct=reshape(&
       & (/ 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1,&
       &    1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0/),(/3,8/))
  integer::i1,j1,k1
  integer::i1min=-1,i1max=+1
  integer::j1min=-1*(ndim/2)
  integer::j1max=+1*(ndim/2)
  integer::k1min=-1*(ndim/3)
  integer::k1max=+1*(ndim/3)
  integer::ind,ipos,idim,ilevel,inbor,igrid

  associate(r=>s%r,m=>s%mdl%m)

  ilevel=hash_key(0)
  hash_father(0)=hash_key(0)-1
  hash_father(1:ndim)=hash_key(1:ndim)/2

  ii(1:ndim)=hash_key(1:ndim)-2*hash_father(1:ndim)
  ind=1
  do idim=1,ndim
     ind=ind+2**(idim-1)*ii(idim)
  end do

  hash_nbor(0)=hash_father(0)
#if NDIM>2
  do k1=0,1
     hash_nbor(3)=hash_father(3)+k1*shift_oct(3,ind)
     ii(3)=start_oct(3,ind)+k1*shift_oct(3,ind)
#endif
#if NDIM>1
     do j1=0,1
        hash_nbor(2)=hash_father(2)+j1*shift_oct(2,ind)
        ii(2)=start_oct(2,ind)+j1*shift_oct(2,ind)
#endif
#if NDIM>0
        do i1=0,1
           hash_nbor(1)=hash_father(1)+i1*shift_oct(1,ind)
           ii(1)=start_oct(1,ind)+i1*shift_oct(1,ind)
#endif
           inbor=1
           do idim=1,ndim
              inbor=inbor+2**(idim-1)*ii(idim)
           end do
           ! Periodic boundary conditions
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(hash_nbor(idim) <m%box_ckey_min(idim,ilevel-1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel-1)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel-1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel-1)
              endif
           enddo
           ! Store lower left neighbor coordinates
           if(inbor==1)hash_ref(1:ndim)=hash_nbor(1:ndim)
           ! Get grid into memory and lock it if remote
           call get_grid(s,hash_nbor,igrid,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
           igrid_twotondim_nbor(inbor)=igrid
#if NDIM>0
        end do
#endif
#if NDIM>1
     end do
#endif
#if NDIM>2
  end do
#endif

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
           ! Periodic boundary conditions
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(hash_nbor(idim) <m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
              endif
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
           ! Periodic boundary conditions
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(ii(idim)<0)ii(idim)=ii(idim)+m%box_ckey_max(idim,ilevel-1)-m%box_ckey_min(idim,ilevel-1)
              endif
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
subroutine get_twondim_nbor_parent_cell(s,hash_key,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim, twondim
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(0:ndim)::hash_key
  integer,dimension(0:twondim)::igrid_nbor
  integer,dimension(0:twondim)::ind_nbor
  logical::flush_cache,fetch_cache
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
  integer::ind,ipos,idim,ilevel,inbor,igrid

  associate(r=>s%r,m=>s%mdl%m)

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
  call get_grid(s,hash_father,igrid,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
  igrid_nbor(0)=igrid
  ind_nbor(0)=ind
  
  ! Deal with neighboring father cells
  do inbor=1,twondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift(1:ndim,inbor)
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
           if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
        endif
     enddo
     hash_father(1:ndim)=hash_nbor(1:ndim)/2
     ii(1:ndim)=hash_nbor(1:ndim)-2*hash_father(1:ndim)
     ind=1
     do idim=1,ndim
        ind=ind+2**(idim-1)*ii(idim)
     end do

     ! Get grid into memory and lock it if remote 
     call get_grid(s,hash_father,igrid,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
     igrid_nbor(inbor)=igrid
     ind_nbor(inbor)=ind
  end do

  end associate
  
end subroutine get_twondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_parent_cell(s,hash_key,igrid,ind,flush_cache,fetch_cache,lock)
  use amr_parameters, only: ndim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::igrid,ind
  logical::flush_cache,fetch_cache
  logical,optional::lock
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
  call get_grid(s,hash_father,igrid,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=lock)

end subroutine get_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache(m,child_grid)
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer::child_grid
  !
  ! This routine locks a cache line because
  ! it will be updated later.
  !
  integer::icache
  if(child_grid>m%ngridmax)then
     icache=child_grid-m%ngridmax
     if(.not.m%locked(icache))then
        m%nlocked=m%nlocked+1
        m%nlocked_max=max(m%nlocked,m%nlocked_max)
     endif
     m%locked(icache)=.true.
  endif
end subroutine lock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine unlock_cache(m,child_grid)
  use amr_commons, only: mesh_t
  use mdl_module
  implicit none
  type(mesh_t)::m
  integer::child_grid
  !
  ! This routine unlocks a cache line because
  ! it has been updated and can be flushed.
  !
  integer::icache
  if(child_grid>m%ngridmax)then
     icache=child_grid-m%ngridmax
     if(m%locked(icache))m%nlocked=m%nlocked-1
     m%locked(icache)=.false.
  endif
end subroutine unlock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_grid(s,hash_key,child,flush_cache,fetch_cache,lock,use_ghost)
  use mdl_module
  use amr_parameters, only: ndim, nhilbert
  use amr_commons, only: oct
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use boundaries, only: get_bound
  use cache_commons
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(0:ndim)::hash_key
  integer::child
  logical::flush_cache,fetch_cache
  logical,optional::lock
  logical,optional::use_ghost
  !
  ! This routine acquires the grid 
  ! corresponding to the input hash key.
  ! Two logicals set the type of acquire.
  ! "fetch" means read-only, "flush" means write-only.
  ! and both used together means read-write.
  !
  integer(kind=4)::full_hash
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_child, hash_ref
  integer::i,ind,idim,ivar,iskip,ilevel,info,ibound
  integer::grid_cpu,ntile_response,icounter,child_grid
  integer::send_request_id,response_id  
  logical::failed_request,in_rank,in_domain,do_ghost
  integer::null_val=-1
  integer::child_ref
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status
#endif

  associate(r=>s%r,g=>s%g,m=>s%mdl%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(mdl,MPI_REQUEST_NULL)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Access hash table
  child=hash_getp(m%grid_dict,hash_key)

  if(present(lock))then
     if(lock)call lock_cache(m,child)
  endif
  if(child==-1)then
     child=0
     return
  endif
  if(child>0)return

  ! Now we know child_grid=0
  child_grid = 0

  ! Check if we do null grid or ghost grid
  do_ghost=.false.
  if(present(use_ghost))then
     if(use_ghost)do_ghost=.true.
  endif

  ! Compute the Hilbert key
  ilevel=hash_key(0)
  ix(1:ndim)=hash_key(1:ndim)

  ! Check if in computational domain
  in_domain = .true.
  do idim = 1, ndim
     in_domain = in_domain .and. ix(idim) .ge. m%box_ckey_min(idim,ilevel) .and. ix(idim) .lt. m%box_ckey_max(idim,ilevel)
  end do

  if(in_domain)then

#ifdef WITHOUTMPI
     child=0
     return
#endif

#ifndef WITHOUTMPI
     ! Compute Hilbert key
     hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

     ! Check if grid sits inside processor boundaries
     in_rank = ge_keys(hk,m%domain(ilevel)%b(1:nhilbert,mdl_self(mdl)-1)).and. &
          &    gt_keys(m%domain(ilevel)%b(1:nhilbert,mdl_self(mdl)),hk)
     if (in_rank)then ! The grid does not exist in the local domain
        child=0
        return
     endif

     ! Determine parent processor
     grid_cpu = m%domain(ilevel)%get_rank(hk)

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
        call check_mail(mdl,response_id)

        ! Wait for ISEND completion to free memory in corresponding MPI buffer
        call MPI_WAIT(send_request_id,send_request_status,info)

        ! Check header for type of response
        iskip=1
        failed_request=(mdl%recv_fetch_array(iskip)==-1)
        iskip=iskip+1

        ! If grid does not exist and we don't do ghost grid, store null in the cache
        ! The output grid index is zero
        if(failed_request)then

           if(.not. do_ghost)then
              ! Delete old null grid if occupied
              if(m%occupied_null(m%free_null))then
                 hash_child(0)=m%lev_null(m%free_null)
                 hash_child(1:ndim)=m%ckey_null(1:ndim,m%free_null)
                 call hash_free(m%grid_dict,hash_child)
              endif
              call hash_setp(m%grid_dict,hash_key,null_val)
              m%occupied_null(m%free_null)=.true.
              m%lev_null(m%free_null)=ilevel
              m%ckey_null(1:ndim,m%free_null)=hash_key(1:ndim)
              ! Go to next free cache line
              m%free_null=m%free_null+1
              m%nnull=m%nnull+1
              if(m%free_null.GT.m%ncachemax)then
                 m%free_null=1
              endif
              if(m%nnull.GT.m%ncachemax)m%nnull=m%ncachemax
           endif

           child_grid = 0

        ! If grid exists, store incoming tile in the cache
        else

           ! Number of tiles in the response buffer
           ntile_response=mdl%recv_fetch_array(iskip)
           iskip=iskip+1

           ! Loop over tiles
           do i=1,ntile_response

              ! Find next available cache line
              if(m%locked(m%free_cache))then
                 icounter=0
                 do while(m%locked(m%free_cache))
                    m%free_cache=m%free_cache+1
                    icounter=icounter+1
                    if(m%free_cache>m%ncachemax)m%free_cache=1
                    if(icounter>m%ncachemax)then
                       write(*,*)'PE ',mdl_self(mdl),'cache entirely locked'
                       stop
                    endif
                 end do
              end if

              ! Next available grid in memory
              child=m%ngridmax+m%free_cache

              ! Get grid coordinates from message header
              hash_child(0)=mdl%recv_fetch_array(iskip)
              hash_child(1:ndim)=mdl%recv_fetch_array(iskip+1:iskip+ndim)
              iskip=iskip+ndim+1

              ! If grid does not already exist, create it in local memory
              if(hash_getp(m%grid_dict,hash_child)<=0)then

                 ! If next cache line is occupied, free it.
                 if(m%occupied(m%free_cache))call destage(mdl,m%ngridmax+m%free_cache)

                 call hash_setp(m%grid_dict,hash_child,child)

                 m%occupied(m%free_cache)=.true.
                 m%parent_cpu(m%free_cache)=grid_cpu
                 m%dirty(m%free_cache)=.false.
                 m%ghost_parent_grid(m%free_cache)=0
                 m%ghost_parent_cell(m%free_cache)=0

                 ! Set the grid index of the requested grid
                 if(same_keys(hash_key,hash_child))then
                    child_grid=child
                 endif

                 ! Unpack response to fetch request
                 call unpack_fetch%proc(m,child,mdl%size_msg_array,mdl%recv_fetch_array(iskip:iskip+mdl%size_msg_array-1),hash_child)

                 ! If we also have also a flush cache...
                 ! This is for combined read-write cache operations
                 if(flush_cache)then
                    m%dirty(m%free_cache)=.true.

                    ! Set initialisation rule for combiner operations
                    call init_flush%proc(m,child,hash_child)

                 endif

                 ! Go to next free cache line
                 m%free_cache=m%free_cache+1
                 m%ncache=m%ncache+1
                 if(m%free_cache.GT.m%ncachemax)then
                    m%free_cache=1
                 endif
                 if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

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
              if(m%free_cache>m%ncachemax)m%free_cache=1
           end do
        end if
        if(m%occupied(m%free_cache))call destage(mdl,m%ngridmax+m%free_cache)

        ! Set grid index to a virtual grid in local memory
        child_grid=m%ngridmax+m%free_cache
        call hash_setp(m%grid_dict,hash_key,child_grid)

        ! Store the grid coordinates
        m%grid(child_grid)%lev=hash_key(0)
        m%grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
        m%occupied(m%free_cache)=.true.
        m%parent_cpu(m%free_cache)=grid_cpu
        m%dirty(m%free_cache)=.true.
        m%ghost_parent_grid(m%free_cache)=0
        m%ghost_parent_cell(m%free_cache)=0

        ! Set initialisation rule for combiner operation
        call init_flush%proc(m,child_grid,hash_key)

        ! Go to next free cache line
        m%free_cache=m%free_cache+1
        m%ncache=m%ncache+1
        if(m%free_cache.GT.m%ncachemax)then
           m%free_cache=1
        endif
        if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

     end if
#endif

  ! Cartesian key is outside conputational domain
  ! We need to create a boundary condition grid in the cache
  else

     call get_bound(s,ix,ilevel,ibound)

     hash_ref(0)=ilevel
     hash_ref(1:ndim)=ix(1:ndim)
     hash_ref(r%bound_dir(ibound))=hash_ref(r%bound_dir(ibound))+r%bound_shift(ibound)

     child_ref=hash_getp(m%grid_dict,hash_ref)

     if (child_ref>0)then

        ! If next cache line is occupied, free it.
        if(m%locked(m%free_cache))then
           do while(m%locked(m%free_cache))
              m%free_cache=m%free_cache+1
              if(m%free_cache>m%ncachemax)m%free_cache=1
           end do
        end if
        if(m%occupied(m%free_cache))call destage(mdl,m%ngridmax+m%free_cache)

        ! Set grid index to a virtual grid in local memory
        child_grid=m%ngridmax+m%free_cache
        call hash_setp(m%grid_dict,hash_key,child_grid)

        ! Store the grid coordinates
        m%grid(child_grid)%lev=hash_key(0)
        m%grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
        m%occupied(m%free_cache)=.true.
        m%parent_cpu(m%free_cache)=0
        m%dirty(m%free_cache)=.false.
        m%ghost_parent_grid(m%free_cache)=0
        m%ghost_parent_cell(m%free_cache)=0

        ! Set initialisation rule for boundary grid using reference grid
        if(associated(init_bound%proc))then
           call init_bound%proc(r,g,m,child_grid,child_ref,ibound)
        endif

        ! Go to next free cache line
        m%free_cache=m%free_cache+1
        m%ncache=m%ncache+1
        if(m%free_cache.GT.m%ncachemax)then
           m%free_cache=1
        endif
        if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

     else

        if(.not. do_ghost)then
           ! Delete old null grid if occupied
           if(m%occupied_null(m%free_null))then
              hash_child(0)=m%lev_null(m%free_null)
              hash_child(1:ndim)=m%ckey_null(1:ndim,m%free_null)
              call hash_free(m%grid_dict,hash_child)
           endif
           call hash_setp(m%grid_dict,hash_key,null_val)
           m%occupied_null(m%free_null)=.true.
           m%lev_null(m%free_null)=ilevel
           m%ckey_null(1:ndim,m%free_null)=hash_key(1:ndim)
           ! Go to next free cache line
           m%free_null=m%free_null+1
           m%nnull=m%nnull+1
           if(m%free_null.GT.m%ncachemax)then
              m%free_null=1
           endif
           if(m%nnull.GT.m%ncachemax)m%nnull=m%ncachemax
        endif

        child_grid = 0

     endif

  endif

  if(child_grid>0)then
     child=child_grid
     if(present(lock))then
        if(lock)call lock_cache(m,child)
     endif
  else
     child=0
  endif

  end associate
end subroutine get_grid
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module nbors_utils

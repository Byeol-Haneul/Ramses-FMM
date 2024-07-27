module nbors_utils
contains
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_threetondim_nbor_parent_cell(s,hash_key,hash_dict,grid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,twotondim,threetondim
  use amr_commons, only: nbor, oct
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(1:threetondim)::ind_nbor
  type(nbor),dimension(1:threetondim)::grid_nbor
  !
  ! This routine computes and acquire the 3**ndim neighboring father cells 
  ! for the input hash_key. The output arrays are the father cells
  ! parent oct indices and their associated cell indices within the oct.
  ! The corresponding data can be accessed using: grid(igrid)%data(ind).
  ! If the grid index is zero, it means that this oct does not exist.
  ! Note that the 2**ndim grids are all locked if remote.
  !
  type(nbor),dimension(1:twotondim)::grid_twotondim_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor,hash_ref,hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer,dimension(1:3,1:8)::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer,dimension(1:3,1:8)::start_oct=reshape(&
       & (/ 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1,&
       &    1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0/),(/3,8/))
  integer::i1,j1,k1
  integer::i1min=-1
  integer::i1max=+1
  integer::j1min=-1*(ndim/2)
  integer::j1max=+1*(ndim/2)
  integer::k1min=-1*(ndim/3)
  integer::k1max=+1*(ndim/3)
  integer::ind,ipos,idim,ilevel,inbor
  type(oct),pointer::gridp

  associate(r=>s%r,g=>s%g,m=>s%m)

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
             if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel-1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel-1)-1
             if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel-1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel-1)
          endif
       enddo
       ! Store lower left neighbor coordinates
       if(inbor==1)hash_ref(1:ndim)=hash_nbor(1:ndim)
       ! Get grid into memory and lock it if remote
       call get_grid(s,hash_nbor,hash_dict,gridp,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
       grid_twotondim_nbor(inbor)%p=>gridp
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
                   if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
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
                   if(ii(idim)<m%box_ckey_min(idim,ilevel-1))ii(idim)=ii(idim)+m%box_ckey_max(idim,ilevel-1)
                endif
             enddo
             ind=1
             do idim=1,ndim
                ind=ind+2**(idim-1)*ii(idim)
             end do
             grid_nbor(inbor)%p=>grid_twotondim_nbor(ind)%p
          end do
       end do
    end do

  end associate

end subroutine get_threetondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_twondim_nbor_parent_cell(s,hash_key,hash_dict,grid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_parameters, only: ndim,twotondim,twondim
  use amr_commons, only: oct, nbor
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(0:twondim)::ind_nbor
  type(nbor),dimension(0:twondim)::grid_nbor
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
  type(oct),pointer::gridp

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
  call get_grid(s,hash_father,hash_dict,gridp,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
  grid_nbor(0)%p=>gridp
  ind_nbor(0)=ind
  
  ! Deal with neighboring father cells
  do inbor=1,twondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift(1:ndim,inbor)
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
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
     call get_grid(s,hash_father,hash_dict,gridp,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=.true.)
     grid_nbor(inbor)%p=>gridp
     ind_nbor(inbor)=ind
  end do

  end associate
  
end subroutine get_twondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_parent_cell(s,hash_key,hash_dict,gridp,ind,flush_cache,fetch_cache,lock)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use hash
  implicit none
  type(ramses_t)::s
  logical::flush_cache,fetch_cache
  logical,optional::lock
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer::ind
  type(oct),pointer::gridp
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
  call get_grid(s,hash_father,hash_dict,gridp,flush_cache=flush_cache,fetch_cache=fetch_cache,lock=lock)

end subroutine get_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine lock_cache(s,child)
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  type(oct),pointer::child
  !
  ! This routine locks a cache line because
  ! it will be updated later.
  !
  integer::child_grid
  integer::icache
  if(.not.associated(child))return
  child_grid=(loc(child)-loc(s%m%grid(1)))/(loc(s%m%grid(2))-loc(s%m%grid(1)))+1
  if(child_grid>s%r%ngridmax)then
     icache=child_grid-s%r%ngridmax
     s%m%locked(icache)=.true.
  endif
end subroutine lock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine unlock_cache(s,child)
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use mdl_module
  implicit none
  type(ramses_t)::s
  type(oct),pointer::child
#ifdef MDL2
  call mdl_release(s%mdl%mdl2,0,child)
#else
  !
  ! This routine unlocks a cache line because
  ! it has been updated and can be flushed.
  !
  integer::child_grid
  integer::icache
  if(.not.associated(child))return
  child_grid=(loc(child)-loc(s%m%grid(1)))/(loc(s%m%grid(2))-loc(s%m%grid(1)))+1
  if(child_grid>s%r%ngridmax)then
     icache=child_grid-s%r%ngridmax
     s%m%locked(icache)=.false.
  endif
#endif
end subroutine unlock_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine get_grid(s,hash_key,hash_dict,child,flush_cache,fetch_cache,lock)
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_F_POINTER, C_ASSOCIATED, C_BOOL
  use mdl_module
  use amr_parameters, only: ndim,nhilbert,twotondim
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
  type(oct),pointer::child
  logical::flush_cache,fetch_cache
  logical,optional::lock
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  !
  ! This routine acquires the grid 
  ! corresponding to the input hash key.
  ! Two logicals set the type of acquire.
  ! "fetch" means read-only, "flush" means write-only.
  ! and both used together means read-write.
  !
  logical::absent
  logical(c_bool)::do_lock,modify,virtual
  integer(kind=4)::full_hash
  integer(kind=8),dimension(1:nhilbert)::hk
  integer(kind=8),dimension(1:ndim)::ix
  integer(kind=8),dimension(0:ndim)::hash_child, hash_ref
  integer::i,ind,idim,ivar,iskip,ichild,ilevel,info,ibound
  integer::grid_cpu,ntile_response,icounter,child_grid
  integer::send_request_id,response_id  
  logical::failed_request,in_rank,in_domain
  type(oct),pointer::child_ref
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status
#endif

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifdef MDL2
  full_hash = hash_func(hash_key)
  do_lock = .false.
  modify = flush_cache
  virtual = .not.fetch_cache
  if (present(lock)) then
    if (lock) do_lock=.true.
  endif
  call c_f_pointer(cache_fetch(mdl%mdl2,0,hash=full_hash,key=hash_key,lock=do_lock,modify=modify,virtual=virtual),child)
#else /* Not MDL2 */

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(s,MPI_REQUEST_NULL,hash_dict)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Access hash table
  call c_f_pointer(hash_getp(hash_dict,hash_key,absent),child)
  if (present(lock).and.associated(child)) then
     if (lock) call lock_cache(s,child)
  endif
  if (associated(child).or.absent) return

  ! Now we know child_grid=0
  child_grid = 0

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
     nullify(child)
     return
#endif

#ifndef WITHOUTMPI
     ! Conpute Hilbert key
     hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

     ! Check if grid sits inside processor boundaries
     !  if (m%domain_hilbert(ilevel)%in_rank(hk)) return
     in_rank = ge_keys(hk,m%domain_hilbert(ilevel)%b(1:nhilbert,mdl_self(mdl)-1)).and. &
          &    gt_keys(m%domain_hilbert(ilevel)%b(1:nhilbert,mdl_self(mdl)),hk)
     if (in_rank)then ! The grid does not exist in the local domain
        nullify(child)
        return
     endif

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
           call hash_setp(hash_dict,hash_key)
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
           child_grid = 0

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
              if(.not.C_ASSOCIATED(hash_getp(hash_dict,hash_child)))then

                 if(m%occupied(m%free_cache))call destage(s,r%ngridmax+m%free_cache,hash_dict)

                 call hash_setp(hash_dict,hash_child,m%grid(ichild))

                 m%occupied(m%free_cache)=.true.
                 m%parent_cpu(m%free_cache)=grid_cpu
                 m%dirty(m%free_cache)=.false.

                 ! Set the grid index of the requested grid
                 if(same_keys(hash_key,hash_child))then
                    child_grid=ichild
                 endif

                 ! Unpack response to fetch request
                 call unpack_fetch%proc(m%grid(ichild),mdl%size_msg_array,mdl%recv_fetch_array(iskip:iskip+mdl%size_msg_array-1),hash_child)

                 ! If we also have also a flush cache...
                 ! This is for combined read-write cache operations
                 if(flush_cache)then
                    m%dirty(m%free_cache)=.true.

                    ! Set initialisation rule for combiner operations
                    call init_flush%proc(m%grid(ichild),hash_child)

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
        call hash_setp(hash_dict,hash_key,m%grid(child_grid))

        ! Store the grid coordinates
        m%grid(child_grid)%lev=hash_key(0)
        m%grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
        m%occupied(m%free_cache)=.true.
        m%parent_cpu(m%free_cache)=grid_cpu
        m%dirty(m%free_cache)=.true.

        ! Set initialisation rule for combiner operation
        call init_flush%proc(m%grid(child_grid),hash_key)

        ! Go to next free cache line
        m%free_cache=m%free_cache+1
        m%ncache=m%ncache+1
        if(m%free_cache.GT.r%ncachemax)then
           m%free_cache=1
        endif
        if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax

     end if
#endif

  ! Cartesian key is outside conputational domain
  ! We need to create a boundary condition grid in the cache
  else

     call get_bound(s,ix,ilevel,ibound)

     hash_ref(0)=ilevel
     hash_ref(1:ndim)=ix(1:ndim)
     hash_ref(r%bound_dir(ibound))=hash_ref(r%bound_dir(ibound))+r%bound_shift(ibound)

     call c_f_pointer(hash_getp(hash_dict,hash_ref,absent),child_ref)

     if (associated(child_ref))then

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
        call hash_setp(hash_dict,hash_key,m%grid(child_grid))

        ! Store the grid coordinates
        m%grid(child_grid)%lev=hash_key(0)
        m%grid(child_grid)%ckey(1:ndim)=hash_key(1:ndim)
        m%occupied(m%free_cache)=.true.
        m%parent_cpu(m%free_cache)=0
        m%dirty(m%free_cache)=.false.

        ! Set initialisation rule for boundary grid using reference grid
        if(associated(init_bound%proc))then
           call init_bound%proc(r,g,m,m%grid(child_grid),child_ref,ibound)
        endif

        ! Go to next free cache line
        m%free_cache=m%free_cache+1
        m%ncache=m%ncache+1
        if(m%free_cache.GT.r%ncachemax)then
           m%free_cache=1
        endif
        if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax

     else

        ! Delete old null grid if occupied
        if(m%occupied_null(m%free_null))then
           hash_child(0)=m%lev_null(m%free_null)
           hash_child(1:ndim)=m%ckey_null(1:ndim,m%free_null)
           call hash_free(hash_dict,hash_child)
        endif
        call hash_setp(hash_dict,hash_key)
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
        child_grid = 0

     endif

  endif

  if(child_grid>0)then
     child => m%grid(child_grid)
     if (present(lock).and.associated(child)) then
       if (lock) call lock_cache(s,child)
     endif
  else
     nullify(child)
  endif

#endif /* MDL2 */

  end associate
end subroutine get_grid
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine interpol_phi(mdl,m,grid_nbor,ind_nbor,ccc,bbb,tfrac,phi_int)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: nbor,oct
  use amr_commons, only: mesh_t
  use mdl_module
  implicit none
  type(mdl_t)::mdl
  type(mesh_t)::m
  integer,dimension(1:threetondim)::ind_nbor
  type(nbor),dimension(1:threetondim)::grid_nbor
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbb
  real(dp)::tfrac
  real(dp),dimension(1:twotondim)::phi_int
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (using CIC) and - if adaptive 
  ! timestepping is on - also in time (using linear extrapolation 
  ! of the change in phi during the last coarse step onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::ind,ind_average,ind_father
  integer::ind_nbr,ind_cen
  real(dp)::coeff,add
  type(oct),pointer::grid_cen,grid_nbr
#ifdef GRAV
  ! Store central cell
  grid_cen=>grid_nbor(threetondim/2+1)%p
  ind_cen=ind_nbor(threetondim/2+1)

  ! Third order phi interpolation
  do ind=1,twotondim
     phi_int(ind)=0d0
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbb(ind_average)
        grid_nbr=>grid_nbor(ind_father)%p
        ind_nbr=ind_nbor(ind_father)
        if (.not.associated(grid_nbr)) then 
           write(*,*)'no all neighbors present in interpol_phi...'
           call mdl_abort(mdl) ! Remove in case it happens
           add=coeff*(grid_cen%phi(ind_cen)+&
                & (grid_cen%phi(ind_cen)-grid_cen%phi_old(ind_cen))*tfrac)
        else
           add=coeff*(grid_nbr%phi(ind_nbr)+&
                & (grid_nbr%phi(ind_nbr)-grid_nbr%phi_old(ind_nbr))*tfrac)
        endif
        phi_int(ind)=phi_int(ind)+add
     end do
  end do
#endif
end subroutine interpol_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_save_phi_old(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SAVE_PHI_OLD,pst%iUpper+1,input_size,0,ilevel)
     call r_save_phi_old(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call save_phi_old(pst%s%m,ilevel)
  endif

end subroutine r_save_phi_old
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine save_phi_old(m,ilevel)
  use amr_parameters, only: ndim,dp,twotondim,threetondim
  use amr_commons, only: mesh_t
  implicit none
  type(mesh_t)::m
  integer ilevel
  ! Save the old potential for time extrapolation in case of subcycling
  integer::ind,igrid

#ifdef GRAV
  ! Loop over level grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Save phi      
        m%grid(igrid)%phi_old(ind)=m%grid(igrid)%phi(ind)
     end do
  end do
#endif

end subroutine save_phi_old
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module nbors_utils

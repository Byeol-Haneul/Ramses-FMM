!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_threetondim_nbor_parent_cell(hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_commons
  use hash
  implicit none
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
  integer,dimension(1:twotondim),save::igrid_twotondim_nbor
  integer(kind=8),dimension(0:ndim)::hash_nbor,hash_ref,hash_father
  integer(kind=8),dimension(1:ndim)::ii
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer::i1,j1,k1
  integer,save::i1min=-1
  integer,save::i1max=+1
  integer,save::j1min=0*(1-ndim/2)-1*(ndim/2)
  integer,save::j1max=0*(1-ndim/2)+1*(ndim/2)
  integer,save::k1min=0*(1-ndim/3)-1*(ndim/3)
  integer,save::k1max=0*(1-ndim/3)+1*(ndim/3)
  integer::ind,ipos,idim,get_grid,ilevel,inbor

  ilevel=hash_key(0)

  hash_father(0)=hash_key(0)-1

  ! Gather twotondim neighboring father grids
  do inbor=1,twotondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift_oct(1:ndim,inbor)
     ! Periodic boundary conditons
     do idim=1,ndim
        if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
        if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
     enddo
     hash_father(1:ndim)=hash_nbor(1:ndim)/2
     ! Store lower left neighbor coordinates 
     if(inbor==1)hash_ref(1:ndim)=hash_father(1:ndim)
     ! Get grid into memory and lock it if remote 
     ipos=get_grid(hash_father,hash_dict,flush_cache,fetch_cache)
     call lock_cache(ipos)
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
              if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
              if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
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
              if(ii(idim)<0)ii(idim)=ii(idim)+ckey_max(ilevel-1)
           enddo
           ind=1
           do idim=1,ndim
              ind=ind+2**(idim-1)*ii(idim)
           end do
           igrid_nbor(inbor)=igrid_twotondim_nbor(ind)
        end do
     end do
  end do

end subroutine get_threetondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_twondim_nbor_parent_cell(hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache)
  use amr_commons
  use hash
  implicit none
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
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::ind,ipos,idim,get_grid,ilevel,inbor

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
  ipos=get_grid(hash_father,hash_dict,flush_cache,fetch_cache)
  call lock_cache(ipos)
  igrid_nbor(0)=ipos
  ind_nbor(0)=ind
  
  ! Deal with neighboring father cells
  do inbor=1,twondim
     hash_nbor(1:ndim)=hash_key(1:ndim)+shift(1:ndim,inbor)

     ! Periodic boundary conditons
     do idim=1,ndim
        if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
        if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
     enddo
     hash_father(1:ndim)=hash_nbor(1:ndim)/2
     ii(1:ndim)=hash_nbor(1:ndim)-2*hash_father(1:ndim)
     ind=1
     do idim=1,ndim
        ind=ind+2**(idim-1)*ii(idim)
     end do

     ! Get grid into memory and lock it if remote 
     ipos=get_grid(hash_father,hash_dict,flush_cache,fetch_cache)
     call lock_cache(ipos)
     igrid_nbor(inbor)=ipos
     ind_nbor(inbor)=ind
  end do

end subroutine get_twondim_nbor_parent_cell
!###############################################################
!###############################################################
!###############################################################
!###############################################################
integer function get_parent_cell(hash_key,hash_dict,flush_cache,fetch_cache) result(parent_cell)
  use amr_commons
  use hash
  implicit none
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  !
  ! This routine acquires the parent cell of the grid 
  ! corresponding to the input hash key.
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
  ipos=get_grid(hash_father,hash_dict,flush_cache,fetch_cache)
  parent_cell=0
  if(ipos>0)parent_cell=(ipos-1)*twotondim+ind
end function get_parent_cell
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
!##############################################################
!##############################################################
!##############################################################
!##############################################################
integer function get_grid(hash_key,hash_dict,flush_cache,fetch_cache) result(child_grid)
  use amr_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(0:ndim)::hash_key
  type(hash_table)::hash_dict
  !
  ! This routine acquires the grid 
  ! corresponding to the input hash key.
  !
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix
  integer(kind=8),dimension(0:ndim)::hash_child
  integer::i,ind,idim,ivar,ichild,ilevel,info,icpu,grid_cpu,ntile_response,icounter
  integer::send_request_id
  type(request),save::send_request
  integer::response_id
  type(int4_msg),save::response_flag
  type(realdp_msg),save::response_hydro
  type(small_realdp_msg),save::response_poisson
  type(large_realdp_msg),save::response_refine
  type(twin_realdp_msg),save::response_mg
  type(three_realdp_msg),save::response_interpol
  logical::failed_request,send_request_completed
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::send_request_status
#endif
 
#ifndef WITHOUTMPI

  ! If counter is good, check on incoming messages and perform actions
  if(mail_counter==32)then
     call check_mail(MPI_REQUEST_NULL,hash_dict)
     mail_counter=0
  endif
  mail_counter=mail_counter+1
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
  ix(1,1:ndim)=hash_key(1:ndim)
  call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

  ! Check if grid sits inside processor boundaries
  if(    ge_keys(hk(1,1:nhilbert),bound_hilbert_key(1:nhilbert,myid-1,ilevel)).AND. &
       & gt_keys(bound_hilbert_key(1:nhilbert,myid,ilevel),hk(1,1:nhilbert)))then
     return
  endif

  ! Determine parent processor
  do icpu=1,ncpu
     if(    ge_keys(hk(1,1:nhilbert),bound_hilbert_key(1:nhilbert,icpu-1,ilevel)).AND. &
          & gt_keys(bound_hilbert_key(1:nhilbert,icpu,ilevel),hk(1,1:nhilbert)))then
        grid_cpu=icpu
        exit
     end if
  end do

  !============================================
  ! We have a fetch and possibly a flush cache
  ! Fetch alone means read-only cache operations.
  ! Both together means read-write cache.
  !============================================
  if(fetch_cache)then

     ! Send a request to the relevant cpu
     send_request%lev=ilevel
     send_request%ckey(1:ndim)=hash_key(1:ndim)
     
     ! Post RECV for the expected response
     if(cache_operation_type.EQ.operation_type_flag)then
        call MPI_IRECV(response_flag,1,new_mpi_int4_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_hydro)then
        call MPI_IRECV(response_hydro,1,new_mpi_realdp_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_poisson)then
        call MPI_IRECV(response_poisson,1,new_mpi_small_realdp_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_refine)then
        call MPI_IRECV(response_refine,1,new_mpi_large_realdp_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_mg)then
        call MPI_IRECV(response_mg,1,new_mpi_twin_realdp_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     if(cache_operation_type.EQ.operation_type_interpol)then
        call MPI_IRECV(response_interpol,1,new_mpi_three_realdp_msg,&
             & grid_cpu-1,msg_tag,MPI_COMM_WORLD,response_id,info)  
     endif
     call MPI_ISEND(send_request,1,new_mpi_request,&
          & grid_cpu-1,request_tag,MPI_COMM_WORLD,send_request_id,info)

     ! While waiting for reply, check on incoming messages and perform actions
     call check_mail(response_id,hash_dict)

     ! Wait for ISEND completion to free memory in corresponding MPI buffer
     call MPI_WAIT(send_request_id,send_request_status,info)
     
     if(cache_operation_type.EQ.operation_type_flag)then
        failed_request=response_flag%type==-1
        ntile_response=response_flag%ntile
     endif
     if(cache_operation_type.EQ.operation_type_hydro)then
        failed_request=response_hydro%type==-1
        ntile_response=response_hydro%ntile
     endif
     if(cache_operation_type.EQ.operation_type_poisson)then
        failed_request=response_poisson%type==-1
        ntile_response=response_poisson%ntile
     endif
     if(cache_operation_type.EQ.operation_type_refine)then
        failed_request=response_refine%type==-1
        ntile_response=response_refine%ntile
     endif
     if(cache_operation_type.EQ.operation_type_mg)then
        failed_request=response_mg%type==-1
        ntile_response=response_mg%ntile
     endif
     if(cache_operation_type.EQ.operation_type_interpol)then
        failed_request=response_interpol%type==-1
        ntile_response=response_interpol%ntile
     endif

     ! If grid does not exist, store -1 in the cache
     ! The output grid index is still zero
     if(failed_request)then

        ! Delete old null grid if occupied
        if(occupied_null(free_null))then
           hash_child(0)=lev_null(free_null)
           hash_child(1:ndim)=ckey_null(1:ndim,free_null)
           call hash_free(hash_dict,hash_child)
        endif
        call hash_set(hash_dict,hash_key,-1)
        occupied_null(free_null)=.true.
        lev_null(free_null)=ilevel
        ckey_null(1:ndim,free_null)=hash_key(1:ndim)

        ! Go to next free cache line
        free_null=free_null+1
        nnull=nnull+1
        if(free_null.GT.ncachemax)then
           free_null=1
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
           if(cache_operation_type.EQ.operation_type_poisson)then
              hash_child(0)=response_poisson%lev(i)
              hash_child(1:ndim)=response_poisson%ckey(1:ndim,i)
           endif
           if(cache_operation_type.EQ.operation_type_refine)then
              hash_child(0)=response_refine%lev(i)
              hash_child(1:ndim)=response_refine%ckey(1:ndim,i)
           endif
           if(cache_operation_type.EQ.operation_type_mg)then
              hash_child(0)=response_mg%lev(i)
              hash_child(1:ndim)=response_mg%ckey(1:ndim,i)
           endif
           if(cache_operation_type.EQ.operation_type_interpol)then
              hash_child(0)=response_interpol%lev(i)
              hash_child(1:ndim)=response_interpol%ckey(1:ndim,i)
           endif
           ichild=ngridmax+free_cache

           if(hash_get(hash_dict,hash_child).EQ.0)then

              if(occupied(free_cache))call destage(ngridmax+free_cache,hash_dict)

              call hash_set(hash_dict,hash_child,ichild)
              
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

              ! Operations of type "flag"
              if(cache_operation_type.EQ.operation_type_flag)then

                 if(cache_operation.EQ.operation_split)then
                    do ind=1,twotondim
                       if(response_flag%int4(ind,i)==1)then
                          grid(ichild)%refined(ind)=.true.
                       else
                          grid(ichild)%refined(ind)=.false.
                       end if
                    end do
                 else
                    grid(ichild)%flag1(1:twotondim)=response_flag%int4(1:twotondim,i)
                 endif

              endif

              ! Operation of type "hydro"
              if(cache_operation_type.EQ.operation_type_hydro)then
                 do ind=1,twotondim
                    if(response_hydro%int4(ind,i)==1)then
                       grid(ichild)%refined(ind)=.true.
                    else
                       grid(ichild)%refined(ind)=.false.
                    end if
#ifdef HYDRO
                    do ivar=1,nvar
                       grid(ichild)%uold(ind,ivar)=response_hydro%realdp(ind,ivar,i)
                    end do
#endif
                 end do
              endif

              ! Operations of type "poisson"
              if(cache_operation_type.EQ.operation_type_poisson)then
#ifdef GRAV
                 if(cache_operation.EQ.operation_restrict_res)then
                    do ind=1,twotondim
                       grid(ichild)%f(ind,3)=response_poisson%realdp(ind,i)
                    end do
                 else if(cache_operation.EQ.operation_scan)then
                    do ind=1,twotondim
                       grid(ichild)%f(ind,3)=response_poisson%realdp(ind,i)
                    end do
                 else if(cache_operation.EQ.operation_cg)then
                    do ind=1,twotondim
                       grid(ichild)%f(ind,2)=response_poisson%realdp(ind,i)
                    end do
                 else
                    do ind=1,twotondim
                       grid(ichild)%phi(ind)=response_poisson%realdp(ind,i)
                    end do
                 endif
#endif
              endif

              ! Operations of type "refine"
              if(cache_operation_type.EQ.operation_type_refine)then
                 do ind=1,twotondim
                    if(response_refine%int4(ind,i)==1)then
                       grid(ichild)%refined(ind)=.true.
                    else
                       grid(ichild)%refined(ind)=.false.
                    end if
#ifdef HYDRO
                    do ivar=1,nvar
                       grid(ichild)%uold(ind,ivar)=response_refine%realdp_hydro(ind,ivar,i)
                    end do
#endif
#ifdef GRAV
                    do idim=1,ndim
                       grid(ichild)%f(ind,idim)=response_refine%realdp_poisson(ind,idim,i)
                    end do
                    grid(ichild)%phi(ind)=response_refine%realdp_poisson(ind,ndim+1,i)
                    grid(ichild)%phi_old(ind)=response_refine%realdp_poisson(ind,ndim+2,i)
#endif
                 end do
              endif
              
              ! Operations of type "multigrid"
              if(cache_operation_type.EQ.operation_type_mg)then
                 do ind=1,twotondim
#ifdef GRAV
                    grid(ichild)%phi(ind)=response_mg%realdp_phi(ind,i)
                    grid(ichild)%f(ind,3)=response_mg%realdp_dis(ind,i)
#endif
                 end do
              endif
              
              ! Operations of type "interpol"
              if(cache_operation_type.EQ.operation_type_interpol)then
                 if(cache_operation.EQ.operation_kick)then
                    do ind=1,twotondim
#ifdef GRAV
                       grid(ichild)%f(ind,1)=response_interpol%realdp_phi(ind,i)
                       grid(ichild)%f(ind,2)=response_interpol%realdp_phi_old(ind,i)
                       grid(ichild)%f(ind,3)=response_interpol%realdp_dis(ind,i)
#endif
                    end do
                 else
                    do ind=1,twotondim
#ifdef GRAV
                       grid(ichild)%phi(ind)=response_interpol%realdp_phi(ind,i)
                       grid(ichild)%phi_old(ind)=response_interpol%realdp_phi_old(ind,i)
                       grid(ichild)%f(ind,3)=response_interpol%realdp_dis(ind,i)
#endif
                    end do
                 endif
              endif
              
              ! If we also have a flush cache...
              ! This is only for combined read-write cache operations
              if(flush_cache)then
                 dirty(free_cache)=.true.
                 
                 !================================================
                 ! Set initialisation rule for combiner operations
                 !================================================
                 
#ifdef HYDRO                 
                 ! Initialisation rule for Godunov update
                 if(cache_operation.EQ.operation_godunov)then           
                    grid(ichild)%unew(1:twotondim,1:nvar)=0.0
                 endif
#endif                 
                 ! Initialisation rule for derefine
                 if(cache_operation.EQ.operation_derefine)then           
                    grid(ichild)%refined(1:twotondim)=.true.
                 endif
#ifdef GRAV                 
                 ! Initialisation rule for restrict multigrid residual
                 if(cache_operation.EQ.operation_restrict_res)then           
                    grid(ichild)%f(1:twotondim,2)=0.0
                 endif
#endif                 
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

     !=================================================
     ! If we have only a flush cache (write-only cache) 
     !=================================================
  else if(flush_cache)then   

     ! If next cache line is occupied, free it.
     if(locked(free_cache))then
        do while(locked(free_cache))
           free_cache=free_cache+1
           if(free_cache>ncachemax)free_cache=1
        end do
     end if
     if(occupied(free_cache))call destage(ngridmax+free_cache,hash_dict)

     ! Set grid index to a virtual grid in local memory
     child_grid=ngridmax+free_cache
     call hash_set(hash_dict,hash_key,child_grid)

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

#ifdef HYDRO
     ! Initialisation rule for hydro upload
     if(cache_operation.EQ.operation_upload)then           
        grid(child_grid)%uold(1:twotondim,1:nvar)=0.0
     endif

     ! Initialisation rule for hydro multipole
     if(cache_operation.EQ.operation_multipole)then           
        grid(child_grid)%unew(1:twotondim,1:ndim+1)=0.0
     endif
#endif

#ifdef GRAV
     ! Initialisation rule for mass deposition
     if(cache_operation.EQ.operation_rho)then           
        grid(child_grid)%rho(1:twotondim)=0.0
     endif

     ! Initialisation rule for MG mask restriction
     if(cache_operation.EQ.operation_restrict_mask)then           
        grid(child_grid)%f(1:twotondim,3)=0.0
     endif
#endif

     ! Go to next free cache line
     free_cache=free_cache+1
     ncache=ncache+1
     if(free_cache.GT.ncachemax)then
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
subroutine check_mail(comm_id,hash_dict)
  use amr_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(hash_table)::hash_dict
  integer::comm_id
  !
  integer::i,ind,ivar,idim,info,ipos,igrid,ichild,grid_cpu,ilevel,itile,ntile_reply
  logical::comm_completed,request_received,flush_received=.false.,reply_sent
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status,comm_status
#endif
  integer(kind=4), dimension(1:nvector),save::dummy_state
  integer(kind=8), dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8), dimension(1:nvector,1:ndim),save::ix
  integer(kind=8),dimension(0:ndim)::hash_key,hash_child
  !
#ifndef WITHOUTMPI
  comm_completed=.false.
  do while (.not. comm_completed)

     ! USE A NESTED DO WHILE LOOP FOR THE 2 ADD. TESTS

     !===========================
     ! Check for incoming request
     !===========================
     request_received=.true.
     do while(request_received)
        call MPI_Test(request_id,request_received,request_status,info)
        if(request_received)then
           
           ! Assemble a reply and send it back
           ilevel=recv_request%lev
           hash_key(0)=ilevel
           hash_key(1:ndim)=recv_request%ckey(1:ndim)
           igrid=hash_get(hash_dict,hash_key)
           grid_cpu=request_status(MPI_SOURCE)+1
           
           ! If grid does not exist, send a null reply
           if(igrid.EQ.0)then
              if(cache_operation_type.EQ.operation_type_flag)then
                 reply_flag(grid_cpu)%type=-1
              endif
              if(cache_operation_type.EQ.operation_type_hydro)then
                 reply_hydro(grid_cpu)%type=-1
              endif
              if(cache_operation_type.EQ.operation_type_poisson)then
                 reply_poisson(grid_cpu)%type=-1
              endif
              if(cache_operation_type.EQ.operation_type_refine)then
                 reply_refine(grid_cpu)%type=-1
              endif
              if(cache_operation_type.EQ.operation_type_mg)then
                 reply_mg(grid_cpu)%type=-1
              endif
              if(cache_operation_type.EQ.operation_type_interpol)then
                 reply_interpol(grid_cpu)%type=-1
              endif
              
              ! Otherwise, assemble a proper reply with a complete tile
           else
              itile=(igrid-head_cache(ilevel))/ntilemax
              ntile_reply=MIN(tail_cache(ilevel)-itile*ntilemax-head_cache(ilevel)+1,ntilemax)
              
              ! Store type of reply and number of entries
              if(cache_operation_type.EQ.operation_type_flag)then
                 reply_flag(grid_cpu)%type=1
                 reply_flag(grid_cpu)%ntile=ntile_reply
              endif
              if(cache_operation_type.EQ.operation_type_hydro)then
                 reply_hydro(grid_cpu)%type=1
                 reply_hydro(grid_cpu)%ntile=ntile_reply
              endif
              if(cache_operation_type.EQ.operation_type_poisson)then
                 reply_poisson(grid_cpu)%type=1
                 reply_poisson(grid_cpu)%ntile=ntile_reply
              endif
              if(cache_operation_type.EQ.operation_type_refine)then
                 reply_refine(grid_cpu)%type=1
                 reply_refine(grid_cpu)%ntile=ntile_reply
              endif
              if(cache_operation_type.EQ.operation_type_mg)then
                 reply_mg(grid_cpu)%type=1
                 reply_mg(grid_cpu)%ntile=ntile_reply
              endif
              if(cache_operation_type.EQ.operation_type_interpol)then
                 reply_interpol(grid_cpu)%type=1
                 reply_interpol(grid_cpu)%ntile=ntile_reply
              endif
              
              ! Store data, depending on reply type
              do i=1,ntile_reply
                 ipos=head_cache(ilevel)+itile*ntilemax+i-1
                 
                 ! Reply of type flag
                 if(cache_operation_type.EQ.operation_type_flag)then
                    reply_flag(grid_cpu)%lev(i)=grid(ipos)%lev
                    reply_flag(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)

                    if(cache_operation.EQ.operation_split)then
                       do ind=1,twotondim
                          if(grid(ipos)%refined(ind))then
                             reply_flag(grid_cpu)%int4(ind,i)=1
                          else
                             reply_flag(grid_cpu)%int4(ind,i)=0
                          endif
                       end do
                    else
                       do ind=1,twotondim
                          reply_flag(grid_cpu)%int4(ind,i)=grid(ipos)%flag1(ind)
                       end do
                    endif
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
#ifdef HYDRO
                       do ivar=1,nvar
                          reply_hydro(grid_cpu)%realdp(ind,ivar,i)=grid(ipos)%uold(ind,ivar)
                       end do
#endif
                    end do
                 endif
                 
                 ! Reply of type poisson
                 if(cache_operation_type.EQ.operation_type_poisson)then
                    reply_poisson(grid_cpu)%lev(i)=grid(ipos)%lev
                    reply_poisson(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
#ifdef GRAV
                    if(cache_operation.EQ.operation_restrict_res)then
                       do ind=1,twotondim
                          reply_poisson(grid_cpu)%realdp(ind,i)=grid(ipos)%f(ind,3)
                       end do
                    else if(cache_operation.EQ.operation_scan)then
                       do ind=1,twotondim
                          reply_poisson(grid_cpu)%realdp(ind,i)=grid(ipos)%f(ind,3)
                       end do
                    else if(cache_operation.EQ.operation_cg)then
                       do ind=1,twotondim
                          reply_poisson(grid_cpu)%realdp(ind,i)=grid(ipos)%f(ind,2)
                       end do
                    else
                       do ind=1,twotondim
                          reply_poisson(grid_cpu)%realdp(ind,i)=grid(ipos)%phi(ind)
                       end do
                    endif
#endif
                 endif
                 
                 ! Reply of type refine
                 if(cache_operation_type.EQ.operation_type_refine)then
                    reply_refine(grid_cpu)%lev(i)=grid(ipos)%lev
                    reply_refine(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
                    do ind=1,twotondim
                       if(grid(ipos)%refined(ind))then
                          reply_refine(grid_cpu)%int4(ind,i)=1
                       else
                          reply_refine(grid_cpu)%int4(ind,i)=0
                       endif
#ifdef HYDRO
                       do ivar=1,nvar
                          reply_refine(grid_cpu)%realdp_hydro(ind,ivar,i)=grid(ipos)%uold(ind,ivar)
                       end do
#endif
#ifdef GRAV
                       do idim=1,ndim
                          reply_refine(grid_cpu)%realdp_poisson(ind,idim,i)=grid(ipos)%f(ind,idim)
                       end do
                       reply_refine(grid_cpu)%realdp_poisson(ind,ndim+1,i)=grid(ipos)%phi(ind)
                       reply_refine(grid_cpu)%realdp_poisson(ind,ndim+2,i)=grid(ipos)%phi_old(ind)
#endif
                    end do
                 endif

                 ! Reply of type multigrid
                 if(cache_operation_type.EQ.operation_type_mg)then
                    reply_mg(grid_cpu)%lev(i)=grid(ipos)%lev
                    reply_mg(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
                    do ind=1,twotondim
#ifdef GRAV
                       reply_mg(grid_cpu)%realdp_dis(ind,i)=grid(ipos)%f(ind,3)
                       reply_mg(grid_cpu)%realdp_phi(ind,i)=grid(ipos)%phi(ind)
#endif
                    end do
                 endif

                 ! Reply of type interpol
                 if(cache_operation_type.EQ.operation_type_interpol)then
                    reply_interpol(grid_cpu)%lev(i)=grid(ipos)%lev
                    reply_interpol(grid_cpu)%ckey(1:ndim,i)=grid(ipos)%ckey(1:ndim)
                    if(cache_operation.EQ.operation_kick)then
                       do ind=1,twotondim
#ifdef GRAV
                          reply_interpol(grid_cpu)%realdp_phi(ind,i)=grid(ipos)%f(ind,1)
                          reply_interpol(grid_cpu)%realdp_phi_old(ind,i)=grid(ipos)%f(ind,2)
                          reply_interpol(grid_cpu)%realdp_dis(ind,i)=grid(ipos)%f(ind,3)
#endif
                       end do
                    else
                       do ind=1,twotondim
#ifdef GRAV
                          reply_interpol(grid_cpu)%realdp_phi(ind,i)=grid(ipos)%phi(ind)
                          reply_interpol(grid_cpu)%realdp_phi_old(ind,i)=grid(ipos)%phi_old(ind)
                          reply_interpol(grid_cpu)%realdp_dis(ind,i)=grid(ipos)%f(ind,3)
#endif
                       end do
                    endif
                 endif

              end do
           endif
           
           ! Wait for the old SEND to free memory in corresponding MPI buffer
           call MPI_WAIT(reply_id(grid_cpu),reply_status,info)
           
           ! Send back the reply
           if(cache_operation_type.EQ.operation_type_flag)then
              call MPI_ISEND(reply_flag(grid_cpu),1,new_mpi_int4_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           if(cache_operation_type.EQ.operation_type_hydro)then
              call MPI_ISEND(reply_hydro(grid_cpu),1,new_mpi_realdp_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           if(cache_operation_type.EQ.operation_type_poisson)then
              call MPI_ISEND(reply_poisson(grid_cpu),1,new_mpi_small_realdp_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           if(cache_operation_type.EQ.operation_type_refine)then
              call MPI_ISEND(reply_refine(grid_cpu),1,new_mpi_large_realdp_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           if(cache_operation_type.EQ.operation_type_mg)then
              call MPI_ISEND(reply_mg(grid_cpu),1,new_mpi_twin_realdp_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           if(cache_operation_type.EQ.operation_type_interpol)then
              call MPI_ISEND(reply_interpol(grid_cpu),1,new_mpi_three_realdp_msg,&
                   & grid_cpu-1,msg_tag,MPI_COMM_WORLD,reply_id(grid_cpu),info)
           endif
           
           ! Post a new RECV for request
           call MPI_IRECV(recv_request,1,new_mpi_request,&
                & MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,request_id,info)
        endif
     end do

     !=========================
     ! Check for incoming flush
     !=========================
     flush_received=.true.
     do while(flush_received)
        call MPI_Test(flush_id,flush_received,flush_status,info)
        if(flush_received)then

        !===========================================================
        ! Combine received data to local memory using combiner rules
        !===========================================================

        !------------------------------------------------------
        ! Combiner rules for initflag
        !------------------------------------------------------
        if(cache_operation.EQ.operation_initflag)then
           do i=1,recv_flush_flag%nflush
              hash_child(0)=recv_flush_flag%lev(i)
              hash_child(1:ndim)=recv_flush_flag%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
              do ind=1,twotondim
                 grid(ichild)%flag1(ind)=MAX(grid(ichild)%flag1(ind),recv_flush_flag%int4(ind,i))
              end do
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for derefine
        !------------------------------------------------------
        if(cache_operation.EQ.operation_derefine)then
           do i=1,recv_flush_flag%nflush
              hash_child(0)=recv_flush_flag%lev(i)
              hash_child(1:ndim)=recv_flush_flag%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
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

        !------------------------------------------------------
        ! Combiner rules for hydro upload
        !------------------------------------------------------
        if(cache_operation.EQ.operation_upload)then
           do i=1,recv_flush_hydro%nflush
              hash_child(0)=recv_flush_hydro%lev(i)
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef HYDRO
              do ivar=1,nvar
                 do ind=1,twotondim
                    if(grid(ichild)%refined(ind))then
                       grid(ichild)%uold(ind,ivar)=grid(ichild)%uold(ind,ivar)&
                            & +recv_flush_hydro%realdp(ind,ivar,i)
                    endif
                 end do
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for hydro multipole
        !------------------------------------------------------
        if(cache_operation.EQ.operation_multipole)then
           do i=1,recv_flush_hydro%nflush
              hash_child(0)=recv_flush_hydro%lev(i)
              hash_child(1:ndim)=recv_flush_hydro%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef HYDRO
              do ivar=1,ndim+1
                 do ind=1,twotondim
                    if(grid(ichild)%refined(ind))then
                       grid(ichild)%unew(ind,ivar)=grid(ichild)%unew(ind,ivar)&
                            & +recv_flush_hydro%realdp(ind,ivar,i)
                    endif
                 end do
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for godunov update
        !------------------------------------------------------
        if(cache_operation.EQ.operation_godunov)then
           do i=1,recv_flush_refine%nflush
              hash_child(0)=recv_flush_refine%lev(i)
              hash_child(1:ndim)=recv_flush_refine%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef HYDRO
              do ivar=1,nvar
                 do ind=1,twotondim
                    grid(ichild)%unew(ind,ivar)=grid(ichild)%unew(ind,ivar)&
                         & +recv_flush_refine%realdp_hydro(ind,ivar,i)
                 end do
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for mass deposition
        !------------------------------------------------------
        if(cache_operation.EQ.operation_rho)then
           do i=1,recv_flush_poisson%nflush
              hash_child(0)=recv_flush_poisson%lev(i)
              hash_child(1:ndim)=recv_flush_poisson%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef GRAV
              do ind=1,twotondim
                 grid(ichild)%rho(ind)=grid(ichild)%rho(ind)&
                      & +recv_flush_poisson%realdp(ind,i)
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for MG mask restriction
        !------------------------------------------------------
        if(cache_operation.EQ.operation_restrict_mask)then
           do i=1,recv_flush_poisson%nflush
              hash_child(0)=recv_flush_poisson%lev(i)
              hash_child(1:ndim)=recv_flush_poisson%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef GRAV
              do ind=1,twotondim
                 grid(ichild)%f(ind,3)=grid(ichild)%f(ind,3)&
                      & +recv_flush_poisson%realdp(ind,i)
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for MG residual restriction
        !------------------------------------------------------
        if(cache_operation.EQ.operation_restrict_res)then
           do i=1,recv_flush_poisson%nflush
              hash_child(0)=recv_flush_poisson%lev(i)
              hash_child(1:ndim)=recv_flush_poisson%ckey(1:ndim,i)
              ichild=hash_get(hash_dict,hash_child)
#ifdef GRAV
              do ind=1,twotondim
                 if(grid(ichild)%f(ind,3)>0.0)then
                    grid(ichild)%f(ind,2)=grid(ichild)%f(ind,2)&
                         & +recv_flush_poisson%realdp(ind,i)
                 endif
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for refinements
        !------------------------------------------------------
        if(cache_operation.EQ.operation_refine)then

           do i=1,recv_flush_refine%nflush
              ilevel=recv_flush_refine%lev(i)
              hash_child(0)=ilevel
              hash_child(1:ndim)=recv_flush_refine%ckey(1:ndim,i)

              ! Compute Hilbert keys of new octs
              ix(1,1:ndim)=hash_child(1:ndim)
              call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

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
              grid(ichild)%hkey(1:nhilbert)=hk(1,1:nhilbert)
              grid(ichild)%refined(1:twotondim)=.false.
              grid(ichild)%flag1(1:twotondim)=0
              grid(ichild)%flag2(1:twotondim)=0
              grid(ichild)%superoct=1
              
              ! Insert new grid in hash table
              call hash_set(hash_dict,hash_child,ichild)

              ! Flush hydro variables
#ifdef HYDRO
              do ind=1,twotondim
                 do ivar=1,nvar
                    grid(ichild)%uold(ind,ivar)=recv_flush_refine%realdp_hydro(ind,ivar,i)
                 end do
              end do
#endif
              ! Flush gravity variables
#ifdef GRAV
              do ind=1,twotondim
                 do idim=1,ndim
                    grid(ichild)%f(ind,idim)=recv_flush_refine%realdp_poisson(ind,idim,i)
                 end do
                 grid(ichild)%phi(ind)=recv_flush_refine%realdp_poisson(ind,ndim+1,i)
                 grid(ichild)%phi_old(ind)=recv_flush_refine%realdp_poisson(ind,ndim+2,i)
              end do
#endif
           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for load balance
        !------------------------------------------------------
        if(cache_operation.EQ.operation_loadbalance)then

           do i=1,recv_flush_refine%nflush
              ilevel=recv_flush_refine%lev(i)
              hash_child(0)=ilevel
              hash_child(1:ndim)=recv_flush_refine%ckey(1:ndim,i)

              ! Compute Hilbert keys of new octs
              ix(1,1:ndim)=hash_child(1:ndim)
              call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

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
              grid(ichild)%hkey(1:nhilbert)=hk(1,1:nhilbert)
              do ind=1,twotondim
                 if(recv_flush_refine%int4(ind,i)==1)then
                    grid(ichild)%refined(ind)=.true.
                 else
                    grid(ichild)%refined(ind)=.false.
                 endif
                 ! Flush hydro variables
#ifdef HYDRO
                 do ivar=1,nvar
                    grid(ichild)%uold(ind,ivar)=recv_flush_refine%realdp_hydro(ind,ivar,i)
                 end do
#endif                 
                 ! Flush gravity variables
#ifdef GRAV
                 do idim=1,ndim
                    grid(ichild)%f(ind,idim)=recv_flush_refine%realdp_poisson(ind,idim,i)
                 end do
                 grid(ichild)%phi(ind)=recv_flush_refine%realdp_poisson(ind,ndim+1,i)
                 grid(ichild)%phi_old(ind)=recv_flush_refine%realdp_poisson(ind,ndim+2,i)
#endif
              end do
              grid(ichild)%flag1(1:twotondim)=0
              grid(ichild)%flag2(1:twotondim)=0
              grid(ichild)%superoct=1

              ! Insert new grid in hash table
              call hash_set(hash_dict,hash_child,ichild)

           end do
        endif

        !------------------------------------------------------
        ! Combiner rules for building multigrid hierarchy
        !------------------------------------------------------
        if(cache_operation.EQ.operation_build_mg)then

           do i=1,recv_flush_poisson%nflush
              ilevel=recv_flush_poisson%lev(i)
              hash_child(0)=ilevel
              hash_child(1:ndim)=recv_flush_poisson%ckey(1:ndim,i)

              ! MG oct could have been activated before,
              ! so we need to check if it exists
              ichild=hash_get(hash_dict,hash_child)

              ! If it doesn't exist then create it
              if(ichild==0)then 
                 
              ! Compute Hilbert keys of new octs
              ix(1,1:ndim)=hash_child(1:ndim)
              call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

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
              grid(ichild)%hkey(1:nhilbert)=hk(1,1:nhilbert)
              grid(ichild)%refined(1:twotondim)=.false.
              grid(ichild)%flag1(1:twotondim)=0
              grid(ichild)%flag2(1:twotondim)=0
              grid(ichild)%superoct=1
#ifdef GRAV
              ! Initialize gravity variables
              do ind=1,twotondim
                 grid(ichild)%f(ind,1:ndim)=0.0
                 grid(ichild)%phi(ind)=0.0
                 grid(ichild)%phi_old(ind)=0.0
              end do
#endif
              ! Insert new grid in hash table
              call hash_set(hash_dict,hash_child,ichild)
              
              endif

           end do
        endif

        !=================================
        ! Post a new RECV for flush
        !=================================
        if(cache_operation_type.EQ.operation_type_flag)then
           call MPI_IRECV(recv_flush_flag,1,new_mpi_int4_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_hydro)then
           call MPI_IRECV(recv_flush_hydro,1,new_mpi_realdp_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_poisson)then
           call MPI_IRECV(recv_flush_poisson,1,new_mpi_small_realdp_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_refine)then
           call MPI_IRECV(recv_flush_refine,1,new_mpi_large_realdp_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_mg)then
           call MPI_IRECV(recv_flush_mg,1,new_mpi_twin_realdp_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        if(cache_operation_type.EQ.operation_type_interpol)then
           call MPI_IRECV(recv_flush_interpol,1,new_mpi_three_realdp_flush,&
                & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
        endif
        endif
     end do

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
subroutine destage(igrid,hash_dict)
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(hash_table)::hash_dict
  integer::igrid
  !
  integer::ind,ivar,idim,ipos,info,icache,iflush,grid_cpu
  integer::send_flush_id
  integer(kind=8),dimension(0:ndim)::hash_key
  !
#ifndef WITHOUTMPI

  hash_key(0)=grid(igrid)%lev
  hash_key(1:ndim)=grid(igrid)%ckey(1:ndim)
  ipos=hash_get(hash_dict,hash_key)

  if(hash_get(hash_dict,hash_key).EQ.0)then
     write(*,*)'PE ',myid,' trying to free non existing grid'
     stop
  endif

  call hash_free(hash_dict,hash_key)

  ! Check if the destage requires a flush
  icache=igrid-ngridmax

  if(dirty(icache))then
     grid_cpu=parent_cpu(icache)
     dirty(icache)=.false.
  
     !------------------------------------------------------
     ! Filling the flush buffer for initflag
     !------------------------------------------------------
     if(cache_operation.EQ.operation_initflag)then
        if(send_flush_flag(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_flag(grid_cpu),1,new_mpi_int4_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_flag(grid_cpu)%nflush=0
        endif
        send_flush_flag(grid_cpu)%nflush=send_flush_flag(grid_cpu)%nflush+1
        iflush=send_flush_flag(grid_cpu)%nflush
        send_flush_flag(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_flag(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        send_flush_flag(grid_cpu)%int4(1:twotondim,iflush)=grid(igrid)%flag1(1:twotondim)
     endif
     
     !------------------------------------------------------
     ! Filling the flush buffer for derefine
     !------------------------------------------------------
     if(cache_operation.EQ.operation_derefine)then
        if(send_flush_flag(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_flag(grid_cpu),1,new_mpi_int4_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_flag(grid_cpu)%nflush=0
        endif
        send_flush_flag(grid_cpu)%nflush=send_flush_flag(grid_cpu)%nflush+1
        iflush=send_flush_flag(grid_cpu)%nflush
        send_flush_flag(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_flag(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        do ind=1,twotondim
           if(grid(igrid)%refined(ind))then
              send_flush_flag(grid_cpu)%int4(ind,iflush)=1
           else
              send_flush_flag(grid_cpu)%int4(ind,iflush)=0
           endif
        end do
     endif
     
     !------------------------------------------------------
     ! Filling the flush buffer for hydro upload
     !------------------------------------------------------
     if(cache_operation.EQ.operation_upload)then
        if(send_flush_hydro(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_hydro(grid_cpu)%nflush=0
        endif
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef HYDRO
        do ind=1,twotondim
           do ivar=1,nvar
              send_flush_hydro(grid_cpu)%realdp(ind,ivar,iflush)=grid(igrid)%uold(ind,ivar)
           end do
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for multipole
     !------------------------------------------------------
     if(cache_operation.EQ.operation_multipole)then
        if(send_flush_hydro(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_hydro(grid_cpu),1,new_mpi_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_hydro(grid_cpu)%nflush=0
        endif
        send_flush_hydro(grid_cpu)%nflush=send_flush_hydro(grid_cpu)%nflush+1
        iflush=send_flush_hydro(grid_cpu)%nflush
        send_flush_hydro(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_hydro(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef HYDRO
        do ind=1,twotondim
           do ivar=1,ndim+1
              send_flush_hydro(grid_cpu)%realdp(ind,ivar,iflush)=grid(igrid)%unew(ind,ivar)
           end do
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for Godunov solver
     !------------------------------------------------------
     if(cache_operation.EQ.operation_godunov)then
        if(send_flush_refine(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_refine(grid_cpu),1,new_mpi_large_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_refine(grid_cpu)%nflush=0
        endif
        send_flush_refine(grid_cpu)%nflush=send_flush_refine(grid_cpu)%nflush+1
        iflush=send_flush_refine(grid_cpu)%nflush
        send_flush_refine(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_refine(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef HYDRO
        do ind=1,twotondim
           do ivar=1,nvar
              send_flush_refine(grid_cpu)%realdp_hydro(ind,ivar,iflush)=grid(igrid)%unew(ind,ivar)
           end do
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for refining the AMR grid
     !------------------------------------------------------
     if(cache_operation.EQ.operation_refine)then
        if(send_flush_refine(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_refine(grid_cpu),1,new_mpi_large_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_refine(grid_cpu)%nflush=0
        endif
        send_flush_refine(grid_cpu)%nflush=send_flush_refine(grid_cpu)%nflush+1
        iflush=send_flush_refine(grid_cpu)%nflush
        send_flush_refine(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_refine(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef HYDRO
        do ind=1,twotondim
           do ivar=1,nvar
              send_flush_refine(grid_cpu)%realdp_hydro(ind,ivar,iflush)=grid(igrid)%uold(ind,ivar)
           end do
        end do
#endif
#ifdef GRAV
        do ind=1,twotondim
           do idim=1,ndim
              send_flush_refine(grid_cpu)%realdp_poisson(ind,idim,iflush)=grid(igrid)%f(ind,idim)
           end do
           send_flush_refine(grid_cpu)%realdp_poisson(ind,ndim+1,iflush)=grid(igrid)%phi(ind)
           send_flush_refine(grid_cpu)%realdp_poisson(ind,ndim+2,iflush)=grid(igrid)%phi_old(ind)
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for load balancing
     !------------------------------------------------------
     if(cache_operation.EQ.operation_loadbalance)then
        if(send_flush_refine(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_refine(grid_cpu),1,new_mpi_large_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_refine(grid_cpu)%nflush=0
        endif
        send_flush_refine(grid_cpu)%nflush=send_flush_refine(grid_cpu)%nflush+1
        iflush=send_flush_refine(grid_cpu)%nflush
        send_flush_refine(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_refine(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
        do ind=1,twotondim
           if(grid(igrid)%refined(ind))then
              send_flush_refine(grid_cpu)%int4(ind,iflush)=1
           else
              send_flush_refine(grid_cpu)%int4(ind,iflush)=0
           endif
#ifdef HYDRO
           do ivar=1,nvar
              send_flush_refine(grid_cpu)%realdp_hydro(ind,ivar,iflush)=grid(igrid)%uold(ind,ivar)
           end do
#endif
#ifdef GRAV
           do idim=1,ndim
              send_flush_refine(grid_cpu)%realdp_poisson(ind,idim,iflush)=grid(igrid)%f(ind,idim)
           end do
           send_flush_refine(grid_cpu)%realdp_poisson(ind,ndim+1,iflush)=grid(igrid)%phi(ind)
           send_flush_refine(grid_cpu)%realdp_poisson(ind,ndim+2,iflush)=grid(igrid)%phi_old(ind)
#endif
        end do
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for gravity mass deposition
     !------------------------------------------------------
     if(cache_operation.EQ.operation_rho)then
        if(send_flush_poisson(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_poisson(grid_cpu),1,new_mpi_small_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_poisson(grid_cpu)%nflush=0
        endif
        send_flush_poisson(grid_cpu)%nflush=send_flush_poisson(grid_cpu)%nflush+1
        iflush=send_flush_poisson(grid_cpu)%nflush
        send_flush_poisson(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_poisson(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef GRAV
        do ind=1,twotondim
           send_flush_poisson(grid_cpu)%realdp(ind,iflush)=grid(igrid)%rho(ind)
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for MG mask restriction
     !------------------------------------------------------
     if(cache_operation.EQ.operation_restrict_mask)then
        if(send_flush_poisson(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_poisson(grid_cpu),1,new_mpi_small_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_poisson(grid_cpu)%nflush=0
        endif
        send_flush_poisson(grid_cpu)%nflush=send_flush_poisson(grid_cpu)%nflush+1
        iflush=send_flush_poisson(grid_cpu)%nflush
        send_flush_poisson(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_poisson(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef GRAV
        do ind=1,twotondim
           send_flush_poisson(grid_cpu)%realdp(ind,iflush)=grid(igrid)%f(ind,3)
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for MG residual restriction
     !------------------------------------------------------
     if(cache_operation.EQ.operation_restrict_res)then
        if(send_flush_poisson(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_poisson(grid_cpu),1,new_mpi_small_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_poisson(grid_cpu)%nflush=0
        endif
        send_flush_poisson(grid_cpu)%nflush=send_flush_poisson(grid_cpu)%nflush+1
        iflush=send_flush_poisson(grid_cpu)%nflush
        send_flush_poisson(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_poisson(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
#ifdef GRAV
        do ind=1,twotondim
           send_flush_poisson(grid_cpu)%realdp(ind,iflush)=grid(igrid)%f(ind,2)
        end do
#endif
     endif

     !------------------------------------------------------
     ! Filling the flush buffer for MG hierarchy build
     !------------------------------------------------------
     if(cache_operation.EQ.operation_build_mg)then
        if(send_flush_poisson(grid_cpu)%nflush==nflushmax)then
           ! Post send
           call MPI_ISSEND(send_flush_poisson(grid_cpu),1,new_mpi_small_realdp_flush,&
                & grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_poisson(grid_cpu)%nflush=0
        endif
        send_flush_poisson(grid_cpu)%nflush=send_flush_poisson(grid_cpu)%nflush+1
        iflush=send_flush_poisson(grid_cpu)%nflush
        send_flush_poisson(grid_cpu)%lev(iflush)=grid(igrid)%lev
        send_flush_poisson(grid_cpu)%ckey(1:ndim,iflush)=grid(igrid)%ckey(1:ndim)
     endif

  endif
#endif
end subroutine destage
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine close_cache(hash_dict)
  use amr_commons
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  type(hash_table)::hash_dict
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
     if(occupied(icache))call destage(igrid,hash_dict)
     occupied(icache)=.false.
     dirty(icache)=.false.
  end do
  free_cache=1
  ncache=0

  do icache=1,nnull
     if(occupied_null(icache))then
        hash_child(0)=lev_null(icache)
        hash_child(1:ndim)=ckey_null(1:ndim,icache)
        call hash_free(hash_dict,hash_child)
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
           call MPI_ISSEND(send_flush_flag(icpu),1,new_mpi_int4_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_flag(icpu)%nflush=0
        endif
     endif     
     if(cache_operation_type.EQ.operation_type_hydro)then
        if(send_flush_hydro(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_hydro(icpu),1,new_mpi_realdp_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_hydro(icpu)%nflush=0
        endif
     endif
     if(cache_operation_type.EQ.operation_type_poisson)then
        if(send_flush_poisson(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_poisson(icpu),1,new_mpi_small_realdp_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_poisson(icpu)%nflush=0
        endif
     endif
     if(cache_operation_type.EQ.operation_type_refine)then
        if(send_flush_refine(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_refine(icpu),1,new_mpi_large_realdp_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_refine(icpu)%nflush=0
        endif
     endif
     if(cache_operation_type.EQ.operation_type_mg)then
        if(send_flush_mg(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_mg(icpu),1,new_mpi_twin_realdp_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_mg(icpu)%nflush=0
        endif
     endif
     if(cache_operation_type.EQ.operation_type_interpol)then
        if(send_flush_interpol(icpu)%nflush>0)then
           ! Post send
           call MPI_ISSEND(send_flush_interpol(icpu),1,new_mpi_three_realdp_flush,&
                & icpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)  
           ! While waiting for completion, check on incoming messages and perform actions
           call check_mail(send_flush_id,hash_dict)
           send_flush_interpol(icpu)%nflush=0
        endif
     endif
  end do
  
    ! CHECK-IN CHECK-OUT
  if(myid.NE.1)then
     call MPI_ISEND(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(close_id,hash_dict)
     call MPI_IRECV(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(close_id,hash_dict)
  else
     do icpu=2,ncpu
        call MPI_IRECV(dummy_int,1,MPI_INTEGER,&
             & MPI_ANY_SOURCE,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(close_id,hash_dict)
     end do
     do icpu=2,ncpu
        call MPI_ISEND(dummy_int,1,MPI_INTEGER,&
             & icpu-1,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(close_id,hash_dict)
     end do
  endif

  ! Barrier to get the last flush message
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  call check_mail(MPI_REQUEST_NULL,hash_dict)

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
subroutine open_cache(cache_operation_init,domain_decompos_init)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::cache_operation_init
  integer::domain_decompos_init
  !
  integer::info,icpu
  !
#ifndef WITHOUTMPI

  do icpu=1,ncpu
     reply_id(icpu)=MPI_REQUEST_NULL
  end do
  mail_counter=0

  do icpu=1,ncpu
     send_flush_interpol(icpu)%nflush=0
     send_flush_mg(icpu)%nflush=0
     send_flush_flag(icpu)%nflush=0
     send_flush_hydro(icpu)%nflush=0
     send_flush_poisson(icpu)%nflush=0
     send_flush_refine(icpu)%nflush=0
  end do

  ! Domain decomposition to use
  if(domain_decompos_init==domain_decompos_amr)then
     bound_hilbert_key=bound_key_level
     head_cache(levelmin:nlevelmax)=head
     tail_cache(levelmin:nlevelmax)=tail
  endif
  if(domain_decompos_init==domain_decompos_mg )then
     bound_hilbert_key=bound_key_mg
     head_cache(1:nlevelmax)=head_mg
     tail_cache(1:nlevelmax)=tail_mg
  endif

  ! Cache operation to perform
  cache_operation=cache_operation_init

  ! Default operation type
  cache_operation_type=operation_type_flag

  ! Operations of type "flag"
  if(cache_operation.EQ.operation_initflag)cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_smooth)cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_derefine)cache_operation_type=operation_type_flag
  if(cache_operation.EQ.operation_split)cache_operation_type=operation_type_flag

  ! Operations of type "hydro"
  if(cache_operation.EQ.operation_hydro)cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_upload)cache_operation_type=operation_type_hydro
  if(cache_operation.EQ.operation_multipole)cache_operation_type=operation_type_hydro

  ! Operations of type "poisson"
  if(cache_operation.EQ.operation_cg)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_phi)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_rho)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_build_mg)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_restrict_mask)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_restrict_res)cache_operation_type=operation_type_poisson
  if(cache_operation.EQ.operation_scan)cache_operation_type=operation_type_poisson

  ! Operations of type "refine"
  if(cache_operation.EQ.operation_godunov)cache_operation_type=operation_type_refine
  if(cache_operation.EQ.operation_refine)cache_operation_type=operation_type_refine
  if(cache_operation.EQ.operation_loadbalance)cache_operation_type=operation_type_refine

  ! Operations of type "mg"
  if(cache_operation.EQ.operation_mg)cache_operation_type=operation_type_mg

  ! Operations of type "interpol"
  if(cache_operation.EQ.operation_interpol)cache_operation_type=operation_type_interpol
  if(cache_operation.EQ.operation_kick)cache_operation_type=operation_type_interpol

  ! Post the first RECV for request
  call MPI_IRECV(recv_request,1,new_mpi_request,&
       & MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,request_id,info)
  
  ! Post the first RECV for flush
  if(cache_operation_type.EQ.operation_type_flag)then
     call MPI_IRECV(recv_flush_flag,1,new_mpi_int4_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_hydro)then
     call MPI_IRECV(recv_flush_hydro,1,new_mpi_realdp_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_poisson)then
     call MPI_IRECV(recv_flush_poisson,1,new_mpi_small_realdp_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_refine)then
     call MPI_IRECV(recv_flush_refine,1,new_mpi_large_realdp_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_mg)then
     call MPI_IRECV(recv_flush_mg,1,new_mpi_twin_realdp_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif
  if(cache_operation_type.EQ.operation_type_interpol)then
     call MPI_IRECV(recv_flush_interpol,1,new_mpi_three_realdp_flush,&
          & MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,flush_id,info)
  endif

#endif

end subroutine open_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################

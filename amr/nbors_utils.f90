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
  use amr_commons, only: oct
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
  type(oct),pointer::child

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
           
           ! If grid does not exist, send a null reply
           if(.not.ASSOCIATED(child))then
              
              ! Store type corresponding to a null reply
              iskip=mdl%size_fetch_array*(grid_cpu-1)+1
              mdl%send_fetch_array(iskip)=-1

           ! Otherwise, assemble a proper reply with a complete tile
           else
              igrid=(loc(child)-loc(s%m%grid(1)))/(loc(s%m%grid(2))-loc(s%m%grid(1)))+1

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

                 ! Get grid from hash table
                 call c_f_pointer(hash_getp(hash_dict,hash_child),child)
                 if(ASSOCIATED(child))then
                    call unpack_flush%proc(child,mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1))
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
                    call hash_setp(hash_dict,hash_child,m%grid(ichild))

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
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_F_POINTER, C_ASSOCIATED
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
  integer::ind,ivar,idim,info,icache,iflush,grid_cpu
  integer::send_flush_id,iskip,nflush
  integer(kind=8),dimension(0:ndim)::hash_key

#ifndef WITHOUTMPI

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  hash_key(0)=m%grid(igrid)%lev
  hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)

  if(.not.C_ASSOCIATED(hash_getp(hash_dict,hash_key)))then
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
  integer::close_tag=7,close_id
  integer,dimension(1)::dummy_int
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
  use flag_utils, only:pack_fetch_flag, unpack_fetch_flag,&
                      init_flush_initflag, pack_flush_initflag, unpack_flush_initflag
  use load_balance_module, only: pack_flush_loadbalance, unpack_flush_loadbalance
  use refine_utils, only: init_flush_derefine,pack_flush_derefine,unpack_flush_derefine,&
                          pack_flush_refine,unpack_flush_refine,&
                          pack_fetch_refine,unpack_fetch_refine
  use upload_module, only: init_flush_upload,pack_flush_upload,unpack_flush_upload
  use rho_fine_module, only: init_flush_multipole,pack_flush_multipole,unpack_flush_multipole,&
                        init_flush_rho,pack_flush_rho,unpack_flush_rho,pack_fetch_split,unpack_fetch_split
  use move_fine_module, only: pack_fetch_kick,unpack_fetch_kick
  use godunov_fine_module, only: init_flush_godunov,pack_flush_godunov,unpack_flush_godunov
  use phi_fine_cg_module, only: pack_fetch_cg,unpack_fetch_cg,pack_fetch_interpol,unpack_fetch_interpol
  use multigrid_fine_commons, only: pack_flush_build_mg,unpack_flush_build_mg
  use multigrid_fine_coarse, only: init_flush_restrict_mask,pack_flush_restrict_mask,unpack_flush_restrict_mask,&
                                  init_flush_restrict_res,pack_flush_restrict_res,unpack_flush_restrict_res,&
                                  pack_fetch_restrict_res,unpack_fetch_restrict_res,&
                                  pack_fetch_scan,pack_fetch_mg,pack_fetch_phi,&
                                  unpack_fetch_scan,unpack_fetch_mg,unpack_fetch_phi
  use hydro_flag_module, only: pack_fetch_hydro,unpack_fetch_hydro
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

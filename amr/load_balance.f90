!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_load_balance(r,g,m,p,mdl,ilevel)
  use amr_parameters, only: nhilbert
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to load balance the AMR grid
  ! for all levels strictly larger than ilevel.
  !--------------------------------------------------------------------
  integer,dimension(1:g%ncpu)::noct
  integer,allocatable,dimension(:)::input_array
  integer,allocatable,dimension(:)::output_array
  integer(kind=8),dimension(1:nhilbert)::zero_key=0
  integer(kind=8),dimension(1:nhilbert,0:g%ncpu)::bound_key
  integer::ilev,icpu,input_size,output_size
  
  if(g%ncpu==1)return
  if(ilevel==r%nlevelmax)return

  if(r%verbose)write(*,111)ilevel
111 format(' Load balancing for all levels greater than ',I2)

  ! Compute the new domain decomposition
  do ilev=ilevel+1,r%nlevelmax

     if(m%noct_tot(ilev)>0)then
        
        ! Collect number of oct in each cpu for current level
        call r_collect_noct(r,g,m,p,mdl,g%ncpu,1,g%ncpu,ilev,noct)

        ! Compute input array
        input_size=g%ncpu+1
        allocate(input_array(1:input_size))
        input_array(1)=ilev
        input_array(2:input_size)=noct(1:g%ncpu)

        ! Allocate output array
        output_size=2*nhilbert*(g%ncpu+1)
        allocate(output_array(1:output_size))
        
        ! Compute and collect new Hilbert key boundaries for the new domain decomposition
        call r_collect_bound_key(r,g,m,p,mdl,g%ncpu,g%ncpu+1,output_size,input_array,output_array)
        bound_key=reshape(transfer(output_array,zero_key),[nhilbert,g%ncpu+1])
        deallocate(input_array,output_array)

        ! Finalize new domain decomposition
        bound_key(1:nhilbert,0)=zero_key
        do icpu=1,g%ncpu
           if(gt_keys(bound_key(1:nhilbert,icpu-1),bound_key(1:nhilbert,icpu)))then
              bound_key(1:nhilbert,icpu)=bound_key(1:nhilbert,icpu-1)
           endif
        end do
        bound_key(1:nhilbert,g%ncpu)=m%hkey_max(1:nhilbert,ilev)

        ! Scatter new domain decomposition to all processors
        input_size=2*nhilbert*(g%ncpu+1)+1
        allocate(input_array(1:input_size))
        input_array(1)=ilev
        input_array(2:input_size)=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),input_array)
        call r_broadcast_bound_key(r,g,m,p,mdl,g%ncpu,input_size,0,input_array)
        deallocate(input_array)

     endif

  end do
  ! End loop over finer levels
  
  ! Redistribute the grid across CPU according to the new domains
  call r_load_balance(r,g,m,p,mdl,g%ncpu,1,0,ilevel)

end subroutine m_load_balance
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_broadcast_bound_key(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array)
  use amr_parameters, only: nhilbert
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array

  integer::next_range,next_cpu
  integer::ilevel
  integer(kind=8),dimension(1)::dummy

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_BROADCAST_BOUND_KEY,next_cpu,next_range,input_size,output_size,input_array)
     call r_broadcast_bound_key(r,g,m,p,mdl,next_range,input_size,output_size,input_array)
  else
     ilevel=input_array(1)
     m%domain(ilevel)%b=reshape(transfer(input_array(2:input_size),dummy),[nhilbert,g%ncpu+1])
  endif

end subroutine r_broadcast_bound_key
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_collect_bound_key(r,g,m,p,mdl,cpu_range,input_size,output_size,input_array,output_array)
  use amr_parameters, only: nhilbert
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  use hilbert
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::next_range,next_cpu
  integer,dimension(:),allocatable::next_output_array
  
  integer::ilevel
  integer,dimension(:),allocatable::noct
  integer(kind=8),dimension(1)::dummy
  integer(kind=8),dimension(:,:),allocatable::bound_key
  integer(kind=8),dimension(:,:),allocatable::next_bound_key

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_COLLECT_BOUND_KEY,next_cpu,next_range,input_size,output_size,input_array)
     call r_collect_bound_key(r,g,m,p,mdl,next_range,input_size,output_size,input_array,output_array)
     allocate(next_output_array(1:output_size))
     call mdl_get_reply(mdl,next_cpu,output_size,next_output_array)
     allocate(bound_key(1:nhilbert,0:g%ncpu))
     allocate(next_bound_key(1:nhilbert,0:g%ncpu))
     bound_key=reshape(transfer(output_array,dummy),[nhilbert,g%ncpu+1])
     next_bound_key=reshape(transfer(next_output_array,dummy),[nhilbert,g%ncpu+1])
     bound_key=bound_key+next_bound_key
     output_array=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),output_array)
     deallocate(bound_key)
     deallocate(next_bound_key)
     deallocate(next_output_array)
  else
     allocate(noct(1:g%ncpu))
     allocate(bound_key(1:nhilbert,0:g%ncpu))
     bound_key=0
     ilevel=input_array(1)
     noct(1:g%ncpu)=input_array(2:input_size)
     call compute_new_bound_key(r,g,m,ilevel,noct,bound_key)
     output_array=transfer(reshape(bound_key,[nhilbert*(g%ncpu+1)]),output_array)
     deallocate(bound_key)
     deallocate(noct)
  endif

end subroutine r_collect_bound_key
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_new_bound_key(r,g,m,ilevel,noct,bound_key_target)
  use amr_parameters, only: ndim,twotondim,nhilbert,dp
  use amr_commons, only: run_t,global_t,mesh_t,oct
  use hilbert
  use hash
  use cache_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
  integer::info
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  integer,dimension(1:g%ncpu)::noct
  integer(kind=8),dimension(1:nhilbert,0:g%ncpu)::bound_key_target  
  !----------------------------------------------------
  ! This routine compute the new Hilbert keys so that
  ! perfect load balancing is enforced.
  !----------------------------------------------------
  integer::idom,ioct
  integer::nleft,nright,ileft,iright,istart,nstart
  integer,dimension(1:g%ncpu)::noct_cum
  integer,dimension(1:g%ncpu)::ntarget_cum
  integer(kind=8),dimension(1:nhilbert)::one_key
  real(dp)::xtarget

  ! Hilbert key corresponding to one units
  one_key=0
  one_key(1)=1

  ! Compute cumulative grid counts
  noct_cum(1)=noct(1)
  do idom=2,g%ncpu
     noct_cum(idom)=noct_cum(idom-1)+noct(idom)
  end do

  ! Perfect load balancing
  xtarget=dble(noct_cum(g%ncpu))/dble(g%ncpu)

  ! Find left and right domains
  ileft=0
  iright=-1
  bound_key_target=0
  do idom=1,g%ncpu
     ntarget_cum(idom)=int(dble(idom)*xtarget)
     if(g%myid>1)then
        nleft=noct_cum(g%myid-1)
     else
        nleft=0
     endif
     nright=noct_cum(g%myid)
     IF(nright.GT.nleft)then
        if(ntarget_cum(idom).GT.nleft.AND.ntarget_cum(idom).LE.nright)then
           if(ileft==0)ileft=idom
           iright=MAX(idom,iright)
        endif
     endif
  end do

  ! Find corresponding Hilbert keys
  if(iright.GE.ileft)then
     if(g%myid.GT.1)then
        nstart=noct_cum(g%myid-1)
     else
        nstart=0
     endif
     istart=ileft
     do ioct=m%head(ilevel),m%tail(ilevel)
        nstart=nstart+1
        if(nstart.GE.ntarget_cum(istart))then
           bound_key_target(1:nhilbert,istart)=m%grid(ioct)%hkey(1:nhilbert)+one_key
           istart=istart+1
        endif
        if(istart.GT.iright)exit
     end do
  endif

end subroutine compute_new_bound_key
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_load_balance(r,g,m,p,mdl,cpu_range,input_size,output_size,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  use mdl_parameters
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::cpu_range,input_size,output_size
  integer::ilevel

  integer::next_range,next_cpu

  next_range=cpu_range/2
  next_cpu=g%myid+next_range

  if(next_range>0)then
     call mdl_send_request(mdl,MDL_LOAD_BALANCE,next_cpu,next_range,input_size,output_size,ilevel)
     call r_load_balance(r,g,m,p,mdl,next_range,input_size,output_size,ilevel)
  else
     call load_balance(r,g,m,ilevel)
  endif

end subroutine r_load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,nhilbert,dp
  use amr_commons, only: run_t,global_t,mesh_t,oct
  use hilbert
  use hash
  use cache_commons
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::i
  integer::grid_cpu,ichild,idom,jdom,lastdom,domains_matched
  integer::nleft,nright,ileft,iright,istart,nstart,noverlaps
  integer::ioct,ilev

  integer::j,ibit,ibucket,inew
  integer::noct_zero,head_zero,indx_zero
  integer::skip_bit,ikey,true_level
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nhilbert)::coarse_key
  integer(kind=8),dimension(1:nhilbert),save::hks
  integer(kind=8),dimension(1:nhilbert,1:r%nlevelmax)::key_ref
  integer,dimension(1:r%nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  type(oct)::grid_tmp

  !-----------------------------------------------------
  ! Step 1: dispatch octs and empty slots according to
  ! the new target Hilbert tick marks
  !-----------------------------------------------------
  m%ifree=m%noct_used+1
  do ilev=ilevel+1,r%nlevelmax

     call open_cache(r,g,m,operation_loadbalance,domain_decompos_amr)

     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)

        ! Check if grid sits outside future processor boundaries
        if (.not. m%domain(ilev)%in_rank(m%grid(ioct)%hkey)) then
           
           ! Determine the future processor
           hks = m%grid(ioct)%hkey(1:nhilbert)
           grid_cpu = m%domain(ilev)%get_rank(hks)

           ! If next cache line is occupied, free it.
           if(m%occupied(m%free_cache))call destage(r,g,m,r%ngridmax+m%free_cache,m%grid_dict)
           ! Set grid index to a virtual grid in local cache memory
           ichild=r%ngridmax+m%free_cache
           m%occupied(m%free_cache)=.true.
           m%parent_cpu(m%free_cache)=grid_cpu
           m%dirty(m%free_cache)=.true.
           ! Go to next free cache line
           m%free_cache=m%free_cache+1
           m%ncache=m%ncache+1
           if(m%free_cache.GT.r%ncachemax)m%free_cache=1
           if(m%ncache.GT.r%ncachemax)m%ncache=r%ncachemax

           ! Copy all data to the cache grid
           m%grid(ichild)=m%grid(ioct)

           ! Set grid level to zero
           m%grid(ioct)%lev=0
           ! Free grid from hash table
           hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
           call hash_free(m%grid_dict,hash_key)
           
           ! Insert new cache grid in hash table
           call hash_set(m%grid_dict,hash_key,ichild)
        
        endif

     end do

     call close_cache(r,g,m,m%grid_dict)

  end do

  !-----------------------------------------------------
  ! Step 2: sort new octs and empty slots according to
  ! their level (using counting sort algorithm).
  !-----------------------------------------------------
  allocate(noct_level(r%levelmin:r%nlevelmax))
  allocate(head_level(r%levelmin:r%nlevelmax))
  allocate(indx_level(r%levelmin:r%nlevelmax))
  ! Count number of octs per bucket
  noct_level=0
  noct_zero=0

  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level>0)then
        noct_level(true_level)=noct_level(true_level)+1
     else
        noct_zero=noct_zero+1
     end if
  end do
  head_level(ilevel+1)=m%tail(ilevel)+1
  do ilev=ilevel+2,r%nlevelmax
     head_level(ilev)=head_level(ilev-1)+noct_level(ilev-1)
  end do
  head_zero=head_level(r%nlevelmax)+noct_level(r%nlevelmax)

  ! Allocate main swap table
  if(m%ifree.GT.head_level(ilevel+1))then
  allocate(swap_table(head_level(ilevel+1):m%ifree-1))

  ! Build index permutation table
  indx_level=head_level
  indx_zero=head_zero
  do ioct=m%tail(ilevel)+1,m%ifree-1
     true_level=m%grid(ioct)%lev
     if(true_level>0)then
        swap_table(indx_level(true_level))=ioct
        indx_level(true_level)=indx_level(true_level)+1
     else
        swap_table(indx_zero)=ioct
        indx_zero=indx_zero+1
     end if
  end do

  !-----------------------------------------------------
  ! Step 3: sort octs level by level according to their
  ! Hilbert key using LSD Radix Sort algorithm.
  !-----------------------------------------------------
  ! Loop over levels
  do ilev=ilevel+1,r%nlevelmax
     if(noct_level(ilev)>0)then

        ! Allocate temporary swap table just for the level
        allocate(swap_tmp(head_level(ilev):head_level(ilev)+noct_level(ilev)-1))

        ! Loop over useful bits at that level
        do ibit=ilev,1,-1

           ! Get bit and key to read from
           skip_bit=ndim*(ilev-ibit)
           ikey = skip_bit / bits_per_int(ndim) + 1
           skip_bit = mod(skip_bit, bits_per_int(ndim))

           ! Count octs in buckets
           bucket_count=0
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim)
              bucket_count(ibucket)=bucket_count(ibucket)+1
           end do

           ! Compute offsets
           bucket_offset(0)=head_level(ilev)
           do ibucket=1,twotondim-1
              bucket_offset(ibucket)=bucket_offset(ibucket-1)+bucket_count(ibucket-1)
           end do

           ! Sort according to Hilbert key
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              ibucket=ibits(m%grid(ioct)%hkey(ikey),skip_bit,ndim)
              swap_tmp(bucket_offset(ibucket))=ioct
              bucket_offset(ibucket)=bucket_offset(ibucket)+1
           end do

           ! Store permutations in swap table
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              swap_table(inew)=swap_tmp(inew)
           end do
        end do

        ! Deallocate tmp swap array
        deallocate(swap_tmp)
     endif
  end do

  !-----------------------------------------------------
  ! Step 4: Apply permutations directly in main memory
  ! Remember: swap_table(inew)=iold means:
  ! New data at position inew COMES FROM
  ! Old data at position iold.
  !-----------------------------------------------------
  ! Perform the swap
  do j=head_level(ilevel+1),m%ifree-1
     if(j.NE.swap_table(j))then
        hash_key(0)=m%grid(j)%lev
        hash_key(1:ndim)=m%grid(j)%ckey(1:ndim)
        if(m%grid(j)%lev>0)then
           call hash_free(m%grid_dict,hash_key)
        endif
        grid_tmp=m%grid(j)
        i=j
        inew=swap_table(j)
        do while(inew.NE.j)
           m%grid(i)=m%grid(inew)
           hash_key(0)=m%grid(inew)%lev
           hash_key(1:ndim)=m%grid(inew)%ckey(1:ndim)
           if(m%grid(inew)%lev>0)then
              call hash_free(m%grid_dict,hash_key)
              call hash_set(m%grid_dict,hash_key,i)
           endif
           swap_table(i)=i
           i=inew
           inew=swap_table(inew)
        end do
        m%grid(i)=grid_tmp
        hash_key(0)=m%grid(i)%lev
        hash_key(1:ndim)=m%grid(i)%ckey(1:ndim)
        if(m%grid(i)%lev>0)then
           call hash_set(m%grid_dict,hash_key,i)
        end if
        swap_table(i)=i
     endif
  end do
  deallocate(swap_table)
  endif

  !-----------------------------------------------------
  ! Step 5: Clean up final AMR structure
  !-----------------------------------------------------
  do ilev=ilevel+1,r%nlevelmax
     m%head(ilev)=head_level(ilev)
     m%tail(ilev)=head_level(ilev)+noct_level(ilev)-1
     m%noct(ilev)=noct_level(ilev)
  end do
  m%noct_used=m%tail(r%nlevelmax)
  deallocate(noct_level,head_level,indx_level)

  !-----------
  ! Super-octs
  !-----------
  do ilev=1,r%nlevelmax
     npatch(ilev)=twotondim**ilev
  end do
  do ilev=ilevel+1,r%nlevelmax
     n_same=0
     key_ref=0
     key_ref(1,1:r%nlevelmax)=-1
     do ioct=m%head(ilev),m%tail(ilev)
        m%grid(ioct)%superoct=1
        coarse_key(1:nhilbert)=m%grid(ioct)%hkey(1:nhilbert)
        do i=1,MIN(ilev-1,r%nsuperoct)
           coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
           if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
              n_same(i)=n_same(i)+1
           else
              n_same(i)=1
              key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
           endif
           if(n_same(i).EQ.npatch(i))then
              m%grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
           endif
        end do
     end do
  end do

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifndef WITHOUTMPI
subroutine balance_part(r,g,m,p,ilevel)
  use amr_parameters, only: nhilbert,nvector,ndim,i8b,dp
  use pm_parameters, only: part_memory
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use domain_m, only: domain_t
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  integer::ilevel
  !---------------------------------------------------------------------
  ! This routine will dispatch particles across processors according to
  ! their Hilbert key and for a given domain decomposition.
  ! It assumes that particles are sorted according to their levels
  ! and within the level, according to their Hilbert key.
  ! It can be used only if routine rho has been called once before.
  !---------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(MPI_STATUS_SIZE,g%ncpu)::statuses
#endif
  integer(kind=8), dimension(1:nvector,1:nhilbert),save::hk_ref
  integer(kind=8), dimension(1:nvector,1:ndim),save::ix_ref
  integer(kind=4), dimension(1:nvector),save::dummy_state

  integer,dimension(1:ndim),save::ix
  integer::i,istart,ipart,jpart,idim,grid_cpu
  integer::ilev,idom,icpu,mydom,ndom,count_loc,recv_cnt_tot,send_cnt_tot
  integer::nbuffer,countrecv,countsend,tag=101
  real(kind=8)::dx_loc

  integer,allocatable,dimension(:)::send_cnt,recv_cnt
  integer,allocatable,dimension(:)::send_oft,recv_oft,offset_cpu
  integer,dimension(g%ncpu)::reqsend,reqrecv

  real(kind=8),dimension(:),allocatable::x_recv_buf,x_send_buf
  integer(i8b),dimension(:),allocatable::l_recv_buf,l_send_buf
  integer,dimension(:),allocatable::i_recv_buf,i_send_buf

  integer(kind=8)::unbalance
  integer(kind=8),dimension(1:nhilbert)::diff_key
  integer(kind=8),dimension(1:nhilbert),save::hks
  type(domain_t),allocatable,dimension(:)::domain_part
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target,bound_key_new
  integer(kind=8),allocatable,dimension(:,:)::bound_key_left,bound_key_right
  integer::npart_lev,npart_lev_tot,iter
  integer,dimension(0:g%ncpu)::npart_cum,npart_cum_tot
  integer,dimension(1:g%ncpu)::npart_dom,npart_dom_tot
  real(dp)::xpart_target,xcum_target

  real(dp),dimension(1:ndim),save::xp_tmp,vp_tmp
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

#ifndef WITHOUTMPI
  if(g%ncpu==1)return
  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel

  !####################################################
  ! Default for particle domains are grid domains
  !####################################################
  allocate(domain_part(ilevel:r%nlevelmax+1))
  do ilev=ilevel,r%nlevelmax+1
     call domain_part(ilev)%copy(m%domain(ilev))
  enddo

  !###############################################
  ! Determine particle domains if needed
  !###############################################
  if(part_memory)then
     
     !#############################
     ! Allocate temporary work space
     !#############################
     ndom=maxval(domain_part%n)
     allocate(bound_key_target(1:nhilbert,0:ndom))
     allocate(bound_key_new(1:nhilbert,0:ndom))
     allocate(bound_key_left(1:nhilbert,0:ndom))
     allocate(bound_key_right(1:nhilbert,0:ndom))
     
     ! Loop over levels
     do ilev=ilevel,r%nlevelmax
        mydom = domain_part(ilev)%r2d(g%myid)
        ndom  = domain_part(ilev)%n

        dx_loc=r%boxlen/2**ilev
     
        ! Compute number of particles
        npart_lev=p%tailp(ilev)-p%headp(ilev)+1
        npart_lev_tot=0
        call MPI_ALLREDUCE(npart_lev,npart_lev_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        if(npart_lev_tot.EQ.0)cycle
        xpart_target=dble(npart_lev_tot)/dble(ndom)

        npart_dom=0
        npart_dom(mydom)=npart_lev
        call MPI_ALLREDUCE(npart_dom,npart_dom_tot,ndom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        npart_dom=npart_dom_tot

        npart_cum=0
        do idom=1,ndom
           npart_cum(idom)=npart_cum(idom-1)+npart_dom(idom)
        end do

        iter=0
        if(g%myid==1.and.r%verbose)write(*,*)"====================================="
        if(g%myid==1.and.r%verbose)write(*,'("Level=",I4," npart=",I10)')ilev,npart_lev_tot

        !#########################################################
        ! Sort particle according to current level Hilbert key
        !#########################################################
        do i=p%headp(ilev),p%tailp(ilev)
           p%sortp(i)=i
        end do
        ix=0
        call sort_hilbert(r,g,p,p%headp(ilev),p%tailp(ilev),ix,0,1,ilev-1)

        ! Compute first guess domain decomposition
        bound_key_target(1:nhilbert,0:ndom)=domain_part(ilev)%b(1:nhilbert,0:ndom)
        bound_key_new=bound_key_target
        bound_key_left=0
        do idom=0,ndom
           bound_key_right(1:nhilbert,idom)=m%hkey_max(1:nhilbert,ilev)
        end do

        unbalance=10
        iter=0

        !#########################################################
        ! Find new Hilbert tick marks by dichotomy
        !#########################################################
        do while (unbalance.GT.1)
           iter=iter+1
           
           ! Compute number of particles above tick marks
           npart_cum=0

           ! Loop over particles in Hilbert order
           do i=p%headp(ilev),p%tailp(ilev)
              ipart=p%sortp(i)

              ! Compute Hilbert key of particle parent grid
              ix_ref(1,1:ndim) = int(p%xp(ipart,1:ndim)/(2*dx_loc))
              call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)
              
              do idom=1,ndom
                 if(gt_keys(bound_key_target(1:nhilbert,idom),hk_ref(1,1:nhilbert)))then
                    npart_cum(idom)=npart_cum(idom)+1
                 end if
              enddo
           end do

           ! Compute global histogram
           call MPI_ALLREDUCE(npart_cum,npart_cum_tot,ndom+1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
           npart_cum=npart_cum_tot

           unbalance=0
           do idom=1,ndom-1
              xcum_target=dble(idom)*xpart_target
              if(npart_cum(idom)>xcum_target)then
                 bound_key_new(1:nhilbert,idom)=average_keys(bound_key_left(1:nhilbert,idom),bound_key_target(1:nhilbert,idom))
                 bound_key_right(1:nhilbert,idom)=bound_key_target(1:nhilbert,idom)
              else
                 bound_key_new(1:nhilbert,idom)=average_keys(bound_key_right(1:nhilbert,idom),bound_key_target(1:nhilbert,idom))
                 bound_key_left(1:nhilbert,idom)=bound_key_target(1:nhilbert,idom)
              endif
              diff_key=difference_keys(bound_key_right(1:nhilbert,idom),bound_key_left(1:nhilbert,idom))
#if NHILBERT==1
              unbalance=MAX(unbalance,ABS(diff_key(1)))
#endif
#if NHILBERT==2
              unbalance=MAX(unbalance,ABS(diff_key(1))+1000*ABS(diff_key(2)))
#endif
#if NHILBERT==3
              unbalance=MAX(unbalance,ABS(diff_key(1))+1000*ABS(diff_key(2))+1000*ABS(diff_key(3)))
#endif
           end do
                            
           bound_key_target=bound_key_new
           
        end do
        if(g%myid==1.and.r%verbose)write(*,'("iter=",I4,1X,17(I10,1X))')iter,npart_cum
        if(g%myid==1.and.r%verbose)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(idom)*xpart_target),idom=0,ndom)
        !#########################################################
        ! Store new Hilbert tick marks after convergence
        !#########################################################
        domain_part(ilev)%b(1:nhilbert,0:ndom)=bound_key_target(1:nhilbert,0:ndom)

     end do
     if(g%myid==1.and.r%verbose)write(*,*)"====================================="

     !#############################
     ! Deallocate work space
     !#############################
     deallocate(bound_key_target)
     deallocate(bound_key_new)
     deallocate(bound_key_left)
     deallocate(bound_key_right)
     
  endif

  !#################################
  ! Balance particles across cpus
  !#################################

  !#############################
  ! Allocate work space
  !#############################
  allocate(send_cnt(1:g%ncpu))
  allocate(recv_cnt(1:g%ncpu))
  allocate(recv_oft(1:g%ncpu))
  allocate(send_oft(1:g%ncpu))
  allocate(offset_cpu(1:g%ncpu))

  !#####################################
  ! Compute number of particles to send
  !#####################################
  send_cnt=0
  count_loc=0

  ! Loop over levels
  do ilev=ilevel,r%nlevelmax
     ix_ref=-1
     dx_loc=r%boxlen/2**ilev

     ! Loop over particles
     do ipart=p%headp(ilev),p%tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int(p%xp(ipart,1:ndim)/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1,1:ndim)))then
           ix_ref(1,1:ndim)=ix(1:ndim)
           grid_cpu = g%myid

           ! Compute Hilbert key of particle parent grid
           call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)

           ! Check if grid sits outside future processor boundaries
           if (.not. domain_part(ilev)%in_rank(hk_ref(1,1:nhilbert))) then
              ! Determine the future processor
              hks = hk_ref(1,1:nhilbert)
              grid_cpu = domain_part(ilev)%get_rank(hks)
           endif
        endif

        if(grid_cpu.EQ.g%myid)then
           ! Update count for local cpu
           count_loc=count_loc+1
        else
           ! Update count for remote cpu
           send_cnt(grid_cpu)=send_cnt(grid_cpu)+1
        end if

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  !#####################################
  ! Compute number of particles to receive
  !#####################################
  call MPI_ALLTOALL(send_cnt(1),1,MPI_INTEGER,recv_cnt(1),1,MPI_INTEGER,MPI_COMM_WORLD,info)
  send_cnt_tot=SUM(send_cnt)
  recv_cnt_tot=SUM(recv_cnt)

  !#####################################
  ! Compute offsets
  !#####################################
  recv_oft=0; send_oft=0
  offset_cpu(1)=p%headp(ilevel)-1+count_loc
  do icpu=2,g%ncpu
     offset_cpu(icpu)=offset_cpu(icpu-1)+send_cnt(icpu-1)
     recv_oft(icpu)=recv_oft(icpu-1)+recv_cnt(icpu-1)
     send_oft(icpu)=send_oft(icpu-1)+send_cnt(icpu-1)
  end do

  !#####################################
  ! Shift to the right particles to send
  !#####################################
  send_cnt=0
  count_loc=0

  ! Loop over levels
  do ilev=ilevel,r%nlevelmax
     ix_ref=-1
     dx_loc=r%boxlen/2**ilev

     ! Loop over particles
     do ipart=p%headp(ilev),p%tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int(p%xp(ipart,1:ndim)/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1,1:ndim)))then
           ix_ref(1,1:ndim)=ix(1:ndim)
           grid_cpu=g%myid

           ! Compute Hilbert key of particle parent grid
           call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)

           ! Check if grid sits outside future processor boundaries
           if (.not. domain_part(ilev)%in_rank(hk_ref(1,1:nhilbert))) then
              ! Determine the future processor
              hks = hk_ref(1,1:nhilbert)
              grid_cpu = domain_part(ilev)%get_rank(hks)
           endif
        endif

        if(grid_cpu.EQ.g%myid)then
           ! Update count for local cpu
           count_loc=count_loc+1
           p%workp(ipart)=p%headp(ilevel)-1+count_loc
        else
           ! Update count for remote cpu
           send_cnt(grid_cpu)=send_cnt(grid_cpu)+1
           p%workp(ipart)=offset_cpu(grid_cpu)+send_cnt(grid_cpu)
        end if

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  !#####################################
  ! Swap particles using new index table
  !#####################################
  do ipart=p%headp(ilevel),p%tailp(r%nlevelmax)
     do while(p%workp(ipart).NE.ipart)
        ! Swap new index
        jpart=p%workp(ipart)
        p%workp(ipart)=p%workp(jpart)
        p%workp(jpart)=jpart
        ! Swap positions
        xp_tmp(1:ndim)=p%xp(ipart,1:ndim)
        p%xp(ipart,1:ndim)=p%xp(jpart,1:ndim)
        p%xp(jpart,1:ndim)=xp_tmp(1:ndim)
        ! Swap velocities
        vp_tmp(1:ndim)=p%vp(ipart,1:ndim)
        p%vp(ipart,1:ndim)=p%vp(jpart,1:ndim)
        p%vp(jpart,1:ndim)=vp_tmp(1:ndim)
        ! Swap masses
        mp_tmp=p%mp(ipart)
        p%mp(ipart)=p%mp(jpart)
        p%mp(jpart)=mp_tmp
        ! Swap ids
        idp_tmp=p%idp(ipart)
        p%idp(ipart)=p%idp(jpart)
        p%idp(jpart)=idp_tmp
        ! Swap levels
        levelp_tmp=p%levelp(ipart)
        p%levelp(ipart)=p%levelp(jpart)
        p%levelp(jpart)=levelp_tmp
     end do
  end do

  !###################################################################
  ! Set new number of particles in local processor
  !###################################################################
  p%npart=p%headp(ilevel)-1+count_loc+recv_cnt_tot
  p%tailp(r%nlevelmax)=p%npart

  !###################################################################
  ! Swap particles positions, velocities and masses between processors
  !###################################################################
  allocate(x_recv_buf(1:recv_cnt_tot))
  allocate(x_send_buf(1:send_cnt_tot))

  !#########################
  ! Swap positions
  !#########################
  do idim=1,ndim

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%xp(ipart,idim)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%xp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !#########################
  ! Swap velocities
  !#########################
  do idim=1,ndim

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%vp(ipart,idim)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%vp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !#########################
  ! Swap masses
  !#########################
  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     x_send_buf(i)=p%mp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(x_send_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%mp(ipart)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  deallocate(x_recv_buf,x_send_buf)

  !#########################
  ! Swap levels
  !#########################
  allocate(i_recv_buf(1:recv_cnt_tot))
  allocate(i_send_buf(1:send_cnt_tot))

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(i_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     i_send_buf(i)=p%levelp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
        call MPI_ISEND(i_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%levelp(ipart)=i_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  deallocate(i_send_buf,i_recv_buf)

  !#########################
  ! Swap ids
  !#########################
  allocate(l_recv_buf(1:recv_cnt_tot))
  allocate(l_send_buf(1:send_cnt_tot))

  countrecv=0
  do icpu=1,g%ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
#ifndef LONGINT
        call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
        call MPI_IRECV(l_recv_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
     endif
  end do

  do i=1,send_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     l_send_buf(i)=p%idp(ipart)
  end do

  countsend=0
  do icpu=1,g%ncpu
     nbuffer=send_cnt(icpu)
     if(nbuffer>0) then
        countsend=countsend+1
        istart=send_oft(icpu)+1
#ifndef LONGINT
        call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(l_send_buf(istart),nbuffer,MPI_INTEGER8,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  do i=1,recv_cnt_tot
     ipart=p%headp(ilevel)-1+count_loc+i
     p%idp(ipart)=l_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  deallocate(l_recv_buf,l_send_buf)

  !#############################
  ! Deallocate work space
  !#############################
  deallocate(send_cnt)
  deallocate(recv_cnt)
  deallocate(recv_oft)
  deallocate(send_oft)
  deallocate(offset_cpu)
  do ilev=ilevel,r%nlevelmax+1
    call domain_part(ilev)%destroy
  end do
  deallocate(domain_part)

  !##################################
  ! Put all particles in level ilevel
  !##################################
  p%tailp(ilevel)=p%npart
  do ilev=ilevel+1,r%nlevelmax
     p%headp(ilev)=p%npart+1
     p%tailp(ilev)=p%npart
  end do

  call MPI_ALLREDUCE(p%npart,p%npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)

#endif

111 format('   Entering balance_part for level',i2)

end subroutine balance_part
#endif
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::ilevel
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::igrid,i,ind,jlevel,info
  integer::icpu,grid_cpu,ichild
  integer::nleft,nright,ileft,iright,istart,nstart
  integer::ilev,ioct
  integer,dimension(:),allocatable::noct_cpu,noct_cum
  integer,dimension(:),allocatable::ntarget_cum
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target_tot
  real(dp)::xtarget

  integer::icell,j,ibit,ibucket,inew,iold,iold_true
  integer::noct_zero,head_zero,indx_zero
  integer::ncreate_tot,nkill_tot
  integer::parent_cell,skip_bit,ikey,true_level
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key,zero_key
  integer(kind=8),dimension(1:nhilbert,1:nlevelmax)::key_ref
  integer,dimension(1:nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok_free,ok_all,ok
  type(oct)::grid_tmp
  type(oct_hydro)::fluid_tmp
  type(oct_grav)::grav_tmp

#ifndef WITHOUTMPI
  if(ncpu==1)return
  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel
  
  ! Constants
  zero_key=0
  one_key=0
  one_key(1)=1

  !-----------------------------------------------------
  ! Step 1: determine the new Hilbert tick marks
  !-----------------------------------------------------
  allocate(noct_cpu(1:ncpu))
  allocate(noct_cum(1:ncpu))
  allocate(ntarget_cum(1:ncpu))
  allocate(bound_key_target(1:nhilbert,0:ncpu))
  allocate(bound_key_target_tot(1:nhilbert,0:ncpu))

  ! Compute new Hilbert tick marks
  do ilev=ilevel+1,nlevelmax

     if(noct_tot(ilev)>0)then
        
        noct_cpu=0
        noct_cpu(myid)=noct(ilev)
        call MPI_ALLREDUCE(noct_cpu,noct_cum,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        noct_cpu=noct_cum
        do icpu=2,ncpu
           noct_cum(icpu)=noct_cum(icpu-1)+noct_cpu(icpu)
        end do
        
        xtarget=dble(noct_cum(ncpu))/dble(ncpu)
        
        ileft=0
        iright=-1
        bound_key_target=0
        do icpu=1,ncpu
           ntarget_cum(icpu)=int(dble(icpu)*xtarget)
           if(myid>1)then
              nleft=noct_cum(myid-1)
           else
              nleft=0
           endif
           nright=noct_cum(myid)
           IF(nright.GT.nleft)then
              if(ntarget_cum(icpu).GT.nleft.AND.ntarget_cum(icpu).LE.nright)then
                 if(ileft==0)ileft=icpu
                 iright=MAX(icpu,iright)
              endif
           endif
        end do
        
        if(iright.GE.ileft)then
           if(myid.GT.1)then
              nstart=noct_cum(myid-1)
           else
              nstart=0
           endif
           istart=ileft
           do ioct=head(ilev),tail(ilev)
              nstart=nstart+1
              if(nstart.GE.ntarget_cum(istart))then
                 bound_key_target(1:nhilbert,istart)=grid(ioct)%hkey(1:nhilbert)+one_key
                 istart=istart+1
              endif
              if(istart.GT.iright)exit
           end do
        endif
        
        call MPI_ALLREDUCE(bound_key_target,bound_key_target_tot,nhilbert*(ncpu+1),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
        bound_key_target=bound_key_target_tot
        
        bound_key_target(1:nhilbert,0)=zero_key
        do icpu=1,ncpu
           if(gt_keys(bound_key_target(1:nhilbert,icpu-1),bound_key_target(1:nhilbert,icpu)))then
              bound_key_target(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu-1)
           endif
        end do
        bound_key_target(1:nhilbert,ncpu)=hkey_max(1:nhilbert,ilev)
        do icpu=0,ncpu
           bound_key_level(1:nhilbert,icpu,ilev)=bound_key_target(1:nhilbert,icpu)
        end do

     endif

  end do
  deallocate(noct_cpu,noct_cum,ntarget_cum)
  deallocate(bound_key_target)
  deallocate(bound_key_target_tot)

  !-----------------------------------------------------
  ! Step 2: dispatch octs and empty slots according to
  ! the new target Hilbert tick marks
  !-----------------------------------------------------
  ifree=noct_used+1
  do ilev=ilevel+1,nlevelmax

     call open_cache(operation_loadbalance,domain_decompos_amr)

     hash_key(0)=ilev
     do ioct=head(ilev),tail(ilev)

        ! Check if grid sits outside future processor boundaries
        if(    gt_keys(bound_key_level(1:nhilbert,myid-1,ilev),grid(ioct)%hkey(1:nhilbert)).OR. &
             & ge_keys(grid(ioct)%hkey(1:nhilbert),bound_key_level(1:nhilbert,myid,ilev)))then
           
           ! Determine the future processor
           do icpu=1,ncpu
              if(    ge_keys(grid(ioct)%hkey(1:nhilbert),bound_key_level(1:nhilbert,icpu-1,ilev)).AND. &
                   & gt_keys(bound_key_level(1:nhilbert,icpu,ilev),grid(ioct)%hkey(1:nhilbert)))then
                 grid_cpu=icpu
              end if
           end do

           ! If next cache line is occupied, free it.
           if(occupied(free_cache))call destage(ngridmax+free_cache,grid_dict)
           ! Set grid index to a virtual grid in local cache memory
           ichild=ngridmax+free_cache
           occupied(free_cache)=.true.
           parent_cpu(free_cache)=grid_cpu
           dirty(free_cache)=.true.
           ! Go to next free cache line
           free_cache=free_cache+1
           ncache=ncache+1
           if(free_cache.GT.ncachemax)free_cache=1
           if(ncache.GT.ncachemax)ncache=ncachemax

           ! Copy all data to the cache grid
           grid(ichild)=grid(ioct)

           ! Set grid level to zero
           grid(ioct)%lev=0
           ! Free grid from hash table
           hash_key(1:ndim)=grid(ioct)%ckey(1:ndim)
           call hash_free(grid_dict,hash_key)
           
           ! Insert new cache grid in hash table
           call hash_set(grid_dict,hash_key,ichild)
        
        endif

     end do

     call close_cache(grid_dict)

  end do

  !-----------------------------------------------------
  ! Step 3: sort new octs and empty slots according to
  ! their level (using counting sort algorithm).
  !-----------------------------------------------------
  allocate(noct_level(levelmin:nlevelmax))
  allocate(head_level(levelmin:nlevelmax))
  allocate(indx_level(levelmin:nlevelmax))
  ! Count number of octs per bucket
  noct_level=0
  noct_zero=0

  do ioct=tail(ilevel)+1,ifree-1
     true_level=grid(ioct)%lev
     if(true_level>0)then
        noct_level(true_level)=noct_level(true_level)+1
     else
        noct_zero=noct_zero+1
     end if
  end do
  head_level(ilevel+1)=tail(ilevel)+1
  do ilev=ilevel+2,nlevelmax
     head_level(ilev)=head_level(ilev-1)+noct_level(ilev-1)
  end do
  head_zero=head_level(nlevelmax)+noct_level(nlevelmax)

  ! Allocate main swap table
  if(ifree.GT.head_level(ilevel+1))then
  allocate(swap_table(head_level(ilevel+1):ifree-1))

  ! Build index permutation table
  indx_level=head_level
  indx_zero=head_zero
  do ioct=tail(ilevel)+1,ifree-1
     true_level=grid(ioct)%lev
     if(true_level>0)then
        swap_table(indx_level(true_level))=ioct
        indx_level(true_level)=indx_level(true_level)+1
     else
        swap_table(indx_zero)=ioct
        indx_zero=indx_zero+1
     end if
  end do

  !-----------------------------------------------------
  ! Step 4: sort octs level by level according to their
  ! Hilbert key using LSD Radix Sort algorithm.
  !-----------------------------------------------------
  ! Loop over levels
  do ilev=ilevel+1,nlevelmax
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
!!$              if(    gt_keys(bound_key_level(1:nhilbert,myid-1,ilev),grid(ioct)%hkey(1:nhilbert)).OR. &
!!$                   & ge_keys(grid(ioct)%hkey(1:nhilbert),bound_key_level(1:nhilbert,myid,ilev)))then
!!$                 write(*,*)'PE ',myid,'######### ',ioct
!!$                 write(*,*)grid(ioct)%hkey
!!$                 write(*,*)grid(ioct)%lev
!!$                 write(*,*)grid(ioct)%ckey
!!$                 write(*,*)bound_key_level(1:nhilbert,myid-1,ilev)
!!$                 write(*,*)bound_key_level(1:nhilbert,myid  ,ilev)
!!$                 stop
!!$              endif
              ibucket=ibits(grid(ioct)%hkey(ikey),skip_bit,ndim)
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
              ibucket=ibits(grid(ioct)%hkey(ikey),skip_bit,ndim)
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
  ! Step 5: Apply permutations directly in main memory
  ! Remember: swap_table(inew)=iold means:
  ! New data at position inew COMES FROM
  ! Old data at position iold.
  !-----------------------------------------------------
  ! Perform the swap
  do j=head_level(ilevel+1),ifree-1
     if(j.NE.swap_table(j))then
        hash_key(0)=grid(j)%lev
        hash_key(1:ndim)=grid(j)%ckey(1:ndim)
        if(grid(j)%lev>0)then
           call hash_free(grid_dict,hash_key)
        endif
        grid_tmp=grid(j)
        if(hydro)then
           fluid_tmp=fluid(j)
        endif
        if(poisson)then
           grav_tmp=grav(j)
        endif
        i=j
        inew=swap_table(j)
        do while(inew.NE.j)
           grid(i)=grid(inew)
           hash_key(0)=grid(inew)%lev
           hash_key(1:ndim)=grid(inew)%ckey(1:ndim)
           if(grid(inew)%lev>0)then
              call hash_free(grid_dict,hash_key)
              call hash_set(grid_dict,hash_key,i)
           endif
           swap_table(i)=i
           i=inew
           inew=swap_table(inew)
        end do
        grid(i)=grid_tmp
        if(hydro)then
           fluid(i)=fluid_tmp
        endif
        if(poisson)then
           grav(i)=grav_tmp
        endif
        hash_key(0)=grid(i)%lev
        hash_key(1:ndim)=grid(i)%ckey(1:ndim)
        if(grid(i)%lev>0)then
           call hash_set(grid_dict,hash_key,i)
        end if
        swap_table(i)=i
     endif
  end do
  deallocate(swap_table)
  endif

  !-----------------------------------------------------
  ! Step 6: Clean up final AMR structure
  !-----------------------------------------------------
  do ilev=ilevel+1,nlevelmax
     head(ilev)=head_level(ilev)
     tail(ilev)=head_level(ilev)+noct_level(ilev)-1
     noct(ilev)=noct_level(ilev)
  end do
  noct_used=tail(nlevelmax)
  deallocate(noct_level,head_level,indx_level)

  !-----------
  ! Super-octs
  !-----------
  do ilev=1,nlevelmax
     npatch(ilev)=twotondim**ilev
  end do
  do ilev=ilevel+1,nlevelmax
     n_same=0
     key_ref=0
     key_ref(1,1:nlevelmax)=-1
     do ioct=head(ilev),tail(ilev)
        grid(ioct)%superoct=1
        coarse_key(1:nhilbert)=grid(ioct)%hkey(1:nhilbert)
        do i=1,MIN(ilev-1,nsuperoct)
           coarse_key(1:nhilbert)=coarsen_key(coarse_key(1:nhilbert),ilev-1) ! ilev-1 used to speed up only
           if(eq_keys(coarse_key(1:nhilbert),key_ref(1:nhilbert,i)))then
              n_same(i)=n_same(i)+1
           else
              n_same(i)=1
              key_ref(1:nhilbert,i)=coarse_key(1:nhilbert)
           endif
           if(n_same(i).EQ.npatch(i))then
              grid(ioct-npatch(i)+1:ioct)%superoct=npatch(i)
           endif
        end do
     end do
  end do

  !---------------------
  ! Total number of octs
  !---------------------
  do ilev=ilevel+1,nlevelmax
     noct_tot(ilev)=noct(ilev)
     noct_min(ilev)=noct(ilev)
     noct_max(ilev)=noct(ilev)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(noct(ilev),noct_tot(ilev),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(noct(ilev),noct_min(ilev),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(noct(ilev),noct_max(ilev),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  end do

  noct_used_max=noct_used
  noct_used_tot=noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(noct_used,noct_used_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(noct_used,noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

#endif

111 format(' Load balancing for all levels greater than ',I2)
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine balance_part(ilevel)
  use amr_commons
  use pm_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  !
  ! This routine will dispatch particles across processors according to
  ! their Hilbert key and for a given domain decomposition.
  ! It assumes that particles are sorted according to their levels
  ! and within the level, according to their Hilbert key.
  ! It can only be called after routine rho has been called.
  !
  integer(kind=8), dimension(1:nvector,1:nhilbert),save::hk_ref
  integer(kind=8), dimension(1:nvector,1:ndim),save::ix_ref
  integer(kind=4), dimension(1:nvector),save::dummy_state

  integer,dimension(1:ndim),save::ix
  integer::i,istart,info,ipart,jpart,idim,grid_cpu
  integer::ilev,icpu,jcpu,count_loc,recv_cnt_tot,send_cnt_tot
  integer::nbuffer,countrecv,countsend,tag=101
  real(kind=8)::dx_loc

  integer,allocatable,dimension(:)::send_cnt,recv_cnt
  integer,allocatable,dimension(:)::send_oft,recv_oft,offset_cpu
  integer,dimension(ncpu)::reqsend,reqrecv

  real(kind=8),dimension(:),allocatable::x_recv_buf,x_send_buf
  integer(i8b),dimension(:),allocatable::l_recv_buf,l_send_buf
  integer,dimension(:),allocatable::i_recv_buf,i_send_buf

  integer(kind=8)::unbalance
  integer(kind=8),allocatable,dimension(:,:,:)::bound_key_part
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target,bound_key_new
  integer(kind=8),allocatable,dimension(:,:)::bound_key_left,bound_key_right
  integer(kind=8),dimension(1:nhilbert)::diff_key,ave_key,one_key
  integer,dimension(1:ncpu)::npart_proc,npart_proc_tot
  integer::npart_lev,npart_lev_tot,iter
  integer,dimension(0:ncpu)::npart_cum,npart_cum_tot
  integer,dimension(1:ncpu)::npart_cpu,npart_cpu_tot
  real(dp)::xpart_target,xcum_target

  real(dp),dimension(1:ndim),save::xp_tmp,vp_tmp
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

#ifndef WITHOUTMPI
  if(ncpu==1)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel
  one_key=0
  one_key(1)=1

  !####################################################
  ! Default for particle domains are grid domains
  !####################################################
  allocate(bound_key_part(1:nhilbert,0:ncpu,1:nlevelmax+1))
  bound_key_part=bound_key_level

  !###############################################
  ! Determine particle domains if needed
  !###############################################
  if(part_memory)then
     
     !#############################
     ! Allocate temporary work space
     !#############################
     allocate(bound_key_target(1:nhilbert,0:ncpu))
     allocate(bound_key_new(1:nhilbert,0:ncpu))
     allocate(bound_key_left(1:nhilbert,0:ncpu))
     allocate(bound_key_right(1:nhilbert,0:ncpu))
     
     ! Loop over levels
     do ilev=ilevel,nlevelmax
        dx_loc=boxlen/2**ilev
     
        ! Compute number of particles
        npart_lev=tailp(ilev)-headp(ilev)+1
        npart_lev_tot=0
        call MPI_ALLREDUCE(npart_lev,npart_lev_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        if(npart_lev_tot.EQ.0)cycle
        xpart_target=dble(npart_lev_tot)/dble(ncpu)

        npart_cpu=0
        npart_cpu(myid)=npart_lev
        call MPI_ALLREDUCE(npart_cpu,npart_cpu_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        npart_cpu=npart_cpu_tot

        npart_cum=0
        do icpu=1,ncpu
           npart_cum(icpu)=npart_cum(icpu-1)+npart_cpu(icpu)
        end do

        iter=0
        if(myid==1)write(*,*)"====================================="
        if(myid==1)write(*,'("Level=",I4," npart=",I10)')ilev,npart_lev_tot
!        if(myid==1)write(*,'(16(I10,1X))')npart_cpu
!        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(icpu)*xpart_target),icpu=0,ncpu)
!        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,npart_cum

        !#########################################################
        ! Sort particle according to current level Hilbert key
        !#########################################################
        do i=headp(ilev),tailp(ilev)
           sortp(i)=i
        end do
        ix=0
        call sort_hilbert(headp(ilev),tailp(ilev),ix,0,1,ilev-1)

        ! Compute first guess domain decomposition
        bound_key_target(1:nhilbert,0:ncpu)=bound_key_part(1:nhilbert,0:ncpu,ilev)
        bound_key_new=bound_key_target
        bound_key_left=0
        do icpu=0,ncpu
           bound_key_right(1:nhilbert,icpu)=hkey_max(1:nhilbert,ilev)
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
           do i=headp(ilev),tailp(ilev)
              ipart=sortp(i)

              ! Compute Hilbert key of particle parent grid
              ix_ref(1,1:ndim) = int(xp(ipart,1:ndim)/(2*dx_loc))
              call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)
              
              ! Determine in which processor the grid sits
              do icpu=1,ncpu
                 if(gt_keys(bound_key_target(1:nhilbert,icpu),hk_ref(1,1:nhilbert)))then
                    npart_cum(icpu)=npart_cum(icpu)+1
                 end if
              end do

           end do

           ! Compute global histogram
           call MPI_ALLREDUCE(npart_cum,npart_cum_tot,ncpu+1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
           npart_cum=npart_cum_tot

           unbalance=0
           do icpu=1,ncpu-1
              xcum_target=dble(icpu)*xpart_target
              if(npart_cum(icpu)>xcum_target)then
                 ave_key=average_keys(bound_key_left(1:nhilbert,icpu),bound_key_target(1:nhilbert,icpu))
                 bound_key_new(1:nhilbert,icpu)=ave_key
                 bound_key_right(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu)
              else
                 ave_key=average_keys(bound_key_right(1:nhilbert,icpu),bound_key_target(1:nhilbert,icpu))
                 bound_key_new(1:nhilbert,icpu)=ave_key
                 bound_key_left(1:nhilbert,icpu)=bound_key_target(1:nhilbert,icpu)
              endif
              diff_key=difference_keys(bound_key_right(1,icpu),bound_key_left(1,icpu))
              unbalance=MAX(unbalance,maxval(diff_key))
           end do
                            
           bound_key_target=bound_key_new
           
        end do
        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,npart_cum
        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(icpu)*xpart_target),icpu=0,ncpu)
        !#########################################################
        ! Store new Hilbert tick marks after convergence
        !#########################################################
        bound_key_part(1:nhilbert,0:ncpu,ilev)=bound_key_target(1:nhilbert,0:ncpu)

     end do
     if(myid==1)write(*,*)"====================================="

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
  allocate(send_cnt(1:ncpu))
  allocate(recv_cnt(1:ncpu))
  allocate(recv_oft(1:ncpu))
  allocate(send_oft(1:ncpu))
  allocate(offset_cpu(1:ncpu))

  !#####################################
  ! Compute number of particles to send
  !#####################################
  send_cnt=0
  count_loc=0

  ! Loop over levels
  do ilev=ilevel,nlevelmax
     ix_ref=-1
     dx_loc=boxlen/2**ilev

     ! Loop over particles
     do ipart=headp(ilev),tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int(xp(ipart,1:ndim)/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1,1:ndim)))then
           ix_ref(1,1:ndim)=ix(1:ndim)
           grid_cpu = myid

           ! Compute Hilbert key of particle parent grid
           call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)

           ! Check if grid sits outside future processor boundaries
           if(    gt_keys(bound_key_part(1:nhilbert,myid-1,ilev),hk_ref(1,1:nhilbert)).OR. &
                & ge_keys(hk_ref(1,1:nhilbert),bound_key_part(1:nhilbert,myid,ilev)))then              
              ! Determine the future processor
              do icpu=1,ncpu
                 if(    ge_keys(hk_ref(1,1:nhilbert),bound_key_part(1:nhilbert,icpu-1,ilev)).AND. &
                      & gt_keys(bound_key_part(1:nhilbert,icpu,ilev),hk_ref(1,1:nhilbert)))then
                    grid_cpu=icpu
                 end if
              end do
           endif
        endif

        if(grid_cpu.EQ.myid)then
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

!  if(myid==1)write(*,*)'Counting done'
!  write(*,'(34(I6,1x))')myid,count_loc,send_cnt

  !#####################################
  ! Compute number of particles to receive
  !#####################################
  call MPI_ALLTOALL(send_cnt(1),1,MPI_INTEGER,recv_cnt(1),1,MPI_INTEGER,MPI_COMM_WORLD,info)
  send_cnt_tot=SUM(send_cnt)
  recv_cnt_tot=SUM(recv_cnt)
!  write(*,'("R ",34(I6,1x))')myid,recv_cnt

  !#####################################
  ! Compute offsets
  !#####################################
  recv_oft=0; send_oft=0
  offset_cpu(1)=headp(ilevel)-1+count_loc
  do icpu=2,ncpu
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
  do ilev=ilevel,nlevelmax
     ix_ref=-1
     dx_loc=boxlen/2**ilev

     ! Loop over particles
     do ipart=headp(ilev),tailp(ilev)

        ! Determine in which cpu particle should sit.
        ix = int(xp(ipart,1:ndim)/(2*dx_loc))
        if(.NOT. ALL(ix.EQ.ix_ref(1,1:ndim)))then
           ix_ref(1,1:ndim)=ix(1:ndim)
           grid_cpu=myid

           ! Compute Hilbert key of particle parent grid
           call hilbert_key(ix_ref,hk_ref,dummy_state,0,ilev-1,1)

           ! Check if grid sits outside future processor boundaries
           if(    gt_keys(bound_key_part(1:nhilbert,myid-1,ilev),hk_ref(1,1:nhilbert)).OR. &
                & ge_keys(hk_ref(1,1:nhilbert),bound_key_part(1:nhilbert,myid,ilev)))then
              ! Determine the future processor
              do icpu=1,ncpu
                 if(    ge_keys(hk_ref(1,1:nhilbert),bound_key_part(1:nhilbert,icpu-1,ilev)).AND. &
                      & gt_keys(bound_key_part(1:nhilbert,icpu,ilev),hk_ref(1,1:nhilbert)))then
                    grid_cpu=icpu
                 end if
              end do
           endif
        endif

        if(grid_cpu.EQ.myid)then
           ! Update count for local cpu
           count_loc=count_loc+1
           workp(ipart)=headp(ilevel)-1+count_loc
        else
           ! Update count for remote cpu
           send_cnt(grid_cpu)=send_cnt(grid_cpu)+1
           workp(ipart)=offset_cpu(grid_cpu)+send_cnt(grid_cpu)
        end if

     end do
     ! End loop over particles
  end do
  ! End loop over levels

  !#####################################
  ! Swap particles using new index table
  !#####################################
  do ipart=headp(ilevel),tailp(nlevelmax)
     do while(workp(ipart).NE.ipart)
        ! Swap new index
        jpart=workp(ipart)
        workp(ipart)=workp(jpart)
        workp(jpart)=jpart
        ! Swap positions
        xp_tmp(1:ndim)=xp(ipart,1:ndim)
        xp(ipart,1:ndim)=xp(jpart,1:ndim)
        xp(jpart,1:ndim)=xp_tmp(1:ndim)
        ! Swap velocities
        vp_tmp(1:ndim)=vp(ipart,1:ndim)
        vp(ipart,1:ndim)=vp(jpart,1:ndim)
        vp(jpart,1:ndim)=vp_tmp(1:ndim)
        ! Swap masses
        mp_tmp=mp(ipart)
        mp(ipart)=mp(jpart)
        mp(jpart)=mp_tmp
        ! Swap ids
        idp_tmp=idp(ipart)
        idp(ipart)=idp(jpart)
        idp(jpart)=idp_tmp
        ! Swap levels
        levelp_tmp=levelp(ipart)
        levelp(ipart)=levelp(jpart)
        levelp(jpart)=levelp_tmp
     end do
  end do

!  if(myid==1)write(*,*)'Swap done'
!  write(*,*)'S ',myid,count_loc,send_cnt_tot,recv_cnt_tot

  !###################################################################
  ! Set new number of particles in local processor
  !###################################################################
  npart=headp(ilevel)-1+count_loc+recv_cnt_tot
  tailp(nlevelmax)=npart

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
  do icpu=1,ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=headp(ilevel)-1+count_loc+i
     x_send_buf(i)=xp(ipart,idim)
  end do

  countsend=0
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     xp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !#########################
  ! Swap velocities
  !#########################
  do idim=1,ndim

  countrecv=0
  do icpu=1,ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=headp(ilevel)-1+count_loc+i
     x_send_buf(i)=vp(ipart,idim)
  end do

  countsend=0
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     vp(ipart,idim)=x_recv_buf(i)
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  end do

  !#########################
  ! Swap masses
  !#########################
  countrecv=0
  do icpu=1,ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(x_recv_buf(istart),nbuffer,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=headp(ilevel)-1+count_loc+i
     x_send_buf(i)=mp(ipart)
  end do

  countsend=0
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     mp(ipart)=x_recv_buf(i)
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
  do icpu=1,ncpu
     nbuffer=recv_cnt(icpu)
     if(nbuffer>0)then
        countrecv=countrecv+1
        istart=recv_oft(icpu)+1
        call MPI_IRECV(i_recv_buf(istart),nbuffer,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     endif
  end do

  do i=1,send_cnt_tot
     ipart=headp(ilevel)-1+count_loc+i
     i_send_buf(i)=levelp(ipart)
  end do

  countsend=0
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     levelp(ipart)=i_recv_buf(i)
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
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     l_send_buf(i)=idp(ipart)
  end do

  countsend=0
  do icpu=1,ncpu
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
     ipart=headp(ilevel)-1+count_loc+i
     idp(ipart)=l_recv_buf(i)
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
  deallocate(bound_key_part)

!  if(myid==1)write(*,*)'Swap done'
!  write(*,*)'SWAP ',myid,headp(ilevel),tailp(nlevelmax)

  !##################################
  ! Put all particles in level ilevel
  !##################################
  tailp(ilevel)=npart
  do ilev=ilevel+1,nlevelmax
     headp(ilev)=npart+1
     tailp(ilev)=npart
  end do

  call MPI_ALLREDUCE(npart,npart_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)

#endif

111 format('   Entering balance_part for level',i2)

end subroutine balance_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance(ilevel)
  use amr_commons
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
  integer(kind=8),allocatable,dimension(:)::bound_key_target
  real(dp)::xtarget

  integer::icell,j,ibit,ibucket,inew,iold,iold_true
  integer::noct_zero,head_zero,indx_zero
  integer::ncreate_tot,nkill_tot
  integer::parent_cell,skip_bit,true_level
  integer::get_parent_cell
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nlevelmax)::key_ref
  integer(kind=8), dimension(1:ndim)::cart_key
  integer(kind=8)::coarse_key
  integer,dimension(1:nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok_free,ok_all,ok
  type(oct)::oct_tmp


#ifndef WITHOUTMPI
  if(ncpu==1)return
  if(myid==1)write(*,111)ilevel
  
!!$  if(verbose)then
!!$     write(*,*)'Input mesh structure'
!!$     do ilev=levelmin,nlevelmax
!!$        if(noct_tot(ilev)>0)write(*,999)ilev,noct_tot(ilev),noct_min(ilev),noct_max(ilev),noct_tot(ilev)/ncpu
!!$     end do
!!$  end if

!!$  if(myid==1)then
!!$     write(*,*)'Old Hilbert tick marks'
!!$     do ilev=ilevel,nlevelmax
!!$        write(*,'(40(I6,1X))')(bound_key_level(icpu,ilev),icpu=0,ncpu)
!!$     end do
!!$  end if

  !-----------------------------------------------------
  ! Step 1: determine the new Hilbert tick marks
  !-----------------------------------------------------
  allocate(noct_cpu(1:ncpu))
  allocate(noct_cum(1:ncpu))
  allocate(ntarget_cum(1:ncpu))
  allocate(bound_key_target(0:ncpu))
  ! Compute new Hilbert tick marks
  do ilev=ilevel,nlevelmax

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
              bound_key_target(istart)=grid(ioct)%hkey+1
              istart=istart+1
           endif
           if(istart.GT.iright)exit
        end do
     endif

     call MPI_ALLREDUCE(bound_key_target,bound_key_target,ncpu+1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)

     bound_key_target(0)=0
     do icpu=1,ncpu
        bound_key_target(icpu)=max(bound_key_target(icpu),bound_key_target(icpu-1))
     end do
     bound_key_target(ncpu)=bound_key_level(ncpu,ilev)

     do icpu=0,ncpu
        bound_key_level(icpu,ilev)=bound_key_target(icpu)
     end do

  end do
  deallocate(noct_cpu,noct_cum,ntarget_cum)
  deallocate(bound_key_target)

!!$  if(myid==1)then
!!$     write(*,*)'New Hilbert tick marks'
!!$     do ilev=ilevel,nlevelmax
!!$        write(*,'(40(I6,1X))')(bound_key_level(icpu,ilev),icpu=0,ncpu)
!!$     end do
!!$  end if

  !-----------------------------------------------------
  ! Step 2: dispatch octs and empty slots according to
  ! the new target Hilbert tick marks
  !-----------------------------------------------------
  ifree=noct_used+1
  do ilev=ilevel,nlevelmax

     cache_operation=operation_loadbalance
     call open_cache

     hash_key(0)=ilev
     do ioct=head(ilev),tail(ilev)

        ! Check if grid sits outside future processor boundaries
        if(    grid(ioct)%hkey.lt.bound_key_level(myid-1,ilev).OR. &
             & grid(ioct)%hkey.ge.bound_key_level(myid  ,ilev))then
           
           ! Determine the future processor
           do icpu=1,ncpu
              if(    grid(ioct)%hkey.ge.bound_key_level(icpu-1,ilev).AND. &
                   & grid(ioct)%hkey.lt.bound_key_level(icpu  ,ilev))then
                 grid_cpu=icpu
              end if
           end do

           ! If next cache line is occupied, free it.
           if(occupied(free_cache))call destage(ngridmax+free_cache)
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

     call close_cache

  end do

!  if(myid==1)write(*,*)'Dispatch completed'

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
     if(true_level.NE.0)then
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
     if(true_level.NE.0)then
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
           skip_bit=ndim*(ilev-ibit) ! Carefull: this works only up to 63 bits !!
           ! Count octs in buckets
           bucket_count=0
           do inew=head_level(ilev),head_level(ilev)+noct_level(ilev)-1
              ioct=swap_table(inew)
              if(    grid(ioct)%hkey.lt.bound_key_level(myid-1,ilev).OR. &
                   & grid(ioct)%hkey.ge.bound_key_level(myid  ,ilev))then
                 write(*,*)'PE ',myid,'######### ',ioct
                 write(*,*)grid(ioct)%hkey
                 write(*,*)grid(ioct)%lev
                 write(*,*)grid(ioct)%ckey
                 write(*,*)bound_key_level(myid-1,ilev)
                 write(*,*)bound_key_level(myid  ,ilev)
                 stop
              endif
              ibucket=ibits(grid(ioct)%hkey,skip_bit,ndim)
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
              ibucket=ibits(grid(ioct)%hkey,skip_bit,ndim)
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
        if(grid(j)%lev>0)call hash_free(grid_dict,hash_key)
        oct_tmp=grid(j)
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
        grid(i)=oct_tmp
        hash_key(0)=grid(i)%lev
        hash_key(1:ndim)=grid(i)%ckey(1:ndim)
        if(grid(i)%lev>0)then
           call hash_set(grid_dict,hash_key,i)
        end if
        swap_table(i)=i
     endif
  end do
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
     key_ref=-1
     do ioct=head(ilev),tail(ilev)
        grid(ioct)%superoct=1
        coarse_key=grid(ioct)%hkey
        do i=1,MIN(ilev-1,nsuperoct)
           coarse_key=coarse_key/twotondim
           if(coarse_key.EQ.key_ref(i))then
              n_same(i)=n_same(i)+1
           else
              n_same(i)=1
              key_ref(i)=coarse_key
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
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(noct_used,noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

!!$  if(verbose)then
!!$     write(*,*)'Output mesh structure'
!!$     do ilev=levelmin,nlevelmax
!!$        if(noct_tot(ilev)>0)write(*,999)ilev,noct_tot(ilev),noct_min(ilev),noct_max(ilev),noct_tot(ilev)/ncpu
!!$     end do
!!$  end if
#endif

111 format(' Load balancing for all levels greater than ',I2)
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance(ilevel)
  use amr_commons
  use hydro_commons
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
  integer::icpu,grid_cpu,ichild,idom,jdom,mydom,ndom,lastdom,domains_matched
  integer::nleft,nright,ileft,iright,istart,nstart,noverlaps
  integer::ilev,ioct
  integer,dimension(:),allocatable::noct_dom,noct_cum
  integer,dimension(:),allocatable::ntarget_cum
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target_tot
  real(dp)::xtarget

  integer::icell,j,ibit,ibucket,inew,iold,iold_true
  integer::noct_zero,head_zero,indx_zero
  integer::ncreate_tot,nkill_tot
  integer::skip_bit,ikey,true_level
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key,zero_key,bleft,bright
  integer(kind=8),dimension(1:nhilbert),save::hks
  integer(kind=8),dimension(1:nhilbert,1:nlevelmax)::key_ref
  integer,dimension(1:nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok_free,ok_all,ok
  type(oct)::grid_tmp
  logical,allocatable,dimension(:,:) :: overlap
  integer,allocatable,dimension(:)   :: overlaps

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
  ndom = maxval(domain%n)
  allocate(noct_dom(1:ndom))
  allocate(noct_cum(1:ndom))
  allocate(ntarget_cum(1:ndom))
  allocate(bound_key_target(1:nhilbert,0:ndom))
  allocate(bound_key_target_tot(1:nhilbert,0:ndom))

  ! Compute new Hilbert tick marks
  do ilev=ilevel+1,nlevelmax

     if(noct_tot(ilev)>0)then
        mydom=domain(ilev)%r2d(myid)
        ndom=domain(ilev)%n
        if (ndom .ne. ncpu) then
           if (myid==1) write(*,*) &
              'Load_balance: Not yet support for multiple domains per MPI rank yet. Please rewrite the following section of code'
           if (myid==1) write(*,*) &
              'Load_balance: You have to scan over the octs in the rank and find out how many belong to each domain'
           stop
        end if
        noct_dom=0
        noct_dom(mydom)=noct(ilev)
        call MPI_ALLREDUCE(noct_dom,noct_cum,ndom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        noct_dom=noct_cum
        do idom=2,ndom
           noct_cum(idom)=noct_cum(idom-1)+noct_dom(idom)
        end do
        
        xtarget=dble(noct_cum(ndom))/dble(ndom)
        
        ileft=0
        iright=-1
        bound_key_target=0
        do idom=1,ndom
           ntarget_cum(idom)=int(dble(idom)*xtarget)
           if(mydom>1)then
              nleft=noct_cum(mydom-1)
           else
              nleft=0
           endif
           nright=noct_cum(mydom)
           IF(nright.GT.nleft)then
              if(ntarget_cum(idom).GT.nleft.AND.ntarget_cum(idom).LE.nright)then
                 if(ileft==0)ileft=idom
                 iright=MAX(idom,iright)
              endif
           endif
        end do
        
        if(iright.GE.ileft)then
           if(mydom.GT.1)then
              nstart=noct_cum(mydom-1)
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

        call MPI_ALLREDUCE(bound_key_target,bound_key_target_tot,nhilbert*(ndom+1),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
        bound_key_target=bound_key_target_tot

        bound_key_target(1:nhilbert,0)=zero_key
        do idom=1,ndom
           if(gt_keys(bound_key_target(1:nhilbert,idom-1),bound_key_target(1:nhilbert,idom)))then
              bound_key_target(1:nhilbert,idom)=bound_key_target(1:nhilbert,idom-1)
           endif
        end do
        bound_key_target(1:nhilbert,ndom)=hkey_max(1:nhilbert,ilev)
        domain(ilev)%b(1:nhilbert,:)=bound_key_target(1:nhilbert,0:ndom)

        !-----------------------------------------------------------------------
        ! Redefine the rank-to-domain mapping to optimize the connectivity to ilev-1
        !-----------------------------------------------------------------------
        allocate(overlap(ndom,ndom),overlaps(ndom))
        domain(ilev)%d2r=0
        domain(ilev)%r2d=0
        !
        do idom=1,ndom
           bleft  = coarsen_key(domain(ilev)%b(1:nhilbert,idom-1),ilev-1)
           bright = coarsen_key(domain(ilev)%b(1:nhilbert,idom),ilev-1)
           ! Test if domain from ilev-1 overlap with this domain
           do jdom=1,ndom
              if ( (ge_keys(bleft,domain(ilev-1)%b(1:nhilbert,jdom-1)) .and. &
                    ge_keys(domain(ilev-1)%b(1:nhilbert,jdom),bleft))  .or. &
                   (ge_keys(bright,domain(ilev-1)%b(1:nhilbert,jdom-1)) .and. &
                    ge_keys(domain(ilev-1)%b(1:nhilbert,jdom),bright))  ) then
                 overlap(idom,jdom) = .true.
              else
                 overlap(idom,jdom) = .false.
              end if
           enddo
        enddo
        !
        ! Select domain to rank matching based on intervals at ilev-1 that only overlap with one interval at ilev
        noverlaps = 1
        domains_matched = 0
        do while(noverlaps .ne. 0 .and. domains_matched < ndom)
          overlaps = count(overlap,dim=1) ! Count how many intervals overlap with each interval at ilev-1
          noverlaps = 0
          do jdom=1,ndom
             if (overlaps(jdom)==1) then
                do idom=1,ndom
                   if (overlap(idom,jdom)) then
                      noverlaps = noverlaps + 1
                      domain(ilev)%d2r(idom) = domain(ilev-1)%d2r(jdom)
                      domain(ilev)%r2d(domain(ilev)%d2r(idom)) = idom
                      domains_matched = domains_matched + 1
                      overlap(idom,:) = .false.
                      overlap(:,jdom) = .false.
                      exit
                   endif
                enddo
             endif
          enddo
        enddo
        ! Match the rest accepting multiple overlaps, taking the first available match
        noverlaps = 1
        do while(noverlaps .gt. 0 .and. domains_matched < ndom)
          overlaps = count(overlap,dim=1) ! Count how many intervals overlap with each interval at ilev-1
          noverlaps = 0
          do jdom=1,ndom
             if (overlaps(jdom) > 0) then
                do idom=1,ndom
                   if (overlap(idom,jdom)) then
                      noverlaps = noverlaps + 1
                      domain(ilev)%d2r(idom) = domain(ilev-1)%d2r(jdom)
                      domain(ilev)%r2d(domain(ilev)%d2r(idom)) = idom
                      domains_matched = domains_matched + 1
                      overlap(idom,:) = .false.
                      overlap(:,jdom) = .false.
                      exit
                   endif
                enddo
             endif
          enddo
        enddo
        ! Hopefully very few are left, and are put linearly according to leftover slots in the domain2rank array
        if (domains_matched < ndom) then
           lastdom = 1
           do idom=1,ndom
             if (domain(ilev)%r2d(idom)==0) then
                do jdom=lastdom,ndom
                   if (domain(ilev)%d2r(jdom)==0) then
                      domain(ilev)%d2r(jdom) = idom
                      domain(ilev)%r2d(idom) = jdom
                      lastdom = jdom+1
                      exit
                   endif
                enddo
             endif
           enddo
        endif

        deallocate(overlaps,overlap)

     endif

  end do
  deallocate(noct_dom,noct_cum,ntarget_cum)
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
        if (.not. domain(ilev)%in_rank(grid(ioct)%hkey)) then
           
           ! Determine the future processor
           hks = grid(ioct)%hkey(1:nhilbert)
           grid_cpu = domain(ilev)%get_rank(hks)

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
subroutine load_balance_2(r,g,m,ilevel)
  use amr_parameters, only: ndim,twotondim,nhilbert,dp
  use amr_commons, only: run_t,global_t,mesh_t,oct
  use hilbert
  use hash
  use cache_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::igrid,i,ind,jlevel,info
  integer::icpu,grid_cpu,ichild,idom,jdom,mydom,ndom,lastdom,domains_matched
  integer::nleft,nright,ileft,iright,istart,nstart,noverlaps
  integer::ilev,ioct
  integer,dimension(:),allocatable::noct_dom,noct_cum
  integer,dimension(:),allocatable::ntarget_cum
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target_tot
  real(dp)::xtarget

  integer::icell,j,ibit,ibucket,inew,iold,iold_true
  integer::noct_zero,head_zero,indx_zero
  integer::ncreate_tot,nkill_tot
  integer::skip_bit,ikey,true_level
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:ndim)::cart_key
  integer(kind=8),dimension(1:nhilbert)::coarse_key,one_key,zero_key,bleft,bright
  integer(kind=8),dimension(1:nhilbert),save::hks
  integer(kind=8),dimension(1:nhilbert,1:r%nlevelmax)::key_ref
  integer,dimension(1:r%nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok_free,ok_all,ok
  type(oct)::grid_tmp
  logical,allocatable,dimension(:,:) :: overlap
  integer,allocatable,dimension(:)   :: overlaps

#ifndef WITHOUTMPI
  if(g%ncpu==1)return
  if(ilevel==r%nlevelmax)return
  if(r%verbose)write(*,111)ilevel
  
  ! Constants
  zero_key=0
  one_key=0
  one_key(1)=1

  !-----------------------------------------------------
  ! Step 1: determine the new Hilbert tick marks
  !-----------------------------------------------------
  ndom = maxval(m%domain%n)
  allocate(noct_dom(1:ndom))
  allocate(noct_cum(1:ndom))
  allocate(ntarget_cum(1:ndom))
  allocate(bound_key_target(1:nhilbert,0:ndom))
  allocate(bound_key_target_tot(1:nhilbert,0:ndom))

  ! Compute new Hilbert tick marks
  do ilev=ilevel+1,r%nlevelmax

     if(m%noct_tot(ilev)>0)then
        mydom=m%domain(ilev)%r2d(g%myid)
        ndom=m%domain(ilev)%n
        if (ndom .ne. g%ncpu) then
           if (g%myid==1) write(*,*) &
              'Load_balance: Not yet support for multiple domains per MPI rank yet. Please rewrite the following section of code'
           if (g%myid==1) write(*,*) &
              'Load_balance: You have to scan over the octs in the rank and find out how many belong to each domain'
           stop
        end if
        noct_dom=0
        noct_dom(mydom)=m%noct(ilev)
        call MPI_ALLREDUCE(noct_dom,noct_cum,ndom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        noct_dom=noct_cum
        do idom=2,ndom
           noct_cum(idom)=noct_cum(idom-1)+noct_dom(idom)
        end do
        
        xtarget=dble(noct_cum(ndom))/dble(ndom)
        
        ileft=0
        iright=-1
        bound_key_target=0
        do idom=1,ndom
           ntarget_cum(idom)=int(dble(idom)*xtarget)
           if(mydom>1)then
              nleft=noct_cum(mydom-1)
           else
              nleft=0
           endif
           nright=noct_cum(mydom)
           IF(nright.GT.nleft)then
              if(ntarget_cum(idom).GT.nleft.AND.ntarget_cum(idom).LE.nright)then
                 if(ileft==0)ileft=idom
                 iright=MAX(idom,iright)
              endif
           endif
        end do
        
        if(iright.GE.ileft)then
           if(mydom.GT.1)then
              nstart=noct_cum(mydom-1)
           else
              nstart=0
           endif
           istart=ileft
           do ioct=m%head(ilev),m%tail(ilev)
              nstart=nstart+1
              if(nstart.GE.ntarget_cum(istart))then
                 bound_key_target(1:nhilbert,istart)=m%grid(ioct)%hkey(1:nhilbert)+one_key
                 istart=istart+1
              endif
              if(istart.GT.iright)exit
           end do
        endif

        call MPI_ALLREDUCE(bound_key_target,bound_key_target_tot,nhilbert*(ndom+1),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
        bound_key_target=bound_key_target_tot

        bound_key_target(1:nhilbert,0)=zero_key
        do idom=1,ndom
           if(gt_keys(bound_key_target(1:nhilbert,idom-1),bound_key_target(1:nhilbert,idom)))then
              bound_key_target(1:nhilbert,idom)=bound_key_target(1:nhilbert,idom-1)
           endif
        end do
        bound_key_target(1:nhilbert,ndom)=m%hkey_max(1:nhilbert,ilev)
        m%domain(ilev)%b(1:nhilbert,:)=bound_key_target(1:nhilbert,0:ndom)

        !-----------------------------------------------------------------------
        ! Redefine the rank-to-domain mapping to optimize the connectivity to ilev-1
        !-----------------------------------------------------------------------
        allocate(overlap(ndom,ndom),overlaps(ndom))
        m%domain(ilev)%d2r=0
        m%domain(ilev)%r2d=0
        !
        do idom=1,ndom
           bleft  = coarsen_key(m%domain(ilev)%b(1:nhilbert,idom-1),ilev-1)
           bright = coarsen_key(m%domain(ilev)%b(1:nhilbert,idom),ilev-1)
           ! Test if domain from ilev-1 overlap with this domain
           do jdom=1,ndom
              if ( (ge_keys(bleft,m%domain(ilev-1)%b(1:nhilbert,jdom-1)) .and. &
                    ge_keys(m%domain(ilev-1)%b(1:nhilbert,jdom),bleft))  .or. &
                   (ge_keys(bright,m%domain(ilev-1)%b(1:nhilbert,jdom-1)) .and. &
                    ge_keys(m%domain(ilev-1)%b(1:nhilbert,jdom),bright))  ) then
                 overlap(idom,jdom) = .true.
              else
                 overlap(idom,jdom) = .false.
              end if
           enddo
        enddo
        !
        ! Select domain to rank matching based on intervals at ilev-1 that only overlap with one interval at ilev
        noverlaps = 1
        domains_matched = 0
        do while(noverlaps .ne. 0 .and. domains_matched < ndom)
          overlaps = count(overlap,dim=1) ! Count how many intervals overlap with each interval at ilev-1
          noverlaps = 0
          do jdom=1,ndom
             if (overlaps(jdom)==1) then
                do idom=1,ndom
                   if (overlap(idom,jdom)) then
                      noverlaps = noverlaps + 1
                      m%domain(ilev)%d2r(idom) = m%domain(ilev-1)%d2r(jdom)
                      m%domain(ilev)%r2d(m%domain(ilev)%d2r(idom)) = idom
                      domains_matched = domains_matched + 1
                      overlap(idom,:) = .false.
                      overlap(:,jdom) = .false.
                      exit
                   endif
                enddo
             endif
          enddo
        enddo
        ! Match the rest accepting multiple overlaps, taking the first available match
        noverlaps = 1
        do while(noverlaps .gt. 0 .and. domains_matched < ndom)
          overlaps = count(overlap,dim=1) ! Count how many intervals overlap with each interval at ilev-1
          noverlaps = 0
          do jdom=1,ndom
             if (overlaps(jdom) > 0) then
                do idom=1,ndom
                   if (overlap(idom,jdom)) then
                      noverlaps = noverlaps + 1
                      m%domain(ilev)%d2r(idom) = m%domain(ilev-1)%d2r(jdom)
                      m%domain(ilev)%r2d(m%domain(ilev)%d2r(idom)) = idom
                      domains_matched = domains_matched + 1
                      overlap(idom,:) = .false.
                      overlap(:,jdom) = .false.
                      exit
                   endif
                enddo
             endif
          enddo
        enddo
        ! Hopefully very few are left, and are put linearly according to leftover slots in the domain2rank array
        if (domains_matched < ndom) then
           lastdom = 1
           do idom=1,ndom
             if (m%domain(ilev)%r2d(idom)==0) then
                do jdom=lastdom,ndom
                   if (m%domain(ilev)%d2r(jdom)==0) then
                      m%domain(ilev)%d2r(jdom) = idom
                      m%domain(ilev)%r2d(idom) = jdom
                      lastdom = jdom+1
                      exit
                   endif
                enddo
             endif
           enddo
        endif

        deallocate(overlaps,overlap)

     endif

  end do
  deallocate(noct_dom,noct_cum,ntarget_cum)
  deallocate(bound_key_target)
  deallocate(bound_key_target_tot)

  !-----------------------------------------------------
  ! Step 2: dispatch octs and empty slots according to
  ! the new target Hilbert tick marks
  !-----------------------------------------------------
  m%ifree=m%noct_used+1
  do ilev=ilevel+1,r%nlevelmax

     call open_cache_2(r,g,m,operation_loadbalance,domain_decompos_amr)

     hash_key(0)=ilev
     do ioct=m%head(ilev),m%tail(ilev)

        ! Check if grid sits outside future processor boundaries
        if (.not. m%domain(ilev)%in_rank(m%grid(ioct)%hkey)) then
           
           ! Determine the future processor
           hks = m%grid(ioct)%hkey(1:nhilbert)
           grid_cpu = m%domain(ilev)%get_rank(hks)

           ! If next cache line is occupied, free it.
           if(m%occupied(m%free_cache))call destage_2(r,g,m,r%ngridmax+m%free_cache,m%grid_dict)
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

     call close_cache_2(r,g,m,m%grid_dict)

  end do

  !-----------------------------------------------------
  ! Step 3: sort new octs and empty slots according to
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
  ! Step 4: sort octs level by level according to their
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
  ! Step 5: Apply permutations directly in main memory
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
  ! Step 6: Clean up final AMR structure
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

  !---------------------
  ! Total number of octs
  !---------------------
  do ilev=ilevel+1,r%nlevelmax
     m%noct_tot(ilev)=m%noct(ilev)
     m%noct_min(ilev)=m%noct(ilev)
     m%noct_max(ilev)=m%noct(ilev)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_tot(ilev),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_min(ilev),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(m%noct(ilev),m%noct_max(ilev),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  end do

  m%noct_used_max=m%noct_used
  m%noct_used_tot=m%noct_used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(m%noct_used,m%noct_used_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m%noct_used,m%noct_used_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif

#endif

111 format(' Load balancing for all levels greater than ',I2)
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine load_balance_2
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
  integer::ilev,idom,icpu,jcpu,mydom,ndom,count_loc,recv_cnt_tot,send_cnt_tot
  integer::nbuffer,countrecv,countsend,tag=101
  real(kind=8)::dx_loc

  integer,allocatable,dimension(:)::send_cnt,recv_cnt
  integer,allocatable,dimension(:)::send_oft,recv_oft,offset_cpu
  integer,dimension(ncpu)::reqsend,reqrecv

  real(kind=8),dimension(:),allocatable::x_recv_buf,x_send_buf
  integer(i8b),dimension(:),allocatable::l_recv_buf,l_send_buf
  integer,dimension(:),allocatable::i_recv_buf,i_send_buf

  integer(kind=8)::unbalance
  integer(kind=8),dimension(1:nhilbert)::diff_key
  integer(kind=8),dimension(1:nhilbert),save::hks
  type(domain_t),allocatable,dimension(:)::domain_part
  integer(kind=8),allocatable,dimension(:,:)::bound_key_target,bound_key_new
  integer(kind=8),allocatable,dimension(:,:)::bound_key_left,bound_key_right
  integer,dimension(1:ncpu)::npart_proc,npart_proc_tot
  integer::npart_lev,npart_lev_tot,iter
  integer,dimension(0:ncpu)::npart_cum,npart_cum_tot
  integer,dimension(1:ncpu)::npart_dom,npart_dom_tot
  real(dp)::xpart_target,xcum_target

  real(dp),dimension(1:ndim),save::xp_tmp,vp_tmp
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

#ifndef WITHOUTMPI
  if(ncpu==1)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !####################################################
  ! Default for particle domains are grid domains
  !####################################################
  allocate(domain_part(ilevel:nlevelmax+1))
  do ilev=ilevel,nlevelmax+1
     call domain_part(ilev)%copy(domain(ilev))
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
     do ilev=ilevel,nlevelmax
        mydom = domain_part(ilev)%r2d(myid)
        ndom  = domain_part(ilev)%n

        dx_loc=boxlen/2**ilev
     
        ! Compute number of particles
        npart_lev=tailp(ilev)-headp(ilev)+1
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
        if(myid==1)write(*,*)"====================================="
        if(myid==1)write(*,'("Level=",I4," npart=",I10)')ilev,npart_lev_tot
!        if(myid==1)write(*,'(16(I10,1X))')npart_dom
!        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(idom)*xpart_target),idom=0,ndom)
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
        bound_key_target(1:nhilbert,0:ndom)=domain_part(ilev)%b(1:nhilbert,0:ndom)
        bound_key_new=bound_key_target
        bound_key_left=0
        do idom=0,ndom
           bound_key_right(1:nhilbert,idom)=hkey_max(1:nhilbert,ilev)
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
        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,npart_cum
        if(myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(idom)*xpart_target),idom=0,ndom)
        !#########################################################
        ! Store new Hilbert tick marks after convergence
        !#########################################################
        domain_part(ilev)%b(1:nhilbert,0:ndom)=bound_key_target(1:nhilbert,0:ndom)

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
           if (.not. domain_part(ilev)%in_rank(hk_ref(1,1:nhilbert))) then
              ! Determine the future processor
              hks = hk_ref(1,1:nhilbert)
              grid_cpu = domain_part(ilev)%get_rank(hks)
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
           if (.not. domain_part(ilev)%in_rank(hk_ref(1,1:nhilbert))) then
              ! Determine the future processor
              hks = hk_ref(1,1:nhilbert)
              grid_cpu = domain_part(ilev)%get_rank(hks)
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
  do ilev=ilevel,nlevelmax+1
    call domain_part(ilev)%destroy
  end do
  deallocate(domain_part)

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
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine balance_part_2(r,g,m,p,ilevel)
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
  !
  ! This routine will dispatch particles across processors according to
  ! their Hilbert key and for a given domain decomposition.
  ! It assumes that particles are sorted according to their levels
  ! and within the level, according to their Hilbert key.
  ! It can only be called after routine rho has been called.
  !
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,g%ncpu)::statuses
#endif
  integer(kind=8), dimension(1:nvector,1:nhilbert),save::hk_ref
  integer(kind=8), dimension(1:nvector,1:ndim),save::ix_ref
  integer(kind=4), dimension(1:nvector),save::dummy_state

  integer,dimension(1:ndim),save::ix
  integer::i,istart,info,ipart,jpart,idim,grid_cpu
  integer::ilev,idom,icpu,jcpu,mydom,ndom,count_loc,recv_cnt_tot,send_cnt_tot
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
  integer,dimension(1:g%ncpu)::npart_proc,npart_proc_tot
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
        if(g%myid==1)write(*,*)"====================================="
        if(g%myid==1)write(*,'("Level=",I4," npart=",I10)')ilev,npart_lev_tot

        !#########################################################
        ! Sort particle according to current level Hilbert key
        !#########################################################
        do i=p%headp(ilev),p%tailp(ilev)
           p%sortp(i)=i
        end do
        ix=0
        call sort_hilbert_2(r,g,p,p%headp(ilev),p%tailp(ilev),ix,0,1,ilev-1)

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
        if(g%myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,npart_cum
        if(g%myid==1)write(*,'("iter=",I4,1X,17(I10,1X))')iter,(int(dble(idom)*xpart_target),idom=0,ndom)
        !#########################################################
        ! Store new Hilbert tick marks after convergence
        !#########################################################
        domain_part(ilev)%b(1:nhilbert,0:ndom)=bound_key_target(1:nhilbert,0:ndom)

     end do
     if(g%myid==1)write(*,*)"====================================="

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

end subroutine balance_part_2
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

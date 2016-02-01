subroutine refine_all
  use amr_commons
  implicit none

  integer::ilevel

  if(verbose)write(*,*)'Entering refine_all' 

  do ilevel=levelmin,nlevelmax-1
     call refine_fine(ilevel)
  end do

  if(verbose)write(*,*)'Complete refine_all'

end subroutine refine_all
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine refine_fine(ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !---------------------------------------------------------
  ! This routine refines cells at level ilevel if cells
  ! are flagged for refinement and are not already refined.
  ! This routine destroys refinements for cells that are 
  ! not flagged for refinement and are already refined.
  ! For single time-stepping, numerical rules are 
  ! automatically satisfied. For adaptive time-stepping,
  ! numerical rules are checked before refining any cell.
  !---------------------------------------------------------
  integer::igrid,icell,i,j,ibit,ibucket,ilev,ind,inew,ioct,iold,iold_true
  integer::noct_zero,head_zero,indx_zero,info
  integer::ncreate_tot,nkill_tot
  integer::parent_cell,skip_bit,true_level
  integer::get_parent_cell
  integer::ind_cell,ind_parent
  integer(kind=8),dimension(0:ndim)::hash_key
  integer(kind=8),dimension(1:nlevelmax)::key_ref
  integer(kind=8)::coarse_key
  integer,dimension(1:nlevelmax)::n_same,npatch
  integer,dimension(:),allocatable::noct_level,head_level,indx_level
  integer,dimension(:),allocatable::swap_table,swap_tmp
  integer,dimension(0:twotondim-1)::bucket_count,bucket_offset
  logical::ok_free,ok_all,ok
  type(oct)::oct_tmp

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !---------------------------------------------------
  ! Step 1: if a cell is flagged for refinement and
  ! if it is not already refined, create its son grid.
  !---------------------------------------------------        
  ncreate=0
  ifree=noct_used+1
  do ilev=ilevel,nlevelmax-1

     cache_operation=operation_refine
     call open_cache

     do ioct=head(ilev),tail(ilev)
        do ind=1,twotondim
           ok   = grid(ioct)%flag1(ind)==1 .and. &
                & .not.grid(ioct)%refined(ind)
           if(ok)then
              ind_parent=ioct
              ind_cell=ind
              call make_new_oct(ind_parent,ind_cell,ilev+1)
              ncreate=ncreate+1
           endif
        end do
     end do

     call close_cache

  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncreate,ncreate_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,112)ncreate_tot
#else
  if(verbose)write(*,112)ncreate
#endif

  !----------------------------------------------------------
  ! Step 2: if the parent cell is not flagged for refinement,
  ! but it is refined, then destroy the child grid.
  !----------------------------------------------------------
  nkill=0  
  do ilev=ilevel+1,nlevelmax

     cache_operation=operation_derefine
     call open_cache

     hash_key(0)=ilev
     do ioct=head(ilev),tail(ilev)
        hash_key(1:ndim)=grid(ioct)%ckey(1:ndim)
        parent_cell=get_parent_cell(hash_key,.true.,.true.)
        igrid=(parent_cell-1)/twotondim+1
        icell=parent_cell-(igrid-1)*twotondim
        ok   = grid(igrid)%flag1(icell)==0 .and. &
             & grid(igrid)%refined(icell)
        if(ok)then
           ! Set grid level to zero
           grid(ioct)%lev=0
           ! Set parent cell to "unrefined" status
           grid(igrid)%refined(icell)=.false.
           ! Free grid from hash table
           call hash_free(grid_dict,hash_key)
           nkill=nkill+1
        end if
     end do

     call close_cache

  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nkill,nkill_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,112)nkill_tot
#else
  if(verbose)write(*,112)nkill
#endif

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
222 continue

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
!!$  if(myid==1)then
!!$     write(*,*)'Temporary mesh structure'
!!$     do ilev=levelmin,nlevelmax
!!$        if(noct_tot(ilev)>0)write(*,999)ilev,noct_tot(ilev),noct_min(ilev),noct_max(ilev),noct_tot(ilev)/ncpu
!!$     end do
!!$  end if
!!$999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

111 format('   Entering refine_fine for level ',I2)
112 format('   ==> Make ',i6,' sub-grids')
113 format('   ==> Kill ',i6,' sub-grids')

end subroutine refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_new_oct(iparent,icell,ilevel)
  use amr_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::iparent,icell
  !--------------------------------------------------------------
  ! This routine creates a children oct at level ilevel.
  ! ilevel is thus the level of the new children oct.
  ! The new oct is labelled using the new index ichild.
  ! The parent cell is labeled with the parent oct index iparent
  ! and the cell index icell (from 1 to 8).
  !--------------------------------------------------------------
  integer::icpu,idim,ichild,ind,nstride,grid_cpu
  integer(kind=4), dimension(1:nvector),save::dummy_state
  integer(kind=8), dimension(1:nvector),save::hk0,hk1,hk2
  integer(kind=8), dimension(1:nvector),save::ix,iy,iz
  integer(kind=8), dimension(1:ndim),save::cart_key
  integer(kind=8), dimension(0:ndim)::hash_key

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mail_counter==32)then
     call check_mail(MPI_REQUEST_NULL)
     mail_counter=0
  endif
  mail_counter=mail_counter+1
#endif

  !=================================
  ! Create new octs into main memory
  !=================================
  ! Compute Cartesian keys of new octs
  do idim=1,ndim
     nstride=2**(idim-1)
     cart_key(idim)=2*grid(iparent)%ckey(idim)+MOD((icell-1)/nstride,2)
  end do

  ! Compute Hilbert keys of new octs
#if NDIM==1
  ix(1)=cart_key(1)
  call hilbert1d(ix,hk0,1)
#endif
#if NDIM==2
  ix(1)=cart_key(1)
  iy(1)=cart_key(2)
  call hilbert2d(ix,iy,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif
#if NDIM==3
  ix(1)=cart_key(1)
  iy(1)=cart_key(2)
  iz(1)=cart_key(3)
  call hilbert3d(ix,iy,iz,hk2,hk1,hk0,dummy_state,0,ilevel-1,1)
#endif

  ! Check if grid sits inside processor boundaries
  if(    hk0(1).ge.bound_key_level(myid-1,ilevel).AND. &
       & hk0(1).lt.bound_key_level(myid  ,ilevel))then


     ! Set grid index to a virtual grid in local main memory
     ichild=ifree

     ! Go to next main memory free line
     ifree=ifree+1
     if(ifree.GT.ngridmax)then
        write(*,*)'No more free memory'
        write(*,*)'Increase ngridmax'
        call clean_abort
     end if

  ! Otherwise, determine parent processor and use the cache
  else
     do icpu=1,ncpu
        if(    hk0(1).ge.bound_key_level(icpu-1,ilevel).AND. &
             & hk0(1).lt.bound_key_level(icpu  ,ilevel))then
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
  endif

  grid(ichild)%lev=ilevel
  grid(ichild)%ckey(1:ndim)=cart_key(1:ndim)
  grid(ichild)%hkey=hk0(1)
  grid(ichild)%refined(1:twotondim)=.false.
  grid(ichild)%flag1(1:twotondim)=0
  grid(ichild)%flag2(1:twotondim)=0
  grid(ichild)%superoct=1

  ! Insert new grid in hash table
  hash_key(0)=ilevel
  hash_key(1:ndim)=cart_key(1:ndim)
  call hash_set(grid_dict,hash_key,ichild)

  ! Set status of parent cell to "refined"
  grid(iparent)%refined(icell)=.true.

  ! Inject parent hydro variables into new children ones
  if(.not.init)then
     ! Interpolate hydro variables
#ifdef SOLVERhydro
     do ind=1,twotondim
        grid(ichild)%uold(ind,1:nvar)=grid(iparent)%uold(icell,1:nvar)
     enddo
#endif
     ! Interpolate gravity variables
#ifdef GRAV
     do ind=1,twotondim
        grid(ichild)%f(ind,1:ndim)=grid(iparent)%f(icell,1:ndim)
        grid(ichild)%phi(ind)=grid(iparent)%phi(icell)
        grid(ichild)%phi_old(ind)=grid(iparent)%phi_old(icell)
     enddo
#endif
  endif

end subroutine make_new_oct
!###############################################################
!###############################################################
!###############################################################
!###############################################################
  

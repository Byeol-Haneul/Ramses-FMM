!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag
  use amr_commons
  implicit none
  integer::ilevel

  if(verbose)write(*,*)'Entering flag'
  do ilevel=nlevelmax-1,levelmin,-1
     call flag_fine(ilevel,2)
  end do
  if(verbose)write(*,*)'Complete flag'

end subroutine flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine flag_fine(ilevel,icount)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !--------------------------------------------------------
  ! This routine builds the refinement map at level ilevel.
  !--------------------------------------------------------
  integer::iexpand,nflag_tot,info
  
  if(ilevel==nlevelmax)return
  if(ilevel<levelmin)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Step 1: initialize refinement map to minimal refinement rules
  call init_flag(ilevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,*) '  ==> end step 1',nflag_tot
#else
  if(verbose)write(*,*) '  ==> end step 1',nflag
#endif
  
  ! Step 2: make one cubic buffer around flagged cells,
  ! in order to enforce numerical rule.
  call build_smooth(ilevel)
  call smooth_fine_fast(ilevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,*) '  ==> end step 2',nflag_tot
#else
  if(verbose)write(*,*) '  ==> end step 2',nflag
#endif

  ! Step 3: if cell satisfies user-defined physical citeria,
  ! then flag cell for refinement.
  call userflag_fine(ilevel)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,*) '  ==> end step 3',nflag_tot
#else
  if(verbose)write(*,*) '  ==> end step 3',nflag
#endif

  ! Step 4: make nexpand cubic buffers around flagged cells.
  do iexpand=1,nexpand(ilevel)
     call smooth_fine_fast(ilevel)
  end do
  call clean_smooth
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,*) '  ==> end step 4',nflag_tot
#else
  if(verbose)write(*,*) '  ==> end step 4',nflag
#endif

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  if(verbose)write(*,112)nflag_tot
#else
  if(verbose)write(*,112)nflag
#endif

  ! In case of adaptive time step ONLY, check for refinement rules.
  if(ilevel>levelmin)then
     if(icount<nsubcycle(ilevel-1))then
        call ensure_ref_rules(ilevel)
     end if
  end if

111 format('   Entering flag_fine for level ',I2)
112 format('   ==> Flag ',i6,' cells')

end subroutine flag_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flag(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  !-------------------------------------------
  ! This routine initialize the refinement map
  ! to a minimal state in order to satisfy the
  ! refinement rules.
  !-------------------------------------------
  integer::igrid,ichild,icell,ind,parent_cell
  integer::get_parent_cell
  logical::ok
  integer(kind=8),dimension(0:ndim)::hash_key

  ! Initialize flag1 to 0 for ilevel grids
  nflag=0
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        grid(igrid)%flag1(ind)=0
     end do
  end do
  !---------------------------------------------------------
  ! Set flag1 to 1 if cell is refined and  contains a 
  ! flagged son or a refined son.
  ! This ensures that refinement rules are satisfied.
  !---------------------------------------------------------
  call open_cache(operation_initflag,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ichild=head(ilevel+1),tail(ilevel+1)
     hash_key(1:ndim)=grid(ichild)%ckey(1:ndim)
     parent_cell=get_parent_cell(hash_key,grid_dict,.true.,.false.)
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim
     ok=.false.
     ! Loop over cells
     do ind=1,twotondim
        ok=(ok.or.(grid(ichild)%refined(ind)))
        ok=(ok.or.(grid(ichild)%flag1(ind)==1))
     end do
     if(ok)then
        grid(igrid)%flag1(icell)=1
        nflag=nflag+1
     endif
  end do

  call close_cache(grid_dict)

end subroutine init_flag
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine userflag_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return
  
  ! Refinement rules for the gravity solver
  if(poisson)call poisson_flag(ilevel)

  ! Refinement rules for the hydro solver
  if(hydro)call hydro_flag(ilevel)

end subroutine userflag_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine smooth_fine(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth,count_nbor,ig,ih,in
  integer::i,iskip,ngrid
  integer::igrid,idim,ind,i_nbor,igrid_nbor,icell_nbor
  integer::parent_cell,get_parent_cell,get_grid
  integer,dimension(1:3),save::n_nbor=(/1,2,2/)
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer,dimension(0:twondim)::igridn

  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,&
       &   0,-1,0,0,1,0,&
       &   0,0,-1,0,0,1/),(/3,6/))
  integer,dimension(1:8,1:6),save::ggg=reshape(&
       & (/1,0,1,0,1,0,1,0,&
       &   0,2,0,2,0,2,0,2,&
       &   3,3,0,0,3,3,0,0,&
       &   0,0,4,4,0,0,4,4,&
       &   5,5,5,5,0,0,0,0,&
       &   0,0,0,0,6,6,6,6/),(/8,6/))
  integer,dimension(1:8,1:6),save::hhh=reshape(&
       & (/2,1,4,3,6,5,8,7,&
       &   2,1,4,3,6,5,8,7,&
       &   3,4,1,2,7,8,5,6,&
       &   3,4,1,2,7,8,5,6,&
       &   5,6,7,8,1,2,3,4,&
       &   5,6,7,8,1,2,3,4/),(/8,6/))

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return

  hash_nbor(0)=ilevel

  ! Loop over steps
  do ismooth=1,ndim

     ! Initialize flag2 to 0
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%flag2(ind)=0
        end do
     end do

     call open_cache(operation_smooth,domain_decompos_amr)

     ! Count neighbors and set flag2 accordingly
     do igrid=head(ilevel),tail(ilevel)

        ! Get neighboring octs
        igridn(0)=igrid
        do i_nbor=1,twondim
           hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,i_nbor)
           ! Periodic boundary conditons
           do idim=1,ndim
              if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
              if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
           enddo
           igridn(i_nbor)=get_grid(hash_nbor,grid_dict,.false.,.true.)
           call lock_cache(igridn(i_nbor))
        end do

        ! Count neighbors and set flag2 accordingly        
        do ind=1,twotondim
           count_nbor=0
           do in=1,twondim
              ig=ggg(ind,in)
              igrid_nbor=igridn(ig)
              icell_nbor=hhh(ind,in)
              if(igrid_nbor>0)then
                 count_nbor=count_nbor+grid(igrid_nbor)%flag1(icell_nbor)
              endif
           end do
           ! flag2 cell if necessary
           if(count_nbor>=n_nbor(ismooth))then
              grid(igrid)%flag2(ind)=1
           endif
        end do

        do i_nbor=1,twondim
           call unlock_cache(igridn(i_nbor))
        end do

     end do
     ! End loop over grids

    call close_cache(grid_dict)

     ! Set flag1=1 for cells with flag2=1
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           if(grid(igrid)%flag1(ind)==1)grid(igrid)%flag2(ind)=0
        end do
        do ind=1,twotondim
           if(grid(igrid)%flag2(ind)==1)then
              grid(igrid)%flag1(ind)=1
              nflag=nflag+1
           endif
        end do
     end do     

  end do
  ! End loop over steps
  
end subroutine smooth_fine
!############################################################
!############################################################
!############################################################
!############################################################
subroutine smooth_fine_fast(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  ! -------------------------------------------------------------------
  ! Dilatation operator.
  ! This routine makes one cell width cubic buffer around flag1 cells 
  ! at level ilevel by following these 3 steps:
  ! step 1: flag1 cells with at least 1 flag1 neighbors (if ndim > 0) 
  ! step 2: flag1 cells with at least 2 flag1 neighbors (if ndim > 1) 
  ! step 3: flag1 cells with at least 2 flag1 neighbors (if ndim > 2) 
  ! Array flag2 is used as temporary workspace.
  ! -------------------------------------------------------------------
  integer::ismooth,count_nbor,ig,id,ih,in
  integer::i,iskip,ngrid
  integer::igrid,idim,ind,inbor,igrid_nbor,icell_nbor
  integer::icpu,nbuffer,istart,info
  integer::parent_cell,get_parent_cell,get_grid
  integer,dimension(1:3),save::n_nbor=(/1,2,2/)
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::igridn
  integer,dimension(1:twotondim,0:twondim),save::flag_nbor
  integer::countrecv,countsend,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,&
       &   0,-1,0,0,1,0,&
       &   0,0,-1,0,0,1/),(/3,6/))
  integer,dimension(1:8,1:6),save::ggg=reshape(&
       & (/1,0,1,0,1,0,1,0,&
       &   0,2,0,2,0,2,0,2,&
       &   3,3,0,0,3,3,0,0,&
       &   0,0,4,4,0,0,4,4,&
       &   5,5,5,5,0,0,0,0,&
       &   0,0,0,0,6,6,6,6/),(/8,6/))
  integer,dimension(1:8,1:6),save::hhh=reshape(&
       & (/2,1,4,3,6,5,8,7,&
       &   2,1,4,3,6,5,8,7,&
       &   3,4,1,2,7,8,5,6,&
       &   3,4,1,2,7,8,5,6,&
       &   5,6,7,8,1,2,3,4,&
       &   5,6,7,8,1,2,3,4/),(/8,6/))

  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel)==0)return

  hash_nbor(0)=ilevel

  ! Loop over steps
  do ismooth=1,ndim

     ! Initialize flag2 to 0
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%flag2(ind)=0
        end do
     end do

#ifndef WITHOUTMPI

     ! Update boundary conditions
     countrecv=0
     do icpu=1,ncpu
        nbuffer=send_cnt(icpu)
        if(nbuffer>0)then
           countrecv=countrecv+1
           istart=send_oft(icpu)*twotondim+1
           call MPI_IRECV(flag_send_buf(istart),nbuffer*twotondim, &
                & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
        endif
     end do
     
     do ind=1,twotondim
        do i=1,recv_tot
           igrid=grid_recv_buf(i)
           istart=(i-1)*twotondim+ind
           flag_recv_buf(istart)=grid(igrid)%flag1(ind)
        end do
     end do
     
     countsend=0
     do icpu=1,ncpu
        nbuffer=recv_cnt(icpu)
        if(nbuffer>0) then
           countsend=countsend+1
           istart=recv_oft(icpu)*twotondim+1
           call MPI_ISEND(flag_recv_buf(istart),nbuffer*twotondim, &
                & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
        end if
     end do
     
     ! Wait for full completion of receives
     call MPI_WAITALL(countrecv,reqrecv,statuses,info)
     
     do ind=1,twotondim
        do i=1,send_tot
           istart=(i-1)*twotondim+ind
           flag_remote(i,ind)=flag_send_buf(istart)
        end do
     end do
     
     ! Wait for full completion of sends
     call MPI_WAITALL(countsend,reqsend,statuses,info)
     
#endif

     ! Count neighbors and set flag2 accordingly
     do igrid=head(ilevel),tail(ilevel)

        ! Get central oct flag1
        do ind=1,twotondim
           flag_nbor(ind,0)=grid(igrid)%flag1(ind)
        end do

        ! Get neighboring octs flag1
        do inbor=1,twondim
           
           ! Get neighbouring grid using read-only cache
           igridn=nbor_indx(igrid,inbor)
           
           ! If grid exists, then copy into array
           if(igridn>0)then
              do ind=1,twotondim
                 flag_nbor(ind,inbor)=grid(igridn)%flag1(ind)
              end do
           else if (igridn==0) then
              do ind=1,twotondim
                 flag_nbor(ind,inbor)=0
              end do
           else
              do ind=1,twotondim
                 flag_nbor(ind,inbor)=flag_remote(-igridn,ind)
              end do
           endif

        end do
        ! End loop over neighboring octs

        ! Count neighbors and set flag2 accordingly        
        do ind=1,twotondim
           count_nbor=0
           do in=1,twondim
              ig=ggg(ind,in)
              id=hhh(ind,in)
              count_nbor=count_nbor+flag_nbor(id,ig)
           end do
           ! flag2 cell if necessary
           if(count_nbor>=n_nbor(ismooth))then
              grid(igrid)%flag2(ind)=1
           endif
        end do

     end do
     ! End loop over grids

     ! Set flag1=1 for cells with flag2=1
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           if(grid(igrid)%flag1(ind)==1)grid(igrid)%flag2(ind)=0
        end do
        do ind=1,twotondim
           if(grid(igrid)%flag2(ind)==1)then
              grid(igrid)%flag1(ind)=1
              nflag=nflag+1
           endif
        end do
     end do     

  end do
  ! End loop over steps
  
end subroutine smooth_fine_fast
!############################################################
!############################################################
!############################################################
!############################################################
subroutine ensure_ref_rules(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------
  ! This routine determines if all grids at level ilevel are 
  ! surrounded by 26 neighboring grids, in order to enforce the 
  ! strict refinement rule. 
  ! Used in case of adaptive time steps only.
  !-----------------------------------------------------------------
  integer::get_grid
  integer::idim,ind,igrid,ichild
  integer::i1,j1,k1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer(kind=8),dimension(0:ndim)::hash_nbor
  logical::ok

  ! Integer constants
  i1min=0; i1max=0; j1min=0; j1max=0; k1min=0; k1max=0
#if NDIM>0
  i1max=2
#endif
#if NDIM>1
  j1max=2
#endif
#if NDIM>2
  k1max=2
#endif

  call open_cache(operation_smooth,domain_decompos_amr)

  hash_nbor(0)=ilevel
  do igrid=head(ilevel),tail(ilevel)
     
     ok=.true.

     ! Loop over 3x3x3 neighboring father cells
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max

              ! Compute neighboring grid Cartesian index
#if NDIM>0
              hash_nbor(1)=grid(igrid)%ckey(1)+i1-1.0
#endif
#if NDIM>1
              hash_nbor(2)=grid(igrid)%ckey(2)+j1-1.0
#endif
#if NDIM>2
              hash_nbor(3)=grid(igrid)%ckey(3)+k1-1.0
#endif
              ! Periodic boundary conditons
              do idim=1,ndim
                 if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
                 if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
              enddo

              ! Get neighboring grid index
              ichild=get_grid(hash_nbor,grid_dict,.false.,.true.)
              ok=ok.and.(ichild>0)

           end do
        end do
     end do

     if(.not. ok)then
        do ind=1,twotondim
           grid(igrid)%flag1(ind)=0
        end do
     end if

  end do

  call close_cache(grid_dict)

end subroutine ensure_ref_rules 
!############################################################
!############################################################
!############################################################
!############################################################

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
  
subroutine build_smooth(ilevel)
  use amr_commons
  use poisson_commons
  use hilbert
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  
  integer,intent(in)::ilevel
  !
  integer::get_grid
  integer::icoarselevel,igrid,inbor,idim,ipos,ichild,icpu,grid_cpu,ind,info
  integer::i,igridn,iremote
  integer(kind=8),dimension(0:ndim)::hash_key,hash_father,hash_nbor
  integer,dimension(1:ndim)::cart_key
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer(kind=4),dimension(1:nvector),save::dummy_state
  integer(kind=8),dimension(1:nvector,1:nhilbert),save::hk
  integer(kind=8),dimension(1:nvector,1:ndim),save::ix
  integer(kind=8),dimension(1:nhilbert),save::hks

#ifndef WITHOUTMPI

  allocate(nremote(1:ncpu))  
  allocate(nbor_indx(head(ilevel):tail(ilevel),1:twondim))

  hash_nbor(0)=ilevel
  
  nremote=0
  call open_cache(operation_smooth,domain_decompos_amr)
  
  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)
     
     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo
        
        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)
        
        ! If grid exists, determine if remote
        if(igridn<=ngridmax)then
           nbor_indx(igrid,inbor)=igridn
        else
           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)
           
           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           hks=hk(1,1:nhilbert)
           grid_cpu = domain(ilevel)%get_rank(hks)
           nremote(grid_cpu)=nremote(grid_cpu)+1
        end if

     end do
     ! End loop over neighbors
     
  end do
  ! End loop over grids
  
  call close_cache(grid_dict)

  ! Build communicator
  allocate(nalltoall(1:ncpu,1:ncpu))
  allocate(nalltoall_tot(1:ncpu,1:ncpu))
  nalltoall=0; nalltoall_tot=0
  do icpu=1,ncpu
     nalltoall(myid,icpu)=nremote(icpu)
  end do
  call MPI_ALLREDUCE(nalltoall,nalltoall_tot,ncpu*ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nalltoall=nalltoall_tot
  deallocate(nalltoall_tot)
  allocate(send_cnt(1:ncpu))
  allocate(send_oft(1:ncpu))
  allocate(recv_cnt(1:ncpu))
  allocate(recv_oft(1:ncpu))

  send_cnt=0; send_oft=0; send_tot=0
  recv_cnt=0; recv_oft=0; recv_tot=0
  do icpu=1,ncpu
     send_cnt(icpu)=nalltoall(myid,icpu)
     recv_cnt(icpu)=nalltoall(icpu,myid)
     send_tot=send_tot+send_cnt(icpu)
     recv_tot=recv_tot+recv_cnt(icpu)
     if(icpu<ncpu)then
        send_oft(icpu+1)=send_oft(icpu)+nalltoall(myid,icpu)
        recv_oft(icpu+1)=recv_oft(icpu)+nalltoall(icpu,myid)
     endif
  end do
  deallocate(nalltoall)

#if NDIM>0
  allocate(x_send_buf(1:send_tot))
#endif
#if NDIM>1
  allocate(y_send_buf(1:send_tot))
#endif
#if NDIM>2
  allocate(z_send_buf(1:send_tot))
#endif

  hash_nbor(0)=ilevel

  nremote=0
  call open_cache(operation_smooth,domain_decompos_amr)

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! Gather twondim neighboring grids
     do inbor=1,twondim

        hash_nbor(1:ndim)=grid(igrid)%ckey(1:ndim)+shift(1:ndim,inbor)

        ! Periodic boundary conditons
        do idim=1,ndim
           if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
           if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
        enddo

        ! Get neighbouring grid using read-only cache
        igridn=get_grid(hash_nbor,grid_dict,.false.,.true.)

        ! If grid is remote
        if(igridn>ngridmax)then

           ! Compute Cartesian keys of new oct
           cart_key(1:ndim)=hash_nbor(1:ndim)

           ! Compute Hilbert keys of new octs
           ix(1,1:ndim)=cart_key(1:ndim)
           call hilbert_key(ix,hk,dummy_state,0,ilevel-1,1)

           ! Determine parent processor and increment counter
           hks=hk(1,1:nhilbert)
           grid_cpu = domain(ilevel)%get_rank(hks)
           nremote(grid_cpu)=nremote(grid_cpu)+1
           iremote=send_oft(grid_cpu)+nremote(grid_cpu)
#if NDIM>0
           x_send_buf(iremote)=cart_key(1)
#endif
#if NDIM>1
           y_send_buf(iremote)=cart_key(2)
#endif
#if NDIM>2
           z_send_buf(iremote)=cart_key(3)
#endif
           ! Store negative index
           nbor_indx(igrid,inbor)=-iremote
        end if

     end do
     ! End loop over neighbors

  end do
  ! End loop over grids

  call close_cache(grid_dict)

  deallocate(nremote)

#if NDIM>0
  allocate(x_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(x_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             x_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(x_send_buf)
#endif
#if NDIM>1
  allocate(y_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(y_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             y_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(y_send_buf)
#endif
#if NDIM>2
  allocate(z_recv_buf(1:recv_tot))
  call MPI_ALLTOALLV(z_send_buf,send_cnt,send_oft,MPI_INTEGER, &
       &             z_recv_buf,recv_cnt,recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  deallocate(z_send_buf)
#endif

  allocate(grid_recv_buf(1:recv_tot))

  hash_key(0)=ilevel
  do i=1,recv_tot
     hash_key(1)=x_recv_buf(i)
     hash_key(2)=y_recv_buf(i)
     hash_key(3)=z_recv_buf(i)
     ipos=hash_get(grid_dict,hash_key)
     grid_recv_buf(i)=ipos
  end do
  deallocate(x_recv_buf)
  deallocate(y_recv_buf)
  deallocate(z_recv_buf)

  allocate(flag_send_buf(1:send_tot*twotondim))
  allocate(flag_recv_buf(1:recv_tot*twotondim))
  allocate(flag_remote(1:send_tot,1:twotondim))
#endif

end subroutine build_smooth

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine clean_smooth
  use amr_commons
  use poisson_commons
  implicit none

#ifndef WITHOUTMPI
  
  deallocate(nbor_indx)
  deallocate(flag_remote)
  deallocate(flag_send_buf)
  deallocate(flag_recv_buf)
  deallocate(grid_recv_buf)
  deallocate(send_cnt)
  deallocate(send_oft)
  deallocate(recv_cnt)
  deallocate(recv_oft)

#endif

end subroutine clean_smooth


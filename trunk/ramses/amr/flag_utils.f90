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
  integer::ilevel,icount
  !--------------------------------------------------------
  ! This routine builds the refinement map at level ilevel.
  !--------------------------------------------------------
  integer::iexpand

  if(ilevel==nlevelmax)return
  if(ilevel<levelmin)return
  if(noct(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Step 1: initialize refinement map to minimal refinement rules
  call init_flag(ilevel)
  if(verbose)write(*,*) '  ==> end step 1',nflag

  ! Step 2: make one cubic buffer around flagged cells,
  ! in order to enforce numerical rule.
  call smooth_fine(ilevel)
  if(verbose)write(*,*) '  ==> end step 2',nflag

  ! Step 3: if cell satisfies user-defined physical citeria,
  ! then flag cell for refinement.
  call userflag_fine(ilevel)    
  if(verbose)write(*,*) '  ==> end step 3',nflag

  ! Step 4: make nexpand cubic buffers around flagged cells.
  do iexpand=1,nexpand(ilevel)
     call smooth_fine(ilevel)
  end do
  if(verbose)write(*,*) '  ==> end step 4',nflag

  if(verbose)write(*,112)nflag

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
  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ichild=head(ilevel+1),tail(ilevel+1)
     hash_key(1:ndim)=grid(ichild)%ckey(1:ndim)
     parent_cell=get_parent_cell(hash_key)
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
  if(noct(ilevel)==0)return
 
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
  integer::i,ncache,iskip,ngrid
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
  if(noct(ilevel)==0)return

  hash_nbor(0)=ilevel

  ! Loop over steps
  do ismooth=1,ndim

     ! Initialize flag2 to 0
     do igrid=head(ilevel),tail(ilevel)
        do ind=1,twotondim
           grid(igrid)%flag2(ind)=0
        end do
     end do

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
           igridn(i_nbor)=get_grid(hash_nbor)
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

end subroutine smooth_fine
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
              ichild=get_grid(hash_nbor)
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

end subroutine ensure_ref_rules 
!############################################################
!############################################################
!############################################################
!############################################################

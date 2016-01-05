subroutine refine_all
  use amr_commons
  implicit none

  integer::ilevel

  if(verbose)write(*,*)'Entering refine_all' 

  do ilevel=levelmin,nlevelmax-1
     call refine_fine(ilevel)
     call build_comm(ilevel+1)
     call make_virtual_fine_int(cpu_map(1),ilevel+1)
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
  integer::ncache,ngrid
  integer::igrid,icell,i
  integer::ind,iskip,info,icpu,ibound
  integer::ncreate_tmp,nkill_tmp
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,1:ndim)::cart_key
  logical::ok_free,ok_all,ok

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !--------------------------
  ! Compute authorization map
  !--------------------------
  call authorize_fine(ilevel)

  !---------------------------------------------------
  ! Step 1: if cell is flagged for refinement and
  ! if it is not already refined, create a son grid.
  !---------------------------------------------------        
  ncreate=0
  ifree=noct_tot+1
  do ilev=ilevel,nlevelmax
     imin=head(ilev)
     ntot=noct(ilev)
     imax=imin+noct-1
     do ioct=imin,imax
        do ind=1,twotondim
           ok   = grid(ioct)%flag2(ind)==1 .and. &
                & grid(ioct)%flag1(ind)==1 .and. &
                & .not.grid(ioct)%refined(ind)
           if(ok)then
              ind_child(1)=ifree
              ind_parent(1)=ioct
              ind_cell(1)=ind
              call make_new_oct(ind_child,ind_parent,ind_cell,ilev+1,1)
              ncreate=ncreate+1
              ifree=ifree+1
              if(ifree.GT.ngridmax)then
                 write(*,*)'No more free memory'
                 write(*,*)'Increase ngridmax'
                 call clean_abort
              end if
        end do
     end do
  end do
  if(verbose)write(*,112)ncreate

  !-----------------------------------------------------
  ! Step 2: if cell is not flagged for refinement, but
  ! it is refined, then destroy its child grid.
  !-----------------------------------------------------
  nkill=0  
  do ilev=ilevel+1,nlevelmax
     imin=head(ilev)
     ntot=noct(ilev)
     imax=imin+noct-1
     do ioct=imin,imax
        do ind=1,twotondim
           parent_cell=get_parent_cell(ioct)
           itile=parent_cell%tile
           igrid=parent_cell%grid
           icell=parent_cell%cell
           ok   = cache(itile)%flag1(igrid,icell)==0 .and. &
                & cache(itile)%refined(igrid,icell)
           if(ok)then
              grid(ioct)%lev=0
              cache(itile)%refined(igrid,icell)=.false.
              cache(itile)%dirty=.true.
              cache(itile)%clean=.false.
              nkill=kill+1
           end if
        end do
     end do
  end do
  if(verbose)write(*,112)nkill

  !-----------------------------------------------------
  ! Step 3: sort new octs and empty slots according to 
  ! the grid level.
  !-----------------------------------------------------


  ! Compute grid number statistics at level ilevel+1
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(1,ilevel+1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(2,ilevel+1),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(3,ilevel+1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#else
  tmp_long=numbl(myid,ilevel+1)
  call MPI_ALLREDUCE(tmp_long,numbtot(1,ilevel+1),1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tmp_long,numbtot(2,ilevel+1),1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tmp_long,numbtot(3,ilevel+1),1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  call MPI_ALLREDUCE(used_mem,used_mem_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbtot(1,ilevel+1)=numbl(myid,ilevel+1)
  numbtot(2,ilevel+1)=numbl(myid,ilevel+1)
  numbtot(3,ilevel+1)=numbl(myid,ilevel+1)
  used_mem_tot=used_mem
#endif
  numbtot(4,ilevel+1)=numbtot(1,ilevel+1)/ncpu

111 format('   Entering refine_fine for level ',I2)
112 format('   ==> Make ',i6,' sub-grids')
113 format('   ==> Kill ',i6,' sub-grids')

end subroutine refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_new_oct(ind_child,ind_parent,ind_cell,ilevel,ngrid)
  use amr_commons
  use hydro_commons, ONLY:uold
  use poisson_commons, ONLY:f, phi,phi_old
  implicit none
  integer::ngrid,ilevel
  integer,dimension(1:nvector)::ind_child,ind_parent,ind_cell
  !--------------------------------------------------------------
  ! This routine creates children octs at level ilevel.
  ! ilevel is thus the level of the new children octs.
  ! The new octs are labelled using the new index ind_child(:).
  ! The parent cells are labeled with their parent oct index ind_parent(:)
  ! and their cell index ind_cell(:) (from 1 to 8).
  !--------------------------------------------------------------
  integer::idim,ind,igrid,iparent,ichild,icell
  integer(kind=4), dimension(1:nvector),save::dummy_state
  integer(kind=8), dimension(1:nvector),save::hk0,hk1,hk2
  integer(kind=8), dimension(1:nvector),save::ix,iy,iz
  integer(kind=8), dimension(1:nvector,1:ndim),save::cart_key
  integer(kind=8), dimension(0:ndim)::hash_key

  !=================================
  ! Create new octs into main memory
  !=================================
  ! Compute Cartesian keys of new octs
  do igrid=1,ngrid
     iparent=ind_parent(igrid)
     icell=ind_cell(igrid)-1
     do idim=1,ndim
        nstride=2**(idim-1)
        cart_key(idim)=2*grid(iparent)%ckey(idim)+MOD(icell/nstride,2)
     end do
  end do

  ! Compute Hilbert keys of new octs
#if NDIM==1
  ix(1:ngrid)=cart_key(1:ngrid,1)
  call hilbert1d(ix,hk0,ngrid)
#endif
#if NDIM==2
  ix(1:ngrid)=cart_key(1:ngrid,1)
  iy(1:ngrid)=cart_key(1:ngrid,2)
  call hilbert2d(ix,iy,hk1,hk0,dummy_state,0,ilevel-1,ngrid)
#endif
#if NDIM==3
  ix(1:ngrid)=cart_key(1:ngrid,1)
  iy(1:ngrid)=cart_key(1:ngrid,2)
  iz(1:ngrid)=cart_key(1:ngrid,3)
  call hilbert3d(ix,iy,iz,hk2,hk1,hk0,dummy_state,0,ilevel-1,ngrid)
#endif

  ! Loop over the new octs
  do igrid=1,ngrid
     ! Insert new oct in main resident memory
     ichild=ind_child(igrid)
     grid(ichild)%ckey(1:ndim)=cart_key(igrid,1:ndim)
     grid(ichild)%hkey=hk0(igrid)
     grid(ichild)%refined(1:twotondim)=.false.
     ! Insert new grid in hash table
     hash_key(0)=ilevel
     hash_key(1:ndim)=cart_key(igrid,1:ndim)
     call hash_set(grid_dict,hash_key,ichild)
  end do

  !=======================================
  ! Set status of parent cell to "refined"
  !=======================================
  ! Loop over the parent octs
  do igrid=1,ngrid
     iparent=ind_parent(igrid)
     icell=ind_cell(igrid)
     grid(iparent)%refined(icell)=.true.
  end do

  !=====================================================
  ! Inject parent hydro variables into new children ones
  !=====================================================
  ! No interpolation yet (to be done later)...
  if(.not.init)then
     !============================
     ! Interpolate hydro variables
     !============================
     if(hydro)then
        do igrid=1,ngrid
           ichild=ind_child(igrid)
           iparent=ind_parent(igrid)
           icell=ind_cell(igrid)
           do ind=1,twotondim
              grid(ichild)%uold(ind,1:nvar)=grid(iparent)%uold(icell,1:nvar)
           enddo
        enddo
     endif
     !==============================
     ! Interpolate gravity variables
     !==============================
     if(poisson)then
        do igrid=1,ngrid
           ichild=ind_child(igrid)
           iparent=ind_parent(igrid)
           icell=ind_cell(igrid)
           do ind=1,twotondim
              grid(ichild)%f(ind,1:ndim)=grid(iparent)%f(icell,1:ndim)
              grid(ichild)%phi(ind)=grid(iparent)%phi(icell)
              grid(ichild)%phi_old(ind)=grid(iparent)%phi_old(icell)
           enddo
        enddo
     endif
  endif

end subroutine make_new_oct
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine kill_oct(ind_grid,ngrid)
  use amr_commons
  use hydro_commons, ONLY:uold
  use poisson_commons, ONLY:f, phi,phi_old
  implicit none
  integer::ngrid
  integer,dimension(1:nvector)::ind_grid
  !--------------------------------------------------------------
  ! This routine creates children octs at level ilevel.
  ! ilevel is thus the level of the new children octs.
  ! The new octs are labelled using the new index ind_child(:).
  ! The parent cells are labeled with their parent oct index ind_parent(:)
  ! and their cell index ind_cell(:) (from 1 to 8).
  !--------------------------------------------------------------
  integer::igrid,ioct

  ! Set level to zero
  do igrid=1,ngrid
     ioct=ind_grid(igrid)
     grid(ioct)%lev=0
  end do

end subroutine kill_oct
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_grid_fine(ind_grid,ind_cell,ind,ilevel,nn,ibound,boundary_region)
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY:f, phi,phi_old
  implicit none
  integer::nn,ind,ilevel,ibound
  logical::boundary_region
  integer,dimension(1:nvector)::ind_grid,ind_cell
  !--------------------------------------------------------------
  ! This routine create new grids at level ilevel (ilevel >= 2)
  ! contained in father cell ind_cell(:).
  ! ind_grid(:) is the number of the grid that contains 
  ! the father cell and ind = 1, 2, 3, 4, 5, 6, 7, or 8.
  ! The actual father cell number is:
  ! ind_cell = ncoarse + (ind-1)*ngridmax + ind_grid
  ! WARNING: USE THIS ROUTINE WITH CARE, SINCE IT ASSUMES THAT
  ! ALL FATHER CELL'S NEIGHBORS DO EXIST !!!
  !--------------------------------------------------------------
  integer ::idim,igrid,iskip,icpu
  integer ::i,j,ix,iy,iz,ivar,nx_loc
  integer ,dimension(1:nvector)          ,save::ind_grid_son
  integer ,dimension(1:nvector,0:twondim),save::ind_fathers
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:twondim),save::indn

  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::xc,skip_loc
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0

  real(dp),dimension(1:nvector,1:ndim),save::xx
  integer ,dimension(1:nvector),save::cc

  logical::error

  ! Mesh spacing in father level
  dx=0.5D0**(ilevel-1)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Get nn new grids from free memory
  do i=1,nn
     igrid=headf
     ind_grid_son(i)=igrid
     headf=next(headf)
     numbf=numbf-1
     used_mem=ngridmax-numbf
  end do
  
  ! Set new grids position
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  if(ndim>0)xc(1)=(dble(ix)-0.5D0)*dx
  if(ndim>1)xc(2)=(dble(iy)-0.5D0)*dx
  if(ndim>2)xc(3)=(dble(iz)-0.5D0)*dx
  do idim=1,ndim
     do i=1,nn
        xg(ind_grid_son(i),idim)=xg(ind_grid(i),idim)+xc(idim)
     end do
  end do

  ! Add grid cells to hash table
  do i = 1, nn
     call add_grid_to_hash_table(ind_grid_son(i), ilevel)
  end do
  
  ! Connect new grids to father cells
  do i=1,nn
     son(ind_cell(i))=ind_grid_son(i)
  end do
  do i=1,nn
     father(ind_grid_son(i))=ind_cell(i)
  end do

  ! Connect news grids to neighboring father cells
  call getnborgrids(ind_grid,igridn,nn)
  call getnborcells(igridn,ind,indn,nn)
  error=.false.
  do j=1,twondim
     do i=1,nn
        nbor(ind_grid_son(i),j)=indn(i,j)
        if(indn(i,j)==0)then
           error=.true.
        end if
     end do
  end do
  if(error)then
     do j=1,twondim
        do i=1,nn
           if(indn(i,j)==0)then
              if(cpu_map(ind_cell(i))==myid)then
                 write(*,*)'Fatal error in make_grid_fine'
                 write(*,*)myid,cpu_map(ind_cell(i))
                 write(*,*)ilevel,j,ibound,boundary_region
                 stop
              endif
           end if
        end do
     end do
  end if
  
  ! Update cpu map
  if(boundary_region)then
     do j=1,twotondim
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,nn
           cpu_map(iskip+ind_grid_son(i))=0
        end do
     end do
  else
     do j=1,twotondim
        iz=(j-1)/4
        iy=(j-1-4*iz)/2
        ix=(j-1-2*iy-4*iz)
        if(ndim>0)xc(1)=(dble(ix)-0.5D0)*dx/2.0d0
        if(ndim>1)xc(2)=(dble(iy)-0.5D0)*dx/2.0d0
        if(ndim>2)xc(3)=(dble(iz)-0.5D0)*dx/2.0d0
        ! Compute cell coordinates
        do idim=1,ndim
           do i=1,nn
              xx(i,idim)=xg(ind_grid_son(i),idim)+xc(idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,nn
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        call cmp_cpumap(xx,cc,nn)
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,nn
           cpu_map(iskip+ind_grid_son(i))=cc(i)
        end do
     end do
  end if

  ! Connect news grids to level ilevel linked list
  if(boundary_region)then
     do i=1,nn
        igrid=ind_grid_son(i)
        if(numbb(ibound,ilevel)>0)then
           next(igrid)=0
           prev(igrid)=tailb(ibound,ilevel)
           next(tailb(ibound,ilevel))=igrid
           tailb(ibound,ilevel)=igrid
           numbb(ibound,ilevel)=numbb(ibound,ilevel)+1
        else
           next(igrid)=0
           prev(igrid)=0
           headb(ibound,ilevel)=igrid
           tailb(ibound,ilevel)=igrid
           numbb(ibound,ilevel)=1
        end if
     end do
  else
     do i=1,nn
        igrid=ind_grid_son(i)
        icpu=cpu_map(ind_cell(i))
        if(numbl(icpu,ilevel)>0)then
           next(igrid)=0
           prev(igrid)=taill(icpu,ilevel)
           next(taill(icpu,ilevel))=igrid
           taill(icpu,ilevel)=igrid
           numbl(icpu,ilevel)=numbl(icpu,ilevel)+1
        else
           next(igrid)=0
           prev(igrid)=0
           headl(icpu,ilevel)=igrid
           taill(icpu,ilevel)=igrid
           numbl(icpu,ilevel)=1
        end if
     end do
  end if

  ! Interpolate parent variables to get new children ones
  if(.not.init .and. .not.balance)then
     ! Get neighboring father cells
     do i=1,nn
        ind_fathers(i,0)=father(ind_grid_son(i))
     end do
     do j=1,twondim
        do i=1,nn
           ind_fathers(i,j)=nbor(ind_grid_son(i),j)
        end do
     end do
     !============================
     ! Interpolate hydro variables
     !============================
     if(hydro)then
        do j=0,twondim
           ! Gather hydro variables
              do ivar=1,nvar
                 do i=1,nn
                    u1(i,j,ivar)=uold(ind_fathers(i,j),ivar)
                 end do
           end do
        end do
        ! Interpolate
        call interpol_hydro(u1,u2,nn)
        ! Scatter to children cells
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
              do ivar=1,nvar
                 do i=1,nn
                    uold(iskip+ind_grid_son(i),ivar)=u2(i,j,ivar)
                 end do
           end do
        enddo
     end if
     !==============================
     ! Interpolate gravity variables
     !==============================
     if(poisson)then
        ! Scatter to children cells
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
           do idim=1,ndim
              do i=1,nn
                 f(iskip+ind_grid_son(i),idim)=f(ind_fathers(i,0),idim)
              end do
           end do
           do i=1,nn
              phi(iskip+ind_grid_son(i))=phi(ind_fathers(i,0))
              phi_old(iskip+ind_grid_son(i))=phi_old(ind_fathers(i,0))
           end do
        end do
     end if
  endif

end subroutine make_grid_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine kill_grid(ind_cell,ilevel,nn,ibound,boundary_region)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use coordinates
  use hash, only: hash_free
  implicit none
  integer::nn,ilevel,ibound
  logical::boundary_region
  integer,dimension(1:nvector)::ind_cell
  !----------------------------------------------------
  ! This routine destroy the grids at level ilevel
  ! contained in father cell ind_cell(:)
  !----------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,j,idim,ind,ivar
  integer,dimension(1:nvector),save::ind_grid_son,ind_cell_son
  integer(kind=8), dimension(1:ndim) :: ix
  
  ! Gather son grids
  do i=1,nn
     ind_grid_son(i)=son(ind_cell(i))
  end do
  
  ! Disconnect son grids from father cells
  do i=1,nn
     son(ind_cell(i))=0
  end do

  ! Disconnect son grids from level ilevel linked list
  if(boundary_region)then
     do i=1,nn
        igrid=ind_grid_son(i)
        if(prev(igrid).ne.0) then
           if(next(igrid).ne.0)then
              next(prev(igrid))=next(igrid)
              prev(next(igrid))=prev(igrid)
           else
              next(prev(igrid))=0
              tailb(ibound,ilevel)=prev(igrid)
           end if
        else
           if(next(igrid).ne.0)then
              prev(next(igrid))=0
              headb(ibound,ilevel)=next(igrid)
           else
              headb(ibound,ilevel)=0
              tailb(ibound,ilevel)=0
           end if
        end if
        numbb(ibound,ilevel)=numbb(ibound,ilevel)-1 
     end do
  else
     do i=1,nn
        igrid=ind_grid_son(i)
        icpu=cpu_map(ind_cell(i))
        if(prev(igrid).ne.0) then
           if(next(igrid).ne.0)then
              next(prev(igrid))=next(igrid)
              prev(next(igrid))=prev(igrid)
           else
              next(prev(igrid))=0
              taill(icpu,ilevel)=prev(igrid)
           end if
        else
           if(next(igrid).ne.0)then
              prev(next(igrid))=0
              headl(icpu,ilevel)=next(igrid)
           else
              headl(icpu,ilevel)=0
              taill(icpu,ilevel)=0
           end if
        end if
        numbl(icpu,ilevel)=numbl(icpu,ilevel)-1 
     end do
  end if

  ! Remove grids from hash table
  do i = 1, nn
     call remove_grid_from_hash_table(ind_grid_son(i), ilevel)
  end do
  
  ! Reset grid variables
  do idim=1,ndim
     do i=1,nn
        xg(ind_grid_son(i),idim)=0.0D0
     end do
  end do
  do i=1,nn
     father(ind_grid_son(i))=0
  end do
  do j=1,twondim
     do i=1,nn
        nbor(ind_grid_son(i),j)=0
     end do
  end do
  if(pic)then
     do i=1,nn
        headp(ind_grid_son(i))=0
        tailp(ind_grid_son(i))=0
        numbp(ind_grid_son(i))=0
     end do
  end if

  ! Reset cell variables
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,nn
        ind_cell_son(i)=iskip+ind_grid_son(i)
     end do
     ! Tree variables
     do i=1,nn
        son     (ind_cell_son(i))=0
        flag1   (ind_cell_son(i))=0
        flag2   (ind_cell_son(i))=0
        cpu_map (ind_cell_son(i))=0
        cpu_map2(ind_cell_son(i))=0
     end do
     ! Gravity variables
     if(poisson)then
        do i=1,nn
           rho(ind_cell_son(i))=0.0D0
           phi(ind_cell_son(i))=0.0D0
           phi_old(ind_cell_son(i))=0.0D0
        end do
        do idim=1,ndim
           do i=1,nn
              f(ind_cell_son(i),idim)=0.0D0
           end do
        end do
     end if
     ! Hydro variables
     if(hydro)then
        do ivar=1,nvar
           do i=1,nn
              uold(ind_cell_son(i),ivar)=0.0D0
              unew(ind_cell_son(i),ivar)=0.0D0
           end do
        end do
     end if
  end do

  ! Put son grids at the tail of the free memory linked list
  do i=1,nn
     igrid=ind_grid_son(i)
     next(tailf)=igrid
     prev(igrid)=tailf
     next(igrid)=0
     tailf=igrid
     numbf=numbf+1
  end do
 
end subroutine kill_grid


subroutine add_grid_to_hash_table(igrid, ilevel)
  use amr_parameters, only: ndim, ngridmax
  use amr_commons,    only: xg, twotondim, ncoarse, cell_dict
  use hash,           only: hash_set
  use coordinates,    only: grid_to_integer
  implicit none

  integer, intent(in) :: igrid, ilevel
  ! Add all cells belonging to igrid to hash table  
  integer :: ind, idim, icell
  integer(kind=8), dimension(0:ndim) :: ix

  do ind = 1, twotondim
     ix(0) = ilevel
     ix(1:ndim) = grid_to_integer(xg(igrid, 1:ndim), ilevel - 1)
     ix(1) = ISHFT(ix(1), 1) + mod(ind - 1, 2)
#if NDIM>1
     ix(2) = ISHFT(ix(2), 1) + mod(ind - 1, 4) / 2
#endif
#if NDIM>2     
     ix(3) = ISHFT(ix(3), 1) + (ind - 1) / 4
#endif
     icell = igrid + ncoarse + ngridmax * (ind - 1)
     call hash_set(cell_dict, ix, icell)
  end do

end subroutine add_grid_to_hash_table

subroutine remove_grid_from_hash_table(igrid, ilevel)
  use amr_parameters, only: ndim, ngridmax
  use amr_commons,    only: xg, twotondim, cell_dict, ncoarse
  use hash,           only: hash_free
  use coordinates,    only: grid_to_integer
  implicit none
  integer, intent(in) :: igrid, ilevel
  ! Remove all cells belonging to igrid to hash table  
  integer :: ind, idim
  integer(kind=8), dimension(0:ndim) :: ix
  
  do ind = 1, twotondim
     ix(0) = ilevel
     ix(1:ndim) = grid_to_integer(xg(igrid, 1:ndim), ilevel - 1)
     ix(1) = ISHFT(ix(1), 1) + mod(ind - 1, 2)
#if NDIM>1
     ix(2) = ISHFT(ix(2), 1) + mod(ind - 1, 4) / 2
#endif
#if NDIM>2 
     ix(3) = ISHFT(ix(3), 1) + (ind - 1) / 4
#endif
     call hash_free(cell_dict, ix)
  end do
end subroutine remove_grid_from_hash_table
  

!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_tree
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------
  ! This subroutine build the particle linked list at the 
  ! coarse level for ALL the particles in the box.
  ! This routine should be used only as initial set up for
  ! the particle tree.
  !------------------------------------------------------
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale

  if(verbose)write(*,*)'  Entering init_tree'

  ! Local constants
  nxny=nx*ny
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !----------------------------------
  ! Initialize particle linked list
  !----------------------------------
  prevp(1)=0; nextp(1)=2
  do ipart=2,npartmax-1
     prevp(ipart)=ipart-1
     nextp(ipart)=ipart+1
  end do
  prevp(npartmax)=npartmax-1; nextp(npartmax)=0
  ! Free memory linked list
  headp_free=npart+1
  tailp_free=npartmax
  numbp_free=tailp_free-headp_free+1
  if(numbp_free>0)then
     prevp(headp_free)=0
  end if
  nextp(tailp_free)=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif

  !--------------
  ! Periodic box
  !--------------
  do idim=1,ndim
     do ipart=1,npart
        if(xp(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
             & xp(ipart,idim)=xp(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
        if(xp(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
             & xp(ipart,idim)=xp(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
     end do
  end do
 
  !----------------------------------
  ! Reset all linked lists at level 1
  !----------------------------------
  do i=1,active(1)%ngrid
     headp(active(1)%igrid(i))=0
     tailp(active(1)%igrid(i))=0
     numbp(active(1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,1)%ngrid
        headp(reception(icpu,1)%igrid(i))=0
        tailp(reception(icpu,1)%igrid(i))=0
        numbp(reception(icpu,1)%igrid(i))=0
     end do
  end do

  !------------------------------------------------
  ! Build linked list at level 1 by vector sweeps
  !------------------------------------------------
  do ipart=1,npart,nvector
     npart1=min(nvector,npart-ipart+1)
     ! Gather particles
     do i=1,npart1
        ind_part(i)=ipart+i-1
     end do
     ! Compute coarse cell
#if NDIM>0
     do i=1,npart1
        ix(i)=int(xp(ind_part(i),1)/scale+skip_loc(1))
     end do
#endif
#if NDIM>1
     do i=1,npart1
        iy(i)=int(xp(ind_part(i),2)/scale+skip_loc(2))
     end do
#endif
#if NDIM>2
     do i=1,npart1
        iz(i)=int(xp(ind_part(i),3)/scale+skip_loc(3))
     end do
#endif
     ! Compute level 1 grid index
     error=.false.
     do i=1,npart1
        ind_grid(i)=son(1+ix(i)+nx*iy(i)+nxny*iz(i))
        if(ind_grid(i)==0)error=.true.
     end do
     if(error)then
        write(*,*)'Error in init_tree'
        write(*,*)'Particles appear in unrefined regions'
        call clean_stop
     end if
     ! Add particle to level 1 linked list
     call add_list(ind_part,ind_grid,ok,npart1)
  end do

  ! Sort particles down to levelmin
  do ilevel=1,levelmin-1
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end do

end subroutine init_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine checks if particles have moved from their parent grid
  ! to one of the 3**ndim neighboring sister grids. The particle is then 
  ! disconnected from the parent grid linked list, and connected to the
  ! corresponding sister grid linked list. If the sister grid does
  ! not exist, the particle is left to its original parent grid.
  ! Particles must not move to a distance greater than direct neighbors
  ! boundaries. Otherwise an error message is issued and the code stops.
  !-----------------------------------------------------------------------
  integer::idim,nx_loc
  real(dp)::dx,scale
  real(dp),dimension(1:3)::xbound
  real(dp),dimension(1:3)::skip_loc
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              ! Gather nvector particles
              if(ip==nvector)then
                 call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering make_tree_fine for level ',I2)

end subroutine make_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by make_tree_fine.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  real(dp),dimension(1:3)::xbound
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_father
  ! Particle-based arrays
  integer,dimension(1:nvector),save::ind_son,igrid_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_father(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_father,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Compute particle position in 3-cube
  error=.false.
  ind_son(1:np)=1
  ok(1:np)=.false.
  do idim=1,ndim
     do j=1,np
        i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
        if(i<0.or.i>2)error=.true.
        ind_son(j)=ind_son(j)+i*3**(idim-1)
        ! Check if particle has escaped from its parent grid
        ok(j)=ok(j).or.i.ne.1 
     end do
  end do
  if(error)then
     write(*,*)'Problem in check_tree at level ',ilevel
     write(*,*)'A particle has moved outside allowed boundaries'
     do idim=1,ndim
        do j=1,np
           i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
           if(i<0.or.i>2)then
              write(*,*)xp(ind_part(j),idim),x0(ind_grid_part(j),idim)
           endif
        end do
     end do
     stop
  end if

  ! Compute neighboring grid index
  do j=1,np
     igrid_son(j)=son(nbors_father_cells(ind_grid_part(j),ind_son(j)))
  end do

  ! If escaped particle sits in unrefined cell, leave it to its parent grid.
  ! For ilevel=levelmin, this should never happen.
  do j=1,np
     if(igrid_son(j)==0)ok(j)=.false.
  end do

  ! Periodic box
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xxx=xp(ind_part(j),idim)/scale+skip_loc(idim)-xg(igrid_son(j),idim)
           if(xxx> xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)-(xbound(idim)-skip_loc(idim))*scale
           endif
           if(xxx<-xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)+(xbound(idim)-skip_loc(idim))*scale
           endif
        endif
     enddo
  enddo

  ! Switch particles linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do
  call remove_list(ind_part,list1,ok,np)
  call add_list(ind_part,list2,ok,np)

end subroutine check_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine sorts particle between ilevel grids and their 
  ! ilevel+1 children grids. Particles are disconnected from their parent 
  ! grid linked list and connected to their corresponding child grid linked 
  ! list. If the  child grid does not exist, the particle is left to its 
  ! original parent grid. 
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel+1)==0)return
  if(verbose)write(*,111)ilevel

  ! Reset all linked lists at level ilevel+1
  do i=1,active(ilevel+1)%ngrid
     headp(active(ilevel+1)%igrid(i))=0
     tailp(active(ilevel+1)%igrid(i))=0
     numbp(active(ilevel+1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,ilevel+1)%ngrid
        headp(reception(icpu,ilevel+1)%igrid(i))=0
        tailp(reception(icpu,ilevel+1)%igrid(i))=0
        numbp(reception(icpu,ilevel+1)%igrid(i))=0
     end do
  end do

  ! Sort particles between ilevel and ilevel+1

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then
                 call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering kill_tree_fine for level ',I2)

end subroutine kill_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine kill_tree_fine.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel   
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute lower left corner of grid
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-dx
     end do
  end do

  ! Select only particles within grid boundaries
  ok(1:np)=.true.
  do idim=1,ndim
     do j=1,np
        xxx=(xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx
        ok(j)=ok(j) .and. (xxx >= 0.d0 .and. xxx < 2.0d0)
     end do
  end do
  
  ! Determines in which son particles sit
  ind_son(1:np)=0
  do idim=1,ndim
     do j=1,np
        i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx)
        ind_son(j)=ind_son(j)+i*2**(idim-1)
     end do
  end do 
  do j=1,np
     ind_son(j)=ncoarse+ind_son(j)*ngridmax+ind_grid(ind_grid_part(j))
  end do

  ! Determine which son cell is refined
  igrid_son(1:np)=0
  do j=1,np
     if(ok(j))igrid_son(j)=son(ind_son(j))
  end do
  do j=1,np
     ok(j)=igrid_son(j)>0
  end do

  ! Compute particle linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do

  ! Remove particles from their original linked lists
  call remove_list(ind_part,list1,ok,np)
  ! Add particles to their new linked lists
  call add_list(ind_part,list2,ok,np)

end subroutine kill_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine disconnects all particles contained in children grids
  ! and connects them to their parent grid linked list.
  !---------------------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,ind,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_grid_son
  logical,dimension(1:nvector),save::ok

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     if(icpu==myid)then
        ncache=active(ilevel)%ngrid
     else
        ncache=reception(icpu,ilevel)%ngrid
     end if
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(icpu==myid)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
        end if
        ! Loop over children grids
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              ind_grid_son(i)=son(ind_cell(i))
           end do
           do i=1,ngrid
              ok(i)=ind_grid_son(i)>0
           end do
           do i=1,ngrid
           if(ok(i))then
           if(numbp(ind_grid_son(i))>0)then
              if(numbp(ind_grid(i))>0)then
                 ! Connect son linked list at the tail of father linked list
                 nextp(tailp(ind_grid(i)))=headp(ind_grid_son(i))
                 prevp(headp(ind_grid_son(i)))=tailp(ind_grid(i))
                 numbp(ind_grid(i))=numbp(ind_grid(i))+numbp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
              else
                 ! Initialize father linked list
                 headp(ind_grid(i))=headp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
                 numbp(ind_grid(i))=numbp(ind_grid_son(i))
              end if

           end if
           end if
           end do
        end do
        ! End loop over children
     end do
     ! End loop over grids
  end do
  ! End loop over cpus
  
111 format('   Entering merge_tree_fine for level ',I2)
  
end subroutine merge_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine move particles across processors boundaries.
  !-----------------------------------------------------------------------
  integer::igrid,ipart,jpart,ncache_tot,next_part
  integer::ip,ipcom,npart1,icpu,ncache
  integer::info,buf_count,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
#endif
  integer,dimension(1:nvector),save::ind_part,ind_list,ind_com
  logical::ok_free,ok_all
  integer::particle_data_width

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifdef WITHOUTMPI
  return
#endif

#ifndef WITHOUTMPI
  ! Count particle sitting in virtual boundaries
  do icpu=1,ncpu
     reception(icpu,ilevel)%npart=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        reception(icpu,ilevel)%npart=reception(icpu,ilevel)%npart+&
             & numbp(reception(icpu,ilevel)%igrid(igrid))
     end do
     sendbuf(icpu)=reception(icpu,ilevel)%npart
  end do

  ! Calculate how many particle properties are being transferred
  particle_data_width = twondim+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  particle_data_width=particle_data_width+1
#endif
  
  ! Allocate communication buffer in emission
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(reception(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(reception(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Gather particle in communication buffer
  do icpu=1,ncpu
     if(reception(icpu,ilevel)%npart>0)then
     ! Gather particles by vector sweeps
     ipcom=0
     ip=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        npart1=numbp(reception(icpu,ilevel)%igrid(igrid))
        ipart =headp(reception(icpu,ilevel)%igrid(igrid))
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <--- Very important !!!
           next_part=nextp(ipart)
           ip=ip+1
           ipcom=ipcom+1
           ind_com (ip)=ipcom
           ind_part(ip)=ipart
           ind_list(ip)=reception(icpu,ilevel)%igrid(igrid)
           reception(icpu,ilevel)%fp(ipcom,1)=igrid
           if(ip==nvector)then
              call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
              ip=0
           end if
           ipart=next_part  ! Go to next particle
        end do
     end do
     if(ip>0)call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
     end if
  end do

  ! Communicate virtual particle number to parent cpu
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate communication buffer in reception
  do icpu=1,ncpu
     emission(icpu,ilevel)%npart=recvbuf(icpu)
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(emission(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(emission(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Receive particles
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countrecv=countrecv+1
#ifndef LONGINT
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif
        buf_count=ncache*particle_data_width
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send particles
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countsend=countsend+1
#ifndef LONGINT
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
        buf_count=ncache*particle_data_width
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Compute total number of newly created particles
  ncache_tot=0
  do icpu=1,ncpu
     ncache_tot=ncache_tot+emission(icpu,ilevel)%npart
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
  ok_free=(numbp_free-ncache_tot)>=0
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
     write(*,*)numbp_free,ncache_tot
     write(*,*)myid
     write(*,*)emission(1:ncpu,ilevel)%npart
     write(*,*)'============================'
     write(*,*)reception(1:ncpu,ilevel)%npart
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  end if

  ! Scatter new particles from communication buffer
  do icpu=1,ncpu
     ! Loop over particles by vector sweeps
     ncache=emission(icpu,ilevel)%npart
     do ipart=1,ncache,nvector
        npart1=min(nvector,ncache-ipart+1)
        do ip=1,npart1
           ind_com(ip)=ipart+ip-1
        end do
        call empty_comm(ind_com,npart1,ilevel,icpu)
     end do
  end do

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(emission(icpu,ilevel)%fp)
        deallocate(emission(icpu,ilevel)%up)
     end if
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(reception(icpu,ilevel)%fp)
        deallocate(reception(icpu,ilevel)%up)
     end if
  end do
#endif

111 format('   Entering virtual_tree_fine for level ',I2)
end subroutine virtual_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fill_comm(ind_part,ind_com,ind_list,np,ilevel,icpu)
  use pm_commons
  use amr_commons
  implicit none
  integer::np,ilevel,icpu
  integer,dimension(1:nvector)::ind_part,ind_com,ind_list
  integer::current_property
  integer::i,idim
  logical,dimension(1:nvector),save::ok=.true.

  ! Gather particle level and identity
  do i=1,np
     reception(icpu,ilevel)%fp(ind_com(i),2)=levelp(ind_part(i))
     reception(icpu,ilevel)%fp(ind_com(i),3)=idp   (ind_part(i))
  end do
  
  ! Gather particle position and velocity
  do idim=1,ndim
     do i=1,np
        reception(icpu,ilevel)%up(ind_com(i),idim     )=xp(ind_part(i),idim)
        reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=vp(ind_part(i),idim)
     end do
  end do
  
  current_property = twondim+1
  ! Gather particle mass
  do i=1,np
     reception(icpu,ilevel)%up(ind_com(i),current_property)=mp(ind_part(i))
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Gather particle potential
  do i=1,np
     reception(icpu,ilevel)%up(ind_com(i),current_property)=ptcl_phi(ind_part(i))
  end do
  current_property = current_property+1
#endif
  
  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_list,ok,np)
  call add_free(ind_part,np)
  
end subroutine fill_comm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine empty_comm(ind_com,np,ilevel,icpu)
  use pm_commons
  use amr_commons
  implicit none
  integer::np,icpu,ilevel
  integer,dimension(1:nvector)::ind_com
  
  integer::i,idim,igrid
  integer,dimension(1:nvector),save::ind_list,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  integer::current_property

  ! Compute parent grid index
  do i=1,np
     igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
     ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
  end do

  ! Add particle to parent linked list
  call remove_free(ind_part,np)
  call add_list(ind_part,ind_list,ok,np)

  ! Scatter particle level and identity
  do i=1,np
     levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
     idp   (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
  end do

  ! Scatter particle position and velocity
  do idim=1,ndim
  do i=1,np
     xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
     vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
  end do
  end do

  current_property = twondim+1

  ! Scatter particle mass
  do i=1,np
     mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Scatter particle phi
  do i=1,np
     ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
  end do
  current_property = current_property+1
#endif

end subroutine empty_comm
!################################################################
!################################################################
!################################################################
!################################################################
! subroutine memory_sort_level(ind_com,np,ilevel,icpu)
!   use pm_commons
!   use amr_commons

!   ! sort all the particles in memory according to their level

!   implicit none
!   integer::np,icpu,ilevel
!   integer,dimension(1:nvector)::ind_com
  
!   integer::i,idim,igrid
!   integer,dimension(1:nvector),save::ind_list,ind_part
!   logical,dimension(1:nvector),save::ok=.true.
!   integer::current_property

!   ! Compute parent grid index
!   do i=1,np
!      igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
!      ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
!   end do

!   ! Add particle to parent linked list
!   call remove_free(ind_part,np)
!   call add_list(ind_part,ind_list,ok,np)

!   ! Scatter particle level and identity
!   do i=1,np
!      levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
!      idp   (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
!   end do

!   ! Scatter particle position and velocity
!   do idim=1,ndim
!   do i=1,np
!      xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
!      vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
!   end do
!   end do

!   current_property = twondim+1

!   ! Scatter particle mass
!   do i=1,np
!      mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1

! #ifdef OUTPUT_PARTICLE_POTENTIAL
!   ! Scatter particle phi
!   do i=1,np
!      ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1
! #endif

! end subroutine memory_sort_level
! !################################################################
! !################################################################
! !################################################################
! !################################################################

! subroutine get_particle_levels

! end subroutine get_particle_levels
! !################################################################
! !################################################################
! !################################################################
! !################################################################
subroutine build_particle_communicator
  use amr_commons,   only:ncpu, myid, bound_key
  use pm_commons
  use pm_parameters, only:npartmax
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  ! ncpu^2 -> ugly, only for a start, replace by point to point communication later
  integer,dimension(1:ncpu,1:ncpu)::npart_alltoall, npart_alltoall_tot
  integer::receive_cpu, ipart, i, info, icpu
  
  npart_alltoall=0
  receive_cpu=1

  do ipart=1,npart_andreas
     do while (part_hkey(ipart,0) > int(bound_key(receive_cpu),kind=8)) 
        receive_cpu=receive_cpu+1
     end do
     npart_alltoall(myid,receive_cpu)=npart_alltoall(myid,receive_cpu)+1
  end do
  npart_alltoall(myid,myid)=0
  call MPI_ALLREDUCE(npart_alltoall,npart_alltoall_tot,ncpu*ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npart_alltoall=npart_alltoall_tot
  
  if(.not. allocated(part_send_cnt))then
     allocate(part_send_cnt(1:ncpu),part_send_oft(1:ncpu))
     allocate(part_recv_cnt(1:ncpu),part_recv_oft(1:ncpu))
  endif
  part_send_cnt=0; part_send_oft=0; part_send_tot=0
  part_recv_cnt=0; part_recv_oft=0; part_recv_tot=0
  do icpu=1,ncpu
     part_send_cnt(icpu)=npart_alltoall(myid,icpu)
     part_recv_cnt(icpu)=npart_alltoall(icpu,myid)
     part_send_tot=part_send_tot+part_send_cnt(icpu)
     part_recv_tot=part_recv_tot+part_recv_cnt(icpu)
     if(icpu<ncpu)then
        part_send_oft(icpu+1)=part_send_oft(icpu)+npart_alltoall(myid,icpu)
        part_recv_oft(icpu+1)=part_recv_oft(icpu)+npart_alltoall(icpu,myid)
     endif
  end do

  !maybe allocate in effective communication routine

  if(allocated(receive_keys))then
     deallocate(receive_keys)
  endif
  if(allocated(send_keys))then
     deallocate(send_keys)
  endif
  allocate(receive_keys(1:part_recv_tot,0:2))
  allocate(send_keys(1:part_send_tot,0:2))
 
  
  
end subroutine build_particle_communicator
!################################################################
!################################################################
!################################################################
!################################################################
subroutine build_histogram_communicator(ilevel)
  use amr_commons,   only: ncpu, myid, bound_key_level, bound_key
  use pm_commons
  use pm_parameters, only: npartmax
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) ::  ilevel
  
  !----------------------------------------------------------------------------
  ! This routine sets up the communication structure for histogrammed particle
  ! quantities. The bins are assumed to be sorted by hilbert key
  !----------------------------------------------------------------------------

    integer::receive_cpu, ibin, i, info, icpu, idest, isource, status
  integer::countrecv, countsend

  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv


  if (nlevelmax>20)then 
     print*, 'problem here with prcision'
  end if


  if(.not. allocated(bin_send_cnt))then
     allocate(bin_send_cnt(1:ncpu),bin_send_oft(1:ncpu))
     allocate(bin_recv_cnt(1:ncpu),bin_recv_oft(1:ncpu))
  endif
  bin_send_cnt=0; bin_send_oft=0; bin_send_tot=0
  bin_recv_cnt=0; bin_recv_oft=0; bin_recv_tot=0

  receive_cpu=1; mybins = 0
  do ibin=1,nbins
     do while (bin_keys(ibin,0) > bound_key_level(receive_cpu, ilevel) ) 
        receive_cpu=receive_cpu+1
     end do
     if (receive_cpu == myid) then
        mybins = mybins +1
     else
        bin_send_cnt(receive_cpu) = bin_send_cnt(receive_cpu)+1
     end if
     if (mybins == 1) mybins_offset = ibin - 1
  end do
 
  countrecv = 0; countsend = 0

  do isource=1,ncpu        
     if (isource .ne. myid) then
        countrecv = countrecv + 1
        call MPI_IRECV(bin_recv_cnt(isource), 1, MPI_INTEGER, isource -1, 1234, MPI_COMM_WORLD, reqrecv(countrecv), info)
     end if
  end do
  do idest=1,ncpu
     if (idest .ne. myid) then
        countsend = countsend + 1
        call MPI_ISEND(bin_send_cnt(idest), 1, MPI_INTEGER, idest -1, 1234, MPI_COMM_WORLD, reqsend(countsend), info)
     end if
  end do

  call MPI_WAITALL(ncpu-1,reqrecv,statuses,info)
  call MPI_WAITALL(ncpu-1,reqsend,statuses,info)

  do icpu=1,ncpu-1
     bin_send_oft(icpu+1)=bin_send_oft(icpu)+bin_send_cnt(icpu)
     bin_recv_oft(icpu+1)=bin_recv_oft(icpu)+bin_recv_cnt(icpu)
  end do

  bin_send_tot=sum(bin_send_cnt); bin_recv_tot=sum(bin_recv_cnt)

  !maybe allocate in effective communication routine

  if(allocated(recv_bin_keys))then
     deallocate(recv_bin_keys)
     deallocate(recv_bin_mass)
     deallocate(send_bin_keys)
     deallocate(send_bin_mass)
  endif

  allocate(send_bin_keys(1:bin_send_tot,0:2))
  allocate(send_bin_mass(1:bin_send_tot))
  allocate(recv_bin_keys(1:bin_recv_tot,0:2))
  allocate(recv_bin_mass(1:bin_recv_tot))

end subroutine build_histogram_communicator
!################################################################
!################################################################
!################################################################
!################################################################
subroutine send_histogram_bins(offset, np, ilevel)
  use amr_commons,   only: ncpu, myid, bound_key_level, son
  use pm_commons
  use pm_parameters, only: npartmax
  use sort,          only: gt_3keys, apply_particle_permutation
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel, offset, np
  integer :: ipos, i, info, ibin, receive_cpu, nb, ipart, ip
  integer :: request, np_small
  integer :: status(MPI_STATUS_SIZE)
  real(dp) :: told
  integer  :: unrefined_pos, refined_pos
  logical :: unrefined

  integer, allocatable, dimension(:) :: bin_index, bin_level
  integer, allocatable, dimension(:) :: recv_bin_index, recv_bin_level
  integer, allocatable, dimension(:) :: send_bin_level


  allocate(bin_index(1:nbins))
  allocate(bin_level(1:nbins))
  allocate(recv_bin_index(1:bin_recv_tot))
  allocate(recv_bin_level(1:bin_recv_tot))
  allocate(send_bin_level(1:bin_send_tot))

  ! Fill communication buffer
  receive_cpu=1
  ipos = 1
  do ibin=1,nbins
     do while (bin_keys(ibin,0) > bound_key_level(receive_cpu, ilevel) ) 
        receive_cpu=receive_cpu+1
     end do
     if (receive_cpu .ne. myid)then
        send_bin_keys(ipos, 0:2) = bin_keys(ibin, 0:2)
        ipos = ipos + 1
     end if
  end do
  
 
  call MPI_IALLTOALLV(send_bin_keys(1:bin_send_tot,0:2),3*bin_send_cnt,3*bin_send_oft,MPI_INTEGER8, &
       &             recv_bin_keys(1:bin_recv_tot,0:2),3*bin_recv_cnt,3*bin_recv_oft,MPI_INTEGER8,MPI_COMM_WORLD,request,info)


  ! Probe domestic cells for refinement 
  do ibin = mybins_offset, mybins_offset + mybins -1, nvector
     nb = min(mybins_offset + mybins - ibin, nvector) 
     call get_cell_index_from_hilbertkey(bin_index(ibin + 1 : ibin + nb), bin_level(ibin + 1 : ibin + nb), &
          bin_keys(ibin + 1: ibin + nb,2), bin_keys(ibin + 1: ibin + nb,1), bin_keys(ibin + 1: ibin + nb,0), nb, ilevel)
  end do

  ! Mark bins corresponding to refined cells
  bin_level=0
  do ibin = mybins_offset + 1, mybins_offset + mybins
     if (son(bin_index(ibin))>0)then
        bin_level(ibin) = 1
     end if
  end do


  ! Finish first communication
  call MPI_WAIT(request,status,info)
  

  ! Probe foreign cells for refinement
  do ibin = 0, bin_recv_tot - 1 , nvector
     nb = min(bin_recv_tot - ibin, nvector) 
     call get_cell_index_from_hilbertkey(recv_bin_index(ibin + 1 : ibin + nb), recv_bin_level(ibin + 1 : ibin + nb), &
          recv_bin_keys(ibin + 1: ibin + nb,2), recv_bin_keys(ibin + 1: ibin + nb,1), recv_bin_keys(ibin + 1: ibin + nb,0), nb, ilevel)
  end do

  recv_bin_level=0
  do ibin = 1, bin_recv_tot
     if (son(recv_bin_index(ibin))>0) recv_bin_level(ibin) = 1
  end do

  ! Send information back
  call MPI_IALLTOALLV(recv_bin_level,bin_recv_cnt,bin_recv_oft,MPI_INTEGER, &
       &             send_bin_level,bin_send_cnt,bin_send_oft,MPI_INTEGER,MPI_COMM_WORLD,request,info)

!  call MPI_ALLTOALLV(recv_bin_level,bin_recv_cnt,bin_recv_oft,MPI_INTEGER, &
!       &             send_bin_level,bin_send_cnt,bin_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)

  

  ! Finish second communicagtion
  call MPI_WAIT(request,status,info)


  ! Read out buffer
  receive_cpu=1
  ipos = 1
  do ibin=1,nbins
     do while (bin_keys(ibin,0) > bound_key_level(receive_cpu, ilevel) ) 
        receive_cpu = receive_cpu + 1
     end do
     if (receive_cpu .ne. myid)then
        bin_level(ibin) = send_bin_level(ipos)
        ipos = ipos + 1
     end if
  end do
  
  ! Find starting indices for refined / unrefined particles

  if (ilevel<nlevelmax)then
     np_small = part_level_offset(ilevel+2) - part_level_offset(ilevel)
     unrefined_pos = offset
     refined_pos = offset
     
     ibin = 1
     unrefined = (bin_level(ibin) == 0)
     do ip = offset + 1, offset + np  
        ipart=part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2)))then
           ibin=ibin+1
           unrefined = (bin_level(ibin) == 0)
        end if
        if (unrefined) then 
           refined_pos = refined_pos + 1
        end if
     end do
     
     part_level_offset(ilevel+1) = refined_pos

     do ip = offset + 1, offset 
        ipart = part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2)))then
           ibin = ibin + 1
           unrefined = (bin_level(ibin) == 1)
        end if
        if (unrefined) then
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos)=ip
        else
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos)=ip
        end if
     end do
    
     part_ind_permutation(offset + 1:offset + np)=part_ind_permutation2(offset +1:offset + np)

     call apply_particle_permutation(offset, np_small, ilevel)
  end if

  deallocate(bin_level, bin_index)
  deallocate(recv_bin_level, recv_bin_index)
  deallocate(send_bin_level)
end subroutine send_histogram_bins
!################################################################
!################################################################
!################################################################
!################################################################
! subroutine part_to_cell_i8(part_array,sortind)
!   use pm_parameters, only: npartmax
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
  
!   integer(kind=8),dimension(1:npartmax)::part_array, sortind
  
!   ! Communication routine that sends a particle-based quantity
!   ! to the MPI domain that owns the respective cell
  
  
  
  

! #ifndef WITHOUTMPI
!   do j=1,part_recv_tot
!      ipart=part_recv_buf(j)-ipart_start(myid)
!      int_part_recv_buf(j)=xx(ipart)
!   end do
!   call MPI_ALLTOALLV(int8_part_recv_buf,part_recv_cnt,part_recv_oft,MPI_INTEGER, &
!        &             int8_part_send_buf,part_send_cnt,part_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  
  
!   deallocate(int8_part_send_buf,int8_part_recv_buf)
! #endif
  
! end subroutine part_to_cell_i8


! !################################################################
! !################################################################
! !################################################################
! !################################################################
!subroutine part_to_cell_dp
  

!  implicit none
  

!end subroutine part_to_cell_dp
! !################################################################
! !################################################################
! !################################################################
! !################################################################

! !################################################################
! !################################################################
! !################################################################
! !################################################################
!subroutine get_cell_index_from_hilbert(cell_index,hkey0,hkey1,hkey2,np,ilevel)
!   use amr_commons
!   ! assuming the last 3 digits are in cartesian order
!   do i=1,np
!      grid_index=hash_get(hkey0(i))
!      ind=mod(hkey0(i),8)
!      iskip=ncoarse+(ind-1)*ngridmax
!      cell_index(i)=iskip+ind_grid
!   end do
! end subroutine get_cell_index

subroutine get_cell_index_from_hilbertkey(cell_index,cell_levl,hilbert_key2,hilbert_key1,hilbert_key0,np,ilevel)
  use amr_commons, only: nlevelmax, nvector, myid
  use hilbert,     only: hilbert3d_reverse
  implicit none
  integer, intent(in)::np,ilevel
  integer(kind=8),dimension(1:nvector)::x,y,z
  integer(kind=8),dimension(1:nvector)::hilbert_key2,hilbert_key1,hilbert_key0
  integer,dimension(1:nvector)::cell_levl, cell_index
  integer::bit_length

  bit_length=ilevel

  call hilbert3d_reverse(x,y,z,hilbert_key2,hilbert_key1,hilbert_key0,bit_length,np)

  call get_cell_index_from_cartesian(cell_index,cell_levl,x,y,z,ilevel,np,bit_length)
  
end subroutine get_cell_index_from_hilbertkey


subroutine get_cell_index_from_cartesian(cell_index,cell_levl,xx,yy,zz,ilevel,n,bit_length)
  use amr_commons
  implicit none

  integer::n,ilevel,bit_length
  integer,dimension(1:nvector)::cell_index,cell_levl
  integer(kind=8),dimension(1:nvector)::xx,yy,zz
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  integer::i,j,ind,iskip,igrid,ind_cell,igrid0
  integer(kind=8)::ii,jj,kk

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     stop
  end if

  if (bit_length>21)then
     print*, 'bit length too big for now'
  end if
  
  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     igrid=igrid0
     do j=1,ilevel 
        ii=ISHFT(xx(i),-bit_length+j)
        jj=ISHFT(yy(i),-bit_length+j)
        kk=ISHFT(zz(i),-bit_length+j)
        ii=mod(ii,2)
        jj=mod(jj,2)
        kk=mod(kk,2)
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index_from_cartesian

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)

     igrid=igrid0
     do j=1,ilevel 
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_particle_histogram(offset, np)
  use pm_commons 
  use amr_commons
  use sort,      only: gt_3keys
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! This routine computes particle histograms. It assumes that particles are 
  ! sorted in memory by hilbert key.  
  !----------------------------------------------------------------------------

  integer,                         save :: ibin, ipart, ip
  integer(kind=8), dimension(0:2), save :: current_bin_key


  ! if there is nothing to do...
  if (.not. np > 0)return
  
  ! Count the number of bins
  nbins=1
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        nbins=nbins+1
        current_bin_key(0:2) = part_hkey(ipart,0:2)
     end if
  end do
  
  ! allocate histograms if necessary
  if(size(bin_mass) < nbins)then
     deallocate(bin_keys)
     deallocate(bin_mass)
     deallocate(bin_count)
  end if
  if (.not. allocated(bin_mass))then
     allocate(bin_keys(nbins,0:2))
     allocate(bin_mass(nbins))
     allocate(bin_count(nbins))
  end if
  bin_mass=0.d0
  bin_count=0
  

  ! label every bin by a key and sum up the particle mass per bin

  ! first particle
  ibin=1
  bin_keys(1,0:2)=part_hkey(part_ind_permutation(offset+1),0:2)  
  bin_mass(1)=mp(part_ind_permutation(offset+1))
  bin_count(1)=1

  ! all other bins
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip=offset+2, offset+np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        ibin=ibin+1
        bin_keys(ibin,0:2)=part_hkey(ipart,0:2)
        current_bin_key(0:2) = part_hkey(ipart,0:2)
     end if
     ! NGP mass assignement --> replace by CIC by using weights instead of mp
     bin_mass(ibin)=bin_mass(ibin)+mp_andreas(ipart)
     bin_count(ibin)=bin_count(ibin)+1
  end do

end subroutine compute_particle_histogram

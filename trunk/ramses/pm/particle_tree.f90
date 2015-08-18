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

subroutine sort_particles(ilevel, use_histograms)
  use pm_commons,  only: npart_andreas, part_level_offset, &
                         nbins, bin_keys, part_hkey, &
                         part_ind_permutation
  use amr_commons, only: dp, myid, levelmin, nlevelmax, ncpu
  use sort,        only: lsd_radix_sort_particles, apply_particle_permutation
  use hilbert,     only: hilbert_for_particle 
  use particle_communication, only: build_communicator
  implicit none

  integer, intent(in) :: ilevel
  logical, intent(in) :: use_histograms
  
  integer :: ndata_remote, ndata_local, local_oft
  integer :: ilev, offset, np, ip
  integer, dimension(1:ncpu, 1:4) :: communicator
  integer, allocatable, dimension(:) :: refined
  
  ! This routine moves all particles that sit in refined cells at level ilevel
  ! to level ilevel + 1. It then sorts the remaining (true) ilevel particles 
  ! by hilbert key.

  offset = part_level_offset(ilevel)
  np = npart_andreas - part_level_offset(ilevel)

  if (ilevel == levelmin)then
     do ilev=levelmin,nlevelmax
        print*, 'nparts on (going in) ',myid, ilev, part_level_offset(ilev + 1) - part_level_offset(ilev)
     end do
  end if
  
  ! Compute hilbert keys (probably move outside of this routine)
  call hilbert_for_particle(offset, npart_andreas - offset, 0, ilevel)
  
  ! Compute a permutation that sorts ALL particles starting from offset
  call lsd_radix_sort_particles(offset, npart_andreas - offset, ilevel, ilevel, .true.)

  if (use_histograms)then
     call compute_particle_histogram(offset, npart_andreas - offset)          
     call build_communicator(communicator, ndata_remote, &
                             nbins, ndata_local, local_oft, &
                             bin_keys(1:nbins, 2), bin_keys(1:nbins, 1), bin_keys(1:nbins, 0), &
                             ilevel)
     allocate(refined(1:nbins))     
     call communicate_refinements(communicator, ndata_remote, &
                                  nbins, ndata_local, local_oft, refined, &
                                  bin_keys(1:nbins, 2), bin_keys(1:nbins, 1), bin_keys(1:nbins, 0), &
                                  ilevel)
     call reshuffle_particles(ilevel, np, nbins, refined, use_histograms)
  else

     ! Need to apply particle permutation here already to have parts sorted
     ! in memory. Using the index insidet build_communicator and
     ! communicate_refinements is possible but will let the code deviate more
     ! from the histogrammed case.
     call apply_particle_permutation(offset, npart_andreas - offset, ilevel) 
     do ip = offset + 1, npart_andreas
        part_ind_permutation(ip) = ip
     end do
     call build_communicator(communicator, ndata_remote, &
                             np, ndata_local, local_oft, &
                             part_hkey(offset + 1: offset + np, 2), &
                             part_hkey(offset + 1: offset + np, 1), &
                             part_hkey(offset + 1: offset + np, 0), ilevel)
     allocate(refined(1:np))
     call communicate_refinements(communicator, ndata_remote, &
                                  np, ndata_local, local_oft, refined, &
                                  part_hkey(offset + 1: offset + np, 2), &
                                  part_hkey(offset + 1: offset + np, 1), &
                                  part_hkey(offset + 1: offset + np, 0), ilevel)

     call reshuffle_particles(ilevel, np, np, refined, use_histograms)
     print*, myid, 'refined total on level ', ilevel, sum(refined)
  end if

  ! Compute NEW number of particles in ilevel
  np = part_level_offset(ilevel + 1) - part_level_offset(ilevel)  

  ! Re-sort remaining (ilevel particles)
  call lsd_radix_sort_particles(offset, np, ilevel, ilevel, .true.)
  call apply_particle_permutation(offset, np, ilevel)
  deallocate(refined)

  if (ilevel == nlevelmax)then
     do ilev=levelmin,nlevelmax
        print*, 'nparts on (going out) ',myid, ilev, part_level_offset(ilev + 1) - part_level_offset(ilev)
     end do
  end if
  
end subroutine sort_particles
!################################################################
!################################################################
!################################################################
!################################################################
subroutine communicate_refinements(communicator, ndata_remote, ndata, ndata_local, ndata_local_oft, &
                                   refined, keys2, keys1, keys0, ilevel)

  use amr_commons,   only: ncpu, myid, bound_key_level, son, nvector, nlevelmax
  use particle_communication, only: part_data_to_domain_i8, part_data_to_domain_i4, &
                                    domain_data_to_part_i4
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) ::  ilevel, ndata
  integer, intent(in) :: ndata_remote, ndata_local, ndata_local_oft
  integer(kind=8), dimension(1:ndata), intent(in) :: keys2, keys1, keys0
  integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
  integer, dimension(1:ndata), intent(inout) :: refined
  
  ! This routine sorts particles between ilevel and ilevel + 1
  ! (formerly known as kill_tree_fine).

  ! The routine can only run after:
  ! call compute_particle_histogram(ilevel)
  ! call build_communicator(...)

  ! The "bins" here are the histogram bins and correspond to 
  ! actual grid cells. "local_... variables" denote properties in the
  ! MPI process which hosts the particles, while "remote_... variables" 
  ! are used for properties in the MPI process which hosts the corresponding
  ! leaf-cell.


  integer  :: idata, ioft, nd
  integer, dimension(1:nvector) :: dummy_int
  
  integer,         allocatable, dimension(:  ) :: remote_refined
  integer(kind=8), allocatable, dimension(:,:) :: remote_keys

  allocate(remote_refined(1:ndata_remote))
  allocate(remote_keys(1:ndata_remote, 0:2))

  call part_data_to_domain_i8(communicator, keys0, remote_keys(:,0))
  call part_data_to_domain_i8(communicator, keys1, remote_keys(:,1))
  call part_data_to_domain_i8(communicator, keys2, remote_keys(:,2))

  ! Probe local cells for refinement (abuse refined to store cell index)
  do ioft = ndata_local_oft, ndata_local_oft + ndata_local -1, nvector
     nd = min(ndata_local_oft + ndata_local - ioft, nvector) 
     call get_cell_index_from_hilbertkey(refined(ioft + 1 : ioft + nd), &
          dummy_int(1 : nd), keys2(ioft + 1: ioft + nd), &
          keys1(ioft + 1: ioft + nd), keys0(ioft + 1: ioft + nd), nd, ilevel)
  end do

  ! Mark data corresponding to refined cells
  do idata = ndata_local_oft + 1, ndata_local_oft + ndata_local
     if (son(refined(idata)) > 0)then
        refined(idata) = 1
     else
        refined(idata) = 0
     end if
  end do
  
  ! Probe remote cells for refinement (abuse remote_refined to store cell index) 
  do ioft = 0, ndata_remote - 1 , nvector
     nd = min(ndata_remote - ioft, nvector) 
     call get_cell_index_from_hilbertkey(remote_refined(ioft + 1 : ioft + nd), &
          dummy_int(1 : nd), remote_keys(ioft + 1: ioft + nd,2), &
          remote_keys(ioft + 1: ioft + nd,1), remote_keys(ioft + 1: ioft + nd,0), nd, ilevel)
  end do

  ! Mark bins corresponding to refined cells
  do idata = 1, ndata_remote
     if (son(remote_refined(idata))>0)then
        remote_refined(idata) = 1
     else
        remote_refined(idata) = 0
     end if
  end do

  ! Send refinement information back
  call domain_data_to_part_i4(communicator, remote_refined, refined)

  deallocate(remote_refined)
  deallocate(remote_keys)
  
end subroutine communicate_refinements
!################################################################
!################################################################
!################################################################
!################################################################
subroutine reshuffle_particles(ilevel, np, ndata, refined, use_histograms)
  use pm_commons,     only: part_level_offset, part_ind_permutation, part_hkey, &
                            bin_keys, part_ind_permutation2
  use sort,           only: gt_3keys, apply_particle_permutation
  use amr_parameters, only: nlevelmax
  implicit none

  
  integer, intent(in) :: ilevel, np, ndata
  integer, dimension(1:ndata), intent(in) :: refined
  logical, intent(in) :: use_histograms
  
  ! Reshuffle particles in memory according to their level
  ! by applying a couting sort on the particles.

  integer  :: offset, ibin, ip, ipart
  integer  :: unrefined_pos, refined_pos
  logical  :: unrefined

  if (ilevel == nlevelmax)return

  offset = part_level_offset(ilevel)

  ! Find starting indices for refined particles
  unrefined_pos = offset; refined_pos = offset          

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)     
     do ip = offset + 1, offset + np  
        ipart = part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2)))then
           ibin=ibin+1
        end if
        if (refined(ibin) == 0) then 
           refined_pos = refined_pos + 1
        end if
     end do
  else
     refined_pos = refined_pos + np - sum(refined)
  end if

  ! set "level boundary" in particle array and rearrange particles
  part_level_offset(ilevel+1) = refined_pos

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (gt_3keys(part_hkey(ipart,0:2), bin_keys(ibin,0:2))) then
           ibin = ibin + 1
        end if
        if (refined(ibin) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  else
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (refined(ipart - offset) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  end if
  
  part_ind_permutation(offset + 1:offset + np) = &
       part_ind_permutation2(offset +1:offset + np)
  call apply_particle_permutation(offset, np, ilevel)

end subroutine reshuffle_particles
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

  call hilbert3d_reverse(x,y,z,hilbert_key2,hilbert_key1,hilbert_key0,ilevel,np)
  call get_cell_index_from_cartesian(cell_index,cell_levl,x,y,z,ilevel,np,ilevel)
  
end subroutine get_cell_index_from_hilbertkey


subroutine get_cell_index_from_cartesian(cell_index,cell_levl,xx,yy,zz,ilevel,n,bit_length)
  use amr_commons
  implicit none

  integer, intent(in)::n,ilevel,bit_length
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
  use sort,        only: gt_3keys
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! This routine computes a particle histogram for np particles in memory, 
  ! starting from offset+1 to offset+np. There must be a precomputed array 
  ! part_ind_permutation which sorts the particles by hilbert key.
  !----------------------------------------------------------------------------

  integer,                         save :: ibin, ipart, ip
  integer(kind=8), dimension(0:2), save :: current_bin_key


  ! if there is nothing to do...
  nbins = 0
  if (.not. np > 0)return
  
  ! Count the number of bins
  nbins = 1
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        nbins=nbins+1
        current_bin_key(0:2) = part_hkey(ipart,0:2)
     end if
  end do
  
  ! Allocate histograms
  if(allocated(bin_count))then
     deallocate(bin_keys)
     deallocate(bin_count)
     deallocate(bin_start_offset)
     deallocate(bin_mass)
  end if
  if (.not. allocated(bin_keys))then
     allocate(bin_keys(nbins,0:2))
     allocate(bin_count(nbins))
     allocate(bin_start_offset(nbins+1))
     allocate(bin_mass(nbins))
  end if
  bin_mass=0; bin_count=0

  ! Label every bin by a key and sum up the particles per bin, store 
  ! the offset of the first particle in each bin in the particle array

  ! First particle in first bin
  ibin=1
  bin_keys(ibin,0:2) = part_hkey(part_ind_permutation(offset+1),0:2)  
  bin_count(ibin) = 1
  bin_start_offset(ibin) = offset

  ! All other particles/bins
  current_bin_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip=offset+2, offset+np
     ipart = part_ind_permutation(ip)
     if (gt_3keys(part_hkey(ipart,0:2), current_bin_key(0:2)))then
        ibin=ibin+1
        bin_start_offset(ibin) = ip - 1 
        bin_keys(ibin,0:2) = part_hkey(ipart,0:2)
        current_bin_key(0:2) = part_hkey(ipart,0:2)
     end if
     bin_count(ibin)=bin_count(ibin)+1
  end do
  bin_start_offset(nbins+1) = offset + np

end subroutine compute_particle_histogram





subroutine count_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  ! ugly routine to count all particles in the simulation level by level
  ! used for debugging
  integer::ilevel
  integer::igrid,jgrid,i,ngrid,ncache
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts

  npts=0  
  do ilevel=1,nlevelmax
     npart2=0
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ig=0
        ip=0
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           npart2=npart2+npart1
           igrid=next(igrid)   ! Go to next grid                                                              
        
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     npart2_tot=npart2
#endif          
     npts(ilevel)=npart2_tot
  end do
  if (myid==1)print*,'total'
  if (myid==1)print*,npts(1:nlevelmax)
  
  

   npts=0  
   do ilevel=1,nlevelmax
      npart2=0        
      
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
         do i=1,ngrid
            npart2=npart2+numbp(ind_grid(i))
         end do
        
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
      npart2_tot=npart2
#endif          
      npts(ilevel)=npart2_tot
   end do
   if (myid==1)print*,'active'
   if (myid==1)print*,npts(1:nlevelmax)
   
   
end subroutine count_parts
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine check_sorted(offset, np)
  use pm_commons,   only : part_hkey, part_ind_permutation
  use sort,         only : ge_3keys
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! Simple helper routine to check whether the np particles starting form offset
  ! are sorted by hilbert key
  !----------------------------------------------------------------------------
  logical,                         save :: ok
  integer,                         save :: ipart, ip
  integer(kind=8), dimension(0:2), save :: current_key

  ok =.true.
  
  current_key(0:2) = part_hkey(part_ind_permutation(offset + 1),0:2)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (.not. ge_3keys(part_hkey(ipart,0:2), current_key(0:2)))then
        ok=.false.
        print*, part_hkey(ipart,0),current_key(0)
        print*, part_hkey(ipart,1),current_key(1)
        print*, part_hkey(ipart,2),current_key(2)
     end if
     current_key(0:2) = part_hkey(ipart,0:2)
  end do
  
  
end subroutine check_sorted

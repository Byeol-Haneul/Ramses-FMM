!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the mass density field to be used as 
  ! source term in the Poisson solver.
  ! The density field is computed for all levels greater than ilevel.
  ! On output, particles are sorted according to the level they sit in
  ! and inside their level, they are sorted in grid Hilbert order.
  !------------------------------------------------------------------
  integer::i,igrid,ind,info
  real(dp)::dx_loc,d_scale,scalar
  real(kind=8),dimension(1:ndim+1)::multipole_in,multipole_out

  if(.not. poisson)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  if(ilevel==levelmin)multipole=0d0

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  do i=nlevelmax,ilevel,-1
     ! Compute mass multipole
     if(hydro)call multipole_fine(i)
     ! Perform CIC using pseudo-particle
     call cic_from_multipole(i)
  end do

  !-------------------------------------------------------
  ! Compute particle contribution to density field
  !-------------------------------------------------------
  if(pic)then
     do i=ilevel,nlevelmax
                               call timer('rho','start')
        call cic_part(i)
                               call timer('particles','start')
        call split_part(i)
                               call timer('rho','start')
     end do
!!$     if(ilevel==levelmin)then
!!$        do i=ilevel,nlevelmax
!!$           if(tailp(i).GE.headp(i))then
!!$              write(*,'("PE ",I4," Level ",I6," npart ",I6)')&
!!$                   & myid,i,tailp(i)-headp(i)+1
!!$           endif
!!$        end do
!!$     endif
  endif

  !--------------------------------------------------------------
  ! Compute multipole contribution from all cpus and set rho_tot
  !--------------------------------------------------------------
#ifndef WITHOUTMPI
  if(ilevel==levelmin)then
     multipole_in=multipole
     call MPI_ALLREDUCE(multipole_in,multipole_out,ndim+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     multipole=multipole_out
  endif
#endif
  rho_tot=multipole(1)/boxlen**ndim
  if(debug)write(*,*)'rho_average=',rho_tot
!!! rho_tot=0d0 ! For non-periodic BC
  
111 format('   Entering rho_fine for level ',I2)
  
end subroutine rho_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multipole_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::igrid,ind,idim,ivar,nstride,ioct,icell
  integer::parent_cell,get_parent_cell
  real(dp),dimension(1:ndim),save::xx
  real(kind=8)::dx_loc,vol_loc,mmm,dd,average
  integer(kind=8),dimension(0:ndim)::hash_key
  logical::leaf_cell

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

#ifdef HYDRO
  ! Initialize multipole fields to zero
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        do idim=1,ndim+1
           grid(igrid)%unew(ind,idim)=0.0D0
        end do
     end do
  end do
#endif

  !-------------------------------------------------------
  ! Compute contribution of leaf cells to mass multipoles
  !-------------------------------------------------------
  do igrid=head(ilevel),tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim

        leaf_cell=grid(igrid)%refined(ind).EQV..FALSE.

        ! For leaf cells only
        if(leaf_cell)then

           ! Cell coordinates
           do idim=1,ndim
              nstride=2**(idim-1)
              xx(idim)=(2*grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx_loc
           end do
#ifdef HYDRO
           ! Add gas mass
           mmm=max(grid(igrid)%uold(ind,1),smallr)*vol_loc
           grid(igrid)%unew(ind,1)=grid(igrid)%unew(ind,1)+mmm
           do idim=1,ndim
              grid(igrid)%unew(ind,idim+1)=grid(igrid)%unew(ind,idim+1)+mmm*xx(idim)
           end do
#endif
           ! Add analytical density profile
           if(gravity_type < 0)then           
              call rho_ana(xx,dd,dx_loc)
              mmm=max(dd,smallr)*vol_loc
#ifdef HYDRO
              grid(igrid)%unew(ind,1)=grid(igrid)%unew(ind,1)+mmm
              do idim=1,ndim
                 grid(igrid)%unew(ind,idim+1)=grid(igrid)%unew(ind,idim+1)+mmm*xx(idim)
              end do
#endif
           end if
        endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !-------------------------------------------------------
  ! Perform octree restriction from level ilevel+1
  !-------------------------------------------------------
  if(ilevel==nlevelmax)return
  if(noct_tot(ilevel+1)==0)return

  call open_cache(operation_multipole,domain_decompos_amr)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=head(ilevel+1),tail(ilevel+1)
     hash_key(1:ndim)=grid(ioct)%ckey(1:ndim)
     ! Get parent cell using a write-only cache
     parent_cell=get_parent_cell(hash_key,grid_dict,.true.,.false.)
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim
#ifdef HYDRO
     ! Average conservative variables
     do ivar=1,ndim+1
        average=0.0d0
        do ind=1,twotondim
           average=average+grid(ioct)%unew(ind,ivar)
        end do
        ! Scatter result to cell
        grid(igrid)%unew(icell,ivar)=average
     end do
#endif
  end do

  call close_cache(grid_dict)

111 format('   Entering multipole_fine for level',i2)

end subroutine multipole_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_from_multipole(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by depositing the gas multipole mass in each cells using CIC.
  ! For pure particle runs, the gas mass deposition is not done
  ! and the routine only set rho to zero.
  !-------------------------------------------------------------------
  integer::igrid,ind

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifdef GRAV
  ! Initialize density field to zero
  do igrid=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        grid(igrid)%rho(ind)=0.0D0
     end do
  end do
#endif  

  if(hydro)call cic_cell(ilevel)

111 format('   Entering cic_from_multipole for level',i2)

end subroutine cic_from_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_cell(ilevel)
  use amr_commons
  use poisson_commons, ONLY:multipole
  implicit none
  integer::ilevel
  !
  !
  real(dp),dimension(1:ndim),save::x,dd,dg
  integer,dimension(1:ndim),save::ig,id
  real(dp),dimension(1:twotondim),save::vol
  integer,dimension(1:ndim,1:twotondim)::ckey
  integer(kind=8),dimension(0:ndim),save::hash_nbor
  integer::inbor,igrid,ind,idim
  integer::ioct,icell,parent_cell,get_parent_cell
  real(kind=8)::dx_loc,vol_loc,mmm
  
  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Use hash table directly for cells (not for grids)
  hash_nbor(0)=ilevel+1

  call open_cache(operation_rho,domain_decompos_amr)

  ! Loop over grids
  do igrid=head(ilevel),tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

#ifdef HYDRO        
        ! Compute pseudo particle mass
        mmm=grid(igrid)%unew(ind,1)

        ! Compute pseudo particle (centre of mass) position
        x(1:ndim)=grid(igrid)%unew(ind,2:ndim+1)/mmm
        
        ! Compute total multipole
        if(ilevel==levelmin)then
           do idim=1,ndim+1
              multipole(idim)=multipole(idim)+grid(igrid)%unew(ind,idim)
           end do
        endif
#endif
        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=x(idim)/dx_loc
        end do
     
        ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
        do idim=1,ndim
           dd(idim)=x(idim)+0.5D0
           id(idim)=dd(idim)
           dd(idim)=dd(idim)-id(idim)
           dg(idim)=1.0D0-dd(idim)
           ig(idim)=id(idim)-1
        end do

        ! Periodic boundary conditions
        do idim=1,ndim
           if(ig(idim)<0)ig(idim)=ckey_max(ilevel+1)-1
           if(id(idim)==ckey_max(ilevel+1))id(idim)=0
        enddo

        ! Compute cloud volumes
#if NDIM==1
        vol(1)=dg(1)
        vol(2)=dd(1)
#endif
#if NDIM==2
        vol(1)=dg(1)*dg(2)
        vol(2)=dd(1)*dg(2)
        vol(3)=dg(1)*dd(2)
        vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
        vol(1)=dg(1)*dg(2)*dg(3)
        vol(2)=dd(1)*dg(2)*dg(3)
        vol(3)=dg(1)*dd(2)*dg(3)
        vol(4)=dd(1)*dd(2)*dg(3)
        vol(5)=dg(1)*dg(2)*dd(3)
        vol(6)=dd(1)*dg(2)*dd(3)
        vol(7)=dg(1)*dd(2)*dd(3)
        vol(8)=dd(1)*dd(2)*dd(3)
#endif

        ! Compute cells Cartesian key
#if NDIM==1
        ckey(1,1)=ig(1)
        ckey(1,2)=id(1)
#endif
#if NDIM==2
        ckey(1:2,1)=(/ig(1),ig(2)/)
        ckey(1:2,2)=(/id(1),ig(2)/)
        ckey(1:2,3)=(/ig(1),id(2)/)
        ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
        ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
        ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
        ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
        ckey(1:3,4)=(/id(1),id(2),ig(3)/)
        ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
        ckey(1:3,6)=(/id(1),ig(2),id(3)/)
        ckey(1:3,7)=(/ig(1),id(2),id(3)/)
        ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif     

#ifdef GRAV
        ! Update mass density
        do inbor=1,twotondim
           hash_nbor(1:ndim)=ckey(1:ndim,inbor)
           ! Get parent cell using write-only cache
           parent_cell=get_parent_cell(hash_nbor,grid_dict,.true.,.false.)
           if(parent_cell>0)then
              ioct=(parent_cell-1)/twotondim+1
              icell=parent_cell-(ioct-1)*twotondim
              grid(ioct)%rho(icell)=grid(ioct)%rho(icell)+mmm*vol(inbor)/vol_loc
           end if
        end do
#endif     
     end do
     ! End loop over cells

  end do
  ! End loop over grids

  call close_cache(grid_dict)

end subroutine cic_cell
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_part(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons, ONLY:multipole
  use hilbert
  implicit none
  integer::ilevel
  !
  !
  real(dp),dimension(1:ndim),save::x,dd,dg
  integer,dimension(1:ndim),save::ig,id,ix
  real(dp),dimension(1:twotondim),save::vol
  integer,dimension(1:ndim,1:twotondim),save::ckey
  integer(kind=8),dimension(0:ndim),save::hash_nbor
  integer::i,ipart,inbor,igrid,ind,idim
  integer::ioct,icell,parent_cell,get_parent_cell
  real(kind=8)::dx_loc,vol_loc,vol2
  
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Compute contribution to multipole
  if(ilevel==levelmin)then
     do i=1,npart
        multipole(1)=multipole(1)+mp(i)
     end do
     do idim=1,ndim
        do i=1,npart
           multipole(idim+1)=multipole(idim+1)+mp(i)*xp(i,idim)
        end do
     end do
  endif

                               call timer('particles','start')
  ! Sort particle according to current level Hilbert key
  do i=headp(ilevel),tailp(nlevelmax)
     sortp(i)=i
  end do
  ix=0
  call sort_hilbert(headp(ilevel),tailp(nlevelmax),ix,0,1,ilevel-1)

                               call timer('rho','start')
  ! Open write-only cache for array rho
  hash_nbor(0)=ilevel+1
  call open_cache(operation_rho,domain_decompos_amr)

  ! Loop over particles in Hilbert order
  do i=headp(ilevel),tailp(nlevelmax)
     ipart=sortp(i)

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        x(idim)=xp(ipart,idim)/dx_loc
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        dd(idim)=x(idim)+0.5D0
        id(idim)=dd(idim)
        dd(idim)=dd(idim)-id(idim)
        dg(idim)=1.0D0-dd(idim)
        ig(idim)=id(idim)-1
     end do
     
     ! Periodic boundary conditions
     do idim=1,ndim
        if(ig(idim)<0)ig(idim)=ckey_max(ilevel+1)-1
        if(id(idim)==ckey_max(ilevel+1))id(idim)=0
     enddo

     ! Compute cloud volumes
#if NDIM==1
     vol(1)=dg(1)
     vol(2)=dd(1)
#endif
#if NDIM==2
     vol(1)=dg(1)*dg(2)
     vol(2)=dd(1)*dg(2)
     vol(3)=dg(1)*dd(2)
     vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
     vol(1)=dg(1)*dg(2)*dg(3)
     vol(2)=dd(1)*dg(2)*dg(3)
     vol(3)=dg(1)*dd(2)*dg(3)
     vol(4)=dd(1)*dd(2)*dg(3)
     vol(5)=dg(1)*dg(2)*dd(3)
     vol(6)=dd(1)*dg(2)*dd(3)
     vol(7)=dg(1)*dd(2)*dd(3)
     vol(8)=dd(1)*dd(2)*dd(3)
#endif

     ! Compute cells Cartesian key
#if NDIM==1
     ckey(1,1)=ig(1)
     ckey(1,2)=id(1)
#endif
#if NDIM==2
     ckey(1:2,1)=(/ig(1),ig(2)/)
     ckey(1:2,2)=(/id(1),ig(2)/)
     ckey(1:2,3)=(/ig(1),id(2)/)
     ckey(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
     ckey(1:3,1)=(/ig(1),ig(2),ig(3)/)
     ckey(1:3,2)=(/id(1),ig(2),ig(3)/)
     ckey(1:3,3)=(/ig(1),id(2),ig(3)/)
     ckey(1:3,4)=(/id(1),id(2),ig(3)/)
     ckey(1:3,5)=(/ig(1),ig(2),id(3)/)
     ckey(1:3,6)=(/id(1),ig(2),id(3)/)
     ckey(1:3,7)=(/ig(1),id(2),id(3)/)
     ckey(1:3,8)=(/id(1),id(2),id(3)/)
#endif

#ifdef GRAV
     ! Update mass density
     do ind=1,twotondim
        hash_nbor(1:ndim)=ckey(1:ndim,ind)
        ! Get parent cell using write-only cache
        parent_cell=get_parent_cell(hash_nbor,grid_dict,.true.,.false.)
        if(parent_cell>0)then
           igrid=(parent_cell-1)/twotondim+1
           icell=parent_cell-(igrid-1)*twotondim
           vol2=mp(ipart)*vol(ind)/vol_loc
           grid(igrid)%rho(icell)=grid(igrid)%rho(icell)+vol2
        endif
     end do
#endif

  end do
  ! End loop over particles
  
  call close_cache(grid_dict)

111 format('   Entering cic_part for level',i2)

end subroutine cic_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine split_part(ilevel)
  use amr_commons
  use pm_commons
  use hilbert
  implicit none
  integer::ilevel
  !
  !
  real(dp),dimension(1:ndim),save::x,xp_tmp,vp_tmp
  integer,dimension(1:ndim),save::ii,ix,ix_ref
  integer(kind=8),dimension(0:ndim),save::hash_key
  integer::i,ipart,jpart,inbor,igrid,ind,idim,ioct,icell,ipos,get_grid
  integer::npart_coarse,npart_fine
  real(kind=8)::dx_loc,vol_loc,vol2
  real(dp)::mp_tmp
  integer::levelp_tmp
  integer(i8b)::idp_tmp

  if(ilevel.GE.nlevelmax)return
  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Open read-only cache for array refined
  hash_key(0)=ilevel
  call open_cache(operation_split,domain_decompos_amr)

  ! Loop over particles
  ix_ref=-1
  npart_coarse=0
  do i=headp(ilevel),tailp(nlevelmax)
     ipart=sortp(i)

     ! Acquire grid using read-only cache
     ix = int(xp(ipart,1:ndim)/(2*dx_loc))
     if(.NOT. ALL(ix.EQ.ix_ref))then
        hash_key(1:ndim)=ix(1:ndim)
        igrid=get_grid(hash_key,grid_dict,.false.,.true.)
        ix_ref=ix
     endif

     ! If particle sits outside current level,
     ! then it is clearly not in a refined cell.
     ! This can happen during second adaptive step
     if(igrid==0)then
        npart_coarse=npart_coarse+1
        levelp(ipart)=-levelp(ipart)
     else
        ! Rescale particle position at level ilevel
        do idim=1,ndim
           x(idim)=xp(ipart,idim)/dx_loc
        end do
        
        ! Shift particle position to to 2x2x2 grid corner
        do idim=1,ndim
           ii(idim)=x(idim)-2*ix_ref(idim)
        end do
        
        ! Compute parent cell
#if NDIM==1
        icell=1+ii(1)
#endif
#if NDIM==2
        icell=1+ii(1)+2*ii(2)
#endif
#if NDIM==3
        icell=1+ii(1)+2*ii(2)+4*ii(3)
#endif
        ! Increase counter if cell is not refined
        if(.NOT.grid(igrid)%refined(icell))then
           npart_coarse=npart_coarse+1
           levelp(ipart)=-levelp(ipart)
        else
           sortp(i)=-sortp(i)
        endif
     endif

  end do
  ! End loop over particles

  call close_cache(grid_dict)
  
  tailp(ilevel)=headp(ilevel)+npart_coarse-1
  headp(ilevel+1)=tailp(ilevel)+1

  ! Loop over fine level particles
  ! This preserves the initial ordering after partioning
  npart_fine=0
  do ipart=headp(ilevel),tailp(nlevelmax)
     if(levelp(ipart)>0)then
        npart_fine=npart_fine+1
        workp(ipart)=headp(ilevel+1)+npart_fine-1
     endif
  end do

  ! Loop over coarse level particles
  ! This enforces Hilbert ordering after partioning
  npart_coarse=0
  do i=headp(ilevel),tailp(nlevelmax)
     ipart=sortp(i)
     if(ipart>0)then
        npart_coarse=npart_coarse+1
        workp(ipart)=headp(ilevel)+npart_coarse-1
        levelp(ipart)=-levelp(ipart)
     endif
  end do

  ! Swap particles using new index table
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

111 format('   Entering split_part for level',i2)

end subroutine split_part
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################

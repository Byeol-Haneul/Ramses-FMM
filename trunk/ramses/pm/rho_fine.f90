!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! ADD NEW DESCRIPTION HERE!
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
  if(ilevel==levelmin.or.icount>1)then
     do i=nlevelmax,ilevel,-1
        ! Compute mass multipole
        if(hydro)call multipole_fine(i)
        ! Perform CIC using pseudo-particle
        call cic_from_multipole(i)
     end do
  endif

  !--------------------------------------------------------------
  ! Compute multipole contribution from all cpus and set rho_tot
  !--------------------------------------------------------------
  if(ilevel==levelmin)then
#ifndef WITHOUTMPI
     multipole_in=multipole
     call MPI_ALLREDUCE(multipole_in,multipole_out,ndim+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     multipole=multipole_out
#endif
     rho_tot=multipole(1)/boxlen**ndim
     if(debug)write(*,*)'rho_average=',rho_tot
!!! rho_tot=0d0 ! For non-periodic BC
  endif
  
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

        ! For leaf cells only
        if(.NOT.grid(igrid)%refined(ind))then

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
  integer ,dimension(1:ndim),save::ig,id
  real(dp),dimension(1:twotondim),save::vol
  integer(kind=8),dimension(0:ndim,1:twotondim)::hash_key
  integer::i_nbor,igrid,ind,idim,ioct,icell,parent_cell,get_parent_cell
  real(kind=8)::dx_loc,vol_loc,mmm
  
  ! Mesh spacing in that level
  dx_loc=boxlen/2**ilevel 
  vol_loc=dx_loc**ndim

  ! Use hash table directly for cells (not for grids)
  hash_key=ilevel+1

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
           if(ig(idim)<0)ig(idim)=ckey_max(ilevel)-1
           if(id(idim)==ckey_max(ilevel))id(idim)=0
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

        ! Compute cells hash key
#if NDIM==1
        hash_key(1,1)=ig(1)
        hash_key(1,2)=id(1)
#endif
#if NDIM==2
        hash_key(1:2,1)=(/ig(1),ig(2)/)
        hash_key(1:2,2)=(/id(1),ig(2)/)
        hash_key(1:2,3)=(/ig(1),id(2)/)
        hash_key(1:2,4)=(/id(1),id(2)/)
#endif
#if NDIM==3
        hash_key(1:3,1)=(/ig(1),ig(2),ig(3)/)
        hash_key(1:3,2)=(/id(1),ig(2),ig(3)/)
        hash_key(1:3,3)=(/ig(1),id(2),ig(3)/)
        hash_key(1:3,4)=(/id(1),id(2),ig(3)/)
        hash_key(1:3,5)=(/ig(1),ig(2),id(3)/)
        hash_key(1:3,6)=(/id(1),ig(2),id(3)/)
        hash_key(1:3,7)=(/ig(1),id(2),id(3)/)
        hash_key(1:3,8)=(/id(1),id(2),id(3)/)
#endif     
        ! Compute parent cell address
#ifdef GRAV
        do i_nbor=1,twotondim
           ! Get parent cell using write-only cache
           parent_cell=get_parent_cell(hash_key(1:ndim,i_nbor),grid_dict,.true.,.false.)
           if(parent_cell>0)then
              ioct=(parent_cell-1)/twotondim+1
              icell=parent_cell-(ioct-1)*twotondim
              grid(ioct)%rho(icell)=grid(ioct)%rho(icell)+mmm*vol(i_nbor)/vol_loc
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

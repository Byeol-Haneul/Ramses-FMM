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
  integer :: ilevel
  !------------------------------------------------------------------
  ! ADD NEW DESCRIPTION HERE!
  !------------------------------------------------------------------
  integer :: particle_level
  integer::iskip,icpu,ind,i,info,nx_loc,ibound,idim,icell, ncell, ilev
  real(dp)::dx,d_scale,scale,dx_loc,scalar
  real(dp)::d0,m_refine_loc,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(kind=8)::total,total_all,total2,total2_all,tms
  real(kind=8),dimension(2)::totals_in,totals_out
  logical::multigrid=.false., ok, first
  real(kind=8),dimension(1:ndim+1)::multipole_in,multipole_out

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  if(ilevel==levelmin)multipole=0d0

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  do i=nlevelmax,ilevel,-1
     ! Compute mass multipole
     if(hydro)call multipole_fine(i)
     ! Perform CIC using pseudo-particle
     call cic_from_multipole(i)
     ! Update boundaries
     call make_virtual_reverse_dp(rho(1),i)
     call make_virtual_fine_dp   (rho(1),i)
  end do

  do ilev = ilevel, nlevelmax
     !--------------------------
     ! Initialize fields to zero
     !--------------------------
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilev)%ngrid
           phi(active(ilev)%igrid(i)+iskip)=0.0D0
        end do
     end do

     !-------------------------------------------------------------------------
     ! Initialize "number density" field to baryon number density in array phi.
     !-------------------------------------------------------------------------
     if(m_refine(ilev)>-1.0d0)then
        d_scale=max(mass_sph/dx_loc**ndim,smallr)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           if(hydro)then
              if(ivar_refine>0)then
                 do i=1,active(ilev)%ngrid
                    scalar=uold(active(ilev)%igrid(i)+iskip,ivar_refine) &
                         & /uold(active(ilev)%igrid(i)+iskip,1)
                    if(scalar>var_cut_refine)then
                       phi(active(ilev)%igrid(i)+iskip)= &
                            & rho(active(ilev)%igrid(i)+iskip)/d_scale
                    endif
                 end do
              else
                 do i=1,active(ilev)%ngrid
                    phi(active(ilev)%igrid(i)+iskip)= &
                         & rho(active(ilev)%igrid(i)+iskip)/d_scale
                 end do
              endif
           endif
        end do
     endif

     !-------------------------------------------------------
     ! Initialize rho and phi to zero in virtual boundaries
     !-------------------------------------------------------
     do icpu=1,ncpu
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,reception(icpu,ilev)%ngrid
              rho(reception(icpu,ilev)%igrid(i)+iskip)=0.0D0
              phi(reception(icpu,ilev)%igrid(i)+iskip)=0.0D0
           end do
        end do
     end do
  end do
  !---------------------------------------------------------
  ! Compute particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles

  if(pic)then
     do particle_level = ilevel, nlevelmax
        call rho_direct_particles(particle_level, ilevel)        
        call rho_histogram_particles(particle_level, ilevel)
     end do
  end if

  do particle_level = ilevel, nlevelmax
     call make_virtual_reverse_dp(rho(1),particle_level)
     call make_virtual_fine_dp   (rho(1),particle_level)
     if(m_refine(particle_level)>-1.0d0)then
        call make_virtual_reverse_dp(phi(1),particle_level)
        call make_virtual_fine_dp   (phi(1),particle_level)
     endif
  end do

  if (ilevel==levelmin) then
     call add_particle_multipole
  end if
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
  if(nboundary==0)then
     rho_tot=multipole(1)/scale**ndim
     if(debug)write(*,*)'rho_average=',rho_tot
  else
     rho_tot=0d0
  endif

  do particle_level = ilevel, nlevelmax
     !----------------------------------------------------
     ! Reset rho and phi in physical boundaries
     !----------------------------------------------------
     do ibound=1,nboundary
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,boundary(ibound,particle_level)%ngrid
              phi(boundary(ibound,particle_level)%igrid(i)+iskip)=0.0
              rho(boundary(ibound,particle_level)%igrid(i)+iskip)=0.0
           end do
        end do
     end do

     !-----------------------------------------
     ! Compute quasi Lagrangian refinement map
     !-----------------------------------------
     ! TODO: FIX QUASI LAGRANGIAN REFINEMENT STRATEGY FOR HISTOGRAM MODE!
     if(m_refine(particle_level)>-1.0d0)then
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(particle_level)%ngrid
              if(phi(active(particle_level)%igrid(i)+iskip)>=m_refine(particle_level))then
                 cpu_map2(active(particle_level)%igrid(i)+iskip)=1
              else
                 cpu_map2(active(particle_level)%igrid(i)+iskip)=0
              end if
           end do
        end do
        ! Update boundaries
        call make_virtual_fine_int(cpu_map2(1),particle_level)
     end if
  end do

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
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by affecting the gas density to leaf cells, and finally
  ! by performing a restriction operation for split cells.
  ! For pure particle runs, the restriction is not necessary and the
  ! routine only set rho to zero. On the other hand, for the Multigrid
  ! solver, the restriction is necessary in any case.
  !-------------------------------------------------------------------
  integer ::ind,i,icpu,ncache,igrid,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,nsplit,ix,iy,iz,iskip_son,ind_son,ind_grid_son,ind_cell_son
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf,ind_split
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector),save::dd
  real(kind=8)::vol,dx,dx_loc,scale,vol_loc,mm
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Initialize fields to zero
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        unew(active(ilevel)%igrid(i)+iskip,1)=0.0D0
     end do
     do idim=1,ndim
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,idim+1)=0.0D0
        end do
     end do
  end do

  ! Compute mass multipoles in each cell
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ! Gather cell indices
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells and compute cell centers
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
              do idim=1,ndim
                 xx(nleaf,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
              end do
           end if
        end do
        
        ! Compute gas multipole for leaf cells only
        if(hydro)then
           do i=1,nleaf
              mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+mm
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do
        endif

        ! Add analytical density profile for leaf cells only
        if(gravity_type < 0)then           
           ! Call user defined routine rho_ana
           call rho_ana(xx,dd,dx_loc,nleaf)
           ! Scatter results to array phi
           do i=1,nleaf
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+dd(i)*vol_loc
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=dd(i)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do           
        end if

        ! Gather split cells
        nsplit=0
        do i=1,ngrid
           if(son(ind_cell(i))>0)then
              nsplit=nsplit+1
              ind_split(nsplit)=ind_cell(i)
           end if
        end do

        ! Add children multipoles
        do ind_son=1,twotondim
           iskip_son=ncoarse+(ind_son-1)*ngridmax
           do i=1,nsplit
              ind_grid_son=son(ind_split(i))
              ind_cell_son=iskip_son+ind_grid_son
              unew(ind_split(i),1)=unew(ind_split(i),1)+unew(ind_cell_son,1)
           end do
           do idim=1,ndim
              do i=1,nsplit
                 ind_grid_son=son(ind_split(i))
                 ind_cell_son=iskip_son+ind_grid_son
                 unew(ind_split(i),idim+1)=unew(ind_split(i),idim+1)+unew(ind_cell_son,idim+1)
              end do
           end do
        end do

     end do
  enddo

  ! Update boundaries
  do idim=1,ndim+1
     call make_virtual_fine_dp(unew(1,idim),ilevel)
  end do

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
  logical::multigrid
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by affecting the gas density to leaf cells, and finally
  ! by performing a restriction operation for split cells.
  ! For pure particle runs, the restriction is not necessary and the
  ! routine only set rho to zero. On the other hand, for the Multigrid
  ! solver, the restriction is necessary in any case.
  !-------------------------------------------------------------------
  integer ::ind,i,j,icpu,ncache,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,ix,iy,iz,igrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Initialize density field to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
  end do
  ! Reset rho in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do
  
  if(hydro)then
     ! Perform a restriction over split cells (ilevel+1)
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call cic_cell(ind_grid,ngrid,ilevel)
     end do
  end if

111 format('   Entering cic_from_multipole for level',i2)

end subroutine cic_from_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_cell(ind_grid,ngrid,ilevel)
  use amr_commons
  use poisson_commons
  use hydro_commons, ONLY: unew
  implicit none
  integer::ngrid,ilevel
  integer,dimension(1:nvector)::ind_grid
  !
  !
  integer::i,j,idim,ind_cell_son,iskip_son,np,ind_son,nx_loc,ind
  integer ,dimension(1:nvector),save::ind_cell,ind_cell_father
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm,ttt
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  logical::error
  
  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  np=ngrid

  ! Compute father cell index
  do i=1,ngrid
     ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather 3x3x3 neighboring parent cells
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Loop over grid cells
  do ind_son=1,twotondim
     iskip_son=ncoarse+(ind_son-1)*ngridmax

     ! Compute pseudo particle (centre of mass) position
     do idim=1,ndim
        do j=1,np
           ind_cell_son=iskip_son+ind_grid(j)
           x(j,idim)=unew(ind_cell_son,idim+1)/unew(ind_cell_son,1)
        end do
     end do
     
     ! Compute total multipole
     if(ilevel==levelmin)then
        do idim=1,ndim+1
           do j=1,np
              ind_cell_son=iskip_son+ind_grid(j)
              multipole(idim)=multipole(idim)+unew(ind_cell_son,idim)
           end do
        end do
     endif

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/scale+skip_loc(idim)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)-(xg(ind_grid(j),idim)-3d0*dx)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/dx
        end do
     end do
     
     ! Gather particle mass
     do j=1,np
        ind_cell_son=iskip_son+ind_grid(j)
        mmm(j)=unew(ind_cell_son,1)
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        do j=1,np
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end do
     end do
     
     ! Check for illegal moves
     error=.false.
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
        end do
     end do
     if(error)then
        write(*,*)'problem in cic'
        do idim=1,ndim
           do j=1,np
              if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
                 write(*,*)x(j,1:ndim)
              endif
           end do
        end do
        stop
     end if

     ! Compute cloud volumes
#if NDIM==1
     do j=1,np
        vol(j,1)=dg(j,1)
        vol(j,2)=dd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)
        vol(j,2)=dd(j,1)*dg(j,2)
        vol(j,3)=dg(j,1)*dd(j,2)
        vol(j,4)=dd(j,1)*dd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
        vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
        vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
        vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
        vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
        vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
        vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
        vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
     end do
#endif
     
     ! Compute parent grids
     do idim=1,ndim
        do j=1,np
           igg(j,idim)=ig(j,idim)/2
           igd(j,idim)=id(j,idim)/2
        end do
     end do
#if NDIM==1
     do j=1,np
        kg(j,1)=1+igg(j,1)
        kg(j,2)=1+igd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
        kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     end do
#endif
     do ind=1,twotondim
        do j=1,np
           igrid(j,ind)=son(nbors_father_cells(j,kg(j,ind)))
        end do
     end do
     
     ! Compute parent cell position
     do idim=1,ndim
        do j=1,np
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end do
     end do
#if NDIM==1
     do j=1,np
        icell(j,1)=1+icg(j,1)
        icell(j,2)=1+icd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end do
#endif
     
     ! Compute parent cell adress
     do ind=1,twotondim
        do j=1,np
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        end do
     end do
     
     ! Update mass density and number density fields
     do ind=1,twotondim
        do j=1,np
           ok(j)=igrid(j,ind)>0
        end do
        do j=1,np
           vol2(j)=mmm(j)*vol(j,ind)/vol_loc
        end do        
        do j=1,np
           if(ok(j))then
              rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
           end if
        end do
     end do
     
  end do
  ! End loop over grid cells

end subroutine cic_cell
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
! subroutine hilbert_allparts(ilevel)
!   use pm_commons
!   use amr_commons
!   use sort, only: quick_sort_keys, qsort_parts_in_mem
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif


!   integer::ilevel
!   integer::igrid,jgrid,i,ngrid,ncache,ipart,jpart
!   integer::ip,npart1,icpu,info
!   integer,dimension(1:nvector)::ind_grid
!   integer,dimension(1:nlevelmax)::npts
! !  real(qdp),dimension(1:nvector)::order
!   real(dp),dimension(1:nvector,1:ndim)::xtest
!   integer(kind=8),dimension(1:nvector)::hkey0,hkey1,hkey2
!   integer::ntot
!   integer,allocatable,dimension(:)::order


!   ntot=0
!   ! Loop over cpus
!   do icpu=1,ncpu
!      igrid=headl(icpu,ilevel)
!      ip=0
!      ! Loop over grids
!      do jgrid=1,numbl(icpu,ilevel)
!         npart1=numbp(igrid)  ! Number of particles in the grid
!         if(npart1>0)then
!            ipart=headp(igrid)
!            ! Loop over particles
!            do jpart=1,npart1
!               ! Save next particle   <--- Very important !!!
!               ip=ip+1
!               xtest(ip,1:ndim)=xp(ipart,1:ndim)
!               if(ip==nvector)then
!                  call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!                  part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!                  part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!                  part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!                  ntot=ntot+ip
!                  ip=0
!               end if
!               ipart=nextp(ipart)  ! Go to next particle
!            end do
!         endif
!         igrid=next(igrid)   ! Go to next grid
!      end do
!      if(ip>0)then 
!         call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!         part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!         part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!         part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!         ntot=ntot+ip
!      end if
!   end do

!   print*,'ntot:',ntot,npart,myid

! !  call qsort_parts_in_mem(ntot,1)

!   allocate(order(1:ntot))
!   call quick_sort_keys(order, ntot)
  
! !  do i=1,ntot
! !     write(*,'(A8,I5,I10,3(I20))'),"myid: ",myid,i,order(i),part_hkey(i,1),part_hkey(i,0)
! !  end do
!   deallocate(order)

! end subroutine hilbert_allparts



! subroutine hilbert_allparts
!   use pm_commons
!   use amr_commons
!   use sort, only: qsort_parts_in_mem
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif


!   integer::ilevel
!   integer::igrid,jgrid,i,ngrid,ncache,jpart
!   integer::ip,npart1,icpu,info
!   integer,dimension(1:nvector)::ind_grid
!   integer,dimension(1:nlevelmax)::npts
! !  real(qdp),dimension(1:nvector)::order
!   real(dp),dimension(1:nvector,1:ndim)::xtest
!   integer(kind=8),dimension(1:nvector)::hkey0,hkey1,hkey2
!   integer::ntot
!   integer,allocatable,dimension(:)::order


!   ntot=0
!   ip=0
!   do while(ntot+ip<npart)
!      ip=ip+1
!      xtest(ip,1:ndim)=xp(ntot+ip,1:ndim)
!      if(ip==nvector)then
!         call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!         part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!         part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!         part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!         ntot=ntot+ip
!         ip=0
!      end if
!   end do
!   if(ip>0)then 
!      call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!      part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!      part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!      part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!      ntot=ntot+ip
!   end if

!   print*,'ntot:',ntot,npart,myid

! !  call qsort_parts_in_mem(ntot,1)

! !   allocate(order(1:ntot))
! !   call quick_sort_keys(order, ntot)
  
! ! !  do i=1,ntot
! ! !     write(*,'(A8,I5,I10,3(I20))'),"myid: ",myid,i,order(i),part_hkey(i,1),part_hkey(i,0)
! ! !  end do
! !   deallocate(order)

! end subroutine hilbert_allparts
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_direct_particles(part_level, min_grid_level)
  use amr_parameters, only: dp, levelmin, nvector, ndim
  use amr_commons,    only: ncpu, myid
  use pm_parameters,  only: n_dump_parts_direct
  use pm_commons,     only: part_level_offset, bin_start_offset, bin_count, &
                            xp, mp, idp, nbins, part_hkey
  use particle_communication, only: build_communicator, part_data_to_domain_dp
  implicit none
  integer, intent(in) :: part_level, min_grid_level
  
  ! This routine deposits all particles that sit at level part_level to 
  ! the grid at level min_grid_level <= grid_level <= part_level.
  
  ! in:           - particle_level
  !               - offset, nparts for local particles 
  !               - mask array which marks particles which are "histogrammed"
  !               - starting offset, number of particles
  !               - current level  
  ! out:          - none  
  ! side effect:  - updates rho field on levels ilevel <= part_level   

  real(dp), dimension(1:nvector, 1:ndim) :: xpart
  real(dp), dimension(1:nvector)         :: mpart
  integer,  dimension(1:ncpu, 1:4)       :: communicator
  integer :: ip, offset, nparts, ibin, ipart, grid_level, npart_direct
  integer :: recv_tot, local_data, local_data_oft

  integer(kind=8), allocatable, dimension(:,:) :: part_hkey_direct
  real(dp), allocatable, dimension(:,:) :: xp_direct, xp_remote
  real(dp), allocatable, dimension(:)   :: mp_direct, mp_remote
  
  offset = part_level_offset(part_level)
  nparts = part_level_offset(part_level+1) - part_level_offset(part_level)

  call compute_particle_histogram(offset, nparts)

  ! Count direct particles
  npart_direct = 0
  do ibin = 1, nbins
     if (bin_count(ibin) < n_dump_parts_direct + 0.5) then
        npart_direct = npart_direct + bin_count(ibin)
     end if
  end do

  
  allocate(part_hkey_direct(1:npart_direct, 0:2))
  allocate(xp_direct(1:npart_direct, 1:3))
  allocate(mp_direct(1:npart_direct))

  ! Fill direct particle arrays
  ip = 0
  ibin = 0
  do ipart = offset + 1, offset + nparts
     if (ipart > bin_start_offset(ibin+1)) ibin = ibin + 1
     ! a bit of a hacky comparison between a float and an integer 
     if (bin_count(ibin) < n_dump_parts_direct + 0.5) then
        ip = ip + 1
        part_hkey_direct(ip, 0:2) = part_hkey(ipart, 0:2)
        xp_direct(ip, 1:ndim)     = xp(ipart, 1:ndim)
        mp_direct(ip)             = mp(ipart)
     end if
  end do

  
  call build_communicator(communicator, recv_tot, npart_direct, local_data, local_data_oft, &
                          part_hkey_direct(1:npart_direct, 2), &
                          part_hkey_direct(1:npart_direct, 1), &
                          part_hkey_direct(1:npart_direct, 0), &
                          part_level)

  deallocate(part_hkey_direct)
  allocate(xp_remote(1:recv_tot, 1:3))
  allocate(mp_remote(1:recv_tot))

  call part_data_to_domain_dp(communicator, xp_direct(:, 1), xp_remote(:, 1))
  call part_data_to_domain_dp(communicator, xp_direct(:, 2), xp_remote(:, 2))
  call part_data_to_domain_dp(communicator, xp_direct(:, 3), xp_remote(:, 3))
  call part_data_to_domain_dp(communicator, mp_direct, mp_remote)

  ! Project local direct particles
  ip = 0
  do ipart = local_data_oft + 1, local_data_oft + local_data 
     ip = ip + 1
     xpart(ip, 1:ndim) = xp_direct(ipart, 1:ndim)
     mpart(ip)         = mp_direct(ipart)
     if (ip == nvector) then
        do grid_level = part_level, min_grid_level, -1
           call cic_amr(xpart, mpart, ip, grid_level)
        end do
        ip = 0
     end if
  end do
  if (ip > 0) then
     do grid_level = part_level, min_grid_level, -1
        call cic_amr(xpart, mpart, ip, grid_level)
     end do
  end if


  ! Project remote direct particles
  ip = 0
  do ipart = 1, recv_tot
     ip = ip + 1
     xpart(ip, 1:ndim) = xp_remote(ipart, 1:ndim)
     mpart(ip)         = mp_remote(ipart)
     if (ip == nvector) then
        do grid_level = part_level, min_grid_level, -1
           call cic_amr(xpart, mpart, ip, grid_level)
        end do
        ip = 0
     end if
  end do
  if (ip > 0) then
     do grid_level = part_level, min_grid_level, -1
        call cic_amr(xpart, mpart, ip, grid_level)
     end do
  end if

  deallocate(xp_remote, mp_remote)
  
contains

  subroutine cic_amr(xpart, mpart, np, grid_level)
    use amr_parameters,  only: static, mass_cut_refine
    use amr_commons,     only: boxlen, icoarse_max, & 
                               icoarse_min, nvector, ndim, nstep_coarse
    use poisson_commons, only: multipole, rho, phi
    use hilbert,         only: hilbert3d
    implicit none
    integer,  intent(in)                               :: np, grid_level
    real(dp), intent(in), dimension(1:nvector)         :: mpart
    real(dp), intent(in), dimension(1:nvector, 1:ndim) :: xpart
    ! This routine deposits nvector particles (local or remote) onto the grid (local)
    ! at level grid_level.

    ! in:           - particle masses
    !               - particle positions
    !               - number of particles
    !               - grid_level 
    
    ! out:          - "corrupted" particle positions -> do not reuse xpart outside of 
    !                 this subroutine  
    
    ! side effect:  - updates rho field on level grid_level

    integer(kind=8), dimension(1:nvector, 0:2),    save :: cloud_hkey
!    integer(kind=8), dimension(1:nvector, 1:ndim), save :: id
    integer(kind=8), dimension(1:nvector, 1:ndim), save :: ix
    integer(kind=4), dimension(1:nvector),         save :: dummy_state
    integer(kind=4), dimension(1:nvector),         save :: parent_cell_level, parent_cell_index
    real(dp),   dimension(1:nvector, 0:1, 1:ndim), save :: cloud_boundary
    real(dp),        dimension(1:nvector),         save :: vol, delta
    real(dp),        dimension(1:nvector, 1:3),    save :: xpart_cart
    logical,         dimension(1:nvector),         save :: ok
    integer,         dimension(1:ndim),            save :: ind
    integer,  save :: idim, nx_loc, ind_cloud, ip
    real(dp), save :: dx, dx_loc, scale, vol_loc, pos_to_cart
    integer(kind=8) :: grid_size

    grid_size = 2**grid_level    
    nx_loc=(icoarse_max-icoarse_min+1)
    scale=boxlen/dble(nx_loc)
    dx = 0.5D0**grid_level
    dx_loc=dx*scale
    vol_loc=dx_loc**ndim
    
    ! Convert particle coordinates in code units
    ! into "cartesian" coordinates at grid_level
    pos_to_cart = 2.0_dp**grid_level / dble(boxlen)
    xpart_cart = xpart * pos_to_cart

    ! compute distances of cloud boundary from nearest "integer coordinate"
    do idim=1,ndim       

       ! upper/right/front boundary of the cloud
       do ip=1,np
          cloud_boundary(ip,1,idim) = xpart_cart(ip, idim) + 0.5D0
       end do
     
       ! upper/rigt/front boundary rel to nearest integer (type conversion here...)
       do ip=1,np
          cloud_boundary(ip,1,idim) = cloud_boundary(ip,1,idim) - floor(cloud_boundary(ip,1,idim), kind=8)
       end do
       
       ! lower/left/back boundary rel to nearest integer
       do ip=1,np
          cloud_boundary(ip,0,idim) = 1._dp - cloud_boundary(ip,1,idim)
       end do
    end do
    
#if NDIM<3
    write(*,*)'add non-3D version of this routine'
    stop
#endif
#if NDIM==1
    ! Loop cloud/cell intersections
    do ind_cloud = 0, 1
       ind(1) = ind_cloud 
       
       ! Compute cloud volume
       do ip=1,np
          vol(ip) = cloud_boundary(ip,ind(1),1) * &
               cloud_boundary(ip,ind(2),2) 
       end do
#endif
#if NDIM==2
    ! Loop cloud/cell intersections
    do ind_cloud = 0, 3
       ind(1) = ind_cloud/2
       ind(2) = mod(ind_cloud,2)
       
       ! Compute cloud volume
       do ip=1,np
          vol(ip) = cloud_boundary(ip,ind(1),1) * &
               cloud_boundary(ip,ind(2),2) 
       end do
#endif
#if NDIM==3
    ! Loop cloud/cell intersections
    do ind_cloud = 0, 7
       ind(1) = ind_cloud/4
       ind(2) = mod(ind_cloud,4)/2
       ind(3) = mod(mod(ind_cloud,4),2)
       
       ! Compute cloud volume
       do ip=1,np
          vol(ip) = cloud_boundary(ip,ind(1),1) * &
               cloud_boundary(ip,ind(2),2) * &
               cloud_boundary(ip,ind(3),3) 
       end do
#endif

       ! Compute cloud corner offset from cloud center
       delta(1:ndim) = ind(1:ndim) - 0.5D0       

       ! Get cell indices which are covered by cloud
       ! (cartesian key -> hilbert key -> cell index)
       ! TODO: Add support for non-periodic boundaries
       ! TODO: Check boundary behaviour - currently particles sitting in cell touching the boundary cause
       ! deviations from the old code.
       do idim = 1, ndim
          do ip = 1, np
             ix(ip,idim) = modulo(floor(xpart_cart(ip,idim) + delta(idim), kind = 8), grid_size)
          end do
       end do

       call hilbert3d(ix(1:np,1), ix(1:np,2), ix(1:np,3), &
            cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), &
            dummy_state, 0, grid_level, np)
       
       call get_cell_index_from_hilbertkey(parent_cell_index(1:np), parent_cell_level(1:np), &
            cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), np, grid_level)

       ! Exclude cloud fraction which lies in coarser level
       do ip = 1, np        
          ok(ip) = (parent_cell_level(ip) == grid_level)
       end do

       ! Add to number density which is stored in phi
       do ip=1,np
          if(ok(ip))then
             phi(parent_cell_index(ip)) = phi(parent_cell_index(ip)) + vol(ip)
          end if
       end do

       ! compute delta rho and add to rho
       do ip = 1, np
          vol(ip) = mpart(ip) * vol(ip) / vol_loc
       end do

       do ip = 1, np
          if (ok(ip)) then
             rho(parent_cell_index(ip)) = rho(parent_cell_index(ip)) + vol(ip)
          end if
       end do

    end do ! end loop over cloud/cell intersections
  end subroutine cic_amr
end subroutine rho_direct_particles

!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_histogram_particles(part_level, min_grid_level)
  use amr_parameters, only: dp, levelmin, icoarse_min, icoarse_max, boxlen
  use pm_parameters,  only: n_dump_parts_direct
  use pm_commons,     only: part_ind_permutation, part_level_offset, bin_start_offset, &
       bin_count, mp
  use amr_commons,    only: myid
  implicit none
  integer, intent(in) :: part_level, min_grid_level
  ! """
  ! This routine deposits the mass of all those part_level particles which have 
  ! not been deposited directly.
  ! General strategy: 'Non-direct' particles are masked and a permutation is
  ! created which allows to access the masked particles in a row.
  ! For the masked particles do cic at all grid_level <= particle_level

  ! in:               - particle_level
  ! 'implicit' input: - part_ind_permutation 
  ! out:              - none                                                       
  ! side effect:      - updates rho field on levels ilevel <= part_level   

  ! highest level subroutine contains:
  ! - subroutine rho_particle_histogram_onelevel(offset, nparts, n_masked, grid_level)
  ! - subroutine cic_histogram(xpart, mpart, bin_nr, np, grid_level, ind_cloud) 
  ! - subroutine dump_histograms(cell_level)
  ! """
  
  
  integer  :: ip, offset, nparts, ibin, ipart, grid_level, n_masked
  offset = part_level_offset(part_level)
  nparts = part_level_offset(part_level+1) - part_level_offset(part_level)

  call compute_particle_histogram(offset, nparts)

  ! Count number of "masked" particles and compute permuation such
  ! that masked particles are accessed first inside the level.
  n_masked = 0
  ibin = 0
  do ipart = offset + 1, offset + nparts
     if (ipart > bin_start_offset(ibin+1)) ibin = ibin + 1
     if (bin_count(ibin) > n_dump_parts_direct + 0.5) then
        n_masked = n_masked + 1
        part_ind_permutation(offset + n_masked) = ipart
     end if
  end do
  
  ! outer loop here over grid levels
!  if (n_masked > 0)then
     do grid_level = min_grid_level, part_level
        call rho_particle_histogram_onelevel(offset, nparts, n_masked, grid_level)
     end do
!  end if

  
contains
  
  subroutine rho_particle_histogram_onelevel(offset, nparts, n_masked, grid_level)
    use amr_parameters,  only: ndim, nvector
    use amr_commons,     only: myid
    use pm_commons,      only: xp, part_ind_permutation, part_ind_permutation2, nbins, bin_start_offset, &
                              bin_count, bin_mass, idp
    use hilbert,         only: hilbert_for_particle
    use sort,            only: lsd_radix_sort_particles
    use poisson_commons, only: rho
    implicit none
    integer, intent(in) :: grid_level, nparts, n_masked, offset

    ! This routine performs the cic for the masked particles at level grid_level.
    ! This is achieved by looping over the 8 "cic-particles", each time
    ! resorting according to the grid-level hilbert key of the cic-particle positions.
    ! The mass per cell/bin is added to rho for each of the 8 histograms individually.
    ! not been deposited directly.
    
    ! in:               - offset, nparts, n_masked, grid_level    
    ! "implicit" input: - part_ind_permutation, 
    ! out:              - none                                                       
    ! side effect:      - updates rho field on levels grid_level
    
    
    integer(kind=4), dimension(1:ndim)   ,          save :: ind
    integer(kind=4), dimension(1:nvector),          save :: bin_nr
    real(dp),        dimension(1:nvector),          save :: delta
    real(dp),        dimension(1:nvector, 1:ndim),  save :: xpart
    real(dp),        dimension(1:nvector),          save :: mpart
    real(dp), save :: dx, dx_loc, scale, vol_loc
    integer,  save :: nx_loc, idim, ind_cloud, ip_sweep

    nx_loc = (icoarse_max - icoarse_min+1)
    scale = boxlen / dble(nx_loc)
    dx = 0.5D0**grid_level
    dx_loc = dx * scale
    vol_loc = dx_loc**ndim

  ! Loop over CIC cloud / cell intersections
    do ind_cloud = 0, 7
       ind(1) = ind_cloud/4
       ind(2) = mod(ind_cloud,4)/2
       ind(3) = mod(mod(ind_cloud,4),2)
       
       ! Compute cloud corner offset from cloud center
       delta(1:ndim) = ind(1:ndim) - 0.5D0       

       ! Relocate particle positions according to index
       do idim = 1, ndim
          do ip = offset + 1, offset + n_masked
             ipart = part_ind_permutation(ip)
             xp(ipart,idim) = xp(ipart,idim) + delta(idim) * dx_loc
             if (xp(ipart, idim) > boxlen) then
                xp(ipart, idim) = xp(ipart, idim) - boxlen
             end if
             if (xp(ipart, idim) < 0.0_dp) then
                xp(ipart, idim) = xp(ipart, idim) + boxlen
             end if
          end do
       end do

       ! Recompute hilbert key for all parts, also the 
       ! "unmasked ones" - > if too slow, try do it only for the masked ones
       call hilbert_for_particle(offset, nparts, 0, grid_level)

       ! Sort hilbert keys
       call lsd_radix_sort_particles(offset, n_masked, grid_level, grid_level, .false.)

       ! Compute "reduced" histogram
       call compute_particle_histogram(offset, n_masked)
       ! Reset bin count as it is computed using cic later
       bin_count = 0.d0; bin_mass = 0.d0;
       
!       print*, 'going in loop', offset, n_masked, nbins, grid_level, part_level, ind_cloud 
       
       ! Loop masked, sorted parts in sweeps and dump the mass for the
       ! given cloud/cell intersection
       ip_sweep = 0
       ibin = 0
       do ip = offset+1,offset+n_masked
          if (ip > bin_start_offset(ibin+1)) ibin = ibin + 1
          ipart = part_ind_permutation(ip)
          ip_sweep = ip_sweep + 1
          xpart(ip_sweep,1:ndim) = xp(ipart,1:ndim)
          mpart(ip_sweep)        = mp(ipart)
          bin_nr(ip_sweep)       = ibin
          
          if (ip_sweep == nvector) then
             call cic_histogram(xpart, mpart, bin_nr, ip_sweep,  grid_level, ind_cloud)
             ip_sweep = 0
          end if
       end do
       if (ip_sweep > 0) then
          call cic_histogram(xpart, mpart, bin_nr, ip_sweep, grid_level, ind_cloud)
       end if

       ! Reset particles to original positions
       do idim = 1, ndim
          do ip = offset + 1, offset + n_masked
             ipart = part_ind_permutation(ip)
             xp(ipart,idim) = xp(ipart,idim) - delta(idim) * dx_loc
             if (xp(ipart, idim) > boxlen) then
                xp(ipart, idim) = xp(ipart, idim) - boxlen
             end if
             if (xp(ipart, idim) < 0.0_dp) then
                xp(ipart, idim) = xp(ipart, idim) + boxlen
             end if
          end do
       end do

       ! Dump bin_mass into rho and bin_count into phi
       call dump_histograms(grid_level)

    end do ! end loop over 8 cic-particles
    ! fix hilber keys for particles!
    call hilbert_for_particle(offset, nparts, 0, grid_level)
    
  end subroutine rho_particle_histogram_onelevel
  
  subroutine cic_histogram(xpart, mpart, bin_nr, np, grid_level, ind_cloud)
    use amr_commons,     only: boxlen, icoarse_max, nvector, ndim
    use pm_commons,      only: bin_keys, bin_mass, nbins, bin_count
    implicit none
    integer,  intent(in)                                  :: np, grid_level, ind_cloud
    integer,  intent(in),    dimension(1:nvector)         :: bin_nr
    real(dp), intent(in),    dimension(1:nvector)         :: mpart
    real(dp), intent(inout), dimension(1:nvector, 1:ndim) :: xpart


    real(dp), dimension(1:nvector),save :: vol, vol_idim
    integer,  dimension(1:ndim),   save :: ind
    integer,  save :: idim, ip
    real(dp), save :: pos_to_cart

#if NDIM<3
    write(*,*)'add non-3D version of this routine'
    stop
#endif

    ! Convert particle coordinates in code units
    ! into "cartesian" coordinates at grid_level
    pos_to_cart = 2.0**grid_level / dble(boxlen)
    xpart = xpart * pos_to_cart
    
    ind(1) = ind_cloud/4
    ind(2) = mod(ind_cloud,4)/2
    ind(3) = mod(mod(ind_cloud,4),2)

    ! compute volume of cloud/cell intersection
    vol(1:np) = 1.d0
    do idim=1,ndim       
       vol_idim(1:np) = xpart(1:np, idim) - floor(xpart(1:np, idim))
       if (ind(idim)==0) vol_idim(1:np) = 1.d0 - vol_idim(1:np)
       vol(1:np) = vol(1:np) * vol_idim(1:np)          
    end do
        
    ! Compute particles per bin
    do ip=1,np
       bin_count(bin_nr(ip)) = bin_count(bin_nr(ip)) + vol(ip)
    end do
    
    ! Compute mass per bin
    do ip = 1, np
       bin_mass(bin_nr(ip)) = bin_mass(bin_nr(ip)) + vol(ip) * mpart(ip)
    end do

  end subroutine cic_histogram

  subroutine dump_histograms(cell_level)
    use amr_parameters,  only: nvector, dp
    use amr_commons,     only: ncpu, myid
    use pm_commons,      only: bin_keys, bin_mass, nbins
    use poisson_commons, only: rho, phi
    use particle_communication, only: build_communicator, part_data_to_domain_dp, &
         part_data_to_domain_i8
    implicit none
    integer, intent(in) :: cell_level


    integer,  dimension(1:ncpu, 1:4)             :: communicator
    integer(kind=8), allocatable, dimension(:,:) :: bin_keys_remote
    real(dp), allocatable, dimension(:)          :: bin_mass_remote
    real(dp), allocatable, dimension(:)          :: bin_count_remote

    integer        , dimension(1:nvector)     , save :: parent_cell_level, parent_cell_index
    integer(kind=8), dimension(1:nvector, 0:2), save :: bkey

    integer,  save :: ib, nb, ibin, recv_tot, local_bins, local_bins_oft, ioft
    real(dp), save :: vol_loc

    vol_loc = (0.5**cell_level * dble(boxlen) )**3    


      call build_communicator(communicator, recv_tot, nbins, local_bins, local_bins_oft, &
                          bin_keys(:, 2), bin_keys(:, 1), bin_keys(:, 0), cell_level)

      allocate(bin_keys_remote(1:recv_tot, 0:2))
      allocate(bin_mass_remote(1:recv_tot))
      allocate(bin_count_remote(1:recv_tot))

      bin_mass_remote = 0.d0
      bin_keys_remote = 0
      bin_count_remote = 0.d0

      call part_data_to_domain_i8(communicator, bin_keys(:, 0), bin_keys_remote(:, 0))
      call part_data_to_domain_i8(communicator, bin_keys(:, 1), bin_keys_remote(:, 1))
      call part_data_to_domain_i8(communicator, bin_keys(:, 2), bin_keys_remote(:, 2))
      call part_data_to_domain_dp(communicator, bin_mass, bin_mass_remote)
      call part_data_to_domain_dp(communicator, bin_count, bin_count_remote)
      
      ! go through bins in sweeps and add mass to corresponding cell
      do ioft = local_bins_oft, local_bins_oft + local_bins - 1, nvector
         nb = min(nvector, local_bins_oft + local_bins - ioft)
         call get_cell_index_from_hilbertkey(parent_cell_index(1:nb), &
              parent_cell_level(1:nb), &
              bin_keys(ioft + 1 : ioft + nb, 2), &
              bin_keys(ioft + 1 : ioft + nb, 1), &
              bin_keys(ioft + 1 : ioft + nb, 0), nb, cell_level)    
         do ib = 1, nb
            ! Don't add mass to coarser levels             
            if (parent_cell_level(ib) == cell_level) then
               rho(parent_cell_index(ib)) = rho(parent_cell_index(ib)) + &
                    bin_mass(ioft + ib) / vol_loc             
               phi(parent_cell_index(ib)) = phi(parent_cell_index(ib)) + &
                    bin_count(ioft + ib)
            end if
         end do
      end do

      ! go through remote bins in sweeps and add mass to corresponding cell
      do ioft = 0, recv_tot - 1, nvector
         nb = min(nvector, recv_tot - ioft)
         call get_cell_index_from_hilbertkey(parent_cell_index(1:nb), &
              parent_cell_level(1:nb), &
              bin_keys_remote(ioft + 1 : ioft + nb, 2), &
              bin_keys_remote(ioft + 1 : ioft + nb, 1), &
              bin_keys_remote(ioft + 1 : ioft + nb, 0), nb, cell_level)    
         do ib = 1, nb
            ! Don't add mass to coarser levels             
            if (parent_cell_level(ib) == cell_level) then
               rho(parent_cell_index(ib)) = rho(parent_cell_index(ib)) + &
                    bin_mass_remote(ioft + ib) / vol_loc
               phi(parent_cell_index(ib)) = phi(parent_cell_index(ib)) + & 
                    bin_mass_remote(ioft + ib) 
            end if
         end do
      end do
      
      
      deallocate(bin_mass_remote, bin_keys_remote, bin_count_remote)
      
    end subroutine dump_histograms

    
end subroutine rho_histogram_particles

subroutine add_particle_multipole
  use amr_parameters, only: ndim
  use amr_commons, only: myid
  use pm_commons, only: xp, mp, npart
  use poisson_commons, only: multipole
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Simple routine to compute the multipole contribution from particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: ipart, idim

  do ipart = 1, npart
     multipole(1) = multipole(1) + mp(ipart)
  end do
  do idim = 1, ndim
     do ipart = 1, npart
        multipole(idim + 1) = multipole(idim + 1) + mp(ipart) * xp(ipart, idim)
     end do
  end do
       
end subroutine add_particle_multipole

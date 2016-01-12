!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(noct(ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  do igrid=head(ilevel),tail(ilevel)
     ind_grid(1)=igrid
     call godfine1(ind_grid,1,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip
  real(dp)::d,u,v,w,e

  if(noct(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do i=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        do ivar=1,nvar
           grid(i)%unew(ind,ivar) = grid(i)%uold(ind,ivar)
        end do
#ifdef DUALENER
        grid(i)%divu(ind) = 0.0
        d=max(grid(i)%uold(ind,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=grid(i)%uold(ind,2)/d
        if(ndim>1)v=grid(i)%uold(ind,3)/d
        if(ndim>2)w=grid(i)%uold(ind,4)/d
        e=grid(i)%uold(ind,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
        do irad=1,nener
           e=e-grid(i)%uold(ind,ndim+2+irad)
        end do
#endif          
        grid(i)%enew(ind) = e
#endif
     end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use hash
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew 
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,irad,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx,fact,d_old

  if(noct(ilevel)==0)return
  if(verbose)write(*,111)ilevel

  dx=0.5d0**ilevel*boxlen

  ! Set uold to unew for myid cells
  do i=head(ilevel),tail(ilevel)
     do ind=1,twotondim
        do ivar=1,nvar
           grid(i)%uold(ind,ivar) = grid(i)%unew(ind,ivar)
        end do
#ifdef DUALENER
        ! Correct total energy if internal energy is too small
        d=max(grid(i)%uold(ind,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=grid(i)%uold(ind,2)/d
        if(ndim>1)v=grid(i)%uold(ind,3)/d
        if(ndim>2)w=grid(i)%uold(ind,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
        do irad=1,nener
           e_kin=e_kin+grid(i)%uold(ind,ndim+2+irad)
        end do
#endif
        e_cons=grid(i)%uold(ind,ndim+2)-e_kin
        e_prim=grid(i)%enew(ind)
        ! Note: here divu=-div.u*dt
        div=abs(grid(i)%divu(ind))*dx/dtnew(ilevel)
        !           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
        if(e_cons<e_trunc)then
           grid(i)%uold(ind,ndim+2)=e_prim+e_kin
        end if
#endif
     end do
  end do

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use hash
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer,ipos
  integer::igrid,icell,ichild,parent_cell,get_parent_cell,indf,parent_cell2
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim
  integer,dimension(1:ndim)::ii
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel*boxlen

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  hash_key(0)=grid(ind_grid(1))%lev
  hash_key(1:ndim)=grid(ind_grid(1))%ckey(1:ndim)
  hash_nbor(0)=hash_key(0)

  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max     

     ! Compute neighboring grid Cartesian index
#if NDIM>0
     hash_nbor(1)=hash_key(1)+i1-1.0
#endif
#if NDIM>1
     hash_nbor(2)=hash_key(2)+j1-1.0
#endif
#if NDIM>2
     hash_nbor(3)=hash_key(3)+k1-1.0
#endif
     ! Periodic boundary conditons
     do idim=1,ndim
        if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
        if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
     enddo

     ! Get neighboring grid index
     ichild=hash_get(grid_dict,hash_nbor)

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! If neighboring grid exists
        if(ichild>0)then
           ! Gather hydro variables
           do ivar=1,nvar
              uloc(i3,j3,k3,ivar)=grid(ichild)%uold(ind_son,ivar)
           end do
           ! Gather refinement flag
           ok(i3,j3,k3)=grid(ichild)%refined(ind_son)
        else
           ! Get parent father cell
           parent_cell=get_parent_cell(hash_nbor)           
           if(parent_cell>0)then
              ind_father=parent_cell
           else
              write(*,*)'parent_cell should exist 1'
              stop
           endif
           igrid=(parent_cell-1)/twotondim+1
           icell=parent_cell-(igrid-1)*twotondim
           do ivar=1,nvar
              uloc(i3,j3,k3,ivar)=grid(igrid)%uold(icell,ivar)
           end do
           ok(i3,j3,k3)=.false.
        end if

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),1)

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           if(ok(i3-i0,j3-j0,k3-k0) .or. ok(i3,j3,k3))then
              flux(i3,j3,k3,ivar,idim)=0.0d0
           end if
        end do
#ifdef DUALENER
        do ivar=1,2
           if(ok(i3-i0,j3-j0,k3-k0) .or. ok(i3,j3,k3))then
              tmp(i3,j3,k3,ivar,idim)=0.0d0
           end if
        end do
#endif
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           grid(ind_grid(1))%unew(ind_son,ivar)=&
                & grid(ind_grid(1))%unew(ind_son,ivar)+ &
                & (flux(i3   ,j3   ,k3   ,ivar,idim) &
                & -flux(i3+i0,j3+j0,k3+k0,ivar,idim))
        end do
#ifdef DUALENER
        ! Update velocity divergence
        grid(ind_grid(1))%divu(ind_son)=&
             & grid(ind_grid(1))%divu(ind_son)+ &
             & (tmp(i3   ,j3   ,k3   ,1,idim) &
             & -tmp(i3+i0,j3+j0,k3+k0,1,idim))
        ! Update internal energy
        grid(ind_grid(1))%enew(ind_son)=&
             & grid(ind_grid(1))%enew(ind_son)+ &
             & (tmp(i3   ,j3   ,k3   ,2,idim) &
             & -tmp(i3+i0,j3+j0,k3+k0,2,idim))
#endif
     end do
     end do
     end do
  end do

  ! If sitting at coarser level, exit. 
  if(ilevel==levelmin)return

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
#if NDIM>0
     hash_nbor(1)=hash_key(1)-i0
#endif
#if NDIM>1
     hash_nbor(2)=hash_key(2)-j0
#endif
#if NDIM>2
     hash_nbor(3)=hash_key(3)-k0
#endif
     ! Periodic boundary conditons
     if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
     if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
     parent_cell=get_parent_cell(hash_nbor)
     if(parent_cell>0)then
        ind_father=parent_cell
     else
        write(*,*)'parent_cell should exist 2'
        write(*,*)parent_cell
        write(*,*)hash_nbor
        write(*,*)hash_key
        stop
     endif
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim

     if(.NOT.grid(igrid)%refined(icell))then
        ! Conservative update of new state variables
        do ivar=1,nvar
           ! Loop over boundary cells
           do k3=k3min,k3max-k0
              do j3=j3min,j3max-j0
                 do i3=i3min,i3max-i0
                    grid(igrid)%unew(icell,ivar)=grid(igrid)%unew(icell,ivar) &
                         & -flux(i3,j3,k3,ivar,idim)*oneontwotondim
                 end do
              end do
           end do
        end do
#ifdef DUELENER
        ! Update velocity divergence
        do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
              do i3=i3min,i3max-i0
                 grid(igrid)%divu(icell)=grid(igrid)%divu(icell) &
                      & -tmp(i3,j3,k3,1,idim)*oneontwotondim
              end do
           end do
        end do
        ! Update internal energy
        do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
              do i3=i3min,i3max-i0
                 grid(igrid)%enew(icell)=grid(igrid)%enew(icell) &
                      & -tmp(i3,j3,k3,2,idim)*oneontwotondim
              end do
           end do
        end do
#endif
     endif

     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
#if NDIM>0
     hash_nbor(1)=hash_key(1)+i0
#endif
#if NDIM>1
     hash_nbor(2)=hash_key(2)+j0
#endif
#if NDIM>2
     hash_nbor(3)=hash_key(3)+k0
#endif
     ! Periodic boundary conditons
     if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
     if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
     parent_cell=get_parent_cell(hash_nbor)
     if(parent_cell>0)then
        ind_father=parent_cell
     else
        write(*,*)'parent_cell should exist 3'
        stop
     endif
     igrid=(parent_cell-1)/twotondim+1
     icell=parent_cell-(igrid-1)*twotondim

     if(.NOT.grid(igrid)%refined(icell))then
        ! Conservative update of new state variables
        do ivar=1,nvar
           ! Loop over boundary cells
           do k3=k3min+k0,k3max
              do j3=j3min+j0,j3max
                 do i3=i3min+i0,i3max
                    grid(igrid)%unew(icell,ivar)=grid(igrid)%unew(icell,ivar) &
                         & +flux(i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                 end do
              end do
           end do
        end do
#ifdef DUALENER
        ! Update velocity divergence
        do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
              do i3=i3min+i0,i3max
                 grid(igrid)%divu(icell)=grid(igrid)%divu(icell) &
                      & +tmp(i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
              end do
           end do
        end do
        ! Update internal energy
        do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
              do i3=i3min+i0,i3max
                 grid(igrid)%enew(icell)=grid(igrid)%enew(incell) &
                      & +tmp(i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
              end do
           end do
        end do
#endif
     end if

  end do
  ! End loop over dimensions

end subroutine godfine1

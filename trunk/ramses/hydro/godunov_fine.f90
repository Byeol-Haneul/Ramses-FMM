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
  integer::i,ivar,igrid

  if(noct_tot(ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  cache_operation=operation_godunov
  call open_cache

  ! Loop over active grids by vector sweeps
  igrid=head(ilevel)
  do while(igrid.LE.tail(ilevel))
     SELECT CASE (grid(igrid)%superoct)
     CASE(1)
        call godfine1(igrid,ilevel,&
             & uloc,gloc,qloc,cloc,&
             & okloc,childloc,parentloc,&
             & flux,tmp,dq,qm,qp,fx,tx,divu,&
             & iu1,iu2,ju1,ju2,ku1,ku2,&
             & io1,io2,jo1,jo2,ko1,ko2,&
             & if1,if2,jf1,jf2,kf1,kf2)
     CASE(2**ndim)
        call godfine1(igrid,ilevel,&
             & uloc_2,gloc_2,qloc_2,cloc_2,&
             & okloc_2,childloc_2,parentloc_2,&
             & flux_2,tmp_2,dq_2,qm_2,qp_2,fx_2,tx_2,divu_2,&
             & iu1_2,iu2_2,ju1_2,ju2_2,ku1_2,ku2_2,&
             & io1_2,io2_2,jo1_2,jo2_2,ko1_2,ko2_2,&
             & if1_2,if2_2,jf1_2,jf2_2,kf1_2,kf2_2)
     CASE(4**ndim)
        call godfine1(igrid,ilevel,&
             & uloc_4,gloc_4,qloc_4,cloc_4,&
             & okloc_4,childloc_4,parentloc_4,&
             & flux_4,tmp_4,dq_4,qm_4,qp_4,fx_4,tx_4,divu_4,&
             & iu1_4,iu2_4,ju1_4,ju2_4,ku1_4,ku2_4,&
             & io1_4,io2_4,jo1_4,jo2_4,ko1_4,ko2_4,&
             & if1_4,if2_4,jf1_4,jf2_4,kf1_4,kf2_4)
     CASE(8**ndim)
        call godfine1(igrid,ilevel,&
             & uloc_8,gloc_8,qloc_8,cloc_8,&
             & okloc_8,childloc_8,parentloc_8,&
             & flux_8,tmp_8,dq_8,qm_8,qp_8,fx_8,tx_8,divu_8,&
             & iu1_8,iu2_8,ju1_8,ju2_8,ku1_8,ku2_8,&
             & io1_8,io2_8,jo1_8,jo2_8,ko1_8,ko2_8,&
             & if1_8,if2_8,jf1_8,jf2_8,kf1_8,kf2_8)
     CASE(16**ndim)
        call godfine1(igrid,ilevel,&
             & uloc_16,gloc_16,qloc_16,cloc_16,&
             & okloc_16,childloc_16,parentloc_16,&
             & flux_16,tmp_16,dq_16,qm_16,qp_16,fx_16,tx_16,divu_16,&
             & iu1_16,iu2_16,ju1_16,ju2_16,ku1_16,ku2_16,&
             & io1_16,io2_16,jo1_16,jo2_16,ko1_16,ko2_16,&
             & if1_16,if2_16,jf1_16,jf2_16,kf1_16,kf2_16)
     CASE(32**ndim)
        call godfine1(igrid,ilevel,&
             & uloc_32,gloc_32,qloc_32,cloc_32,&
             & okloc_32,childloc_32,parentloc_32,&
             & flux_32,tmp_32,dq_32,qm_32,qp_32,fx_32,tx_32,divu_32,&
             & iu1_32,iu2_32,ju1_32,ju2_32,ku1_32,ku2_32,&
             & io1_32,io2_32,jo1_32,jo2_32,ko1_32,ko2_32,&
             & if1_32,if2_32,jf1_32,jf2_32,kf1_32,kf2_32)
     END SELECT
     igrid=igrid+grid(igrid)%superoct
  end do

  call close_cache

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

  if(noct_tot(ilevel)==0)return
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

  if(noct_tot(ilevel)==0)return
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
subroutine godfine1(ind_grid,ilevel,&
     & uloc,gloc,qloc,cloc,&
     & okloc,childloc,parentloc,&
     & flux,tmp,dq,qm,qp,fx,tx,divu,&
     & iu1,iu2,ju1,ju2,ku1,ku2,&
     & io1,io2,jo1,jo2,ko1,ko2,&
     & if1,if2,jf1,jf2,kf1,kf2)
  use amr_commons
  use hash
  implicit none
  integer::ind_grid,ilevel
  integer::iu1,iu2,ju1,ju2,ku1,ku2
  integer::io1,io2,jo1,jo2,ko1,ko2
  integer::if1,if2,jf1,jf2,kf1,kf2
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uloc
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gloc
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::qloc
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2)::cloc
  logical ,dimension(iu1:iu2,ju1:ju2,ku1:ku2)::okloc
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp
  real(dp),dimension(if1:if2,jf1:jf2,kf1:kf2)::divu
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::fx
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2,1:2   )::tx
  integer ,dimension(io1:io2,jo1:jo2,ko1:ko2)::childloc
  integer ,dimension(io1:io2,jo1:jo2,ko1:ko2)::parentloc
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
  integer::get_grid,get_parent_cell
  integer::i,j,ivar,idim,ind_son,iskip,nbuffer,ibuffer,ipos,ind_oct
  integer::igrid,icell,ichild,parent_cell,indf,parent_cell2
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::ii0,jj0,kk0,ii1,jj1,kk1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::ii1min,ii1max,jj1min,jj1max,kk1min,kk1max
  integer::i2min=0,i2max=0,j2min=0,j2max=0,k2min=0,k2max=0
  integer::i3min=1,i3max=1,j3min=1,j3max=1,k3min=1,k3max=1
  real(dp)::dx,scale,oneontwotondim
  integer,dimension(1:ndim)::ii
  integer,dimension(1:ndim)::ckey_corner,ckey
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  logical::okx=.true.,oky=.true.,okz=.true.

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel*boxlen

  ! Integer constants
  i1min=io1; i1max=io2; j1min=jo1; j1max=jo2; k1min=ko1; k1max=ko2
#if NDIM>0
  i2max=1; i3min=iu1+2; i3max=iu2-2
#endif
#if NDIM>1
  j2max=1; j3min=ju1+2; j3max=ju2-2
#endif
#if NDIM>2
  k2max=1; k3min=ku1+2; k3max=ku2-2
#endif

  gloc=0.0

  !---------------------
  ! Gather hydro stencil
  !---------------------
  hash_key(0)=grid(ind_grid)%lev
  hash_nbor(0)=grid(ind_grid)%lev
  ckey_corner(1:ndim)=(grid(ind_grid)%ckey(1:ndim)/(i1max-1))*(i1max-1)
  ind_oct=ind_grid

  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
#if NDIM>2
     okz=(k1>k1min.and.k1<k1max)
#endif
     do j1=j1min,j1max
#if NDIM>1
        oky=(j1>j1min.and.j1<j1max)
#endif
        do i1=i1min,i1max     
#if NDIM>0
           okx=(i1>i1min.and.i1<i1max)
#endif
           ! For inner octs only
           if(okx.and.oky.and.okz)then

              ! Compute relative Cartesian key
              ckey(1:ndim)=grid(ind_oct)%ckey(1:ndim)-ckey_corner(1:ndim)

              ! Store grid index
              ii1=1; jj1=1; kk1=1
#if NDIM>0
              ii1=ckey(1)+1
#endif
#if NDIM>1
              jj1=ckey(2)+1
#endif
#if NDIM>2
              kk1=ckey(3)+1
#endif
              childloc(ii1,jj1,kk1)=ind_oct
              parentloc(ii1,jj1,kk1)=0

              ! Loop over 2x2x2 cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
#if NDIM>0
                       i3=1+2*(ii1-1)+i2
#endif
#if NDIM>1
                       j3=1+2*(jj1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(kk1-1)+k2
#endif             
                       ! Gather hydro variables
                       do ivar=1,nvar
                          uloc(i3,j3,k3,ivar)=grid(ind_oct)%uold(ind_son,ivar)
                       end do
                       ! Gather refinement flag
                       okloc(i3,j3,k3)=grid(ind_oct)%refined(ind_son)
                    end do
                 end do
              end do

              ! Go to next inner oct
              ind_oct=ind_oct+1

           ! For boundary octs only
           else
              ! Compute neighboring grid Cartesian index
#if NDIM>0
              hash_nbor(1)=ckey_corner(1)+i1-1.0
#endif
#if NDIM>1
              hash_nbor(2)=ckey_corner(2)+j1-1.0
#endif
#if NDIM>2
              hash_nbor(3)=ckey_corner(3)+k1-1.0
#endif
              ! Periodic boundary conditons
              do idim=1,ndim
                 if(hash_nbor(idim)<0)hash_nbor(idim)=ckey_max(ilevel)-1
                 if(hash_nbor(idim)==ckey_max(ilevel))hash_nbor(idim)=0
              enddo
              
              ! Get neighboring grid index
              ichild=get_grid(hash_nbor,.false.,.true.)
              parent_cell=0
              if(ichild==0)then

                 ! Get parent father cell
                 parent_cell=get_parent_cell(hash_nbor,.true.,.true.)
                 if(parent_cell==0)then
                    write(*,*)'GODUNOV: parent_cell should exist'
                    write(*,*)'PE ',myid,hash_nbor
                    stop
                 endif
                 igrid=(parent_cell-1)/twotondim+1
                 icell=parent_cell-(igrid-1)*twotondim
                 call lock_cache(igrid)
              else
                 call lock_cache(ichild)
              endif

              ! Store grid index
              childloc(i1,j1,k1)=ichild
              parentloc(i1,j1,k1)=parent_cell

              ! Loop over 2x2x2 cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max                       
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
#if NDIM>0
                       i3=1+2*(i1-1)+i2
#endif
#if NDIM>1
                       j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(k1-1)+k2
#endif             
                       ! If neighboring grid exists
                       if(ichild>0)then
                          ! Gather hydro variables
                          do ivar=1,nvar
                             uloc(i3,j3,k3,ivar)=grid(ichild)%uold(ind_son,ivar)
                          end do
                          ! Gather refinement flag
                          okloc(i3,j3,k3)=grid(ichild)%refined(ind_son)
                       else
                          ! Gather hydro variables
                          do ivar=1,nvar
                             uloc(i3,j3,k3,ivar)=grid(igrid)%uold(icell,ivar)
                          end do
                          ! Gather refinement flag
                          okloc(i3,j3,k3)=.false.
                       end if

                    end do
                 end do
              end do
              ! End loop over 2x2x2 cells
           endif
        end do
     end do
  end do
  ! End over octs

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,qloc,cloc,&
       & flux,tmp,dq,qm,qp,fx,tx,divu,&
       & dx,dx,dx,dtnew(ilevel),&
       & iu1,iu2,ju1,ju2,ku1,ku2,&
       & if1,if2,jf1,jf2,kf1,kf2)
  
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
                 if(okloc(i3-i0,j3-j0,k3-k0) .or. okloc(i3,j3,k3))then
                    flux(i3,j3,k3,ivar,idim)=0.0d0
                 end if
              end do
#ifdef DUALENER
              do ivar=1,2
                 if(okloc(i3-i0,j3-j0,k3-k0) .or. okloc(i3,j3,k3))then
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
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     ii0=0; jj0=0; kk0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
#if NDIM>0
     ii0=1
#endif
#if NDIM>1
     jj0=1
#endif
#if NDIM>2
     kk0=1
#endif
     ! Loop over inner octs
     do k1=k1min+kk0,k1max-kk0
        do j1=j1min+jj0,j1max-jj0
           do i1=i1min+ii0,i1max-ii0
              ! Get oct index
              ind_oct=childloc(i1,j1,k1)
              ! Loop over cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
#if NDIM>0
                       i3=1+2*(i1-1)+i2
#endif
#if NDIM>1
                       j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(k1-1)+k2
#endif
                       ! Update conservative variables new state vector
                       do ivar=1,nvar
                          grid(ind_oct)%unew(ind_son,ivar)=&
                               & grid(ind_oct)%unew(ind_son,ivar)+ &
                               & (flux(i3   ,j3   ,k3   ,ivar,idim) &
                               & -flux(i3+i0,j3+j0,k3+k0,ivar,idim))
                       end do
#ifdef DUALENER
                       ! Update velocity divergence
                       grid(ind_oct)%divu(ind_son)=&
                            & grid(ind_oct)%divu(ind_son)+ &
                            & (tmp(i3   ,j3   ,k3   ,1,idim) &
                            & -tmp(i3+i0,j3+j0,k3+k0,1,idim))
                       ! Update internal energy
                       grid(ind_oct)%enew(ind_son)=&
                            & grid(ind_oct)%enew(ind_son)+ &
                            & (tmp(i3   ,j3   ,k3   ,2,idim) &
                            & -tmp(i3+i0,j3+j0,k3+k0,2,idim))
#endif
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! If sitting in coarsest level, exit. 
  if(ilevel>levelmin)then

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     ii0=0; jj0=0; kk0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
#if NDIM>0
     ii0=1
#endif
#if NDIM>1
     jj0=1
#endif
#if NDIM>2
     kk0=1
#endif
     ii1min=i1min+ii0; ii1max=i1max-ii0
     jj1min=j1min+jj0; jj1max=j1max-jj0
     kk1min=k1min+kk0; kk1max=k1max-kk0
     !----------------------
     ! Left flux at boundary
     !----------------------     
     if(idim==1)then
        ii1min=i1min; ii1max=i1min
     endif
     if(idim==2)then
        jj1min=j1min; jj1max=j1min
     endif
     if(idim==3)then
        kk1min=k1min; kk1max=k1min
     endif
     ! Loop over outer octs on left face
     do k1=kk1min,kk1max
        do j1=jj1min,jj1max
           do i1=ii1min,ii1max
              ! Get oct index
              ind_oct=childloc(i1,j1,k1)
              ! Check that parent cell is not refined
              if(ind_oct==0)then
                 ! Get parent cell index
                 parent_cell=parentloc(i1,j1,k1)
                 igrid=(parent_cell-1)/twotondim+1
                 icell=parent_cell-(igrid-1)*twotondim
                 ! Loop over inner cell left faces
                 do k2=k2min,k2max-k0
                    do j2=j2min,j2max-j0
                       do i2=i2min,i2max-i0
                          i3=1; j3=1; k3=1
#if NDIM>0
                          i3=1+2*(i1+i0-1)+i2
#endif
#if NDIM>1
                          j3=1+2*(j1+j0-1)+j2
#endif
#if NDIM>2
                          k3=1+2*(k1+k0-1)+k2
#endif
                          ! Conservative update of new state variables
                          do ivar=1,nvar
                             grid(igrid)%unew(icell,ivar)=grid(igrid)%unew(icell,ivar) &
                                  & -flux(i3,j3,k3,ivar,idim)*oneontwotondim
                          end do
#ifdef DUALENER
                          ! Update velocity divergence
                          grid(igrid)%divu(icell)=grid(igrid)%divu(icell) &
                               & -tmp(i3,j3,k3,1,idim)*oneontwotondim
                          ! Update internal energy
                          grid(igrid)%enew(icell)=grid(igrid)%enew(icell) &
                               & -tmp(i3,j3,k3,2,idim)*oneontwotondim
#endif
                       end do
                    end do
                 end do
              endif
           end do
        end do
     end do
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     if(idim==1)then
        ii1min=i1max; ii1max=i1max
     endif
     if(idim==2)then
        jj1min=j1max; jj1max=j1max
     endif
     if(idim==3)then
        kk1min=k1max; kk1max=k1max
     endif
     ! Loop over outer octs on right face
     do k1=kk1min,kk1max
        do j1=jj1min,jj1max
           do i1=ii1min,ii1max
              ! Get oct index
              ind_oct=childloc(i1,j1,k1)
              ! Check that parent cell is not refined
              if(ind_oct==0)then
                 ! Get parent cell index
                 parent_cell=parentloc(i1,j1,k1)
                 igrid=(parent_cell-1)/twotondim+1
                 icell=parent_cell-(igrid-1)*twotondim
                 ! Loop over inner cell right faces
                 do k2=k2min+k0,k2max
                    do j2=j2min+j0,j2max
                       do i2=i2min+i0,i2max
                          i3=1; j3=1; k3=1
#if NDIM>0
                          i3=1+2*(i1-i0-1)+i2
#endif
#if NDIM>1
                          j3=1+2*(j1-j0-1)+j2
#endif
#if NDIM>2
                          k3=1+2*(k1-k0-1)+k2
#endif
                          ! Conservative update of new state variables
                          do ivar=1,nvar
                             grid(igrid)%unew(icell,ivar)=grid(igrid)%unew(icell,ivar) &
                                  & +flux(i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                          end do
#ifdef DUALENER
                          ! Update velocity divergence
                          grid(igrid)%divu(icell)=grid(igrid)%divu(icell) &
                               & +tmp(i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
                          ! Update internal energy
                          grid(igrid)%enew(icell)=grid(igrid)%enew(incell) &
                               & +tmp(i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
#endif
                       end do
                    end do
                 end do
                 ! End loop over faces
              endif
           end do
        end do
     end do
     ! End loop over boundary octs
  end do
  ! End loop over dimensions
  endif

  ! Unlock all octs
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max     
           ! Get oct index
           ind_oct=childloc(i1,j1,k1)
           ! Check that parent cell is not refined
           if(ind_oct>0)then
              call unlock_cache(ind_oct)
           else
              ! Get parent cell index
              parent_cell=parentloc(i1,j1,k1)
              igrid=(parent_cell-1)/twotondim+1
              icell=parent_cell-(igrid-1)*twotondim
              call unlock_cache(igrid)
           endif
        end do
     end do
  end do

end subroutine godfine1

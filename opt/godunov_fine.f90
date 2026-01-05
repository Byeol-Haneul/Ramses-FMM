module godunov_fine_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_godunov_fine(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_GODUNOV_FINE,pst%iUpper+1,input_size,0,ilevel)
     call r_godunov_fine(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call godunov_fine(pst%s,ilevel)
  endif

end subroutine r_godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(s,ilevel)
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine, unpack_fetch_refine
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  type(msg_large_realdp)::dummy_large_realdp
  integer::igrid

  if(s%r%verbose.and.s%g%myid==1)write(*,'("   Entering godunov_fine for level ",I2)')ilevel

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  call open_cache(mdl, m, pack_size=storage_size(dummy_large_realdp)/32, &
       pack=pack_fetch_refine, unpack=unpack_fetch_refine, &
       init=init_flush_godunov, flush=pack_flush_godunov, &
       combine=unpack_flush_godunov, bound=init_bound_refine)

  ! Collect or create all boundary (ghost) octs.
  ! These could be octs at processor domain boundaries,
  ! octs at the physical domain boundaries, 
  ! or octs at coarse-fine boundaries.

  call make_boundaries(s,ilevel)

  ! Loop over active grids by vector sweeps
  igrid=m%head(ilevel)
  do while(igrid.LE.m%tail(ilevel))
     SELECT CASE (m%grid(igrid)%superoct)
     CASE(1)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_1)
     CASE(2**ndim)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_2)
     CASE(4**ndim)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_4)
     CASE(8**ndim)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_8)
     CASE(16**ndim)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_16)
     CASE(32**ndim)
        call godfine1(s,igrid,ilevel,m%hydro_w%kernel_32)
     END SELECT
     igrid=igrid+m%grid(igrid)%superoct
  end do

  call close_cache(mdl)

  end associate
end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_unew(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_UNEW,pst%iUpper+1,input_size,0,ilevel)
     call r_set_unew(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call set_unew(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i

#ifdef HYDRO
  ! Set unew to uold for myid cells
  do i = m%head(ilevel),m%tail(ilevel)
     m%unew(:,:,i) = m%uold(:,:,i)
  end do
#endif
#ifdef MHD
  ! Set bnew to bold for myid cells
  do i = m%head(ilevel),m%tail(ilevel)
     m%bnew(:,:,i) = m%bold(:,:,i)
  end do
#endif

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_uold(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_UOLD,pst%iUpper+1,input_size,0,ilevel)
     call r_set_uold(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call set_uold(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew 
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i

#ifdef HYDRO
  ! Set uold to unew
  do i = m%head(ilevel),m%tail(ilevel)
     m%uold(:,:,i) = m%unew(:,:,i)
  end do
#endif
#ifdef MHD
  ! Set bold to bnew
  do i = m%head(ilevel),m%tail(ilevel)
     m%bold(:,:,i) = m%bnew(:,:,i)
  end do
#endif

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(s,ind_grid,ilevel,h)
  use mdl_module
  use amr_parameters, only: ndim, twondim, twotondim
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_commons
  use hash
  implicit none
  type(ramses_t)::s
  integer::ind_grid,ilevel
  type(hydro_kernel_t)::h
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
  integer::ivar,idim,ind_son,ind_oct
  integer::icell,inbor,icache,ichild
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::ii0,jj0,kk0,ii1,jj1,kk1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::ii1min,ii1max,jj1min,jj1max,kk1min,kk1max
#ifdef MHD
  real(kind=8)::dflux,dflux_x,dflux_y,dflux_z,weight
  integer::igrid1,igrid2,igrid3
  integer::icell1,icell2,icell3
  logical::ok1,ok2,ok3
#endif
  integer,dimension(1:ndim)::ckey_corner,ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor
  real(kind=8)::dx,oneontwotondim
  integer::igrid,ichild

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  oneontwotondim=1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=r%boxlen/2**ilevel

  ! Integer constants
  i1min=h%io1; i1max=h%io2; j1min=h%jo1; j1max=h%jo2; k1min=h%ko1; k1max=h%ko2
  i2min=0; i2max=0; j2min=0; j2max=0; k2min=0; k2max=0
  i3min=1; i3max=1; j3min=1; j3max=1; k3min=1; k3max=1
#if NDIM>0
  i2max=1; i3min=h%iu1+2; i3max=h%iu2-2
#endif
#if NDIM>1
  j2max=1; j3min=h%ju1+2; j3max=h%ju2-2
#endif
#if NDIM>2
  k2max=1; k3min=h%ku1+2; k3max=h%ku2-2
#endif

  !---------------------
  ! Gather hydro stencil
  !---------------------
  hash_nbor(0)=m%grid(ind_grid)%lev
  ckey_corner(1:ndim)=(m%grid(ind_grid)%ckey(1:ndim)/(i1max-1))*(i1max-1)

  ! Loop over 3x3x3 octs
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max

           ! Compute Cartesian index
           hash_nbor(1)=ckey_corner(1)+i1-1
#if NDIM>1
           hash_nbor(2)=ckey_corner(2)+j1-1
#endif
#if NDIM>2
           hash_nbor(3)=ckey_corner(3)+k1-1
#endif
           ! Periodic boundary conditions
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
              endif
           enddo

           ! Get neighboring grid index using hash table
           ichild=hash_getp(m%grid_dict,hash_nbor)
           h%childloc(i1,j1,k1)=ichild

        end do
     end do
  end do
  ! End loop over oct

  !-----------------------
  ! Gather hydro variables
  !-----------------------
  ! Loop over 3x3x3 octs
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max
           ! Get oct index
           ichild=h%childloc(i1,j1,k1)
           ! Loop over 2x2x2 cells
           do k2=k2min,k2max
              do j2=j2min,j2max
                 do i2=i2min,i2max                       
                    ind_son=1+i2+2*j2+4*k2
                    i3=1; j3=1; k3=1
                    i3=1+2*(i1-1)+i2
#if NDIM>1
                    j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                    k3=1+2*(k1-1)+k2
#endif             
                    ! Gather hydro variables
                    do ivar=1,nvar
                       h%uloc(i3,j3,k3,ivar)=m%uold(ind_son,ivar,ichild)
                    end do
#ifdef MHD
                    ! Gather MHD variables
                    do ivar=1,6
                       h%bloc(i3,j3,k3,ivar)=m%bold(ind_son,ivar,ichild)
                    end do
#endif
#ifdef GRAV
                    ! Gather self-gravitational acceleration
                    do idim=1,ndim
                       h%gloc(i3,j3,k3,idim)=m%f(ind_son,idim,ichild)
                    end do
#else
                    ! Gather constant gravitational acceleration
                    do idim=1,ndim
                       h%gloc(i3,j3,k3,idim)=r%constant_gravity(idim)
                    end do
#endif
                    ! Gather refinement flag
                    h%okloc(i3,j3,k3)=m%grid(ichild)%refined(ind_son)

                 end do
              end do
           end do
           ! End loop over 2x2x2 cells
        end do
     end do
  end do
  ! End over octs

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(h%uloc,h%gloc,h%qloc,h%cloc,&
       & h%flux,h%tmp,h%dq,h%qm,h%qp,h%fx,h%tx,h%divu,&
#ifdef MHD
       & h%bloc,h%emfx,h%emfy,h%emfz,h%bf,h%dbf,h%Ex,h%Ey,h%Ez,h%qRT,h%qRB,h%qLT,h%qLB,&
#endif
       & dx,g%dtnew(ilevel),&
       & h%iu1,h%iu2,h%ju1,h%ju2,h%ku1,h%ku2,&
       & h%if1,h%if2,h%jf1,h%jf2,h%kf1,h%kf2,&
       & s%h_params)

  ! If finest level, skip
  if(ilevel < r%nlevelmax)then

  !-------------------------------------------------
  ! Reset flux along direction at refined interfaces
  !-------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
        do j3=j3min,j3max+j0
           do i3=i3min,i3max+i0
              if(h%okloc(i3-i0,j3-j0,k3-k0) .or. h%okloc(i3,j3,k3))then
                 do ivar=1,nprim
                    h%flux(i3,j3,k3,ivar,idim)=0.0d0
                 end do
              end if
           end do
        end do
     end do
  end do
#ifdef MHD
  ! In case of pure induction, set Euler fluxes to zero
  if(r%induction)h%flux=0.0d0
#if NDIM>1
  !---------------------------------------------------
  ! Reset electromotive force along z at refined edges
  !---------------------------------------------------
  do k3=k3min,k3max
     do j3=j3min,j3max+1
        do i3=i3min,i3max+1
           if(    h%okloc(i3-1,j3-1,k3) .or. h%okloc(i3  ,j3-1,k3) .or.  &
                & h%okloc(i3-1,j3  ,k3) .or. h%okloc(i3  ,j3  ,k3))then
              h%emfz(i3,j3,k3)=0.0d0
           end if
        end do
     end do
  end do
#if NDIM==3
  !---------------------------------------------------
  ! Reset electromotive force along y at refined edges
  !---------------------------------------------------
  do j3=j3min,j3max
     do i3=i3min,i3max+1
        do k3=k3min,k3max+1
           if(    h%okloc(i3-1,j3,k3-1) .or. h%okloc(i3  ,j3,k3-1) .or.  &
                & h%okloc(i3-1,j3,k3  ) .or. h%okloc(i3  ,j3,k3  ))then
              h%emfy(i3,j3,k3)=0.0d0
           end if
        end do
     end do
  end do
  !---------------------------------------------------
  ! Reset electromotive force along x at refined edges
  !---------------------------------------------------
  do i3=i3min,i3max
     do j3=j3min,j3max+1
        do k3=k3min,k3max+1
           if(    h%okloc(i3,j3-1,k3-1) .or. h%okloc(i3,j3  ,k3-1) .or.  &
                & h%okloc(i3,j3-1,k3  ) .or. h%okloc(i3,j3  ,k3  ))then
              h%emfx(i3,j3,k3)=0.0d0
           end if
        end do
     end do
  end do
#endif
#endif
#endif

  endif

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
     ii0=1
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
              ichild=h%childloc(i1,j1,k1)
              ! Loop over cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
                       i3=1+2*(i1-1)+i2
#if NDIM>1
                       j3=1+2*(j1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(k1-1)+k2
#endif
                       ! Update conservative variables new state vector
                       do ivar=1,5
                          m%unew(ind_son,ivar,ichild)=m%unew(ind_son,ivar,ichild)+ &
                               & (h%flux(i3   ,j3   ,k3   ,ivar,idim) &
                               & -h%flux(i3+i0,j3+j0,k3+k0,ivar,idim))
                       end do
#if NVAR>5
                       do ivar=6,nvar
                          m%unew(ind_son,ivar,ichild)=m%unew(ind_son,ivar,ichild)+ &
                               & (h%flux(i3   ,j3   ,k3   ,ivar+ie-5,idim) &
                               & -h%flux(i3+i0,j3+j0,k3+k0,ivar+ie-5,idim))
                       end do
#endif
#ifdef MHD
#if NDIM<3
                       m%bnew(ind_son,3,ichild)=m%bnew(ind_son,3,ichild)+ &
                            & (h%flux(i3   ,j3   ,k3   ,8,idim) &
                            & -h%flux(i3+i0,j3+j0,k3+k0,8,idim))
                       m%bnew(ind_son,6,ichild)=m%bnew(ind_son,6,ichild)+ &
                            & (h%flux(i3   ,j3   ,k3   ,8,idim) &
                            & -h%flux(i3+i0,j3+j0,k3+k0,8,idim))
#if NDIM==1
                       m%bnew(ind_son,2,ichild)=m%bnew(ind_son,2,ichild)+ &
                            & (h%flux(i3   ,j3   ,k3   ,7,idim) &
                            & -h%flux(i3+i0,j3+j0,k3+k0,7,idim))
                       m%bnew(ind_son,5,ichild)=m%bnew(ind_son,5,ichild)+ &
                            & (h%flux(i3   ,j3   ,k3   ,7,idim) &
                            & -h%flux(i3+i0,j3+j0,k3+k0,7,idim))
#endif
#endif
#endif
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
#ifdef MHD
#if NDIM>1

  !----------------------------------------
  ! Induction system update at level ilevel
  !----------------------------------------
  ii0=1; jj0=1; kk0=0
#if NDIM==3
  kk0=1
#endif
  ! Loop over inner octs
  do k1=k1min+kk0,k1max-kk0
     do j1=j1min+jj0,j1max-jj0
        do i1=i1min+ii0,i1max-ii0
           ! Get oct index
           ichild=h%childloc(i1,j1,k1)
           ! Loop over cells
           do k2=k2min,k2max
              do j2=j2min,j2max
                 do i2=i2min,i2max
                    ind_son=1+i2+2*j2+4*k2
                    i3=1+2*(i1-1)+i2
                    j3=1+2*(j1-1)+j2
                    k3=1
#if NDIM==3
                    k3=1+2*(k1-1)+k2
#endif
                    ! Update Bx using constrained transport
#if NDIM==3
                    dflux_x=( h%emfy(i3,j3,k3)-h%emfy(i3,j3,k3+1) ) &
                         & -( h%emfz(i3,j3,k3)-h%emfz(i3,j3+1,k3) )
#else
                    dflux_x=-(h%emfz(i3,j3,k3)-h%emfz(i3,j3+1,k3))
#endif
                    m%bnew(ind_son,1,ichild)=m%bnew(ind_son,1,ichild)+dflux_x
#if NDIM==3
                    dflux_x=( h%emfy(i3+1,j3,k3)-h%emfy(i3+1,j3,k3+1) ) &
                         & -( h%emfz(i3+1,j3,k3)-h%emfz(i3+1,j3+1,k3) )
#else
                    dflux_x=-(h%emfz(i3+1,j3,k3)-h%emfz(i3+1,j3+1,k3))
#endif
                    m%bnew(ind_son,4,ichild)=m%bnew(ind_son,4,ichild)+dflux_x

                    ! Update By using constrained transport
#if NDIM==3
                    dflux_y=( h%emfz(i3,j3,k3)-h%emfz(i3+1,j3,k3) ) &
                         & -( h%emfx(i3,j3,k3)-h%emfx(i3,j3,k3+1) )
#else
                    dflux_y=(h%emfz(i3,j3,k3)-h%emfz(i3+1,j3,k3))
#endif
                    m%bnew(ind_son,2,ichild)=m%bnew(ind_son,2,ichild)+dflux_y
#if NDIM==3
                    dflux_y=( h%emfz(i3,j3+1,k3)-h%emfz(i3+1,j3+1,k3) ) &
                         & -( h%emfx(i3,j3+1,k3)-h%emfx(i3,j3+1,k3+1) )
#else
                    dflux_y=(h%emfz(i3,j3+1,k3)-h%emfz(i3+1,j3+1,k3))
#endif
                    m%bnew(ind_son,5,ichild)=m%bnew(ind_son,5,ichild)+dflux_y
#if NDIM==3
                    ! Update Bz using constrained transport
                    dflux_z=( h%emfx(i3,j3,k3)-h%emfx(i3,j3+1,k3) ) &
                         & -( h%emfy(i3,j3,k3)-h%emfy(i3+1,j3,k3) )
                    m%bnew(ind_son,3,ichild)=m%bnew(ind_son,3,ichild)+dflux_z
                    dflux_z=( h%emfx(i3,j3,k3+1)-h%emfx(i3,j3+1,k3+1) ) &
                         & -( h%emfy(i3,j3,k3+1)-h%emfy(i3+1,j3,k3+1) )
                    m%bnew(ind_son,6,ichild)=m%bnew(ind_son,6,ichild)+dflux_z
#endif
                 end do
              end do
           end do
        end do
     end do
  end do
#endif
#endif

  ! If coarsest level, skip.
  if(ilevel>r%levelmin)then

  !----------------------------
  ! Gather parent grid and cell
  !----------------------------
  ! Loop over 3x3x3 octs
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max
           ! Set parent grid and cell index to 0
           h%gridloc(i1,j1,k1)=0
           h%cellloc(i1,j1,k1)=0
           ! Get child grid index
           ichild=h%childloc(i1,j1,k1)
           ! Gather parent grid and parent cell only if child is a ghost grid
           if(ichild>m%ngridmax)then
              icache=ichild-m%ngridmax
              if(m%ghost_parent_grid(icache)>0)then
                 h%gridloc(i1,j1,k1)=m%ghost_parent_grid(icache)
                 h%cellloc(i1,j1,k1)=m%ghost_parent_cell(icache)
              endif
           endif
        end do
     end do
  end do

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
     ii0=1
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
              ! Check that grid is a ghost
              igrid=h%gridloc(i1,j1,k1)
              if(igrid>0)then
                 ! Get parent cell index
                 icell=h%cellloc(i1,j1,k1)
                 ! Loop over inner cell left faces
                 do k2=k2min,k2max-k0
                    do j2=j2min,j2max-j0
                       do i2=i2min,i2max-i0
                          i3=1; j3=1; k3=1
                          i3=1+2*(i1+i0-1)+i2
#if NDIM>1
                          j3=1+2*(j1+j0-1)+j2
#endif
#if NDIM>2
                          k3=1+2*(k1+k0-1)+k2
#endif
                          ! Conservative update of new state variables
                          do ivar=1,5
                             m%unew(icell,ivar,igrid)=m%unew(icell,ivar,igrid) &
                                  & -h%flux(i3,j3,k3,ivar,idim)*oneontwotondim
                          end do
#if NVAR>5
                          do ivar=6,nvar
                             m%unew(icell,ivar,igrid)=m%unew(icell,ivar,igrid) &
                                  & -h%flux(i3,j3,k3,ivar+ie-5,idim)*oneontwotondim
                          end do
#endif
#ifdef MHD
#if NDIM<3
                          m%bnew(icell,3,igrid)=m%bnew(icell,3,igrid) &
                               & -h%flux(i3,j3,k3,8,idim)*oneontwotondim
                          m%bnew(icell,6,igrid)=m%bnew(icell,6,igrid) &
                               & -h%flux(i3,j3,k3,8,idim)*oneontwotondim
#if NDIM==1
                          m%bnew(icell,2,igrid)=m%bnew(icell,2,igrid) &
                               & -h%flux(i3,j3,k3,7,idim)*oneontwotondim
                          m%bnew(icell,5,igrid)=m%bnew(icell,5,igrid) &
                               & -h%flux(i3,j3,k3,7,idim)*oneontwotondim
#endif
#endif
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
              ! Check that grid is a ghost
              igrid=h%gridloc(i1,j1,k1)
              if(igrid>0)then
                 ! Get parent cell index
                 icell=h%cellloc(i1,j1,k1)
                 ! Loop over inner cell right faces
                 do k2=k2min+k0,k2max
                    do j2=j2min+j0,j2max
                       do i2=i2min+i0,i2max
                          i3=1; j3=1; k3=1
                          i3=1+2*(i1-i0-1)+i2
#if NDIM>1
                          j3=1+2*(j1-j0-1)+j2
#endif
#if NDIM>2
                          k3=1+2*(k1-k0-1)+k2
#endif
                          ! Conservative update of new state variables
                          do ivar=1,5
                             m%unew(icell,ivar,igrid)=m%unew(icell,ivar,igrid) &
                                  & +h%flux(i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                          end do
#if NVAR>5
                          do ivar=6,nvar
                             m%unew(icell,ivar,igrid)=m%unew(icell,ivar,igrid) &
                                  & +h%flux(i3+i0,j3+j0,k3+k0,ivar+ie-5,idim)*oneontwotondim
                          end do
#endif
#ifdef MHD
#if NDIM<3
                          m%bnew(icell,3,igrid)=m%bnew(icell,3,igrid) &
                               & +h%flux(i3+i0,j3+j0,k3+k0,8,idim)*oneontwotondim
                          m%bnew(icell,6,igrid)=m%bnew(icell,6,igrid) &
                               & +h%flux(i3+i0,j3+j0,k3+k0,8,idim)*oneontwotondim
#if NDIM==1
                          m%bnew(icell,2,igrid)=m%bnew(icell,2,igrid) &
                               & +h%flux(i3+i0,j3+j0,k3+k0,7,idim)*oneontwotondim
                          m%bnew(icell,5,igrid)=m%bnew(icell,5,igrid) &
                               & +h%flux(i3+i0,j3+j0,k3+k0,7,idim)*oneontwotondim
#endif
#endif
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

#ifdef MHD
#if NDIM>1
  !------------------------------------------
  ! Induction system update at level ilevel-1
  !------------------------------------------
  ii0=1; jj0=1
  kk0=0
#if NDIM==3
  kk0=1
#endif
  ii1min=i1min+ii0; ii1max=i1max-ii0

  jj1min=j1min+jj0; jj1max=j1max-jj0
  kk1min=k1min+kk0; kk1max=k1max-kk0
  ! Loop over inner octs
  do k1=kk1min,kk1max
     do j1=jj1min,jj1max
        do i1=ii1min,ii1max

           i3=1+2*(i1-1)
           j3=1+2*(j1-1)
           k3=1
#if NDIM==3
           k3=1+2*(k1-1)
#endif
           !--------------------------------------
           ! Deal with 4 EMFz edges
           !--------------------------------------

           ! Update coarse Bx and By using fine EMFz on X=0 and Y=0 grid edge
           igrid1=h%gridloc(i1  ,j1-1,k1); icell1=h%cellloc(i1  ,j1-1,k1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1-1,j1-1,k1); icell2=h%cellloc(i1-1,j1-1,k1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1-1,j1  ,k1); icell3=h%cellloc(i1-1,j1  ,k1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
#if NDIM==3
              dflux=(h%emfz(i3,j3,k3)+h%emfz(i3,j3,k3+1))*0.25*weight
#else
              dflux=(h%emfz(i3,j3,k3))*0.5*weight
#endif
              if(ok1)grid1%bnew(icell1,5)=grid1%bnew(icell1,5)+dflux
              if(ok1)grid1%bnew(icell1,1)=grid1%bnew(icell1,1)+dflux
              if(ok2)grid2%bnew(icell2,4)=grid2%bnew(icell2,4)+dflux
              if(ok2)grid2%bnew(icell2,5)=grid2%bnew(icell2,5)-dflux
              if(ok3)grid3%bnew(icell3,2)=grid3%bnew(icell3,2)-dflux
              if(ok3)grid3%bnew(icell3,4)=grid3%bnew(icell3,4)-dflux
           endif

           ! Update coarse Bx and By using fine EMFz on X=0 and Y=1 grid edge
           igrid1=h%gridloc(i1-1,j1  ,k1); icell1=h%cellloc(i1-1,j1  ,k1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1-1,j1+1,k1); icell2=h%cellloc(i1-1,j1+1,k1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1  ,j1+1,k1); icell3=h%cellloc(i1  ,j1+1,k1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
#if NDIM==3
              dflux=(h%emfz(i3,j3+2,k3)+h%emfz(i3,j3+2,k3+1))*0.25*weight
#else
              dflux=(h%emfz(i3,j3+2,k3))*0.5*weight
#endif
              if(ok1)grid1%bnew(icell1,4)=grid1%bnew(icell1,4)+dflux
              if(ok1)grid1%bnew(icell1,5)=grid1%bnew(icell1,5)-dflux
              if(ok2)grid2%bnew(icell2,2)=grid2%bnew(icell2,2)-dflux
              if(ok2)grid2%bnew(icell2,4)=grid2%bnew(icell2,4)-dflux
              if(ok3)grid3%bnew(icell3,1)=grid3%bnew(icell3,1)-dflux
              if(ok3)grid3%bnew(icell3,2)=grid3%bnew(icell3,2)+dflux
           endif

           ! Update coarse Bx and By using fine EMFz on X=1 and Y=1 grid edge
           igrid1=h%gridloc(i1  ,j1+1,k1); icell1=h%cellloc(i1  ,j1+1,k1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1+1,j1+1,k1); icell2=h%cellloc(i1+1,j1+1,k1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1+1,j1  ,k1); icell3=h%cellloc(i1+1,j1  ,k1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
#if NDIM==3
              dflux=(h%emfz(i3+2,j3+2,k3)+h%emfz(i3+2,j3+2,k3+1))*0.25*weight
#else
              dflux=(h%emfz(i3+2,j3+2,k3))*0.5*weight
#endif
              if(ok1)grid1%bnew(icell1,2)=grid1%bnew(icell1,2)-dflux
              if(ok1)grid1%bnew(icell1,4)=grid1%bnew(icell1,4)-dflux
              if(ok2)grid2%bnew(icell2,1)=grid2%bnew(icell2,1)-dflux
              if(ok2)grid2%bnew(icell2,2)=grid2%bnew(icell2,2)+dflux
              if(ok3)grid3%bnew(icell3,5)=grid3%bnew(icell3,5)+dflux
              if(ok3)grid3%bnew(icell3,1)=grid3%bnew(icell3,1)+dflux
           endif

           ! Update coarse Bx and By using fine EMFz on X=1 and Y=0 grid edge
           igrid1=h%gridloc(i1+1,j1  ,k1); icell1=h%cellloc(i1+1,j1  ,k1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1+1,j1-1,k1); icell2=h%cellloc(i1+1,j1-1,k1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1  ,j1-1,k1); icell3=h%cellloc(i1  ,j1-1,k1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
#if NDIM==3
              dflux=(h%emfz(i3+2,j3,k3)+h%emfz(i3+2,j3,k3+1))*0.25*weight
#else
              dflux=(h%emfz(i3+2,j3,k3))*0.5*weight
#endif
              if(ok1)grid1%bnew(icell1,1)=grid1%bnew(icell1,1)-dflux
              if(ok1)grid1%bnew(icell1,2)=grid1%bnew(icell1,2)+dflux
              if(ok2)grid2%bnew(icell2,5)=grid2%bnew(icell2,5)+dflux
              if(ok2)grid2%bnew(icell2,1)=grid2%bnew(icell2,1)+dflux
              if(ok3)grid3%bnew(icell3,4)=grid3%bnew(icell3,4)+dflux
              if(ok3)grid3%bnew(icell3,5)=grid3%bnew(icell3,5)-dflux
           endif
#if NDIM==3
           !--------------------------------------
           ! Deal with 4 EMFx edges
           !--------------------------------------

           ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
           igrid1=h%gridloc(i1,j1  ,k1-1); icell1=h%cellloc(i1,j1  ,k1-1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1,j1-1,k1-1); icell2=h%cellloc(i1,j1-1,k1-1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1,j1-1,k1  ); icell3=h%cellloc(i1,j1-1,k1  ); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfx(i3,j3,k3)+h%emfx(i3+1,j3,k3))*0.25*weight
              if(ok1)grid1%bnew(icell1,6)=grid1%bnew(icell1,6)+dflux
              if(ok1)grid1%bnew(icell1,2)=grid1%bnew(icell1,2)+dflux
              if(ok2)grid2%bnew(icell2,5)=grid2%bnew(icell2,5)+dflux
              if(ok2)grid2%bnew(icell2,6)=grid2%bnew(icell2,6)-dflux
              if(ok3)grid3%bnew(icell3,3)=grid3%bnew(icell3,3)-dflux
              if(ok3)grid3%bnew(icell3,5)=grid3%bnew(icell3,5)-dflux
           endif

           ! Update coarse By and Bz using fine EMFx on Y=0 and Z=1 grid edge
           igrid1=h%gridloc(i1,j1-1,k1  ); icell1=h%cellloc(i1,j1-1,k1  ); ok1=(igrid1>0)
           igrid2=h%gridloc(i1,j1-1,k1+1); icell2=h%cellloc(i1,j1-1,k1+1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1,j1  ,k1+1); icell3=h%cellloc(i1,j1  ,k1+1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfx(i3,j3,k3+2)+h%emfx(i3+1,j3,k3+2))*0.25*weight
              if(ok1)grid1%bnew(icell1,5)=grid1%bnew(icell1,5)+dflux
              if(ok1)grid1%bnew(icell1,6)=grid1%bnew(icell1,6)-dflux
              if(ok2)grid2%bnew(icell2,3)=grid2%bnew(icell2,3)-dflux
              if(ok2)grid2%bnew(icell2,5)=grid2%bnew(icell2,5)-dflux
              if(ok3)grid3%bnew(icell3,2)=grid3%bnew(icell3,2)-dflux
              if(ok3)grid3%bnew(icell3,3)=grid3%bnew(icell3,3)+dflux
           endif

           ! Update coarse By and Bz using fine EMFx on Y=1 and Z=1 grid edge
           igrid1=h%gridloc(i1,j1  ,k1+1); icell1=h%cellloc(i1,j1  ,k1+1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1,j1+1,k1+1); icell2=h%cellloc(i1,j1+1,k1+1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1,j1+1,k1  ); icell3=h%cellloc(i1,j1+1,k1  ); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfx(i3,j3+2,k3+2)+h%emfx(i3+1,j3+2,k3+2))*0.25*weight
              if(ok1)grid1%bnew(icell1,3)=grid1%bnew(icell1,3)-dflux
              if(ok1)grid1%bnew(icell1,5)=grid1%bnew(icell1,5)-dflux
              if(ok2)grid2%bnew(icell2,2)=grid2%bnew(icell2,2)-dflux
              if(ok2)grid2%bnew(icell2,3)=grid2%bnew(icell2,3)+dflux
              if(ok3)grid3%bnew(icell3,6)=grid3%bnew(icell3,6)+dflux
              if(ok3)grid3%bnew(icell3,2)=grid3%bnew(icell3,2)+dflux
           endif

           ! Update coarse By and Bz using fine EMFx on Y=1 and Z=0 grid edge
           igrid1=h%gridloc(i1,j1+1,k1  ); icell1=h%cellloc(i1,j1+1,k1  ); ok1=(igrid1>0)
           igrid2=h%gridloc(i1,j1+1,k1-1); icell2=h%cellloc(i1,j1+1,k1-1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1,j1  ,k1-1); icell3=h%cellloc(i1,j1  ,k1-1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfx(i3,j3+2,k3)+h%emfx(i3+1,j3+2,k3))*0.25*weight
              if(ok1)grid1%bnew(icell1,2)=grid1%bnew(icell1,2)-dflux
              if(ok1)grid1%bnew(icell1,3)=grid1%bnew(icell1,3)+dflux
              if(ok2)grid2%bnew(icell2,6)=grid2%bnew(icell2,6)+dflux
              if(ok2)grid2%bnew(icell2,2)=grid2%bnew(icell2,2)+dflux
              if(ok3)grid3%bnew(icell3,5)=grid3%bnew(icell3,5)+dflux
              if(ok3)grid3%bnew(icell3,6)=grid3%bnew(icell3,6)-dflux
           endif

           !--------------------------------------
           ! Deal with 4 EMFy edges
           !--------------------------------------

           ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=0 grid edge
           igrid1=h%gridloc(i1  ,j1,k1-1); icell1=h%cellloc(i1  ,j1,k1-1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1-1,j1,k1-1); icell2=h%cellloc(i1-1,j1,k1-1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1-1,j1,k1  ); icell3=h%cellloc(i1-1,j1,k1  ); ok3=(igrid3>0)

           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfy(i3,j3,k3)+h%emfy(i3,j3+1,k3))*0.25*weight
              if(ok1)grid1%bnew(icell1,6)=grid1%bnew(icell1,6)-dflux
              if(ok1)grid1%bnew(icell1,1)=grid1%bnew(icell1,1)-dflux
              if(ok2)grid2%bnew(icell2,4)=grid2%bnew(icell2,4)-dflux
              if(ok2)grid2%bnew(icell2,6)=grid2%bnew(icell2,6)+dflux
              if(ok3)grid3%bnew(icell3,3)=grid3%bnew(icell3,3)+dflux
              if(ok3)grid3%bnew(icell3,4)=grid3%bnew(icell3,4)+dflux
           endif

           ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=1 grid edge
           igrid1=h%gridloc(i1-1,j1,k1  ); icell1=h%cellloc(i1-1,j1,k1  ); ok1=(igrid1>0)
           igrid2=h%gridloc(i1-1,j1,k1+1); icell2=h%cellloc(i1-1,j1,k1+1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1  ,j1,k1+1); icell3=h%cellloc(i1  ,j1,k1+1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfy(i3,j3,k3+2)+h%emfy(i3,j3+1,k3+2))*0.25*weight
              if(ok1)grid1%bnew(icell1,4)=grid1%bnew(icell1,4)-dflux
              if(ok1)grid1%bnew(icell1,6)=grid1%bnew(icell1,6)+dflux
              if(ok2)grid2%bnew(icell2,3)=grid2%bnew(icell2,3)+dflux
              if(ok2)grid2%bnew(icell2,4)=grid2%bnew(icell2,4)+dflux
              if(ok3)grid3%bnew(icell3,1)=grid3%bnew(icell3,1)+dflux
              if(ok3)grid3%bnew(icell3,3)=grid3%bnew(icell3,3)-dflux
           endif

           ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=1 grid edge
           igrid1=h%gridloc(i1  ,j1,k1+1); icell1=h%cellloc(i1  ,j1,k1+1); ok1=(igrid1>0)
           igrid2=h%gridloc(i1+1,j1,k1+1); icell2=h%cellloc(i1+1,j1,k1+1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1+1,j1,k1  ); icell3=h%cellloc(i1+1,j1,k1  ); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfy(i3+2,j3,k3+2)+h%emfy(i3+2,j3+1,k3+2))*0.25*weight
              if(ok1)grid1%bnew(icell1,3)=grid1%bnew(icell1,3)+dflux
              if(ok1)grid1%bnew(icell1,4)=grid1%bnew(icell1,4)+dflux
              if(ok2)grid2%bnew(icell2,1)=grid2%bnew(icell2,1)+dflux
              if(ok2)grid2%bnew(icell2,3)=grid2%bnew(icell2,3)-dflux
              if(ok3)grid3%bnew(icell3,6)=grid3%bnew(icell3,6)-dflux
              if(ok3)grid3%bnew(icell3,1)=grid3%bnew(icell3,1)-dflux
           endif

           ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=0 grid edge
           igrid1=h%gridloc(i1+1,j1,k1  ); icell1=h%cellloc(i1+1,j1,k1  ); ok1=(igrid1>0)
           igrid2=h%gridloc(i1+1,j1,k1-1); icell2=h%cellloc(i1+1,j1,k1-1); ok2=(igrid2>0)
           igrid3=h%gridloc(i1  ,j1,k1-1); icell3=h%cellloc(i1  ,j1,k1-1); ok3=(igrid3>0)
           if(ok1 .or. ok3)then
              weight=1.0
              if(.not.ok1 .or. .not.ok2 .or. .not.ok3)weight=0.5
              dflux=(h%emfy(i3+2,j3,k3)+h%emfy(i3+2,j3+1,k3))*0.25*weight
              if(ok1)grid1%bnew(icell1,1)=grid1%bnew(icell1,1)+dflux
              if(ok1)grid1%bnew(icell1,3)=grid1%bnew(icell1,3)-dflux
              if(ok2)grid2%bnew(icell2,6)=grid2%bnew(icell2,6)-dflux
              if(ok2)grid2%bnew(icell2,1)=grid2%bnew(icell2,1)-dflux
              if(ok3)grid3%bnew(icell3,4)=grid3%bnew(icell3,4)-dflux
              if(ok3)grid3%bnew(icell3,6)=grid3%bnew(icell3,6)+dflux
           endif
#endif
        end do
     end do
  end do
#endif
#endif

  endif

  end associate

#endif

end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_godunov(mesh,igrid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  type(mesh_t)::mesh
  integer::igrid
  integer(kind=8),dimension(0:ndim)::hash_key

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)

#ifdef HYDRO
  mesh%unew(:,:,igrid)=0.0d0
#endif

#ifdef MHD
  mesh%bnew(:,:,igrid)=0.0d0
#endif

end subroutine init_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_godunov(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_large_realdp)::msg

#ifdef HYDRO
  msg%realdp_hydro=mesh%unew(:,:,igrid)
#endif

#ifdef MHD
  msg%realdp_mhd=mesh%bnew(:,:,igrid)
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_godunov(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  type(mesh_t)::mesh
  integer::igrid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  mesh%unew(:,:,igrid)=mesh%unew(:,:,igrid)+msg%realdp_hydro
#endif

#ifdef MHD
  mesh%bnew(:,:,igrid)=mesh%bnew(:,:,igrid)+msg%realdp_mhd
#endif

end subroutine unpack_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_boundaries(s,ilevel)
  use mdl_module
  use amr_parameters, only: ndim, twondim, twotondim
  use amr_commons, only: nbor, oct
  use hydro_parameters, only: nvar
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine collect remote information and create ghost grids
  ! in the current processor cache memory. It loops over only dirty
  ! octs and fill the entire cache space with neighboring octs.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::dummy_int,close_tag=7,close_id,icpu,info
#endif
  integer::idim,ipass,igrid,ichild,ind_grid,i1,j1,k1
  integer,dimension(1:ndim)::ckey_corner
  integer(kind=8),dimension(0:ndim)::hash_nbor

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  !--------------------
  ! Gather grid stencil
  !--------------------

  ! Loop over dirty octs
  do igrid=1,m%noct_dirty(ilevel)

  ind_grid=m%indx_dirty(m%head_dirty(ilevel)+igrid-1)
  hash_nbor(0)=m%grid(ind_grid)%lev
  ckey_corner(1:ndim)=m%grid(ind_grid)%ckey(1:ndim)

  ! Loop over 3x3x3 neighboring father cells using 27 passes
  i1=1; j1=1; k1=1 ! ipass = 1 skipped
  do ipass = 2, threetondim

     if(MOD(ipass,3) == 2)then ! ipass = 2, 5, 8, 11, 14, 17, 20, 23, 26
        i1=0
     endif
     if(MOD(ipass,3) == 0)then ! ipass = 3, 6, 9, 12, 15, 18, 21, 24, 27
        i1=2
     endif
     if(MOD(ipass,9) == 4)then ! ipass = 4, 13, 22
        i1=1; j1=0
     endif
     if(MOD(ipass,9) == 7)then ! ipass = 7, 16, 25
        i1=1; j1=2
     endif
     if(ipass == 10)then
        i1=1; j1=1; k1=0
     endif
     if(ipass == 19)then
        i1=1; j1=1; k1=2
     endif

     ! Compute neighboring grid Cartesian index
     hash_nbor(1)=ckey_corner(1)+i1-1
#if NDIM>1
     hash_nbor(2)=ckey_corner(2)+j1-1
#endif
#if NDIM>2
     hash_nbor(3)=ckey_corner(3)+k1-1
#endif
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(hash_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
           if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
        endif
     enddo

     ! Get neighboring grid index with read-only cache
     call get_grid(s,hash_nbor,ichild,flush_cache=.false.,fetch_cache=.true.,lock=.true.,use_ghost=.true.)

     ! If grid does not exist, create new ghost grid in cache memory
     if(ichild==0)call make_grid_ghost(s,hash_nbor,ilevel)

  end do
  ! End loop over pass

  end do
  ! End loop over dirty octs

#ifndef WITHOUTMPI
  ! CHECK-IN CHECK-OUT
  if(g%myid.NE.1)then
     call MPI_ISEND(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(mdl,close_id)
     call MPI_IRECV(dummy_int,1,MPI_INTEGER,0,close_tag,MPI_COMM_WORLD,close_id,info)
     call check_mail(mdl,close_id)
  else
     do icpu=2,g%ncpu
        call MPI_IRECV(dummy_int,1,MPI_INTEGER,MPI_ANY_SOURCE,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(mdl,close_id)
     end do
     do icpu=2,g%ncpu
        call MPI_ISEND(dummy_int,1,MPI_INTEGER,icpu-1,close_tag,MPI_COMM_WORLD,close_id,info)
        call check_mail(mdl,close_id)
     end do
  endif
#endif

  end associate

#endif

end subroutine make_boundaries
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_grid_ghost(s,hash_nbor,ilevel)
  use mdl_module
  use amr_parameters, only: ndim, twotondim, twondim
  use hydro_parameters, only: nvar
  use rt_parameters, only: nrtvar, smallnp
  use ramses_commons, only: ramses_t
  use cache_commons
  use nbors_utils
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(0:ndim)::hash_nbor
  integer::ilevel
  !--------------------------------------------------------------
  ! This routine creates a ghost oct at level ilevel.
  ! The new oct is stored in the cache.
  ! Input argument hash_nbor is the hash key of the new oct.
  !--------------------------------------------------------------
  integer::ivar,idim,ind,icell,igrid,inbor,ichild
  integer,dimension(0:twondim)::ind_nbor
  integer,dimension(0:twondim)::igrid_nbor
  real(kind=8),dimension(0:twondim  ,1:nvar)::u1
  real(kind=8),dimension(1:twotondim,1:nvar)::u2
  logical::oknbor
#ifdef MHD
  integer::igridn
  integer,dimension(1:twondim)::igrid_son_nbor
  real(kind=8),dimension(0:twondim  ,1:6)::b1
  real(kind=8),dimension(1:twotondim,1:6)::b2
  real(kind=8),dimension(1:twondim,1:twotondim,1:6)::b3
  integer(kind=8),dimension(0:ndim)::hash_son_nbor
  logical,dimension(1:twondim)::refined
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
#endif

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

#ifndef WITHOUTMPI
  ! If counter is good, check on incoming messages and perform actions
  if(mdl%mail_counter==32)then
     call check_mail(mdl,MPI_REQUEST_NULL)
     mdl%mail_counter=0
  endif
  mdl%mail_counter=mdl%mail_counter+1
#endif

  ! Get parent father cell with read-write cache
  call get_parent_cell(s,hash_nbor,igrid,icell,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
  if(igrid==0)then
     write(*,*)'GODUNOV: parent_cell should exist'
     write(*,*)'PE ',g%myid,hash_nbor
     call mdl_abort(mdl)
  endif

  ! Get 2ndim neighboring father cells with read-write cache
  ! Note that cache grids are locked inside this routine
  call get_twondim_nbor_parent_cell(s,hash_nbor,igrid_nbor,ind_nbor,flush_cache=.true.,fetch_cache=.true.)
  oknbor=.true.
  do inbor=0,twondim
     oknbor=oknbor.and.(igrid_nbor(inbor)>0)
  end do
  if(.not. oknbor)then
     write(*,*)"GODUNOV: parent neighbors should exist"
     write(*,*)'PE ',g%myid,hash_nbor
     write(*,*)igrid_nbor(0)
     do idim=1,ndim
        write(*,*)igrid_nbor(2*idim-1)
        write(*,*)igrid_nbor(2*idim)
     end do
     call mdl_abort(mdl)
  endif

  ! Gather hydro variables
  do inbor=0,twondim
     do ivar=1,nvar
        u1(inbor,ivar)=m%uold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
#ifdef MHD
  ! Gather MHD variables
  do inbor=0,twondim
     do ivar=1,6
        b1(inbor,ivar)=m%bold(ind_nbor(inbor),ivar,igrid_nbor(inbor))
     end do
  end do
  ! Get neighboring children grids
  hash_son_nbor(0)=ilevel
  do inbor=1,twondim
     hash_son_nbor(1:ndim)=hash_nbor(1:ndim)+shift(1:ndim,inbor)
     ! Periodic boundary conditions
     do idim=1,ndim
        if(r%periodic(idim))then
           if(hash_son_nbor(idim)< m%box_ckey_min(idim,ilevel))hash_son_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
           if(hash_son_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_son_nbor(idim)=m%box_ckey_min(idim,ilevel)
        endif
     enddo
     call get_grid(s,hash_son_nbor,igridn,flush_cache=.false.,fetch_cache=.true.,lock=.true.,use_ghost=.true.)
     igrid_son_nbor(inbor)=igridn
     refined(inbor)=(igridn>0)
     if(refined(inbor))then
        do ind=1,twotondim
           do ivar=1,6
              b3(inbor,ind,ivar)=m%bold(ind,ivar,igridn)
           end do
        end do
     endif
  end do
  ! Interpolate using MHD variables
  call interpol_mhd(u1,u2,b1,b2,b3,refined,r%interpol_var,r%interpol_type,r%smallr)
#else
  ! Interpolate using hydro variables
  call interpol_hydro(u1,u2,r%interpol_var,r%interpol_type,r%smallr)
#endif

  ! If next cache line is occupied, free it.
  if(m%locked(m%free_cache))then
     do while(m%locked(m%free_cache))
        m%free_cache=m%free_cache+1
        if(m%free_cache>m%ncachemax)m%free_cache=1
     end do
  end if
  if(m%occupied(m%free_cache))call destage(m,m%ngridmax+m%free_cache)

  ! Set grid index to a virtual grid in local memory
  ichild=m%ngridmax+m%free_cache
  call hash_setp(m%grid_dict,hash_nbor,ichild)

  ! Store grid coordinates
  m%grid(ichild)%lev=hash_nbor(0)
  m%grid(ichild)%ckey(1:ndim)=hash_nbor(1:ndim)
  m%occupied(m%free_cache)=.true.
  m%parent_cpu(m%free_cache)=0
  m%dirty(m%free_cache)=.false.

  ! Store parent cell coordinates
  m%ghost_parent_grid(m%free_cache)=igrid
  m%ghost_parent_cell(m%free_cache)=icell

  ! Set refined to false
  m%grid(ichild)%refined(1:twotondim)=.false.

  ! Store children cell hydro variables
  do ivar=1,nvar
     do ind=1,twotondim
        m%uold(ind,ivar,ichild)=u2(ind,ivar)
     enddo
  end do
#ifdef MHD
  ! Store children cell MHD variables
  do ivar=1,6
     do ind=1,twotondim
        m%bold(ind,ivar,ichild)=b2(ind,ivar)
     enddo
  end do
#endif
#ifdef GRAV
  ! Store children cell gravity variables using straight injection
  do ind=1,twotondim
     m%f(ind,1:ndim,ichild)=m%f(icell,1:ndim,igrid)
     m%phi(ind,ichild)=m%phi(icell,igrid)
     m%phi_old(ind,ichild)=m%phi_old(icell,igrid)
  enddo
#endif

  ! Go to next free cache line
  m%free_cache=m%free_cache+1
  m%ncache=m%ncache+1
  if(m%free_cache.GT.m%ncachemax)then
     m%free_cache=1
  endif
  if(m%ncache.GT.m%ncachemax)m%ncache=m%ncachemax

  end associate

end subroutine make_grid_ghost
!###############################################################
!###############################################################
!###############################################################
!###############################################################
end module godunov_fine_module

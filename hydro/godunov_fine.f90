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
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  type(msg_large_realdp)::dummy_large_realdp
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::igrid

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                hilbert=m%domain,pack_size=storage_size(dummy_large_realdp)/32,&
                pack=pack_fetch_refine,unpack=unpack_fetch_refine,&
                init=init_flush_godunov, flush=pack_flush_godunov,&
                combine=unpack_flush_godunov, bound=init_bound_refine)

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

  call close_cache(s,m%grid_dict)

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
  use amr_parameters, only: ndim,twotondim,dp
  use amr_commons, only: run_t,global_t,mesh_t
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
  do i=m%head(ilevel),m%tail(ilevel)
     m%grid(i)%unew = m%grid(i)%uold
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
  use amr_parameters, only: dp,ndim,twotondim
  use amr_commons, only: run_t,global_t,mesh_t
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
  do i=m%head(ilevel),m%tail(ilevel)
     m%grid(i)%uold=m%grid(i)%unew
  end do
#endif

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(s,ind_grid,ilevel,h)
  use mdl_module
  use amr_parameters, only: ndim,twondim,twotondim,dp
  use amr_commons, only: nbor,oct
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
  integer::igrid,icell,inbor,ipass
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::ii0,jj0,kk0,ii1,jj1,kk1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::ii1min,ii1max,jj1min,jj1max,kk1min,kk1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer,dimension(1:ndim)::ckey_corner,ckey
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  integer,dimension(0:twondim)::ind_nbor
  type(nbor),dimension(0:twondim)::grid_nbor
  real(dp)::dx,oneontwotondim
  real(dp),dimension(0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  logical::okx,oky,okz
  type(oct),pointer::gridp,childp

  i2min=0; i2max=0; j2min=0; j2max=0; k2min=0; k2max=0
  i3min=1; i3max=1; j3min=1; j3max=1; k3min=1; k3max=1
  okx=.true.; oky=.true.; okz=.true.

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)
  
  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=r%boxlen/2**ilevel

  ! Integer constants
  i1min=h%io1; i1max=h%io2; j1min=h%jo1; j1max=h%jo2; k1min=h%ko1; k1max=h%ko2
#if NDIM>0
  i2max=1; i3min=h%iu1+2; i3max=h%iu2-2
#endif
#if NDIM>1
  j2max=1; j3min=h%ju1+2; j3max=h%ju2-2
#endif
#if NDIM>2
  k2max=1; k3min=h%ku1+2; k3max=h%ku2-2
#endif

  ! Reset gravitational acceleration
  h%gloc=0.0d0

  !---------------------
  ! Gather hydro stencil
  !---------------------
  hash_key(0)=m%grid(ind_grid)%lev
  hash_nbor(0)=m%grid(ind_grid)%lev
  ckey_corner(1:ndim)=(m%grid(ind_grid)%ckey(1:ndim)/(i1max-1))*(i1max-1)
  ind_oct=ind_grid

  ! Loop over 3x3x3 neighboring father cells using 27 passes
  do ipass = 1, threetondim

  if(ipass == 1)then
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min; jj1max = j1max
#if NDIM>1
     jj1min = jj1min+1; jj1max = jj1max-1
#endif
     kk1min = k1min; kk1max = k1max;
#if NDIM>2
     kk1min = kk1min+1; kk1max = kk1max-1
#endif
  endif
  if(MOD(ipass,3) == 2)then ! ipass = 2, 5, 8, 11, 14, 17, 20, 23, 26
     ii1min = i1min; ii1max = i1min
  endif
  if(MOD(ipass,3) == 0)then ! ipass = 3, 6, 9, 12, 15, 18, 21, 24, 27
     ii1min = i1max; ii1max = i1max
  endif
  if(MOD(ipass,9) == 4)then ! ipass = 4, 13, 22
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min; jj1max = j1min
  endif
  if(MOD(ipass,9) == 7)then ! ipass = 7, 16, 25
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1max; jj1max = j1max
  endif
  if(ipass == 10)then
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min+1; jj1max = j1max-1
     kk1min = k1min; kk1max = k1min
  endif
  if(ipass == 19)then
     ii1min = i1min+1; ii1max = i1max-1
     jj1min = j1min+1; jj1max = j1max-1
     kk1min = k1max; kk1max = k1max
  endif

  do k1=kk1min,kk1max
#if NDIM>2
     okz=(k1>k1min.and.k1<k1max)
#endif
     do j1=jj1min,jj1max
#if NDIM>1
        oky=(j1>j1min.and.j1<j1max)
#endif
        do i1=ii1min,ii1max
           okx=(i1>i1min.and.i1<i1max)
           ! For inner octs only
           if(okx.and.oky.and.okz)then

              ! Compute relative Cartesian key
              ckey(1:ndim)=m%grid(ind_oct)%ckey(1:ndim)-ckey_corner(1:ndim)

              ! Store grid index
              ii1=1; jj1=1; kk1=1
              ii1=ckey(1)+1
#if NDIM>1
              jj1=ckey(2)+1
#endif
#if NDIM>2
              kk1=ckey(3)+1
#endif
              h%childloc(ii1,jj1,kk1)%p=>m%grid(ind_oct)
              nullify(h%gridloc(ii1,jj1,kk1)%p)
              h%cellloc(ii1,jj1,kk1)=0

              ! Loop over 2x2x2 cells
              do k2=k2min,k2max
                 do j2=j2min,j2max
                    do i2=i2min,i2max
                       ind_son=1+i2+2*j2+4*k2
                       i3=1; j3=1; k3=1
                       i3=1+2*(ii1-1)+i2
#if NDIM>1
                       j3=1+2*(jj1-1)+j2
#endif
#if NDIM>2
                       k3=1+2*(kk1-1)+k2
#endif             
                       ! Gather hydro variables
                       do ivar=1,nvar
                          h%uloc(i3,j3,k3,ivar)=m%grid(ind_oct)%uold(ind_son,ivar)
                       end do
#ifdef GRAV
                       ! Gather gravitational acceleration
                       do idim=1,ndim
                          h%gloc(i3,j3,k3,idim)=m%grid(ind_oct)%f(ind_son,idim)
                       end do
#endif
                       ! Gather refinement flag
                       h%okloc(i3,j3,k3)=m%grid(ind_oct)%refined(ind_son)
                    end do
                 end do
              end do

              ! Go to next inner oct
              ind_oct=ind_oct+1

           ! For boundary octs only
           else
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
                    if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel))hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
                    if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel))hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
                 endif
              enddo
              
              ! Set pointers to null
              icell=0
              nullify(gridp)
              do inbor=0,twondim
                 nullify(grid_nbor(inbor)%p)
              end do

              ! Get neighboring grid index with read-only cache
              call get_grid(s,hash_nbor,m%grid_dict,childp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              if(.not.associated(childp))then

                 ! Get parent father cell with read-write cache
                 call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icell,flush_cache=.true.,fetch_cache=.true.,lock=.true.)
                 if(.not.associated(gridp))then
                    write(*,*)'GODUNOV: parent_cell should exist'
                    write(*,*)'PE ',g%myid,hash_nbor
                    call mdl_abort(mdl)
                 endif

                 ! In case one wants to interpolate using high-order schemes
                 if(r%interpol_type>0)then

                    ! Get 2ndim neighboring father cells with read-write cache
                    ! Note that possible cache grids are locked inside the routine
                    call get_twondim_nbor_parent_cell(s,hash_nbor,m%grid_dict,grid_nbor,ind_nbor,flush_cache=.true.,fetch_cache=.true.)
                    do inbor=0,twondim
                       do ivar=1,nvar
                          u1(inbor,ivar)=grid_nbor(inbor)%p%uold(ind_nbor(inbor),ivar)
                       end do
                    end do

                    ! Interpolate
                    call interpol_hydro(u1,u2,r%interpol_var,r%interpol_type,r%smallr)

                 endif

              endif

              ! Store grid index
              h%childloc(i1,j1,k1)%p=>childp
              h%gridloc (i1,j1,k1)%p=>gridp
              h%cellloc (i1,j1,k1)=icell
              if(r%interpol_type>0)then
                 do inbor=1,twondim
                    h%nborloc(i1,j1,k1,inbor)%p=>grid_nbor(inbor)%p
                 end do
              endif

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
                       ! If neighboring grid exists
                       if(associated(childp))then

                          ! Gather hydro variables
                          do ivar=1,nvar
                             h%uloc(i3,j3,k3,ivar)=childp%uold(ind_son,ivar)
                          end do

#ifdef GRAV
                          ! Gather gravitational acceleration
                          do idim=1,ndim
                             h%gloc(i3,j3,k3,idim)=childp%f(ind_son,idim)
                          end do
#else
                          ! Gather gravitational acceleration
                          do idim=1,ndim
                             h%gloc(i3,j3,k3,idim)=r%constant_gravity(idim)
                          end do
#endif
                          ! Gather refinement flag
                          h%okloc(i3,j3,k3)=childp%refined(ind_son)

                       ! If neighboring grid doesn't exist, interpolate
                       else

                          ! Gather hydro variables
                          do ivar=1,nvar
                             h%uloc(i3,j3,k3,ivar)=gridp%uold(icell,ivar)
                          end do

                          ! Gather interpolated hydro variables
                          if(r%interpol_type>0)then
                             do ivar=1,nvar
                                h%uloc(i3,j3,k3,ivar)=u2(ind_son,ivar)
                             end do
                          endif

#ifdef GRAV
                          ! Gather gravitational acceleration
                          do idim=1,ndim
                             h%gloc(i3,j3,k3,idim)=gridp%f(icell,idim)
                          end do
#else
                          ! Gather gravitational acceleration
                          do idim=1,ndim
                             h%gloc(i3,j3,k3,idim)=r%constant_gravity(idim)
                          end do
#endif
                          ! Gather refinement flag
                          h%okloc(i3,j3,k3)=.false.
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
  end do

  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(h%uloc,h%gloc,h%qloc,h%cloc,&
       & h%flux,h%tmp,h%dq,h%qm,h%qp,h%fx,h%tx,h%divu,&
       & dx,dx,dx,g%dtnew(ilevel),&
       & h%iu1,h%iu2,h%ju1,h%ju2,h%ku1,h%ku2,&
       & h%if1,h%if2,h%jf1,h%jf2,h%kf1,h%kf2,&
       & r%gamma,r%gamma_rad,r%smallr,r%smallc,r%slope_type,r%riemann,r%difmag)
  
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
                 if(h%okloc(i3-i0,j3-j0,k3-k0) .or. h%okloc(i3,j3,k3))then
                    h%flux(i3,j3,k3,ivar,idim)=0.0d0
                 end if
              end do
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
              childp=>h%childloc(i1,j1,k1)%p
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
                       do ivar=1,nvar
                          childp%unew(ind_son,ivar)=childp%unew(ind_son,ivar)+ &
                               & (h%flux(i3   ,j3   ,k3   ,ivar,idim) &
                               & -h%flux(i3+i0,j3+j0,k3+k0,ivar,idim))
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  ! If sitting in coarsest level, skip.
  if(ilevel>r%levelmin)then

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
              ! Get oct index
              childp=>h%childloc(i1,j1,k1)%p
              ! Check that parent cell is not refined
              if(.not.associated(childp))then
                 ! Get parent cell index
                 gridp=>h%gridloc(i1,j1,k1)%p
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
                          do ivar=1,nvar
                             gridp%unew(icell,ivar)=gridp%unew(icell,ivar) &
                                  & -h%flux(i3,j3,k3,ivar,idim)*oneontwotondim
                          end do
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
              childp=>h%childloc(i1,j1,k1)%p
              ! Check that parent cell is not refined
              if(.not.associated(childp))then
                 ! Get parent cell index
                 gridp=>h%gridloc(i1,j1,k1)%p
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
                          do ivar=1,nvar
                             gridp%unew(icell,ivar)=gridp%unew(icell,ivar) &
                                  & +h%flux(i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                          end do
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
           childp=>h%childloc(i1,j1,k1)%p
           ! Check that parent cell is not refined
           if(associated(childp))then
              call unlock_cache(s,childp)
           else
              ! Get parent cell index
              gridp=>h%gridloc(i1,j1,k1)%p
              call unlock_cache(s,gridp)
              ! Get neighbouring parent oct index
              if(r%interpol_type>0)then
                 do inbor=1,twondim
                    gridp=>h%nborloc(i1,j1,k1,inbor)%p
                    call unlock_cache(s,gridp)
                 end do
              endif
           endif
        end do
     end do
  end do

  end associate

#endif

end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_godunov(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%unew(ind,ivar)=0.0d0
     enddo
  enddo
#endif

end subroutine init_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_godunov(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp_hydro(ind,ivar)=grid%unew(ind,ivar)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_godunov(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_large_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%unew(ind,ivar)=grid%unew(ind,ivar)+msg%realdp_hydro(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module godunov_fine_module

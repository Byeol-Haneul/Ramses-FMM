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
  use amr_commons
  use cache_commons
  use cache
  use marshal, only: pack_fetch_refine,unpack_fetch_refine
  use boundaries, only: init_bound_refine
  use gpu_runner, only: gpu_integrator
  use nvtx
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

  if(s%r%verbose.and.s%g%myid==1)write(*,'("   Entering godunov_fine for level ",I2)')ilevel

  call nvtxStartRange("GPU integrator", color=6)!teal
  call gpu_integrator(s, ilevel)
  call nvtxEndRange()
 
  ! Loop over active grids by vector sweeps
!   call nvtxStartRange("CPU Integrator", color=9)!white
!   do igrid=m%head(ilevel),m%tail(ilevel)
!      call godfine1(s,igrid,ilevel,m%hydro_w%kernel_1)
!   end do
!   call nvtxEndRange()

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
!$OMP PARALLEL DO
  do i = m%head(ilevel),m%tail(ilevel)
     m%grid(i)%unew = m%grid(i)%uold
  end do
!$OMP END PARALLEL DO
#endif
#ifdef MHD
  ! Set bnew to bold for myid cells
!$OMP PARALLEL DO
  do i = m%head(ilevel),m%tail(ilevel)
     m%grid(i)%bnew = m%grid(i)%bold
  end do
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO
  do i = m%head(ilevel),m%tail(ilevel)
     m%grid(i)%uold = m%grid(i)%unew
  end do
!$OMP END PARALLEL DO
#endif
#ifdef MHD
  ! Set bold to bnew
!$OMP PARALLEL DO
  do i = m%head(ilevel),m%tail(ilevel)
     m%grid(i)%bold = m%grid(i)%bnew
  end do
!$OMP END PARALLEL DO
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
  integer::ivar,idim,ind,ind_son,ind_oct
  integer::igrid,icell,inbor,ipass
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
  integer::ii0,jj0,kk0,ii1,jj1,kk1
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::ii1min,ii1max,jj1min,jj1max,kk1min,kk1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer,dimension(1:ndim)::ckey_corner,ckey
  integer(kind=8),dimension(0:ndim)::hash_nbor,hash_son_nbor
  integer,dimension(0:twondim)::ind_nbor
  type(nbor),dimension(0:twondim)::grid_nbor
  real(dp)::dx,oneontwotondim
  real(dp),dimension(0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:twotondim,1:nvar)::u2
  logical::okx,oky,okz,oknbor
  logical::ok1,ok2,ok3
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
  hash_nbor(0)=m%grid(ind_grid)%lev
  ckey_corner(1:ndim)=(m%grid(ind_grid)%ckey(1:ndim)/(i1max-1))*(i1max-1)
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
           okx=(i1>i1min.and.i1<i1max)

           !--------------------
           ! For inner octs only
           !--------------------
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
                       ! Gather self-gravitational acceleration
                       do idim=1,ndim
                          h%gloc(i3,j3,k3,idim)=m%grid(ind_oct)%f(ind_son,idim)
                       end do
#else
                       ! Gather constant gravitational acceleration
                       do idim=1,ndim
                          h%gloc(i3,j3,k3,idim)=r%constant_gravity(idim)
                       end do
#endif
                       ! Gather refinement flag
                       h%okloc(i3,j3,k3)=m%grid(ind_oct)%refined(ind_son)
                    end do
                 end do
              end do

              ! Go to next inner oct
              ind_oct=ind_oct+1

           !-----------------------
           ! For boundary octs only
           !-----------------------
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

              ! Store grid index
              h%childloc(i1,j1,k1)%p=>childp
              h%gridloc (i1,j1,k1)%p=>gridp
              h%cellloc (i1,j1,k1)=icell
              do inbor=1,twondim
                 h%nborloc(i1,j1,k1,inbor)%p=>grid_nbor(inbor)%p
              end do

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
                          h%uloc(i3,j3,k3,ivar)=childp%uold(ind_son,ivar)
                       end do
#ifdef GRAV
                       ! Gather self-gravitational acceleration
                       do idim=1,ndim
                          h%gloc(i3,j3,k3,idim)=childp%f(ind_son,idim)
                       end do
#else
                       ! Gather constant gravitational acceleration
                       do idim=1,ndim
                          h%gloc(i3,j3,k3,idim)=r%constant_gravity(idim)
                       end do
#endif
                       ! Gather refinement flag
                       h%okloc(i3,j3,k3)=childp%refined(ind_son)

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
  call unsplit(h%uloc,h%gloc,h%qloc,h%cloc,&
       & h%flux,h%tmp,h%dq,h%qm,h%qp,h%fx,h%tx,h%divu,&
       & dx,g%dtnew(ilevel),&
       & h%iu1,h%iu2,h%ju1,h%ju2,h%ku1,h%ku2,&
       & h%if1,h%if2,h%jf1,h%jf2,h%kf1,h%kf2,&
       & s%h_params)

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
                       do ivar=1,5
                          childp%unew(ind_son,ivar)=childp%unew(ind_son,ivar)+ &
                               & (h%flux(i3   ,j3   ,k3   ,ivar,idim) &
                               & -h%flux(i3+i0,j3+j0,k3+k0,ivar,idim))
                       end do
#if NVAR>5
                       do ivar=6,nvar
                          childp%unew(ind_son,ivar)=childp%unew(ind_son,ivar)+ &
                               & (h%flux(i3   ,j3   ,k3   ,ivar+ie-5,idim) &
                               & -h%flux(i3+i0,j3+j0,k3+k0,ivar+ie-5,idim))
                       end do
#endif
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  !----------------
  ! Unlock all octs
  !----------------
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max     
           ! Get oct index
           childp=>h%childloc(i1,j1,k1)%p
           call unlock_cache(s,childp)
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

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
#ifdef HYDRO
  grid%unew=0.0d0
#endif

#ifdef MHD
  grid%bnew=0.0d0
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

  type(msg_large_realdp)::msg

#ifdef HYDRO
  msg%realdp_hydro=grid%unew
#endif

#ifdef MHD
  msg%realdp_mhd=grid%bnew
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

  type(msg_large_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef HYDRO
  grid%unew=grid%unew+msg%realdp_hydro
#endif

#ifdef MHD
  grid%bnew=grid%bnew+msg%realdp_mhd
#endif

end subroutine unpack_flush_godunov
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module godunov_fine_module

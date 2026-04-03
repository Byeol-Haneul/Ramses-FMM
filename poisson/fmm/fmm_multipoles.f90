module fmm_multipoles
contains
#ifdef GRAV
!###############################################
!###############################################
!###############################################
!###############################################
subroutine m_fmm_multipoles(pst,ilevel)
  use ramses_commons, only: pst_t
  use init_fmm_module, only: fmm_level_t, FMM_MULTIPOLE_STANDARD
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !------------------------------------------------------------------
  ! This master routine computes the mass density field to be used
  ! as source term in the Poisson solver.
  ! The density field is computed for all levels greater than ilevel.
  ! On output, particles are sorted according to their grid level of
  ! refinement, and inside their level, they are sorted according to
  ! their grid Hilbert order.
  !------------------------------------------------------------------
  type(fmm_level_t)::fmm_levels
  integer::i,input_size
  associate(r=>pst%s%r)

  if(.not. r%poisson)return
  if(r%verbose)write(*,'(" Entering fmm_multipoles for level ",I2)')ilevel

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------

  fmm_levels%ilev = ilevel
  fmm_levels%mode = FMM_MULTIPOLE_STANDARD
  input_size = storage_size(fmm_levels)/32

  ! Initialize both AMR and FMM grids. 
  do i = r%bound_levelmin, r%nlevelmax, 1
    fmm_levels%flev=i
    call r_reset_multipoles_taylor(pst, fmm_levels, input_size)
  end do

  ! Add multipoles from AMR grids
  if(r%verbose) print *, "[P2M] LEVEL: ", ilevel
  call r_fmm_multipole_amr2fmm(pst, ilevel, 1)
  if(r%verbose) print *, "      <AMR->FMM> : ", ilevel

  if(r%verbose) print *, "[M2M] LEVEL: ", ilevel
  ! Add multipoles to FMM grids. 
  do i=ilevel-r%level_fmm_to_amr-1,r%bound_levelmin,-1
     fmm_levels%flev=i
     if(r%verbose)write(*,'("     <ACCUMULATION> TREE for AMR LEVEL: ",I2,", TREE LEVEL: ",I2)')ilevel, i
     call r_fmm_multipole_fmm2fmm(pst,fmm_levels,input_size)
  end do

  call update_fmm_local_multipole_level(pst, ilevel)

  end associate

end subroutine m_fmm_multipoles
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_accumulate_fmm_global_multipole(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer,VALUE::input_size
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_ACCUM_GLOBAL_MULTIPOLE_FMM,pst%iUpper+1,input_size,0,ilevel)
     call r_accumulate_fmm_global_multipole(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call accumulate_fmm_global_multipole(pst%s,ilevel)
  endif

end subroutine r_accumulate_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accumulate_fmm_global_multipole(s,ilevel)
  use amr_parameters, only: twotondim, multipole_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t) :: s
  integer, intent(in) :: ilevel

  integer :: ioct, icell
  type(mesh_t), pointer :: m_fmm

  associate(r=>s%r, g=>s%g)
  m_fmm => s%m_fmm_list(ilevel)
#ifdef FMM
  g%multipole%q(1:multipole_size) = 0.0d0
#endif
  if (m_fmm%tail(r%bound_levelmin) < m_fmm%head(r%bound_levelmin)) return

#ifdef FMM
  do ioct=m_fmm%head(r%bound_levelmin),m_fmm%tail(r%bound_levelmin)
     do icell=1,twotondim
        g%multipole%q(1:multipole_size) = g%multipole%q(1:multipole_size) + &
             & m_fmm%multipole(icell,1:multipole_size,ioct)
     end do
  end do
#endif
  end associate
end subroutine accumulate_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine update_fmm_local_multipole_level(pst,ilevel)
  use amr_parameters, only: multipole_size
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t) :: pst
  integer, intent(in) :: ilevel

  associate(g=>pst%s%g)
  call r_accumulate_fmm_global_multipole(pst, ilevel, 1)
  g%multipole_fmm_level(ilevel)%q(1:multipole_size) = g%multipole%q(1:multipole_size)
  end associate
end subroutine update_fmm_local_multipole_level
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sync_fmm_global_multipole(pst)
  use amr_parameters, only: multipole_size
  use amr_commons, only: multipole_t
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t) :: pst

  type(multipole_t) :: multipole_tot
  integer :: input_size, ilevel

  associate(g=>pst%s%g, r=>pst%s%r)
  multipole_tot%q(1:multipole_size) = 0.0d0
  input_size = storage_size(multipole_tot)/32

  do ilevel = r%levelmin, r%nlevelmax
     multipole_tot%q(1:multipole_size) = multipole_tot%q(1:multipole_size) + &
          & g%multipole_fmm_level(ilevel)%q(1:multipole_size)
  end do

  g%multipole_fmm_raw%q(1:multipole_size) = multipole_tot%q(1:multipole_size)
  multipole_tot%q(1:multipole_size) = 0.0d0
  call r_collect_fmm_global_multipole(pst, r%levelmin, input_size, multipole_tot, input_size)
  g%multipole_fmm_raw%q(1:multipole_size) = multipole_tot%q(1:multipole_size)
  call center_fmm_global_multipole(multipole_tot)

  call r_broadcast_fmm_global_multipole(pst, multipole_tot, input_size)
  end associate
end subroutine sync_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine center_fmm_global_multipole(multipole)
  use amr_parameters, only: ndim
  use amr_commons, only: multipole_t
  implicit none
  type(multipole_t), intent(inout) :: multipole

  real(kind=8) :: mass
  real(kind=8), dimension(3) :: center

  mass = multipole%q(1)
  if (mass <= 0.0d0) then
     multipole%q = 0.0d0
     return
  end if

  center = 0.0d0
  center(1:ndim) = multipole%q(2:ndim+1)/mass
#ifdef FMM
#if NDIM==1
  multipole%q(3) = multipole%q(3) - mass*center(1)*center(1)
#endif
#if NDIM==2
  multipole%q(4) = multipole%q(4) - mass*center(1)*center(1)
  multipole%q(5) = multipole%q(5) - mass*center(1)*center(2)
  multipole%q(6) = multipole%q(6) - mass*center(2)*center(2)
#endif
#if NDIM==3
  multipole%q(5)  = multipole%q(5)  - mass*center(1)*center(1)
  multipole%q(6)  = multipole%q(6)  - mass*center(1)*center(2)
  multipole%q(7)  = multipole%q(7)  - mass*center(1)*center(3)
  multipole%q(8)  = multipole%q(8)  - mass*center(2)*center(2)
  multipole%q(9)  = multipole%q(9)  - mass*center(2)*center(3)
  multipole%q(10) = multipole%q(10) - mass*center(3)*center(3)
#endif
#endif
end subroutine center_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_collect_fmm_global_multipole(pst,ilevel,input_size,multipole,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer::ilevel
  type(multipole_t)::multipole,next_multipole

  integer::rID

  multipole%q = 0.0d0
  next_multipole%q = 0.0d0

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COLLECT_GLOBAL_MULTIPOLE_FMM,pst%iUpper+1,input_size,output_size,ilevel)
     call r_collect_fmm_global_multipole(pst%pLower,ilevel,input_size,multipole,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_multipole)
     multipole%q = multipole%q + next_multipole%q
  else
     multipole%q = pst%s%g%multipole_fmm_raw%q
  endif
end subroutine r_collect_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_broadcast_fmm_global_multipole(pst,multipole,input_size)
  use mdl_module
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(multipole_t)::multipole

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_GLOBAL_MULTIPOLE_FMM,pst%iUpper+1,input_size,0,multipole)
     call r_broadcast_fmm_global_multipole(pst%pLower,multipole,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%multipole = multipole
     pst%s%g%rho_tot = pst%s%g%multipole%q(1)/PRODUCT(pst%s%r%box_size(1:ndim))
  endif
end subroutine r_broadcast_fmm_global_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_multipoles(pst,active_levelmin)
  use ramses_commons, only: pst_t
  use init_fmm_module, only: fmm_level_t, FMM_MULTIPOLE_MERGED
  implicit none
  type(pst_t)::pst
  integer, intent(in) :: active_levelmin
  type(fmm_level_t)::fmm_levels
  integer::input_size

  associate(r=>pst%s%r)
  if(.not. r%poisson)return
  if(r%verbose) print *, "[MERGE MULTIPOLES] Building merged multipoles"

  fmm_levels%ilev = active_levelmin
  fmm_levels%flev = r%bound_levelmin
  fmm_levels%mode = FMM_MULTIPOLE_MERGED
  input_size = storage_size(fmm_levels)/32
  call r_fmm_multipole_fmm2fmm(pst,fmm_levels,input_size)
  end associate
end subroutine merge_multipoles
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_multipole_amr2fmm(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer,VALUE::input_size
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_AMR2FMM,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_multipole_amr2fmm(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_multipole_amr2fmm(pst%s,ilevel)
  endif

end subroutine r_fmm_multipole_amr2fmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_multipole_amr2fmm(s,ilevel)
  use amr_parameters, only: ndim, twotondim, multipole_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::ind,idim,ioct,icell,nstride,igrid_fmm
  integer(kind=8),dimension(0:ndim)::hash_key_amr, hash_key_fmm
  integer(kind=8),dimension(1:ndim)::ii
  logical::leaf_cell
  type(msg_large_realdp)::dummy_realdp

  integer :: nq
  real(kind=8), dimension(1:multipole_size) :: multipole
  real(kind=8), dimension(ndim) :: xx
  real(kind=8) :: dx_loc, vol_loc, mmm

  ! Multipole arrays (static)
  real(kind=8) :: monopole
  real(kind=8), dimension(1:ndim) :: dipole
  real(kind=8), dimension(1:int(ndim*(ndim+1)/2)) :: quadrupole

  associate(r=>s%r,m=>s%m,mdl=>s%mdl)

  !---------------------------------------------------
  ! Initialize constants
  !---------------------------------------------------
  nq = int(ndim*(ndim+1)/2)

  ! Mesh spacing for this level
  dx_loc = r%boxlen / 2.0D0**ilevel
  vol_loc = dx_loc**ndim

  call open_cache(mdl,s%m_fmm_list(ilevel),pack_size=storage_size(dummy_realdp)/32,&
                     init=init_flush_multipole, flush=pack_flush_multipole, combine=unpack_flush_multipole)

  ! Loop over levelmin grids.
  hash_key_fmm(0)=ilevel - r%level_fmm_to_amr
  hash_key_amr(0) = ilevel
  do ioct=m%head(ilevel),m%tail(ilevel)
     ! Get fmm grid above level_fmm_to_amr
     hash_key_amr(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     hash_key_fmm(1:ndim)= hash_key_amr(1:ndim)/(2**r%level_fmm_to_amr)
     ii(1:ndim)=hash_key_amr(1:ndim)-(2**r%level_fmm_to_amr)*hash_key_fmm(1:ndim) ! 0 to 2^(level_fmm_to_amr)-1
     ii(1:ndim)=ii(1:ndim)/(2**(r%level_fmm_to_amr-1)) ! 0 or 1
     icell=1
     do idim=1,ndim
       icell=icell+2**(idim-1)*ii(idim) ! 1 to twotondim
     end do
     ! Get fmm grid using a write-only cache
     call get_grid(s,hash_key_fmm,igrid_fmm,flush_cache=.true.,fetch_cache=.false.)

     if (igrid_fmm <= 0) cycle

     multipole = 0.0D0

    ! Loop over cells
     do ind = 1, twotondim
        leaf_cell=m%grid(ioct)%refined(ind).EQV..FALSE.
        ! Reset multipoles for this grid
        monopole   = 0.0D0
        dipole     = 0.0D0
        quadrupole = 0.0D0

        if (leaf_cell) then
          ! Compute cell center coordinates
           do idim = 1, ndim
              nstride = 2**(idim-1)
              xx(idim) = (2*m%grid(ioct)%ckey(idim) + MOD((ind-1)/nstride, 2) + 0.5D0) * dx_loc - m%skip(idim)
           end do

           ! Gas mass contribution
           mmm = m%rho(ind,ioct) * vol_loc
           monopole = monopole + mmm
           dipole   = dipole   + mmm * xx

           ! Quadrupole contribution
#if NDIM==1
           quadrupole(1) = quadrupole(1) + mmm * xx(1)**2       ! quadrupole_xx
#endif
#if NDIM==2
           quadrupole(1) = quadrupole(1) + mmm * xx(1)**2       ! quadrupole_xx
           quadrupole(2) = quadrupole(3) + mmm * xx(1)*xx(2)    ! quadrupole_xy
           quadrupole(3) = quadrupole(3) + mmm * xx(2)**2       ! quadrupole_yy
#endif
#if NDIM==3
           quadrupole(1) = quadrupole(1) + mmm * xx(1)**2        ! quadrupole_xx
           quadrupole(2) = quadrupole(2) + mmm * xx(1)*xx(2)     ! quadrupole_xy
           quadrupole(3) = quadrupole(3) + mmm * xx(1)*xx(3)     ! quadrupole_xz
           quadrupole(4) = quadrupole(4) + mmm * xx(2)**2        ! quadrupole_yy
           quadrupole(5) = quadrupole(5) + mmm * xx(2)*xx(3)     ! quadrupole_yz
           quadrupole(6) = quadrupole(6) + mmm * xx(3)**2        ! quadrupole_zz
#endif
        end if
        multipole(1) = multipole(1) + monopole
        multipole(2:1+ndim) = multipole(2:1+ndim) + dipole
        multipole(2+ndim:1+ndim+nq) = multipole(2+ndim:1+ndim+nq) + quadrupole
     end do  ! cell loop
#ifdef FMM
      ! need to fix so that we loop over all ilevel above
     s%m_fmm_list(ilevel)%multipole(icell,:,igrid_fmm) = s%m_fmm_list(ilevel)%multipole(icell,:,igrid_fmm) + multipole
#endif
  end do
  call close_cache(mdl)
  end associate
end subroutine fmm_multipole_amr2fmm
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_multipole_fmm2fmm(pst,fmm_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_fmm_module, only: fmm_level_t, FMM_MULTIPOLE_MERGED
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  type(fmm_level_t)::fmm_levels
  integer,VALUE::input_size

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_FMM2FMM,pst%iUpper+1,input_size,0,fmm_levels)
     call r_fmm_multipole_fmm2fmm(pst%pLower,fmm_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if (fmm_levels%mode == FMM_MULTIPOLE_MERGED) then
        call fmm_merge_multipoles_all(pst%s,fmm_levels%ilev)
     else
        call fmm_multipole_fmm2fmm(pst%s,pst%s%m_fmm_list(fmm_levels%ilev),fmm_levels%flev)
     end if
  endif

end subroutine r_fmm_multipole_fmm2fmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_multipole_fmm2fmm(s,m_fmm,flev)
  use amr_parameters, only: ndim, twotondim, multipole_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  integer::flev
  !-------------------------------------------------------------------
  ! This routine compute the monopole and dipole of the gas mass and
  ! the analytical profile (if any) within each cell.
  ! For pure particle runs, this is not necessary and the
  ! routine is not even called.
  !-------------------------------------------------------------------
  integer::ind,ioct,icell,igrid
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_large_realdp)::dummy_realdp

  real(kind=8), dimension(1:multipole_size) :: multipole
  type(mesh_t)::m_fmm

  associate(mdl=>s%mdl)
  
  call open_cache(mdl,m_fmm,pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_multipole,unpack=unpack_fetch_multipole,&
                     init=init_flush_multipole, flush=pack_flush_multipole, combine=unpack_flush_multipole)

  ! Loop over finer level grids
  hash_key(0)=flev+1
  do ioct=m_fmm%head(flev+1),m_fmm%tail(flev+1)
     hash_key(1:ndim)=m_fmm%grid(ioct)%ckey(1:ndim)
     ! Get parent cell using a write-only cache
     call get_parent_cell(s,hash_key,igrid,icell,flush_cache=.true.,fetch_cache=.false.)
     if (igrid <= 0) cycle
     multipole = 0.0D0
#ifdef FMM
     do ind=1,twotondim
       multipole = multipole + m_fmm%multipole(ind,:,ioct)
     end do
     m_fmm%multipole(icell,:,igrid) = m_fmm%multipole(icell,:,igrid) + multipole

#endif
  end do
  call close_cache(mdl)
  end associate
end subroutine fmm_multipole_fmm2fmm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fmm_merge_multipoles_all(s,active_levelmin)
  use amr_parameters, only: ndim, twotondim, multipole_size, taylor_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  use hash
  implicit none

  type(ramses_t)::s
  integer, intent(in) :: active_levelmin
  type(mesh_t), pointer :: m_merged
  integer::ilevel, flev, ioct, igrid_merged
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_large_realdp)::dummy_realdp

  associate(r=>s%r, mdl=>s%mdl)
  if (.not. associated(s%m_fmm_merged)) return

  m_merged => s%m_fmm_merged

  do flev=r%bound_levelmin,r%nlevelmax
     if (m_merged%tail(flev) < m_merged%head(flev)) cycle
     do ioct=m_merged%head(flev),m_merged%tail(flev)
#ifdef FMM
        m_merged%multipole(1:twotondim,1:multipole_size,ioct)=0.0D0
        m_merged%taylor_coeff(1:twotondim,1:taylor_size,ioct)=0.0D0
#endif
     end do
  end do

  call open_cache(mdl, m_merged, pack_size=storage_size(dummy_realdp)/32, &
       init=init_flush_multipole, flush=pack_flush_multipole, combine=unpack_flush_multipole)

  do ilevel=active_levelmin,r%nlevelmax
     do flev=r%bound_levelmin,ilevel-r%level_fmm_to_amr
        if (s%m_fmm_list(ilevel)%tail(flev) < s%m_fmm_list(ilevel)%head(flev)) cycle
        hash_key(0)=flev

        do ioct=s%m_fmm_list(ilevel)%head(flev),s%m_fmm_list(ilevel)%tail(flev)
           hash_key(1:ndim)=s%m_fmm_list(ilevel)%grid(ioct)%ckey(1:ndim)
           ! Use the write-only cache so off-rank merged grids receive all
           ! multipole contributions instead of silently dropping remote ones.
           call get_grid(s, hash_key, igrid_merged, flush_cache=.true., fetch_cache=.false.)
           if(igrid_merged<=0) cycle
#ifdef FMM
           m_merged%multipole(:,:,igrid_merged)=m_merged%multipole(:,:,igrid_merged)+s%m_fmm_list(ilevel)%multipole(:,:,ioct)
#endif
        end do
     end do
  end do

  call close_cache(mdl)

  if(r%verbose) print *, "      <MERGED MULTIPOLES> accumulation done"
  end associate
end subroutine fmm_merge_multipoles_all
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_multipole_shift_downward(pst,fmm_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_fmm_module, only: fmm_level_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(fmm_level_t)::fmm_levels
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_SHIFT_DOWNWARD,pst%iUpper+1,input_size,0,fmm_levels)
     call r_fmm_multipole_shift_downward(pst%pLower,fmm_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_multipole_shift_downward(pst%s,pst%s%m_fmm_list(fmm_levels%ilev),fmm_levels%flev)
  endif

end subroutine r_fmm_multipole_shift_downward
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_multipole_shift_downward(s,m_fmm,flev)
  use amr_parameters, only: ndim, twotondim, multipole_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  type(mesh_t)::m_fmm
  integer::flev
  integer::idim,ioct,icell, nstride
  real(kind=8)::mass
  integer(kind=8),dimension(0:ndim)::hash_key
  real(kind=8) :: dx_loc
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  integer(kind=8), dimension(ndim) :: cc_icell! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_icell ! box unit real coordinate

  associate(r=>s%r,m=>s%m)

  dx_loc = r%boxlen / 2.0D0**flev
  hash_key(0)=flev
  mass = 0.0D0
  do ioct=m_fmm%head(flev),m_fmm%tail(flev)
     hash_key(1:ndim)=m_fmm%grid(ioct)%ckey(1:ndim)
     multipole = 0.0D0
     do icell = 1, twotondim
      do idim = 1, ndim
        nstride = 2**(idim-1)
        cc_icell(idim) = 2*hash_key(idim) + MOD((icell-1)/nstride, 2)
        xx_icell(idim) = (cc_icell(idim) + 0.5d0) * dx_loc - m%skip(idim)
      end do
#ifdef FMM
      multipole = m_fmm%multipole(icell, :, ioct)
      call shift_multipole(multipole, xx_icell, multipole_shifted)
      m_fmm%multipole(icell, :, ioct) = multipole_shifted
#endif
    end do
  end do
  end associate
end subroutine fmm_multipole_shift_downward
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_multipole(mesh,igrid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  integer::igrid
  type(mesh_t)::mesh
  integer(kind=8),dimension(0:ndim)::hash_key
  integer :: ind

#ifdef FMM
  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  do ind=1,twotondim
     mesh%grid(igrid)%refined(ind)=.true.
  end do
  mesh%multipole(:,:,igrid)=0.0
#endif
end subroutine init_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_multipole(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim,multipole_size
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
#ifdef FMM
  do ivar=1,multipole_size
     do ind=1,twotondim
        msg%realdp_fmm_multipole(ind,ivar)=mesh%multipole(ind,ivar,igrid)
     end do
  end do
  msg_array=transfer(msg,msg_array)
#endif
end subroutine pack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_multipole(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_large_realdp)::msg
#ifdef FMM
  do ind=1,twotondim
     if(mesh%grid(igrid)%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  msg%realdp_fmm_multipole=mesh%multipole(:,:,igrid)
#endif
  msg_array=transfer(msg,msg_array)
end subroutine pack_fetch_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_multipole(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
#ifdef FMM
  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     endif
  end do
  mesh%multipole(:,:,igrid)=msg%realdp_fmm_multipole
#endif
end subroutine unpack_fetch_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_multipole(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim,multipole_size
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
#ifdef FMM  
  do ivar=1,multipole_size
     do ind=1,twotondim
        mesh%multipole(ind,ivar,igrid)=mesh%multipole(ind,ivar,igrid)+msg%realdp_fmm_multipole(ind,ivar)
     end do
  end do
#endif
end subroutine unpack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_rho(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4_small_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_int4_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=mesh%rho(ind, igrid)
  end do
#endif
  do ind=1,twotondim
     msg%flg(ind)=0
     if (mesh%grid(igrid)%refined(ind)) then
        msg%ref(ind)=1
     else
        msg%ref(ind)=0
     end if
  end do

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_rho
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_rho(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_int4_small_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_int4_small_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     mesh%rho(ind, igrid)=msg%realdp(ind)
  end do
#endif
  do ind=1,twotondim
     if (msg%ref(ind) == 1) then
        mesh%grid(igrid)%refined(ind)=.true.
     else
        mesh%grid(igrid)%refined(ind)=.false.
     end if
  end do

end subroutine unpack_fetch_rho
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_reset_multipoles_taylor(pst,fmm_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_fmm_module, only: fmm_level_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  type(fmm_level_t)::fmm_levels
  integer,VALUE::input_size
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_MULTIPOLES,pst%iUpper+1,input_size,0,fmm_levels)
     call r_reset_multipoles_taylor(pst%pLower,fmm_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if (fmm_levels%flev <= fmm_levels%ilev-pst%s%r%level_fmm_to_amr) then
        call reset_multipoles_taylor(pst%s%r,pst%s%g,pst%s%m_fmm_list(fmm_levels%ilev),fmm_levels%flev)
     else 
        return
     end if
  endif

end subroutine r_reset_multipoles_taylor
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reset_multipoles_taylor(r,g,m,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)   :: r
  type(global_t):: g
  type(mesh_t)  :: m
  integer       :: ilevel
  integer :: igrid
  integer :: first, last
  first = m%head(ilevel)
  last  = m%tail(ilevel)
#ifdef FMM
  do igrid = first, last
    m%multipole(:,:,igrid) = 0.0D0
    m%taylor_coeff(:,:,igrid) = 0.0D0
  end do
#endif
end subroutine reset_multipoles_taylor
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine shift_multipole(multipole_in, a, multipole_out)
  use amr_parameters, only: ndim, multipole_size
  implicit none

  real(kind=8), intent(in)  :: multipole_in(1:multipole_size)
  real(kind=8), intent(in)  :: a(1:ndim)
  real(kind=8), intent(out) :: multipole_out(1:multipole_size)

  ! Locals
  real(kind=8) :: mp_in, mp_out
  real(kind=8) :: dp_in(3), dp_out(3)
  real(kind=8) :: qp_in(3,3), qp_out(3,3)
  integer :: i, j

  !----------------------------------------
  ! Monopole
  !----------------------------------------
  mp_in  = multipole_in(1)
  mp_out = mp_in
  multipole_out(1) = mp_out

  !----------------------------------------
  ! Dipole (extract from multipole_in)
  !----------------------------------------
  dp_in(1:ndim) = multipole_in(2:1+ndim)

  ! Shifted dipole: dp_out = dp_in - a*mp_in
  do i = 1, ndim
     dp_out(i) = dp_in(i) - a(i)*mp_in
     multipole_out(1+i) = dp_out(i)
  end do

  !----------------------------------------
  ! Quadrupole (extract from multipole_in, packed order)
  ! Order: xx, xy, xz, yy, yz, zz
  !----------------------------------------
  qp_in = 0.0d0
#if NDIM==1
  qp_in(1,1) = multipole_in(3)   ! xx
#endif
#if NDIM==2
  qp_in(1,1) = multipole_in(4)   ! xx
  qp_in(1,2) = multipole_in(5); qp_in(2,1) = qp_in(1,2) ! xy
  qp_in(2,2) = multipole_in(6)   ! yy
#endif
#if NDIM==3
  qp_in(1,1) = multipole_in(5)   ! xx
  qp_in(1,2) = multipole_in(6); qp_in(2,1) = qp_in(1,2) ! xy
  qp_in(1,3) = multipole_in(7); qp_in(3,1) = qp_in(1,3) ! xz
  qp_in(2,2) = multipole_in(8)   ! yy
  qp_in(2,3) = multipole_in(9); qp_in(3,2) = qp_in(2,3) ! yz
  qp_in(3,3) = multipole_in(10)  ! zz
#endif

  ! Shifted quadrupole:
  ! qp_out = qp_in - (a dp_in^T + dp_in a^T) + mp_in * (a a^T)
  qp_out = 0.0d0
  do i = 1, ndim
     do j = 1, ndim
        qp_out(i,j) = qp_in(i,j) - (a(i)*dp_in(j) + dp_in(i)*a(j)) + mp_in*a(i)*a(j)
     end do
  end do

  !----------------------------------------
  ! Pack quadrupole back into multipole_out
  !----------------------------------------
#if NDIM==1
  multipole_out(3) = qp_out(1,1)
#endif
#if NDIM==2
  multipole_out(4) = qp_out(1,1)   ! xx
  multipole_out(5) = qp_out(1,2)   ! xy
  multipole_out(6) = qp_out(2,2)   ! yy
#endif
#if NDIM==3
  multipole_out(5)  = qp_out(1,1)  ! xx
  multipole_out(6)  = qp_out(1,2)  ! xy
  multipole_out(7)  = qp_out(1,3)  ! xz
  multipole_out(8)  = qp_out(2,2)  ! yy
  multipole_out(9)  = qp_out(2,3)  ! yz
  multipole_out(10) = qp_out(3,3)  ! zz
#endif
end subroutine shift_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_dump_multipole(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_DUMP_MULTIPOLE,pst%iUpper+1,input_size,0,ilevel)
     call r_dump_multipole(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call dump_multipole(pst%s%r,pst%s%g,pst%s%m_fmm_list(ilevel),ilevel)
  endif

end subroutine r_dump_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dump_multipole(r, g, m, ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: oct, run_t, global_t, mesh_t
  implicit none
  type(run_t)    :: r
  type(global_t) :: g
  type(mesh_t)   :: m
  integer, intent(in)      :: ilevel

  integer :: ioct, icell, idim, nstride
  integer :: unit_debug
  integer(kind=8), dimension(ndim) :: cc_icell
  real(kind=8), dimension(ndim) :: xx_icell
  real(kind=8) :: dx_loc
  character(len=50) :: filename
  ! open debug file
  unit_debug = 99
#ifdef FMM
    write(filename, '(A,I0,A,I0,A)') "out/mult_fmm_", ilevel, "mpi_", g%myid,".out"
#else
    write(filename, '(A,I0,A)') "out/mult_fmm_", ilevel, ".out"
#endif
  open(unit_debug, file=filename, status="replace")
  dx_loc = r%boxlen / 2.0D0**ilevel
  do ioct = m%head(ilevel), m%tail(ilevel)
     do icell = 1, twotondim
        do idim = 1, ndim
          nstride = 2**(idim-1)
          cc_icell(idim) = 2*m%grid(ioct)%ckey(idim) + MOD((icell-1)/nstride, 2)
          xx_icell(idim) = (cc_icell(idim) + 0.5D0) * dx_loc - m%skip(idim)
        end do
#ifdef FMM
        write(unit_debug, '(3I6, E20.4)') cc_icell, m%multipole(icell, 1, ioct)
#endif
     end do
  end do
  close(unit_debug)
end subroutine dump_multipole
#endif
end module fmm_multipoles

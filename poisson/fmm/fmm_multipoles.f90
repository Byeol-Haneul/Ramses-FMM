module fmm_multipoles
contains
#ifdef GRAV
!###############################################
!###############################################
!###############################################
!###############################################
subroutine m_fmm_multipoles(pst,ilevel)
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use amr_commons, only: multipole_t
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
  type(multipole_t)::multipole_tot
  integer::i,input_size
  integer,dimension(1:2)::input_array
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  if(.not. r%poisson)return
  if(r%verbose)write(*,'(" Entering fmm_multipoles for level ",I2)')ilevel

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------

  ! Initialize both AMR and FMM grids. 
  do i = r%bound_levelmin, r%nlevelmax, 1
    call r_reset_multipoles_taylor(pst, i, 1)
  end do

  ! Add multipoles from AMR grids
  do i=r%nlevelmax,r%levelmin,-1
      if(r%verbose)write(*,'(" [M2M] AMR to FMM AMR LEVEL", I2)')i
      call r_fmm_multipole_amr2fmm(pst, i, 1)
  end do

  ! Add multipoles to FMM grids. 
  do i=r%levelmin-r%level_fmm_to_amr-1,r%bound_levelmin,-1
     if(i<1) cycle
     if(r%verbose)write(*,'(" [M2M] Compute multipoles for FMM level ",I2)')i
     call r_fmm_multipole_fmm2fmm(pst,i,1)
  end do

  do i=r%bound_levelmin,r%levelmin-r%level_fmm_to_amr
    if(r%verbose)write(*,'(" [M2M] Downward pass shifting multipoles for FMM level ",I2)')i
    call r_fmm_multipole_shift_downward(pst,i,1)
 end do

  !do i=r%levelmin-r%level_fmm_to_amr,r%bound_levelmin,-1
  !  write(*,'(" [M2M] DUMPING FOR MULT ",I2)')i
  !  call dump_multipole(r, g, m, i)
  !end do
  end associate

end subroutine m_fmm_multipoles
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
  integer,VALUE::input_size
  integer::ilevel

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
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_flag_module, only: pack_fetch_hydro, unpack_fetch_hydro
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
  integer::ind,idim,ivar,ioct,icell,nstride
  real(kind=8)::average
  integer(kind=8),dimension(0:ndim)::hash_key_amr, hash_key_fmm
  integer(kind=8),dimension(1:ndim)::ii
  logical::leaf_cell
  type(oct),pointer::grid_fmm
  type(msg_large_realdp)::dummy_realdp

  integer :: nm, nd, nq
  real(kind=8), dimension(1:multipole_size) :: multipole
  real(kind=8), dimension(ndim) :: xx
  real(kind=8) :: dx_loc, vol_loc, mmm, dd

  ! Multipole arrays (static)
  real(kind=8) :: monopole
  real(kind=8), dimension(1:ndim) :: dipole
  real(kind=8), dimension(1:int(ndim*(ndim+1)/2)) :: quadrupole

  associate(r=>s%r,g=>s%g,m=>s%m)

  !---------------------------------------------------
  ! Initialize constants
  !---------------------------------------------------
  nm = 1
  nd = ndim
  nq = int(ndim*(ndim+1)/2)

  ! Mesh spacing for this level
  dx_loc = r%boxlen / 2.0D0**ilevel
  vol_loc = dx_loc**ndim

  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
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
     call get_grid(s,hash_key_fmm,m%mg_dict,grid_fmm,flush_cache=.true.,fetch_cache=.false.)
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
           mmm = m%grid(ioct)%rho(ind) * vol_loc
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
     grid_fmm%multipole(icell,:) = grid_fmm%multipole(icell,:) + multipole
#endif
  end do
  call close_cache(s,m%mg_dict)
  end associate
end subroutine fmm_multipole_amr2fmm
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_multipole_fmm2fmm(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_FMM2FMM,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_multipole_fmm2fmm(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_multipole_fmm2fmm(pst%s,ilevel)
  endif

end subroutine r_fmm_multipole_fmm2fmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_multipole_fmm2fmm(s,ilevel)
  use amr_parameters, only: ndim, twotondim, multipole_size
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_flag_module, only: pack_fetch_hydro, unpack_fetch_hydro
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
  integer::ind,idim,ivar,ioct,icell
  real(kind=8)::average
  integer(kind=8),dimension(0:ndim)::hash_key
  logical::leaf_cell
  type(oct),pointer::gridp
  type(msg_large_realdp)::dummy_realdp

  integer :: nm, nd, nq
  real(kind=8), dimension(1:multipole_size) :: multipole

  associate(r=>s%r,g=>s%g,m=>s%m)
  
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
                     init=init_flush_multipole, flush=pack_flush_multipole, combine=unpack_flush_multipole)

  ! Loop over finer level grids
  hash_key(0)=ilevel+1
  do ioct=m%head_mg(ilevel+1),m%tail_mg(ilevel+1)
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     ! Get parent cell using a write-only cache
     call get_parent_cell(s,hash_key,m%mg_dict,gridp,icell,flush_cache=.true.,fetch_cache=.false.,lock=.true.)
     multipole = 0.0D0
#ifdef FMM
     do ind=1,twotondim
       multipole = multipole + m%grid(ioct)%multipole(ind,:)
     end do
     gridp%multipole(icell,:) = gridp%multipole(icell,:) + multipole
#endif
  end do
  call close_cache(s,m%mg_dict)
  end associate
end subroutine fmm_multipole_fmm2fmm
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_multipole_shift_downward(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_MULTIPOLE_SHIFT_DOWNWARD,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_multipole_shift_downward(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_multipole_shift_downward(pst%s,ilevel)
  endif

end subroutine r_fmm_multipole_shift_downward
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_multipole_shift_downward(s,ilevel)
  use amr_parameters, only: ndim, twotondim, multipole_size
  use amr_commons, only: oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use hydro_flag_module, only: pack_fetch_hydro, unpack_fetch_hydro
  use cache_commons
  use cache
  implicit none
  type(ramses_t)::s
  integer::ilevel
  integer::idim,ioct,icell, nstride
  real(kind=8)::average
  integer(kind=8),dimension(0:ndim)::hash_key
  real(kind=8) :: dx_loc
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  integer(kind=8), dimension(ndim) :: cc_icell! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_icell ! box unit real coordinate

  associate(r=>s%r,g=>s%g,m=>s%m)

  dx_loc = r%boxlen / 2.0D0**ilevel
  hash_key(0)=ilevel
  do ioct=m%head_mg(ilevel),m%tail_mg(ilevel)
     hash_key(1:ndim)=m%grid(ioct)%ckey(1:ndim)
     multipole = 0.0D0
     do icell = 1, twotondim
      do idim = 1, ndim
        nstride = 2**(idim-1)
        cc_icell(idim) = 2*hash_key(idim) + MOD((icell-1)/nstride, 2)
        xx_icell(idim) = (cc_icell(idim) + 0.5d0) * dx_loc - m%skip(idim)
      end do
#ifdef FMM
      multipole = m%grid(ioct)%multipole(icell, :)
      call shift_multipole(multipole, xx_icell, multipole_shifted)
      m%grid(ioct)%multipole(icell, :) = multipole_shifted
#endif
    end do
  end do
  end associate
end subroutine fmm_multipole_shift_downward
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_multipole(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
#ifdef FMM
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  grid%multipole=0.0D0
#endif
end subroutine init_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_multipole(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim,multipole_size
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
#ifdef FMM
  do ivar=1,multipole_size
     do ind=1,twotondim
        msg%realdp_fmm_multipole(ind,ivar)=grid%multipole(ind,ivar)
     end do
  end do
  msg_array=transfer(msg,msg_array)
#endif
end subroutine pack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_multipole(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim,multipole_size
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
#ifdef FMM  
  do ivar=1,multipole_size
     do ind=1,twotondim
        if(grid%refined(ind))then
           grid%multipole(ind,ivar)=grid%multipole(ind,ivar)+msg%realdp_fmm_multipole(ind,ivar)
        endif
     end do
  end do
#endif
end subroutine unpack_flush_multipole
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_reset_multipoles_taylor(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_RESET_MULTIPOLES,pst%iUpper+1,input_size,0,ilevel)
     call r_reset_multipoles_taylor(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call reset_multipoles_taylor(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_reset_multipoles_taylor
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reset_multipoles_taylor(r,g,m,ilevel)
  use amr_parameters, only: twotondim
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)   :: r
  type(global_t):: g
  type(mesh_t)  :: m
  integer       :: ilevel
  integer :: igrid, ind
  integer :: first, last

  if (ilevel <= r%levelmin-r%level_fmm_to_amr) then
     first = m%head_mg(ilevel)
     last  = m%tail_mg(ilevel)
  else if (ilevel > r%nlevelmax .and. ilevel < r%levelmin) then
     return
  else
     first = m%head(ilevel)
     last  = m%tail(ilevel)
  end if
#ifdef FMM
  do igrid = first, last
    m%grid(igrid)%multipole = 0.0D0
    m%grid(igrid)%taylor_coeff = 0.0D0
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
subroutine dump_multipole(r, g, m, ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: nbor, oct, run_t, global_t, mesh_t
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
    write(filename, '(A,I0,A)') "out/mult_fmm_", ilevel, ".out"
#else
    write(filename, '(A,I0,A)') "out/mult_mg_", ilevel, ".out"
#endif
  open(unit_debug, file=filename, status="replace")
  dx_loc = r%boxlen / 2.0D0**ilevel
  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
     do icell = 1, twotondim
        do idim = 1, ndim
          nstride = 2**(idim-1)
          cc_icell(idim) = 2*m%grid(ioct)%ckey(idim) + MOD((icell-1)/nstride, 2)
          xx_icell(idim) = (cc_icell(idim) + 0.5D0) * dx_loc - m%skip(idim)
        end do
#ifdef FMM
        write(unit_debug, '(3I6, E20.4)') cc_icell, m%grid(ioct)%multipole(icell, 1)
#endif
     end do
  end do
  close(unit_debug)
end subroutine dump_multipole
#endif
end module fmm_multipoles

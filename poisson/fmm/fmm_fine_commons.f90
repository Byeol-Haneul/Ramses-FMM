module fmm_fine_commons
contains
! ------------------------------------------------------------------------
! FMM Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine fmm routines, such as
!   * fmm iterations @ MG levels
!   * MG workspace building
!
! Used variables:
#ifdef FMM
subroutine fmm(pst,ilevel,icount)
  use amr_parameters, only: twotondim
  use poisson_parameters, only: ngs_fine, ngs_coarse, ncycles_coarse_safe
  use ramses_commons, only: pst_t
  use phi_fine_cg_module, only: r_make_initial_phi, in_make_initial_phi_t
  use init_fmm_module, only: m_init_fmm
  use fmm_multipoles, only: m_fmm_multipoles
  use cleanup_fmm_module, only: r_cleanup_fmm
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ilevel,icount
  
  integer :: igrid, ifine, i, iter, allmasked, ilev
  integer,dimension(1:4) :: output_array
  type(in_make_initial_phi_t)::in_make_initial_phi
  
  if(pst%s%r%gravity_type>0)return
  if(pst%s%m%noct_tot(ilevel)==0)return
  
  if(pst%s%r%verbose) print '(A,I2)','Entering fmm at level ',ilevel

  ! ---------------------------------------------------------------------
  ! Build FMM hierarchy in memory
  ! ---------------------------------------------------------------------
  if(ilevel==pst%s%r%levelmin) then
    call m_init_fmm(pst)
    if(pst%s%r%verbose) print '(A)','FMM init done ' 
  endif

  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------
   call m_fmm_multipoles(pst, ilevel) ! do upward pass !

  ! Downward pass for fmm grids. 
   do ilev = 2, pst%s%r%levelmin-pst%s%g%level_fmm_to_amr
     call r_fmm_downward(pst, ilev, 1)
     if(pst%s%r%verbose) print '(A,I2)','[M2L & L2L] Downpass for FMM grids at level done', ilev
   end do

   ! Call direct force calculation
   call r_fmm_amr_intermediate(pst, pst%s%r%levelmin, 1)
   if(pst%s%r%verbose) print '(A,I2)','AMR Intermediate Calculation done', pst%s%r%levelmin
   call r_fmm_amr_direct(pst, pst%s%r%levelmin, 1)
   if(pst%s%r%verbose) print '(A,I2)','Direct Force Calculation done', pst%s%r%levelmin

   do ilev = 1, pst%s%r%levelmin - pst%s%g%level_fmm_to_amr
     call dump_taylor(pst%s%r, pst%s%m, ilev)
   end do 
    
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
   if(ilevel==pst%s%r%levelmin) then 
     call r_cleanup_fmm(pst)
   if(pst%s%r%verbose) print '(A)','FMM cleanup done '
  endif
end subroutine fmm

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive fmm routine for coarse MG levels
! ------------------------------------------------------------------------

recursive subroutine recursive_fmm(pst,ifinelevel)
  use amr_parameters, only: twotondim
  use poisson_parameters, only: ngs_fine, ngs_coarse, ncycles_coarse_safe
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ifinelevel

  integer :: i, igrid, icycle, ncycle
  
  if(ifinelevel<=pst%s%r%levelmin - pst%s%g%level_fmm_to_amr) then
     ! Solve 'directly' :
     return
  end if
     
   ! FMM-solve the upper level
   call recursive_fmm(pst,ifinelevel-1)

   ! Interpolate coarse solution and correct back into fine solution
   !call r_interpolate_and_correct(pst,ifinelevel,1)  
end subroutine recursive_fmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_fmm_downward(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_DOWNWARD,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_downward(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_downward(pst%s,ilevel)
  endif

end subroutine r_fmm_downward
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_downward(s, ilevel)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size, taylor_size
  use amr_commons, only: nbor, oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache

  use fmm_multipoles, only: shift_multipole
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilevel

  integer :: ioct, idim, pcell, icell, inbor, jcell, nstride, counter
  integer(kind=8), dimension(ndim) :: cc_grid, cc_icell, cc_jcell, cc_jcell_periodic! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_igrid, xx_icell, xx_jcell, xx_jcell_periodic, xx_pgrid, dx, diff, offset ! box unit real coordinate
  real(kind=8) :: dx_loc
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor, hash_nbor_periodic, hash_parent
  type(nbor), dimension(1:threetondim) :: grid_nbor
  integer, dimension(1:twotondim) :: ind_nbor_cells

  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp)::dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor
  logical::cycle_flag

  associate(r=>s%r, g=>s%g, m=>s%m)

  ! Open cache for multipoles
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,& 
            hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,& 
            pack=pack_fetch_taylor,unpack=unpack_fetch_taylor,& 
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_nbor(0) = ilevel - 1
  hash_nbor_periodic(0) = ilevel - 1
  hash_parent(0) = ilevel - 1
  dx_loc = r%boxlen / 2.0D0**ilevel

  ! Loop over octs at this level
  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
     temp_taylor(:) = 0.0D0
     hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)

    ! Multipole Shifting
    do icell = 1, twotondim
      call get_cell_pos(hash_key, icell, r%boxlen, xx_icell, cc_icell)
      multipole = m%grid(ioct)%multipole(icell, :)
      call shift_multipole(multipole, xx_icell, multipole_shifted)
      m%grid(ioct)%multipole(icell, :) = multipole_shifted
    end do

     ! Far Field Calculation
     call get_parent_cell(s, hash_key, m%mg_dict, gridp_parent, pcell, flush_cache=.false., fetch_cache=.true.)
     parent_taylor = gridp_parent%taylor_coeff
     hash_parent(1:ndim) = gridp_parent%ckey(1:ndim)
     call get_grid_pos(hash_key, r%boxlen, xx_igrid)
     call get_grid_pos(hash_parent, r%boxlen, xx_pgrid)

     call get_displacement(xx_pgrid, xx_igrid, r%boxlen, dx)
     call shift_taylor(parent_taylor, -dx, temp_taylor)

     ! Get neighboring parent grids (returns 3^n parent level grids)
     call get_intermediate_nbor_grid(s, hash_key, m%mg_dict, grid_nbor, flush_cache=.false., fetch_cache=.true.)
     do inbor = 1, threetondim
       ! Get offset for neighboring parent grids (can reach off bounds)
       do idim=1,ndim
         offset(idim) = MOD((inbor-1)/3**(idim-1), 3) - 1
       end do

       ! get actual positions
       gridp_nbor => grid_nbor(inbor)%p
       hash_nbor(1:ndim) = gridp_nbor%ckey(1:ndim)

       ! get periodic positions for the nbor
       hash_nbor_periodic(1:ndim) = hash_parent(1:ndim) + offset

       do jcell = 1, twotondim
          call get_cell_pos(hash_nbor, jcell, r%boxlen, xx_jcell, cc_jcell)
          cycle_flag = .false.
          do idim=1,ndim
            nstride = 2**(idim-1)
            cc_jcell_periodic(idim) = 2*hash_nbor_periodic(idim) + MOD((jcell-1)/nstride, 2)
            xx_jcell_periodic(idim) = (cc_jcell_periodic(idim) + 0.5d0) * (dx_loc*2)
            if ((cc_jcell_periodic(idim)<m%box_ckey_min(idim, ilevel) .or. cc_jcell_periodic(idim)>=m%box_ckey_max(idim, ilevel))) then
              cycle_flag = .true.
            end if 
          end do 
          ! skip direct neighbors
          if (is_direct_neighbor(hash_key(1:ndim), cc_jcell_periodic, ilevel - 1) .or. cycle_flag) cycle
          ! Shift multipole from origin -> source center (Need to use grid position)
          multipole = gridp_nbor%multipole(jcell,:)
          ! Get taylor coeffs from local
          dx = (hash_key(1:ndim) - cc_jcell_periodic) * r%boxlen / 2**(ilevel-1)
          call calc_taylor_from_multipole(dx, multipole, temp_taylor)
       end do ! over neighboring grid's cells 2^n
     end do ! over neighboring grids 3^n 
     ! Add taylor coefficients from intermediate fields
     m%grid(ioct)%taylor_coeff = m%grid(ioct)%taylor_coeff + temp_taylor
     ! Unlock neighbor octs
     do inbor = 1, threetondim
        call unlock_cache(s, grid_nbor(inbor)%p)
     end do
  end do
  call close_cache(s, m%mg_dict)
  end associate
end subroutine fmm_downward
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_amr_intermediate(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_AMR_INTERMEDIATE,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_amr_intermediate(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_amr_intermediate(pst%s,ilevel)
  endif

end subroutine r_fmm_amr_intermediate
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_amr_intermediate(s, ilevel)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size, taylor_size
  use amr_commons, only: nbor, oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_multipoles, only: shift_multipole
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilevel

  integer :: ioct, idim, ind, icell, jcell, nstride, nfine
  real(kind=8) :: phi, phi_out, fourpi, dx_loc
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_jcell_periodic, offset
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, xx_pgrid, xx_jcell_periodic, diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_nbor, hash_nbor_periodic, prev_hash_fmm_grid

  ! --- Precomputed lists to avoid repeated get_cell_pos calls ---
  real(kind=8), dimension(twotondim, ndim) :: xx_icell_list, xx_jcell_list
  integer(kind=8), dimension(twotondim, ndim) :: cc_icell_list, cc_jcell_list

  type(nbor), dimension(1:threetondim) :: grid_nbor
  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp) :: dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor
  logical :: cycle_flag, neighbors_cached

  associate(r=>s%r, g=>s%g, m=>s%m)
  fourpi = 4.D0*ACOS(-1.0D0)
  if(r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  ! Open cache for multipoles
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,&
                  hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,&
                  pack=pack_fetch_taylor, unpack=unpack_fetch_taylor,&
                  init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_fmm_grid(0) = ilevel - g%level_fmm_to_amr
  prev_hash_fmm_grid(0) = ilevel - g%level_fmm_to_amr
  hash_fmm_cell(0) = ilevel - g%level_fmm_to_amr + 1
  prev_hash_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilevel
  nfine = 2**g%level_fmm_to_amr
  neighbors_cached = .false.

  ! Loop over octs at this level
  do ioct = m%head(ilevel), m%tail(ilevel)
    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / nfine
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim) / (nfine/2)

    ! Check if we need to fetch neighbors & parent Taylor
    if (.not. all(hash_fmm_grid == prev_hash_fmm_grid)) then
      if (neighbors_cached) then
        do ind = 1, threetondim
          call unlock_cache(s, grid_nbor(ind)%p)
        end do
      end if

      call get_grid_pos(hash_fmm_grid, r%boxlen, xx_pgrid)
      call get_grid(s, hash_fmm_grid, m%mg_dict, gridp_parent, flush_cache=.false., fetch_cache=.true.)
      parent_taylor = gridp_parent%taylor_coeff

      call get_intermediate_nbor_grid(s, hash_fmm_cell, m%mg_dict, grid_nbor, flush_cache=.false., fetch_cache=.true.)
      neighbors_cached = .true.
      prev_hash_fmm_grid = hash_fmm_grid
    end if

    ! --- Precompute AMR cell positions ---
    do icell = 1, twotondim
      call get_cell_pos(hash_key, icell, r%boxlen, xx_icell_list(icell, :), cc_icell_list(icell, :))
    end do

    ! Loop over AMR cells
    do icell = 1, twotondim
      xx_icell = xx_icell_list(icell, :)
      cc_icell = cc_icell_list(icell, :)

      ! Far field
      phi = 0.0D0
      call calc_phi(parent_taylor, xx_icell - xx_pgrid, phi)

      ! Intermediate field
      temp_taylor = 0.0D0
      hash_nbor(0) = ilevel - g%level_fmm_to_amr
      do ind = 1, threetondim
        gridp_nbor => grid_nbor(ind)%p
        hash_nbor(1:ndim) = gridp_nbor%ckey(1:ndim)

        ! Precompute neighbor cell positions
        do jcell = 1, twotondim
          call get_cell_pos(hash_nbor, jcell, r%boxlen, xx_jcell_list(jcell, :), cc_jcell_list(jcell, :))
        end do

        do jcell = 1, twotondim
          xx_jcell = xx_jcell_list(jcell, :)
          cc_jcell = cc_jcell_list(jcell, :)

          ! Skip direct neighbor
          if (is_direct_neighbor(hash_fmm_cell(1:ndim), cc_jcell, ilevel - g%level_fmm_to_amr)) cycle

          ! Wrap-around offsets
          do idim = 1, ndim
            offset(idim) = MOD((ind-1)/3**(idim-1), 3) - 1
          end do

          cycle_flag = .false.
          hash_nbor_periodic(1:ndim) = hash_fmm_grid(1:ndim) + offset
          do idim = 1, ndim
            nstride = 2**(idim-1)
            cc_jcell_periodic(idim) = 2*hash_nbor_periodic(idim) + MOD((jcell-1)/nstride, 2)
            xx_jcell_periodic(idim) = (cc_jcell_periodic(idim) + 0.5d0) * (dx_loc*nfine)
#ifdef PERIODIC
            if ((cc_jcell_periodic(idim) < m%box_ckey_min(idim, ilevel)) .or. &
                (cc_jcell_periodic(idim) >= m%box_ckey_max(idim, ilevel))) then
              cycle_flag = .true.
            end if
#endif
          end do
          if (cycle_flag) cycle

          multipole = gridp_nbor%multipole(jcell, 1:multipole_size)
          call get_displacement(xx_icell, xx_jcell_periodic, r%boxlen, diff)
          call calc_phi_from_multipole(diff, multipole, phi_out)
          phi = phi + phi_out
        end do
      end do
      m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + phi
    end do
  end do

  call close_cache(s, m%mg_dict)
  end associate
end subroutine fmm_amr_intermediate
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_amr_direct(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_AMR_DIRECT,pst%iUpper+1,input_size,0,ilevel)
     call r_fmm_amr_direct(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_amr_direct(pst%s,ilevel)
  endif

end subroutine r_fmm_amr_direct
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_amr_direct(s, ilevel)
  use amr_parameters, only: ndim, twotondim, threetondim
  use amr_commons, only: nbor, oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilevel

  integer :: ioct, idim, ind, icell, jcell, jcell_amr, nstride, nfine
  integer :: i, j, k
  real(kind=8) :: phi, dist2, fourpi, dx_loc, dxn
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_fmm_cell
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_direct, hash_prev_fmm_grid, prev_hash_fmm_cell

  type(nbor), dimension(1:threetondim) :: direct_grid_nbor
  integer, dimension(1:threetondim) :: direct_ind_nbor
  type(oct), pointer :: gridp_nbor
  type(msg_large_realdp) :: dummy_realdp
  logical :: cycle_flag, initialized

  ! --- New locals for optimizations ---
  real(kind=8), dimension(twotondim, ndim) :: xx_icell_list, xx_jcell_list
  integer(kind=8), dimension(twotondim, ndim) :: cc_icell_list, cc_jcell_list
  logical :: neighbors_cached

  associate(r=>s%r, g=>s%g, m=>s%m)

  fourpi = 4.D0*ACOS(-1.0D0)
  if (r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  ! Open cache for multipoles (unchanged)
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,&
            hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,&
            pack=pack_fetch_taylor,unpack=unpack_fetch_taylor,&
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_fmm_grid(0) = ilevel - g%level_fmm_to_amr
  hash_prev_fmm_grid(0) = ilevel - g%level_fmm_to_amr
  hash_fmm_cell(0) = ilevel - g%level_fmm_to_amr + 1
  hash_direct(0) = ilevel

  hash_prev_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilevel
  dxn = dx_loc**ndim                      ! precompute (used repeatedly)
  nfine = 2**g%level_fmm_to_amr
  initialized = .false.

  prev_hash_fmm_cell = -huge(0_8)         ! impossible initial key
  neighbors_cached = .false.

  ! Loop over octs at this level
  do ioct = m%head(ilevel), m%tail(ilevel)
    ! set keys for this amr grid
    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / nfine
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim) / (nfine/2)

    ! If parent fmm cell changed, fetch (and unlock previous) neighbor info
    if (.not. neighbors_cached .or. .not. all(hash_fmm_cell == prev_hash_fmm_cell)) then
      if (neighbors_cached) then
        do ind = 1, threetondim
          call unlock_cache(s, direct_grid_nbor(ind)%p)
        end do
      end if
      call get_threetondim_nbor_parent_cell(s, hash_fmm_cell, m%mg_dict, direct_grid_nbor, &
           direct_ind_nbor, flush_cache = .false., fetch_cache = .true.)
      prev_hash_fmm_cell = hash_fmm_cell
      neighbors_cached = .true.
    end if

    ! Precompute all target cell positions for this AMR grid (once)
    do icell = 1, twotondim
      call get_cell_pos(hash_key, icell, r%boxlen, xx_icell_list(icell, :), cc_icell_list(icell, :))
    end do

    ! For each neighbor fmm cell (direct neighbors), compute direct AMR child grids
    do ind = 1, threetondim
      ! Get fmm cell info
      jcell = direct_ind_nbor(ind)
      gridp_nbor => direct_grid_nbor(ind)%p
      do idim = 1, ndim
        nstride = 2**(idim-1)
        cc_fmm_cell(idim) = 2*gridp_nbor%ckey(idim) + MOD((jcell-1)/nstride, 2)
      end do

      ! Find nfine/2 amr grids (loops over i/j/k)
#if NDIM>2
      do k = 1, nfine/2
        hash_direct(3) = (nfine/2) * cc_fmm_cell(3) + k - 1
#endif
#if NDIM>1
      do j = 1, nfine/2
        hash_direct(2) = (nfine/2) * cc_fmm_cell(2) + j - 1
#endif
#if NDIM>0
      do i = 1, nfine/2
        hash_direct(1) = (nfine/2) * cc_fmm_cell(1) + i - 1
#endif
          ! periodic boundary conditions & skipping out-of-box grids
          cycle_flag = .false.
#ifdef PERIODIC
          do idim = 1, ndim
            if (r%periodic(idim)) then
              if (hash_direct(idim) < m%box_ckey_min(idim, ilevel)) then
                hash_direct(idim) = m%box_ckey_max(idim, ilevel) - 1
              end if
              if (hash_direct(idim) >= m%box_ckey_max(idim, ilevel)) then
                hash_direct(idim) = m%box_ckey_min(idim, ilevel)
              end if
            end if
            if (hash_direct(idim) < m%box_ckey_min(idim, ilevel) .OR. &
                hash_direct(idim) >= m%box_ckey_max(idim, ilevel)) then
              cycle_flag = .true.
            end if
          end do
#endif
          if (cycle_flag) cycle

          ! get the neighbor AMR grid (source) and precompute its cell positions (once)
          call get_grid(s, hash_direct, m%grid_dict, gridp_nbor, flush_cache = .false., fetch_cache = .true.)
          do jcell_amr = 1, twotondim
            call get_cell_pos(hash_direct, jcell_amr, r%boxlen, xx_jcell_list(jcell_amr, :), &
                              cc_jcell_list(jcell_amr, :))
          end do

          ! Now loop over target cells and accumulate direct interactions using vector ops
          do icell = 1, twotondim
            phi = 0.0D0
            ! use precomputed target pos
            xx_icell(:) = xx_icell_list(icell, :)

            do jcell_amr = 1, twotondim
              ! Skip self same-grid same-cell interaction
              if (all(hash_direct(1:ndim) == hash_key(1:ndim)) .and. icell == jcell_amr) cycle

              ! vector displacement with periodic wrap
              xx_jcell(:) = xx_jcell_list(jcell_amr, :)
              call get_displacement(xx_jcell, xx_icell, r%boxlen, diff)
              dist2 = sum(diff(:)**2)

              phi = phi - (gridp_nbor%rho(jcell_amr) - g%rho_tot) * dxn / sqrt(dist2)
            end do

            m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + phi
            m%grid(ioct)%f(icell,2) = m%grid(ioct)%phi(icell)
          end do  ! target cells loop

#if NDIM>0
      end do
#endif
#if NDIM>1
      end do
#endif
#if NDIM>2
      end do
#endif
    end do ! end direct force calculation for given amr grid
  end do ! end over all amr grids @ given ilevel

  ! Unlock final cached neighbor grids (if any)
  if (neighbors_cached) then
    do ind = 1, threetondim
      call unlock_cache(s, direct_grid_nbor(ind)%p)
    end do
  end if

  call close_cache(s, m%mg_dict)
  end associate
end subroutine fmm_amr_direct
!################################################################
!################################################################
!################################################################
!################################################################
logical function is_direct_neighbor(cc_icell, cc_jcell, ilevel)
  use amr_parameters, only: ndim
  implicit none
  integer(kind=8), intent(in) :: cc_icell(ndim), cc_jcell(ndim)
  integer, intent(in) :: ilevel
  integer :: d, n, diff
  
  is_direct_neighbor = .true.
  do d = 1, ndim
     diff = abs(cc_icell(d) - cc_jcell(d))
     !if (ilevel > 1) then
      !diff = min(diff, 2**ilevel - diff)
     !end if
     if (diff > 1) then
        is_direct_neighbor = .false.
        return
     end if
  end do
end function is_direct_neighbor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_grid_pos(hash_key, boxlen, pos)
  use amr_parameters, only: ndim
  implicit none
  integer(kind=8), intent(in)  :: hash_key(0:ndim)   ! (0)=ilevel, (1:ndim)=spatial keys
  real(kind=8),    intent(in)  :: boxlen             ! size of domain
  real(kind=8),    intent(out) :: pos(ndim)          ! grid center position in physical units

  integer :: idim, ilevel
  real(kind=8) :: dx_loc

  ilevel = hash_key(0)
  dx_loc = boxlen / 2.0d0**(ilevel-1)
  do idim = 1, ndim
     pos(idim) = (hash_key(idim) + 0.5d0) * dx_loc
  end do
end subroutine get_grid_pos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_cell_pos(hash_key, icell, boxlen, pos, cc_icell)
  use amr_parameters, only: ndim
  implicit none
  integer(kind=8), intent(in)  :: hash_key(0:ndim)   ! (0)=ilevel, (1:ndim)=spatial keys
  integer,        intent(in)   :: icell              ! local cell index (1..2^ndim)
  real(kind=8),   intent(in)   :: boxlen             ! size of domain
  real(kind=8),   intent(out)  :: pos(:)          ! cell center position in physical units
  integer(kind=8),   intent(out)  :: cc_icell(:)          ! cartesian coordinate

  integer :: idim, ilevel, nstride
  real(kind=8) :: dx_loc

  ilevel = hash_key(0)
  dx_loc = boxlen / 2.0d0**ilevel
  do idim = 1, ndim
     nstride = 2**(idim-1)
     cc_icell(idim) = 2*hash_key(idim) + MOD((icell-1)/nstride, 2)
     pos(idim)      = (cc_icell(idim) + 0.5d0) * dx_loc
  end do
end subroutine get_cell_pos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_taylor(grid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  type(oct)::grid
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  grid%multipole=0.0D0
  grid%taylor_coeff=0.0D0
end subroutine init_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_taylor(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
  do ivar=0,ndim+ int(ndim * (ndim+1)/2)
    msg%realdp_fmm_taylor(ivar)=grid%taylor_coeff(ivar)
  end do
  msg_array=transfer(msg,msg_array)
end subroutine pack_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_taylor(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim,taylor_size
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
  
  do ivar=1,taylor_size
    grid%taylor_coeff(ivar)=grid%taylor_coeff(ivar)+msg%realdp_fmm_taylor(ivar)
  end do
end subroutine unpack_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_taylor(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim,multipole_size
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
  msg%realdp_fmm_multipole=grid%multipole
  msg%realdp_fmm_taylor=grid%taylor_coeff
  msg_array=transfer(msg,msg_array)
end subroutine pack_fetch_taylor
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_taylor(grid,msg_size,msg_array,hash_key)
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

  grid%multipole=msg%realdp_fmm_multipole
  grid%taylor_coeff=msg%realdp_fmm_taylor
end subroutine unpack_fetch_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine dump_taylor(r, m, ilevel)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, mesh_t
  implicit none
  type(run_t) :: r
  type(mesh_t) :: m
  integer, intent(in) :: ilevel

  integer :: ioct, icell, unit_debug
  real(kind=8) :: dx_loc
  character(len=256) :: filename

  ! Construct filename based on level
  write(filename, '(A,I0,A)') "./out_fmm/taylor_level", ilevel, ".out"

  unit_debug = 999
  open(unit_debug, file=filename, status="replace")

  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
    write(unit_debug, '(3I6, 20E20.5)') m%grid(ioct)%ckey(1:ndim), m%grid(ioct)%taylor_coeff
  end do

  close(unit_debug)
end subroutine dump_taylor
#endif
end module fmm_fine_commons

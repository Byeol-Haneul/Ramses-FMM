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
#ifdef GRAV
subroutine fmm(pst,ilevel,icount)
  use amr_parameters, only: twotondim, nhilbert
  use poisson_parameters, only: ngs_fine, ngs_coarse, ncycles_coarse_safe
  use ramses_commons, only: pst_t
  use phi_fine_cg_module, only: r_make_initial_phi, in_make_initial_phi_t
  use init_fmm_module, only: r_init_fmm
  use fmm_multipoles!, only: m_fmm_multipoles
  use cleanup_fmm_module, only: r_cleanup_fmm
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ilevel,icount
  
  integer :: igrid, ifine, i, ierr, allmasked, ilev
  integer,dimension(1:4) :: output_array
  type(in_make_initial_phi_t)::in_make_initial_phi
  
  if(pst%s%r%gravity_type>0)return
  if(pst%s%m%noct_tot(ilevel)==0)return
  
  if(pst%s%r%verbose) print '(A,I2)','Entering fmm at level ',ilevel

  ! ---------------------------------------------------------------------
  ! Build FMM hierarchy in memory
  ! ---------------------------------------------------------------------

  if(ilevel==pst%s%r%levelmin) then
    call r_init_fmm(pst, ilevel, 1)
    if(pst%s%r%verbose) print '(A)','FMM init done ' 
  endif

  if(pst%s%r%verbose) print '(A)','FMM init done ' 

  ! ---------------------------------------------------------------------
  ! Initiate solve at fine level
  ! ---------------------------------------------------------------------
   do i = pst%s%r%bound_levelmin, pst%s%r%levelmin-pst%s%r%level_fmm_to_amr, 1
    call r_reset_multipoles_taylor(pst, i, 1)
   end do

   call m_timer(pst,'fmm: multipole upward','start')
   call m_fmm_multipoles(pst, ilevel) ! do upward pass !

  ! Downward pass for fmm grids. 
   call m_timer(pst,'fmm: downward for fmm','start')
   do ilev = pst%s%r%bound_levelmin+1, pst%s%r%levelmin-pst%s%r%level_fmm_to_amr
     call r_fmm_downward(pst, ilev, 1)
     if(pst%s%r%verbose) print '(A,I2)','[M2L & L2L] Downpass for FMM grids at level done', ilev
   end do

   ! Call direct force calculation
   call m_timer(pst,'fmm: amr intermediate force','start')
   call r_fmm_amr_intermediate(pst, pst%s%r%levelmin, 1)
   if(pst%s%r%verbose) print '(A,I2)','AMR Intermediate Calculation done', pst%s%r%levelmin

   call m_timer(pst,'fmm: direct force','start')
   call r_fmm_amr_direct(pst, pst%s%r%levelmin, 1)
   if(pst%s%r%verbose) print '(A,I2)','Direct Force Calculation done', pst%s%r%levelmin

   !do ilev = 1, pst%s%r%levelmin - pst%s%r%level_fmm_to_amr
   !  call dump_taylor(pst%s%r, pst%s%m, ilev)
   !end do 
    
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
   call m_timer(pst,'fmm: cleanup','start')
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
  
  if(ifinelevel<=pst%s%r%levelmin - pst%s%r%level_fmm_to_amr) then
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

  integer :: ioct, idim, pcell, icell, inbor, jcell, nstride
  integer(kind=8), dimension(ndim) :: cc_grid, cc_icell, cc_jcell, cc_jcell_periodic, offset! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_igrid, xx_icell, xx_jcell, xx_jcell_periodic, xx_pgrid, dx, diff ! box unit real coordinate
  real(kind=8) :: dx_loc, dist, D0, D1, D2, D3
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor, hash_nbor_periodic, hash_parent
  type(nbor), dimension(1:threetondim) :: grid_nbor
  integer, dimension(1:twotondim) :: ind_nbor_cells

  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp)::dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor
  real(kind=8), dimension(twotondim, taylor_size) :: accum_taylor
  logical::cycle_flag

  integer(kind=8), dimension(threetondim, ndim) :: offset_list
  real(kind=8), dimension(threetondim, twotondim, twotondim, twotondim) :: D0_list, D1_list, D2_list, D3_list
  real(kind=8), dimension(threetondim, twotondim, twotondim, twotondim, ndim) :: intermediate_diff_list
  logical, dimension(threetondim, twotondim, twotondim) :: direct_neighbor_list

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )
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

  if(m%noct_mg(ilevel)<1) return

  ! jcell to icell
  do inbor = 1, threetondim
    ! calculate offsets
    do idim = 1, ndim
      offset_list(inbor, idim) = MOD((inbor-1)/3**(idim-1), 3) - 1 ! offset by how many fmm parent grids
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(inbor,:) + displacement_list(jcell,:) ! respect to parent grid left corner / fmm grid unit
      do pcell = 1,twotondim
        cycle_flag = .true.
        do idim=1, ndim
          if (abs(cc_jcell(idim) - displacement_list(pcell, idim)) > 1) cycle_flag = .false.
        end do
        direct_neighbor_list(inbor, jcell, pcell) = cycle_flag
        do icell = 1, twotondim
          cc_icell = 2 * displacement_list(pcell, :) + displacement_list(icell, :)
          diff = (cc_icell - 0.5 - 2 * cc_jcell) * dx_loc
          intermediate_diff_list(inbor, jcell, pcell, icell, :) = diff
          dist = sqrt(sum(diff(:)**2))
          !print *, diff/dx_loc, cc_jcell, cc_icell
          D0_list(inbor, jcell, pcell, icell) = 1.0D0 / dist
          D1_list(inbor, jcell, pcell, icell) = -1.0D0 / dist**3
          D2_list(inbor, jcell, pcell, icell) = 3.0D0 / dist**5
          D3_list(inbor, jcell, pcell, icell) = -15.0D0 / dist**7
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
     accum_taylor(:, :) = 0.0D0
     hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)

    call get_parent_cell(s, hash_key, m%mg_dict, gridp_parent, pcell, flush_cache=.false., fetch_cache=.true.)
#ifdef FMM
    parent_taylor = gridp_parent%taylor_coeff(pcell, :)
#endif
    hash_parent(1:ndim) = gridp_parent%ckey(1:ndim)

    do icell = 1, twotondim
      ! Multipole Shifting
      call get_cell_pos(hash_key, icell, r%boxlen, xx_icell, cc_icell)
#ifdef FMM
      multipole = m%grid(ioct)%multipole(icell, :)
      call shift_multipole(multipole, xx_icell, multipole_shifted)
      m%grid(ioct)%multipole(icell, :) = multipole_shifted
#endif
      ! Far Field Calculation
      do idim =1,ndim
        dx(idim) = (displacement_list(icell, idim) - 0.5) * dx_loc
      end do 
      call shift_taylor(parent_taylor, dx, temp_taylor)
      accum_taylor(icell,:) = accum_taylor(icell,:) + temp_taylor
    end do

     ! Get neighboring parent grids (returns 3^n parent level grids)
     call get_intermediate_nbor_grid(s, hash_key, m%mg_dict, grid_nbor, flush_cache=.false., fetch_cache=.true.)
     do inbor = 1, threetondim
       ! Get offset for neighboring parent grids (can reach off bounds)
       do idim=1,ndim
         offset(idim) = MOD((inbor-1)/3**(idim-1), 3) - 1
       end do

       gridp_nbor => grid_nbor(inbor)%p
       hash_nbor_periodic(1:ndim) = hash_parent(1:ndim) + offset

       do jcell = 1, twotondim
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
          if (direct_neighbor_list(inbor, jcell, pcell) .or. cycle_flag) cycle
          ! Shift multipole from origin -> source center (Need to use grid position)
#ifdef FMM
          multipole = gridp_nbor%multipole(jcell,:)
#endif
          ! Get taylor coeffs from local
          do icell=1, twotondim
            dx = intermediate_diff_list(inbor, jcell, pcell, icell, :)
            D0 = D0_list(inbor, jcell, pcell, icell)
            D1 = D1_list(inbor, jcell, pcell, icell)
            D2 = D2_list(inbor, jcell, pcell, icell)
            D3 = D3_list(inbor, jcell, pcell, icell)
            call calc_taylor_from_multipole(dx, D0, D1, D2, D3, multipole, temp_taylor)
            accum_taylor(icell, :) = accum_taylor(icell, :) + temp_taylor
          end do
       end do ! over neighboring grid's cells 2^n
     end do ! over neighboring grids 3^n 
     ! Add taylor coefficients from intermediate fields
#ifdef FMM
     m%grid(ioct)%taylor_coeff = m%grid(ioct)%taylor_coeff + accum_taylor
#endif
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

  integer :: ioct, idim, ind, icell, jcell, nstride, nfine, igrid, nbox, pcell
  real(kind=8) :: phi, phi_out, fourpi, dx_loc
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_jcell_periodic, offset
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, xx_pgrid, xx_jcell_periodic, diff, diff2
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_nbor, hash_nbor_periodic, prev_hash_fmm_grid, prev_hash_fmm_cell

  real(kind=8), dimension(twotondim, ndim) :: xx_icell_list, xx_jcell_list
  integer(kind=8), dimension(twotondim, ndim) :: cc_icell_list, cc_jcell_list

  type(nbor), dimension(1:threetondim) :: grid_nbor
  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp) :: dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole, multipole_shifted
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor
  logical :: cycle_flag, neighbors_cached

  real(kind=8) :: dist, D0, D1, D2
  real(kind=8), dimension(threetondim, ndim) :: offset_list

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )

  ! =====================================================
  ! Runtime-allocated arrays depending on r%level_fmm_to_amr
  ! =====================================================
  real(kind=8), allocatable :: D0_list(:,:,:,:), D1_list(:,:,:,:), D2_list(:,:,:,:)
  real(kind=8), allocatable :: intermediate_diff_list(:,:,:,:,:), cell_diff_list(:,:,:,:)
  real(kind=8), allocatable :: far_diff_list(:,:,:)
  integer, allocatable :: fmm_grid_center_offset(:,:), fmm_cell_center_offset(:,:)
  logical, allocatable :: direct_neighbor_list(:,:,:)

  associate(r=>s%r, g=>s%g, m=>s%m)
  fourpi = 4.D0*ACOS(-1.0D0)
  if(r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  ! Open cache for multipoles
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,&
                  hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,&
                  pack=pack_fetch_taylor, unpack=unpack_fetch_taylor,&
                  init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_fmm_grid(0) = ilevel - r%level_fmm_to_amr
  prev_hash_fmm_grid(0) = ilevel - r%level_fmm_to_amr
  hash_fmm_cell(0) = ilevel - r%level_fmm_to_amr + 1
  prev_hash_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilevel
  nfine = 2**r%level_fmm_to_amr
  neighbors_cached = .false.

  ! nbox is the number of cells at the target AMR level
  nbox = nfine ** ndim

  allocate(D0_list(threetondim, twotondim, nbox, twotondim))
  allocate(D1_list(threetondim, twotondim, nbox, twotondim))
  allocate(D2_list(threetondim, twotondim, nbox, twotondim))

  allocate(intermediate_diff_list(threetondim, twotondim, nbox, twotondim, ndim))
  allocate(far_diff_list(nbox, twotondim, ndim))
  allocate(fmm_grid_center_offset(nbox, ndim))
  allocate(fmm_cell_center_offset(nbox, ndim))
  allocate(direct_neighbor_list(threetondim, twotondim, nbox)) ! we can reduce this if we really need to
  allocate(cell_diff_list(threetondim, twotondim, nbox, ndim))

  ! Precalculate differences
  do igrid=1, nbox
    do idim = 1,ndim
      nstride = nfine**(idim-1)
      fmm_grid_center_offset(igrid, idim) = MOD((igrid-1)/nstride, nfine) - (nfine/2) ! offset by how many amr octs from fmm grid center
      nstride = (nfine/2)**(idim-1)
      fmm_cell_center_offset(igrid, idim) = 2 * MOD(fmm_grid_center_offset(igrid, idim) + (nfine/2), nfine/2) - (nfine/2) ! offset by how many amr cells from fmm cell center
    end do 
    do icell = 1, twotondim
      far_diff_list(igrid, icell, :) = (fmm_cell_center_offset(igrid, :) + displacement_list(icell,:) + 0.5) * dx_loc
    end do
  end do 

  ! jcell to icell
  do ind = 1, threetondim
    ! calculate offsets
    do idim = 1, ndim
      offset_list(ind, idim) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(ind,:) + displacement_list(jcell,:) ! respect to grid left corner / fmm cell unit
      offset = (cc_jcell+ 0.5) * nfine ! offset by how many amr cells
      do igrid=1,nbox
        cc_icell = (fmm_grid_center_offset(igrid, :) + nfine/2)/(nfine/2) ! respect to grid left corner / fmm cell unit
        cycle_flag = .true.
        cell_diff_list(ind, jcell, igrid, :) = cc_icell - cc_jcell
        do idim=1, ndim
          if (abs(cc_icell(idim) - cc_jcell(idim)) > 1) cycle_flag = .false.
        end do
        direct_neighbor_list(ind, jcell, igrid) = cycle_flag
        do icell=1, twotondim
          diff = ((2 * fmm_grid_center_offset(igrid, :) + displacement_list(icell,:) + 0.5) + (- offset(:) + nfine)) * dx_loc
          intermediate_diff_list(ind, jcell, igrid, icell, :) = diff
          dist = sqrt(sum(diff(:)**2))
          D0_list(ind, jcell, igrid, icell) = 1.0D0 / dist
          D1_list(ind, jcell, igrid, icell) = -1.0D0 / dist**3
          D2_list(ind, jcell, igrid, icell) = 3.0D0 / dist**5
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilevel), m%tail(ilevel)
    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / nfine
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim) / (nfine/2)

    igrid = 1
    do idim=1,ndim
      nstride = nfine**(idim-1)
      igrid = igrid + nstride * MOD(hash_key(idim), nfine)
    end do

    ! Check if we need to fetch neighbors & parent Taylor
    if (.not. all(hash_fmm_grid == prev_hash_fmm_grid)) then
      if (neighbors_cached) then
        do ind = 1, threetondim
          call unlock_cache(s, grid_nbor(ind)%p)
        end do
      end if

      call get_grid(s, hash_fmm_grid, m%mg_dict, gridp_parent, flush_cache=.false., fetch_cache=.true.)

      call get_intermediate_nbor_grid(s, hash_fmm_cell, m%mg_dict, grid_nbor, flush_cache=.false., fetch_cache=.true.)
      neighbors_cached = .true.
      prev_hash_fmm_grid = hash_fmm_grid
    end if

    pcell = 1
    do idim=1,ndim
      nstride = 2**(idim-1)
      pcell = pcell + nstride * MOD(hash_fmm_cell(idim), 2)
    end do
#ifdef FMM
    parent_taylor = gridp_parent%taylor_coeff(pcell, :)
#endif
    ! Far field
    do icell = 1, twotondim
      diff = far_diff_list(igrid, icell, :)
      call calc_phi(parent_taylor, diff, phi)
#ifdef FMM
      m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + phi
#endif
    end do

    ! Intermediate field
    do ind = 1, threetondim
      gridp_nbor => grid_nbor(ind)%p
      do jcell = 1, twotondim
        cycle_flag = .false.
        cc_jcell_periodic = hash_fmm_cell(1:ndim) - cell_diff_list(ind, jcell, igrid,:)
        do idim = 1, ndim
          if ((cc_jcell_periodic(idim) < m%box_ckey_min(idim, ilevel - r%level_fmm_to_amr + 1)) .or. &
              (cc_jcell_periodic(idim) >= m%box_ckey_max(idim, ilevel - r%level_fmm_to_amr + 1))) then
            cycle_flag = .true.
          end if
        end do
        if (cycle_flag .or. direct_neighbor_list(ind, jcell, igrid)) cycle
#ifdef FMM
        multipole = gridp_nbor%multipole(jcell, 1:multipole_size)
#endif
        do icell=1, twotondim
          diff  = intermediate_diff_list(ind, jcell, igrid, icell, :)
          D0 = D0_list(ind, jcell, igrid, icell)
          D1 = D1_list(ind, jcell, igrid, icell)
          D2 = D2_list(ind, jcell, igrid, icell)
          call calc_phi_from_multipole(diff, D0, D1, D2, multipole, phi_out)
          m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + phi_out
        end do
      end do
    end do
  end do

  deallocate(D0_list, D1_list, D2_list, intermediate_diff_list, far_diff_list, direct_neighbor_list, cell_diff_list, fmm_grid_center_offset, fmm_cell_center_offset)
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
  use amr_parameters, only: ndim, twotondim, threetondim, nhilbert
  use amr_commons, only: nbor, oct
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilevel

  integer :: ioct, idim, ind, icell, jcell, jcell_amr, nstride, nfine, nbox, jgrid, igrid
  integer :: i, j, k, grid_idx, total_grids, cell_idx
  real(kind=8) :: phi, dist2, fourpi, dx_loc, dxn, dist
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_igrid, cc_jgrid, cc_fmm_cell, offset
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_direct, hash_prev_fmm_grid, prev_hash_fmm_cell
  real(kind=8), dimension(threetondim, ndim) :: offset_list
  type(oct), pointer :: gridp_nbor
  type(msg_small_realdp) :: dummy_realdp
  logical :: cycle_flag, initialized
  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
    [ &
      0, 1, 0, 1, 0, 1, 0, 1,  &
      0, 0, 1, 1, 0, 0, 1, 1,  &
      0, 0, 0, 0, 1, 1, 1, 1   &
    ], [twotondim, ndim] )

  ! Arrays sized for all source cells: threetondim * (nfine/2)^ndim * twotondim
  real(kind=8), dimension(:,:,:), allocatable      :: mm_jcell_list
  real(kind=8), dimension(:,:,:,:,:), allocatable  :: inv_dist

  associate(r=>s%r, g=>s%g, m=>s%m)

  fourpi = 4.D0*ACOS(-1.0D0)
  if (r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  ! Open cache for multipoles (unchanged)
  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
            hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
            pack=pack_fetch_rho, unpack=unpack_fetch_rho,&
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_fmm_grid(0) = ilevel - r%level_fmm_to_amr
  hash_prev_fmm_grid(0) = ilevel - r%level_fmm_to_amr
  hash_fmm_cell(0) = ilevel - r%level_fmm_to_amr + 1
  hash_direct(0) = ilevel

  hash_prev_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilevel
  dxn    = dx_loc**ndim
  nfine  = 2**r%level_fmm_to_amr
  nbox   = (nfine/2)**ndim
  initialized = .false.

  ! Allocate arrays for all possible source cells
  allocate(mm_jcell_list(threetondim, nbox, twotondim))
  allocate(inv_dist(nbox, twotondim, threetondim, nbox, twotondim))

  prev_hash_fmm_cell = -huge(0_8)

  ! Get offset lists
  do ind = 1, threetondim
    do idim = 1, ndim
      offset_list(ind, idim) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids
    end do
  end do

  ! target cell
  do igrid = 1, nbox
    do idim = 1, ndim
      cc_igrid(idim) = MOD((igrid-1)/(nfine/2)**(idim-1), nfine/2)
    end do
    do icell = 1, twotondim
      cc_icell = 2 * cc_igrid + displacement_list(icell, :)
      do ind = 1, threetondim
        do jgrid = 1, nbox
          do idim = 1, ndim
            cc_jgrid(idim) = offset_list(ind, idim) * (nfine/2) + MOD((jgrid-1)/(nfine/2)**(idim-1), nfine/2)
          end do
          do jcell = 1, twotondim
            cc_jcell = 2 * cc_jgrid + displacement_list(jcell, :)
            if (all(cc_icell(1:ndim) == cc_jcell(1:ndim))) then
              inv_dist(igrid, icell, ind, jgrid, jcell) = 1 / dx_loc !! cap it to dx_loc instead of skipping
            else
              diff = (cc_icell - cc_jcell) * dx_loc
              inv_dist(igrid, icell, ind, jgrid, jcell) = 1.d0 / sqrt(sum(diff(:)**2))
            end if
          end do
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilevel), m%tail(ilevel)
    ! set keys for this amr grid
    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / nfine
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim) / (nfine/2)

    igrid = 1
    do idim = 1, ndim
      nstride = (nfine/2)**(idim-1)
      igrid = igrid + nstride * MOD(hash_key(idim), nfine/2)
    end do

    ! If parent fmm cell changed, fetch (and unlock previous) neighbor info
    if (.not. all(hash_fmm_cell == prev_hash_fmm_cell)) then
      prev_hash_fmm_cell = hash_fmm_cell
      
      do ind = 1, threetondim
        do idim = 1, ndim
          offset(idim) = MOD((ind-1)/3**(idim-1), 3) - 1
        end do
        
        ! calculate neighboring fmm_cell's cartesian coordinate
        cc_fmm_cell = hash_fmm_cell(1:ndim) + offset
        i = 1
        j = 1
        k = 1
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
          jgrid = 1 + (k-1)*((nfine/2)**2) + (j-1)*(nfine/2) + (i-1)
          hash_direct(1) = (nfine/2) * cc_fmm_cell(1) + i - 1
#endif
          ! periodic boundary conditions & skipping out-of-box grids
          cycle_flag = .false.
          do idim = 1, ndim
#ifdef PERIODIC
            if (r%periodic(idim)) then
              if (hash_direct(idim) < m%box_ckey_min(idim, ilevel)) then
                hash_direct(idim) = m%box_ckey_max(idim, ilevel) - 1
              end if
              if (hash_direct(idim) >= m%box_ckey_max(idim, ilevel)) then
                hash_direct(idim) = m%box_ckey_min(idim, ilevel)
              end if
            end if
#endif
            if (hash_direct(idim) < m%box_ckey_min(idim, ilevel) .OR. &
                hash_direct(idim) >= m%box_ckey_max(idim, ilevel)) then
              cycle_flag = .true.
            end if
          end do

          if (cycle_flag) then
            mm_jcell_list(ind, jgrid, :) = 0.0d0
          else
            call get_grid(s, hash_direct, m%grid_dict, gridp_nbor, flush_cache = .false., fetch_cache = .true.)
            do jcell = 1, twotondim
#ifdef FMM
              mm_jcell_list(ind, jgrid, jcell) = gridp_nbor%rho(jcell)*dxn
#endif
            end do
          end if
#if NDIM>0
        end do
#endif
#if NDIM>1
        end do
#endif
#if NDIM>2
        end do
#endif
      end do
    end if

    ! Compute interactions for all cells in this AMR grid
    do icell = 1, twotondim
      phi = 0.0D0
      do ind = 1, threetondim
        do jgrid = 1, nbox
          do jcell = 1, twotondim
            phi = phi - mm_jcell_list(ind, jgrid, jcell) * inv_dist(igrid, icell, ind, jgrid, jcell)
          end do 
        end do 
      end do
      m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + phi
    end do
  end do ! end over all amr grids @ given ilevel
  deallocate(mm_jcell_list, inv_dist)
  
  call close_cache(s, m%grid_dict)
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
#ifdef FMM  
  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  grid%multipole=0.0D0
  grid%taylor_coeff=0.0D0
#endif
end subroutine init_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_taylor(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim,taylor_size
  use amr_commons, only: oct
  use cache_commons, only: msg_large_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
#ifdef FMM
  do ind=1,twotondim
    do ivar=1,taylor_size
      msg%realdp_fmm_taylor(ind, ivar)=grid%taylor_coeff(ind, ivar)
    end do
  end do
#endif
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
#ifdef FMM 
  do ind=1,twotondim
    do ivar=1,taylor_size
      grid%taylor_coeff(ind,ivar)=grid%taylor_coeff(ind,ivar)+msg%realdp_fmm_taylor(ind,ivar)
    end do
  end do
#endif
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
#ifdef FMM
  msg%realdp_fmm_multipole=grid%multipole
  msg%realdp_fmm_taylor=grid%taylor_coeff
#endif
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
#ifdef FMM
  grid%multipole=msg%realdp_fmm_multipole
  grid%taylor_coeff=msg%realdp_fmm_taylor
#endif
end subroutine unpack_fetch_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_rho(grid,msg_size,msg_array)
  use amr_parameters, only: twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind
  type(msg_small_realdp)::msg

#ifdef GRAV
  do ind=1,twotondim
     msg%realdp(ind)=grid%rho(ind)
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_rho
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_rho(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: oct
  use cache_commons, only: msg_small_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind
  type(msg_small_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

#ifdef GRAV
  do ind=1,twotondim
     grid%rho(ind)=msg%realdp(ind)
  end do
#endif

end subroutine unpack_fetch_rho
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
#ifdef FMM
  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
    write(unit_debug, '(3I6, 20E20.5)') m%grid(ioct)%ckey(1:ndim), m%grid(ioct)%taylor_coeff
  end do
#endif
  close(unit_debug)
end subroutine dump_taylor
#endif
end module fmm_fine_commons

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
subroutine fmm(pst,ilev,icount)
  use ramses_commons, only: pst_t
  use init_fmm_module, only: r_init_fmm, r_build_fmm, double_level_t, downward_level_t, &
       & FMM_BUILD_STANDARD, FMM_BUILD_MERGED, FMM_TREE_SOURCE_STANDARD, FMM_TREE_SOURCE_MERGED
  use fmm_multipoles, only: m_fmm_multipoles, merge_multipoles
  use cleanup_fmm_module, only: r_cleanup_fmm
  implicit none
  type(pst_t)::pst
  integer,intent(in) :: ilev,icount
  
  integer :: ifine, jlev, flev, input_size, jlev_max_amr_direct
  logical :: use_merged
  type(double_level_t)::double_level
  type(downward_level_t)::downward_levels

  if(pst%s%r%gravity_type>0)return
  if(pst%s%m%noct_tot(ilev)==0)return
  
  if(pst%s%r%verbose) print '(A,I2)','Entering fmm at AMR level ',ilev

  do jlev=ilev,pst%s%r%nlevelmax
    if(pst%s%r%verbose) print '(A,I2)','[Build FMM] ', jlev
    call r_cleanup_fmm(pst, jlev)
    call r_init_fmm(pst, jlev, 1)
    double_level%ilevel=jlev
    double_level%mode=FMM_BUILD_STANDARD
    do ifine=jlev,pst%s%r%bound_levelmin+1,-1
      double_level%ifine=ifine
      call r_build_fmm(pst,double_level,storage_size(double_level)/32)
    end do
  end do
  double_level%ilevel=ilev+1
  double_level%ifine=ilev
  double_level%mode=FMM_BUILD_MERGED
  call r_build_fmm(pst,double_level,storage_size(double_level)/32)

  if(pst%s%r%verbose) print '(A)','FMM Hierarchy done '

  do jlev=ilev,pst%s%r%nlevelmax
    call m_fmm_multipoles(pst, jlev) ! do upward pass !
  end do
  call merge_multipoles(pst,ilev+1)
  use_merged = associated(pst%s%m_fmm_merged) .and. (pst%s%m_fmm_merged%noct_used > 0) .and. &
       &       (ilev < pst%s%r%nlevelmax)

   ! Downward pass for fmm grids. 
   input_size = storage_size(downward_levels)/32
   downward_levels%ilev=ilev
   
   if(pst%s%r%verbose) print *, "[M2L & L2L] LEVEL: ", ilev
   downward_levels%mode = FMM_TREE_SOURCE_STANDARD
   do flev = pst%s%r%bound_levelmin+1, ilev-pst%s%r%level_fmm_to_amr
      downward_levels%flev=flev
      do jlev = max(pst%s%r%levelmin, flev+pst%s%r%level_fmm_to_amr-2), min(ilev, pst%s%r%nlevelmax)
         downward_levels%jlev=jlev
         call r_fmm_downward(pst, downward_levels, input_size)
         if(pst%s%r%verbose) print *,'     <Downpass> (ilev, jlev, flev): ', ilev, jlev, flev
      end do
      if (use_merged) then
         downward_levels%jlev=ilev+1
         downward_levels%mode = FMM_TREE_SOURCE_MERGED
         call r_fmm_downward(pst, downward_levels, input_size)
         if(pst%s%r%verbose) print *,'     <Downpass merged> (ilev, flev): ', ilev, flev
         downward_levels%mode = FMM_TREE_SOURCE_STANDARD
      end if
   end do

   ! Call direct force calculation
   downward_levels%flev=ilev-pst%s%r%level_fmm_to_amr

   !! L2P and M2P from ilev - 1 is done through combined_direct force. 
   !! ilev-2 should also be done via a similar function as combined_direct force 2. 
   if(pst%s%r%verbose) print *, "[L2P & M2P] LEVEL: ", ilev
   downward_levels%mode=FMM_TREE_SOURCE_STANDARD
   downward_levels%jlev=ilev
   call r_fmm_amr_intermediate(pst, downward_levels, input_size)
   if(pst%s%r%verbose) print *,'     <AMR Intermediate> (ilev, jlev)', ilev, ilev
   if (use_merged) then
      downward_levels%jlev=ilev+1
      downward_levels%mode=FMM_TREE_SOURCE_MERGED
      call r_fmm_amr_intermediate(pst, downward_levels, input_size)
      if(pst%s%r%verbose) print *,'     <AMR Intermediate merged> (ilev)', ilev
   else
      do jlev = max(pst%s%r%levelmin, ilev+1), pst%s%r%nlevelmax
         downward_levels%jlev=jlev
         downward_levels%mode=FMM_TREE_SOURCE_STANDARD
         call r_fmm_amr_intermediate(pst, downward_levels, input_size)
         if(pst%s%r%verbose) print *,'     <AMR Intermediate> (ilev, jlev)', ilev, jlev
      end do
   end if

   if(pst%s%r%verbose) print *, "[P2P] LEVEL: ", ilev
   downward_levels%mode=FMM_TREE_SOURCE_STANDARD
   jlev_max_amr_direct = min(pst%s%r%nlevelmax, ilev)
   do jlev = max(pst%s%r%levelmin, ilev-2), jlev_max_amr_direct
      downward_levels%jlev=jlev
      call r_fmm_amr_direct(pst, downward_levels, input_size)
      if(pst%s%r%verbose) print *,'     <Direct Force AMR> (ilev, jlev)', ilev, jlev
   end do
   if (use_merged) then
      if (ilev+1 <= pst%s%r%nlevelmax) then
         downward_levels%jlev=ilev+1
         downward_levels%mode=FMM_TREE_SOURCE_MERGED
         call r_fmm_amr_direct(pst, downward_levels, input_size)
         if(pst%s%r%verbose) print *,'     <Direct Force merged taylor> (ilev)', ilev
      end if
   else
      downward_levels%mode=FMM_TREE_SOURCE_STANDARD
      do jlev = max(pst%s%r%levelmin, ilev+1), pst%s%r%nlevelmax
         downward_levels%jlev=jlev
         call r_fmm_amr_direct(pst, downward_levels, input_size)
         if(pst%s%r%verbose) print *,'     <Direct Force> (ilev, jlev)', ilev, jlev
      end do
   end if

   !do i = 1, pst%s%r%levelmin - pst%s%r%level_fmm_to_amr
   !  call dump_taylor(pst%s%r, pst%s%m_fmm, i)
   !end do 
    
  ! ---------------------------------------------------------------------
  ! Cleanup MG levels after solve complete
  ! ---------------------------------------------------------------------
end subroutine fmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_fmm_downward(pst,downward_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_fmm_module, only: downward_level_t, FMM_TREE_SOURCE_MERGED
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  type(downward_level_t)::downward_levels
  integer,VALUE::input_size
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_DOWNWARD,pst%iUpper+1,input_size,0,downward_levels)
     call r_fmm_downward(pst%pLower,downward_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if (downward_levels%jlev == downward_levels%flev - 1) then
       call fmm_downward_coarse(pst%s, downward_levels%ilev, downward_levels%jlev, downward_levels%flev, &
            use_merged=(downward_levels%mode == FMM_TREE_SOURCE_MERGED))
     else
       call fmm_downward(pst%s, downward_levels%ilev, downward_levels%jlev, downward_levels%flev, &
            use_merged=(downward_levels%mode == FMM_TREE_SOURCE_MERGED))
     end if
  endif

end subroutine r_fmm_downward
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_downward(s, ilev, jlev, flev, use_merged)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size, taylor_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_taylor
  use fmm_multipoles, only: pack_fetch_multipole, unpack_fetch_multipole
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev, flev
  logical, optional, intent(in) :: use_merged

  integer :: ioct, idim, pcell, icell, inbor, jcell, nstride, iact
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_jcell_periodic, offset, ii
  real(kind=8), dimension(ndim) :: dx, diff
  real(kind=8) :: dx_loc, dist, D0, D1, D2, D3
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor_periodic, hash_parent, prev_hash_parent
  integer, dimension(1:threetondim) :: grid_nbors

  integer :: igrid_nbor, igrid_parent
  type(msg_large_realdp)::dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor
  real(kind=8), dimension(twotondim, taylor_size) :: accum_taylor
  logical::cycle_flag, neighbors_cached
  logical :: use_merged_local
  type(mesh_t), pointer :: m_target, m_source

  integer(kind=8), dimension(ndim, threetondim) :: offset_list
  real(kind=8), dimension(twotondim, twotondim, twotondim, threetondim) :: D0_list, D1_list, D2_list, D3_list
  real(kind=8), dimension(ndim, twotondim, twotondim, twotondim, threetondim) :: intermediate_diff_list
  real(kind=8), dimension(multipole_size, twotondim, threetondim) :: multipole_jcell_list
  logical, dimension(threetondim, twotondim, twotondim) :: direct_neighbor_list
  integer, dimension(threetondim) :: source_active_count
  integer, dimension(twotondim, threetondim) :: source_active_idx

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )
  use_merged_local = .false.
  if (present(use_merged)) use_merged_local = use_merged
  m_target => s%m_fmm_list(ilev)
  if (use_merged_local .and. associated(s%m_fmm_merged)) then
    m_source => s%m_fmm_merged
  else
    m_source => s%m_fmm_list(jlev)
  end if
  associate(r=>s%r, m=>s%m, mdl=>s%mdl)

  !if(m%noct_fmm(flev)<1) return
  
  ! Open cache for multipoles
  call open_cache(mdl, m_source, pack_size=storage_size(dummy_realdp)/32,& 
            pack=pack_fetch_taylor,unpack=unpack_fetch_taylor,& 
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  ! (debug logging removed)

  hash_key(0) = flev
  hash_nbor_periodic(0) = flev - 1
  hash_parent(0) = flev - 1
  prev_hash_parent(0) = flev - 1
  prev_hash_parent(1:ndim) = -1
  dx_loc = r%boxlen / 2.0D0**flev
  neighbors_cached = .false.

  ! jcell to icell
  do inbor = 1, threetondim
    ! calculate offsets
    do idim = 1, ndim
      offset_list(idim, inbor) = MOD((inbor-1)/3**(idim-1), 3) - 1 ! offset by how many fmm parent grids
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(:,inbor) + displacement_list(jcell,:) ! respect to parent grid left corner / fmm grid unit
      do pcell = 1,twotondim
        cycle_flag = .true.
        do idim=1, ndim
          if (abs(cc_jcell(idim) - displacement_list(pcell, idim)) > 1) cycle_flag = .false.
        end do
        direct_neighbor_list(inbor, jcell, pcell) = cycle_flag
        do icell = 1, twotondim
          cc_icell = 2 * displacement_list(pcell, :) + displacement_list(icell, :)
          diff = (cc_icell - 0.5 - 2 * cc_jcell) * dx_loc
          intermediate_diff_list(:, icell, pcell, jcell, inbor) = diff
          dist = sqrt(sum(diff(:)**2))
          D0_list(icell, pcell, jcell, inbor) = 1.0D0 / dist
          D1_list(icell, pcell, jcell, inbor) = -1.0D0 / dist**3
          D2_list(icell, pcell, jcell, inbor) = 3.0D0 / dist**5
          D3_list(icell, pcell, jcell, inbor) = -15.0D0 / dist**7
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m_target%head(flev), m_target%tail(flev)
     accum_taylor(:, :) = 0.0D0
     hash_key(1:ndim) = m_target%grid(ioct)%ckey(1:ndim)
     hash_parent(1:ndim) = hash_key(1:ndim)/2 !!! check

    ! Only do L2L if source and target trees are the same.
    if (ilev == jlev) then
      call get_parent_cell(s, hash_key, igrid_parent, pcell, flush_cache=.false., fetch_cache=.true.)
    
#ifdef FMM
      if (igrid_parent > 0) then
        parent_taylor = m_source%taylor_coeff(pcell, :, igrid_parent)
      else
        parent_taylor = 0.0
      end if
#endif

      ! Far Field Calculation
      do icell = 1, twotondim
        do idim =1,ndim
          dx(idim) = (displacement_list(icell, idim) - 0.5) * dx_loc
        end do 
        call shift_taylor(parent_taylor, dx, temp_taylor)
        accum_taylor(icell,:) = accum_taylor(icell,:) + temp_taylor
      end do
    else
      ii(1:ndim)=hash_key(1:ndim)-2*hash_parent(1:ndim)
      pcell=1
      do idim=1,ndim
        pcell=pcell+2**(idim-1)*ii(idim)
      end do
    end if

     if (.not. all(hash_parent == prev_hash_parent)) then
       if (neighbors_cached) then
         do inbor = 1, threetondim
           if (grid_nbors(inbor) > 0) call unlock_cache(m_source, grid_nbors(inbor))
         end do
       end if

       ! Get neighboring parent grids (returns 3^n parent level grids)
       call get_intermediate_nbor_grid(s, hash_key, grid_nbors, flush_cache=.false., fetch_cache=.true.)
       source_active_count(:) = 0
       do inbor = 1, threetondim
         igrid_nbor = grid_nbors(inbor)
         if (igrid_nbor <= 0) cycle
         do jcell = 1, twotondim
#ifdef FMM
           multipole_jcell_list(:, jcell, inbor) = m_source%multipole(jcell,:,igrid_nbor)
           if (any(multipole_jcell_list(:, jcell, inbor) /= 0.0d0)) then
             source_active_count(inbor) = source_active_count(inbor) + 1
             source_active_idx(source_active_count(inbor), inbor) = jcell
           end if
#endif
         end do
       end do
       neighbors_cached = .true.
       prev_hash_parent = hash_parent
     end if

     do inbor = 1, threetondim
       ! Get offset for neighboring parent grids (can reach off bounds)
       do idim=1,ndim
         offset(idim) = MOD((inbor-1)/3**(idim-1), 3) - 1
       end do

       igrid_nbor = grid_nbors(inbor)
       hash_nbor_periodic(1:ndim) = hash_parent(1:ndim) + offset

       if (igrid_nbor<=0) cycle

       do iact = 1, source_active_count(inbor)
          jcell = source_active_idx(iact, inbor)
          cycle_flag = .false.
          do idim=1,ndim
            nstride = 2**(idim-1)
            cc_jcell_periodic(idim) = 2*hash_nbor_periodic(idim) + MOD((jcell-1)/nstride, 2)
            if (cc_jcell_periodic(idim) < m%box_ckey_min(idim, flev) .or. &
                cc_jcell_periodic(idim) >= m%box_ckey_max(idim, flev)) then
              cycle_flag = .true.
            end if
          end do
          if (direct_neighbor_list(inbor, jcell, pcell) .or. cycle_flag) cycle
          ! Shift multipole from origin -> source center (Need to use grid position)
#ifdef FMM
          multipole = multipole_jcell_list(:, jcell, inbor)
#endif
          ! Get taylor coeffs from local
          do icell=1, twotondim
            dx = intermediate_diff_list(:, icell, pcell, jcell, inbor)
            D0 = D0_list(icell, pcell, jcell, inbor)
            D1 = D1_list(icell, pcell, jcell, inbor)
            D2 = D2_list(icell, pcell, jcell, inbor)
            D3 = D3_list(icell, pcell, jcell, inbor)
            call calc_taylor_from_multipole(dx, D0, D1, D2, D3, multipole, temp_taylor)
            accum_taylor(icell, :) = accum_taylor(icell, :) + temp_taylor
          end do
       end do ! over neighboring grid's cells 2^n
	     end do ! over neighboring grids 3^n 
	     ! Add taylor coefficients from intermediate fields
#ifdef FMM
		     m_target%taylor_coeff(:,:,ioct) = m_target%taylor_coeff(:,:,ioct) + accum_taylor
#endif
  end do
  if (neighbors_cached) then
    do inbor = 1, threetondim
      if (grid_nbors(inbor) > 0) call unlock_cache(m_source, grid_nbors(inbor))
    end do
  end if
  call close_cache(mdl)
  end associate
end subroutine fmm_downward
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fmm_downward_coarse(s, ilev, jlev, flev, use_merged)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size, taylor_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use fmm_multipoles, only: pack_fetch_rho, unpack_fetch_rho
  use cache_commons
  use cache
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev, flev
  logical, optional, intent(in) :: use_merged

  integer :: ioct, idim, pcell, icell, inbor, jcell, nstride
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_jcell_periodic, offset, ii
  real(kind=8), dimension(ndim) :: dx, diff
  real(kind=8) :: dx_loc, dist, D0, vol
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor_periodic, hash_parent
  integer, dimension(1:threetondim) :: grid_nbors

  integer :: igrid_nbor
  type(msg_int4_small_realdp)::dummy_rho
  real(kind=8), dimension(taylor_size) :: temp_taylor
  real(kind=8), dimension(twotondim, taylor_size) :: accum_taylor
  logical::cycle_flag
  logical :: use_merged_local
  type(mesh_t), pointer :: m_target, m_source

  integer(kind=8), dimension(ndim, threetondim) :: offset_list
  real(kind=8), dimension(twotondim, twotondim, twotondim, threetondim) :: D0_list
  real(kind=8), dimension(ndim, twotondim, twotondim, twotondim, threetondim) :: intermediate_diff_list
  logical, dimension(threetondim, twotondim, twotondim) :: direct_neighbor_list

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )
  use_merged_local = .false.
  if (present(use_merged)) use_merged_local = use_merged
  m_target => s%m_fmm_list(ilev)
  m_source => s%m
  associate(r=>s%r, m=>s%m, mdl=>s%mdl)
  
  ! Open cache for multipoles
  call open_cache(mdl, m_source, pack_size=storage_size(dummy_rho)/32, pack=pack_fetch_rho, unpack=unpack_fetch_rho)

  hash_key(0) = flev
  hash_nbor_periodic(0) = flev - 1
  hash_parent(0) = flev - 1
  dx_loc = r%boxlen / 2.0D0**flev
  vol = dx_loc ** ndim
  temp_taylor = 0.0D0

  ! jcell to icell
  do inbor = 1, threetondim
    ! calculate offsets
    do idim = 1, ndim
      offset_list(idim, inbor) = MOD((inbor-1)/3**(idim-1), 3) - 1 ! offset by how many fmm parent grids
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(:,inbor) + displacement_list(jcell,:) ! respect to parent grid left corner / fmm grid unit
      do pcell = 1,twotondim
        cycle_flag = .true.
        do idim=1, ndim
          if (abs(cc_jcell(idim) - displacement_list(pcell, idim)) > 1) cycle_flag = .false.
        end do
        direct_neighbor_list(inbor, jcell, pcell) = cycle_flag
        do icell = 1, twotondim
          cc_icell = 2 * displacement_list(pcell, :) + displacement_list(icell, :)
          diff = (cc_icell - 0.5 - 2 * cc_jcell) * dx_loc
          intermediate_diff_list(:, icell, pcell, jcell, inbor) = diff
          dist = sqrt(sum(diff(:)**2))
          D0_list(icell, pcell, jcell, inbor) = 1.0D0 / dist
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m_target%head(flev), m_target%tail(flev)
     accum_taylor(:, :) = 0.0D0
     hash_key(1:ndim) = m_target%grid(ioct)%ckey(1:ndim)
     hash_parent(1:ndim) = hash_key(1:ndim)/2
     ii(1:ndim)=hash_key(1:ndim)-2*hash_parent(1:ndim)
     pcell=1
     do idim=1,ndim
       pcell=pcell+2**(idim-1)*ii(idim)
     end do

     ! Get neighboring parent grids (returns 3^n parent level grids)
     call get_intermediate_nbor_grid(s, hash_key, grid_nbors, flush_cache=.false., fetch_cache=.true.)
     do inbor = 1, threetondim
       ! Get offset for neighboring parent grids (can reach off bounds)
       do idim=1,ndim
         offset(idim) = MOD((inbor-1)/3**(idim-1), 3) - 1
       end do

       igrid_nbor = grid_nbors(inbor)
       hash_nbor_periodic(1:ndim) = hash_parent(1:ndim) + offset

       if (igrid_nbor<=0) cycle
       ! If any has a refined cell, by definition there exists an FMM grid of size of the grid we called.
       if (any(m_source%grid(igrid_nbor)%refined(1:twotondim))) cycle

       do jcell = 1, twotondim
          cycle_flag = .false.
          do idim=1,ndim
            nstride = 2**(idim-1)
            cc_jcell_periodic(idim) = 2*hash_nbor_periodic(idim) + MOD((jcell-1)/nstride, 2)
            if (cc_jcell_periodic(idim) < m%box_ckey_min(idim, flev) .or. &
                cc_jcell_periodic(idim) >= m%box_ckey_max(idim, flev)) then
              cycle_flag = .true.
            end if
          end do
          if (cycle_flag) cycle
          ! skip direct neighbors
          if (direct_neighbor_list(inbor, jcell, pcell)) cycle
          ! Get taylor coeffs from local
          do icell=1, twotondim
            dx = intermediate_diff_list(:, icell, pcell, jcell, inbor)
            D0 = D0_list(icell, pcell, jcell, inbor)
            temp_taylor(1) = D0 * (m_source%rho(jcell,igrid_nbor) * vol)
            accum_taylor(icell, :) = accum_taylor(icell, :) + temp_taylor
          end do
       end do ! over neighboring grid's cells 2^n
	     end do ! over neighboring grids 3^n 
	     ! Add taylor coefficients from intermediate fields
#ifdef FMM
		     m_target%taylor_coeff(:,:,ioct) = m_target%taylor_coeff(:,:,ioct) + accum_taylor
#endif
	     ! Unlock neighbor octs
	     do inbor = 1, threetondim
	        call unlock_cache(m_source, grid_nbors(inbor))
	     end do
  end do
  call close_cache(mdl)
  end associate
end subroutine fmm_downward_coarse
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_amr_intermediate(pst,downward_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use init_fmm_module, only: downward_level_t, FMM_TREE_SOURCE_MERGED
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(downward_level_t)::downward_levels
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_AMR_INTERMEDIATE,pst%iUpper+1,input_size,0,downward_levels)
     call r_fmm_amr_intermediate(pst%pLower,downward_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call fmm_amr_intermediate(pst%s,downward_levels%ilev,downward_levels%jlev, &
          use_merged=(downward_levels%mode == FMM_TREE_SOURCE_MERGED))
  endif

end subroutine r_fmm_amr_intermediate
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_amr_intermediate(s, ilev, jlev, use_merged)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size, taylor_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_taylor
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev
  logical, optional, intent(in) :: use_merged

  integer :: ioct, idim, ind, icell, jcell, nstride, nfine, igrid, nbox, pcell, iact
  real(kind=8) :: phi, phi_out, dx_loc
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_jcell_periodic, offset
  real(kind=8), dimension(ndim) :: diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        prev_hash_fmm_grid

  integer, dimension(1:threetondim) :: grid_nbors
  integer :: igrid_nbor, igrid_parent
  type(msg_large_realdp) :: dummy_realdp
  real(kind=8), dimension(1:multipole_size) :: multipole
  real(kind=8), dimension(taylor_size) :: parent_taylor
  logical :: cycle_flag, neighbors_cached

  real(kind=8) :: dist, D0, D1, D2
  integer(kind=8), dimension(ndim, threetondim) :: offset_list
  logical :: use_merged_local
  type(mesh_t), pointer :: m_fmm

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
  real(kind=8), allocatable :: multipole_jcell_list(:,:,:)
  integer, allocatable :: fmm_grid_center_offset(:,:), fmm_cell_center_offset(:,:)
  integer, allocatable :: source_active_count(:), source_active_idx(:,:)
  logical, allocatable :: direct_neighbor_list(:,:,:)

  use_merged_local = .false.
  if (present(use_merged)) use_merged_local = use_merged
  if (use_merged_local .and. associated(s%m_fmm_merged)) then
    m_fmm => s%m_fmm_merged
  else
    m_fmm => s%m_fmm_list(jlev)
  end if

  associate(r=>s%r, m=>s%m, mdl=>s%mdl)


  ! Open cache for multipoles
  call open_cache(mdl, m_fmm, pack_size=storage_size(dummy_realdp)/32,&
                  pack=pack_fetch_taylor, unpack=unpack_fetch_taylor,&
                  init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilev
  hash_fmm_grid(0) = ilev - r%level_fmm_to_amr
  prev_hash_fmm_grid(0) = ilev - r%level_fmm_to_amr
  hash_fmm_cell(0) = ilev - r%level_fmm_to_amr + 1
  prev_hash_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilev
  nfine = 2**r%level_fmm_to_amr
  neighbors_cached = .false.

  ! nbox is the number of cells at the target AMR level
  nbox = nfine ** ndim

  !ind  : index of source fmm_grid within the neighboring 3^n fmm_grid
  !jcell: index of source fmm_cell within the neighboring fmm_grid
  !igrid: index of oct within target fmm_grid
  !icell: index of target cell within target fmm_grid

  allocate(D0_list(twotondim, nbox, twotondim, threetondim))
  allocate(D1_list(twotondim, nbox, twotondim, threetondim))
  allocate(D2_list(twotondim, nbox, twotondim, threetondim))

  allocate(intermediate_diff_list(ndim, twotondim, nbox, twotondim, threetondim))
  allocate(far_diff_list(ndim, twotondim, nbox))
  allocate(multipole_jcell_list(multipole_size, twotondim, threetondim))
  allocate(fmm_grid_center_offset(ndim, nbox))
  allocate(fmm_cell_center_offset(ndim, nbox))
  allocate(source_active_count(threetondim))
  allocate(source_active_idx(twotondim, threetondim))
  allocate(direct_neighbor_list(twotondim, nbox, threetondim))
  allocate(cell_diff_list(ndim, twotondim, nbox, threetondim))

  ! Precalculate differences
  do igrid=1, nbox
    do idim = 1,ndim
      nstride = nfine**(idim-1)
      fmm_grid_center_offset(idim, igrid) = MOD((igrid-1)/nstride, nfine) - (nfine/2) ! offset by how many amr octs from fmm grid center
      nstride = (nfine/2)**(idim-1)
      fmm_cell_center_offset(idim, igrid) = 2 * MOD(fmm_grid_center_offset(idim, igrid) + (nfine/2), nfine/2) - (nfine/2) ! offset by how many amr cells from fmm cell center
    end do 
    do icell = 1, twotondim
      far_diff_list(:, icell, igrid) = (fmm_cell_center_offset(:, igrid) + displacement_list(icell,:) + 0.5) * dx_loc
    end do
  end do 

  ! jcell to icell
  do ind = 1, threetondim
    ! calculate offsets
    do idim = 1, ndim
      offset_list(idim, ind) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(:,ind) + displacement_list(jcell,:) ! respect to grid left corner / fmm cell unit
      offset = (cc_jcell+ 0.5) * nfine ! offset by how many amr cells
      do igrid=1,nbox
        cc_icell = (fmm_grid_center_offset(:, igrid) + nfine/2)/(nfine/2) ! respect to grid left corner / fmm cell unit
        cycle_flag = .true.
        cell_diff_list(:, jcell, igrid, ind) = cc_icell - cc_jcell
        do idim=1, ndim
          if (abs(cc_icell(idim) - cc_jcell(idim)) > 1) cycle_flag = .false.
        end do
        direct_neighbor_list(jcell, igrid, ind) = cycle_flag
        do icell=1, twotondim
          diff = ((2 * fmm_grid_center_offset(:, igrid) + displacement_list(icell,:) + 0.5) + (- offset(:) + nfine)) * dx_loc
          intermediate_diff_list(:, icell, igrid, jcell, ind) = diff
          dist = sqrt(sum(diff(:)**2))
          D0_list(icell, igrid, jcell, ind) = 1.0D0 / dist
          D1_list(icell, igrid, jcell, ind) = -1.0D0 / dist**3
          D2_list(icell, igrid, jcell, ind) = 3.0D0 / dist**5
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilev), m%tail(ilev)

    if (all(m%grid(ioct)%refined(1:twotondim))) cycle

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
          call unlock_cache(m_fmm, grid_nbors(ind))
        end do
      end if

      if(ilev == jlev) call get_grid(s, hash_fmm_grid, igrid_parent, flush_cache=.false., fetch_cache=.true.)

      call get_intermediate_nbor_grid(s, hash_fmm_cell, grid_nbors, flush_cache=.false., fetch_cache=.true.)
      source_active_count(:) = 0
      do ind = 1, threetondim
        igrid_nbor = grid_nbors(ind)
        if (igrid_nbor <= 0) cycle
        do jcell = 1, twotondim
#ifdef FMM
          multipole_jcell_list(:, jcell, ind) = m_fmm%multipole(jcell, 1:multipole_size, igrid_nbor)
          if (any(multipole_jcell_list(:, jcell, ind) /= 0.0d0)) then
            source_active_count(ind) = source_active_count(ind) + 1
            source_active_idx(source_active_count(ind), ind) = jcell
          end if
#endif
        end do
      end do
      neighbors_cached = .true.
      prev_hash_fmm_grid = hash_fmm_grid
    end if

    !! FAR-FIELD
    if(ilev == jlev) then
      pcell = 1
      do idim=1,ndim
        nstride = 2**(idim-1)
        pcell = pcell + nstride * MOD(hash_fmm_cell(idim), 2)
      end do
#ifdef FMM
      if (igrid_parent > 0) then
        parent_taylor = m_fmm%taylor_coeff(pcell, :, igrid_parent)
      else
        parent_taylor = 0.0
      end if
#endif
      ! Far field
      do icell = 1, twotondim
        diff = far_diff_list(:, icell, igrid)
        call calc_phi(parent_taylor, diff, phi)
#ifdef FMM
        m%phi(icell, ioct) = m%phi(icell, ioct) + phi
#endif
      end do
    end if

    ! Intermediate field
    do ind = 1, threetondim
      igrid_nbor = grid_nbors(ind)
      if (igrid_nbor<=0) cycle

      do iact = 1, source_active_count(ind)
        jcell = source_active_idx(iact, ind)
        cycle_flag = .false.
        cc_jcell_periodic = hash_fmm_cell(1:ndim) - cell_diff_list(:, jcell, igrid, ind)
        do idim = 1, ndim
          if ((cc_jcell_periodic(idim) < m%box_ckey_min(idim, ilev - r%level_fmm_to_amr + 1)) .or. &
              (cc_jcell_periodic(idim) >= m%box_ckey_max(idim, ilev - r%level_fmm_to_amr + 1))) then
            cycle_flag = .true.
          end if
        end do
        if (direct_neighbor_list(jcell, igrid, ind) .or. cycle_flag) cycle
#ifdef FMM
        multipole = multipole_jcell_list(:, jcell, ind)
#endif
        do icell=1, twotondim
          if (m%grid(ioct)%refined(icell)) cycle
          diff  = intermediate_diff_list(:, icell, igrid, jcell, ind)
          D0 = D0_list(icell, igrid, jcell, ind)
          D1 = D1_list(icell, igrid, jcell, ind)
          D2 = D2_list(icell, igrid, jcell, ind)
          call calc_phi_from_multipole(diff, D0, D1, D2, multipole, phi_out)
          m%phi(icell, ioct) = m%phi(icell, ioct) + phi_out
        end do
      end do
    end do
  end do

  deallocate(D0_list, D1_list, D2_list, intermediate_diff_list, far_diff_list, multipole_jcell_list, &
       &     direct_neighbor_list, cell_diff_list, fmm_grid_center_offset, fmm_cell_center_offset, &
       &     source_active_count, source_active_idx)
  call close_cache(mdl)
  end associate
end subroutine fmm_amr_intermediate
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_fmm_amr_direct(pst,downward_levels,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use init_fmm_module, only: downward_level_t, FMM_TREE_SOURCE_MERGED
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(downward_level_t)::downward_levels
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_FMM_AMR_DIRECT,pst%iUpper+1,input_size,0,downward_levels)
     call r_fmm_amr_direct(pst%pLower,downward_levels,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     if (downward_levels%ilev-2==downward_levels%jlev) then
       call fmm_direct_coarsest(pst%s,downward_levels%ilev,downward_levels%jlev)
     else if (downward_levels%ilev-1==downward_levels%jlev) then
       !! HERE WE DO NEAR FIELD + MID FIELD TOGETHER WITH 6^n cells
       call fmm_combined_direct(pst%s,downward_levels%ilev,downward_levels%jlev)
     else if (downward_levels%ilev==downward_levels%jlev) then
       call fmm_amr_direct(pst%s,downward_levels%ilev,downward_levels%jlev)
     else
       call fmm_amr_direct_taylor(pst%s,downward_levels%ilev,downward_levels%jlev, &
            use_merged=(downward_levels%mode == FMM_TREE_SOURCE_MERGED)) ! important to exclude nearest neighbors.
     end if
  endif

end subroutine r_fmm_amr_direct
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_direct_coarsest(s, ilev, jlev)
  use amr_parameters, only: ndim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  use nbors_utils
  use fmm_multipoles, only: pack_fetch_rho, unpack_fetch_rho
  use cache_commons
  use cache
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev

  integer :: ioct, idim, ind, icell, ipcell, jcell, ifinecell, nstride
  real(kind=8) :: dx_loc, vol
  integer(kind=8), dimension(ndim) :: cc_jcell, cc_ifinecell, cc_jcell_periodic
  real(kind=8), dimension(ndim) :: diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, hash_fmm_pgrid, &
                                        prev_hash_fmm_grid

  integer, dimension(1:threetondim) :: grid_nbors, ind_nbors
  integer :: igrid_nbor
  type(msg_int4_small_realdp) :: dummy_rho
  logical :: cycle_flag, neighbors_cached

  real(kind=8) :: dist, D0
  integer(kind=8), dimension(ndim, threetondim) :: offset_list

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )

  ! =====================================================
  ! Runtime-allocated arrays depending on r%level_fmm_to_amr
  ! =====================================================
  real(kind=8), allocatable :: D0_list(:,:,:,:,:)
  real(kind=8), allocatable :: cell_diff_list(:,:,:,:,:,:)

  associate(r=>s%r, m=>s%m, mdl=>s%mdl)

  ! (debug logging removed)

  ! Open cache for multipoles
  call open_cache(mdl, m, pack_size=storage_size(dummy_rho)/32,&
                  pack=pack_fetch_rho, unpack=unpack_fetch_rho)

  hash_key(0) = ilev
  hash_fmm_grid(0)  = ilev - 1
  hash_fmm_pgrid(0) = ilev - 2
  prev_hash_fmm_grid(0) = ilev - 1
  hash_fmm_cell(0) = ilev
  prev_hash_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilev
  vol    = (4 * dx_loc ** ndim) * twotondim ! two-level difference

  neighbors_cached = .false.
  !ind  : index of source fmm_grid within the neighboring 3^n fmm_grid
  !jcell: index of source cell within the neighboring fmm_grid
  !igrid: index of oct within target fmm_grid
  !icell: index of target cell within target fmm_grid

  allocate(D0_list(twotondim, twotondim, twotondim, twotondim, threetondim))
  allocate(cell_diff_list(ndim, twotondim, twotondim, twotondim, twotondim, threetondim))

  ! jcell to icell
  do ind = 1, threetondim
    do idim = 1, ndim
      offset_list(idim, ind) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm pgrids (-1, 0, 1)
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(:,ind) + displacement_list(jcell,:) ! respect to pgrid left corner / fmm grid unit (-2, -1, 0, 1, 2, 3)
      cc_jcell = 4 * cc_jcell ! offset by how many amr cells respect to grid left corner (-6, -2, 2, 6, 10, 14)
      do ipcell=1,twotondim
        do icell=1,twotondim
          do ifinecell=1, twotondim
            cc_ifinecell = 2 * (2 * displacement_list(ipcell, :) + displacement_list(icell, :)) + displacement_list(ifinecell,:)
            diff = (cc_jcell - cc_ifinecell) * dx_loc
            cell_diff_list(:, ifinecell, icell, ipcell, jcell, ind) = diff
            dist = sqrt(sum(diff(:)**2))
            if (dist == 0.0D0) then
              D0_list(ifinecell, icell, ipcell, jcell, ind) = 0.0D0
            else
              D0_list(ifinecell, icell, ipcell, jcell, ind) = 1.0D0 / dist
            end if
          end do
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilev), m%tail(ilev)

    if (all(m%grid(ioct)%refined(1:twotondim))) cycle

    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / 2
    hash_fmm_pgrid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / 4
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim)

    ipcell = 1
    do idim=1,ndim
      nstride = 2**(idim-1)
      ipcell = ipcell + nstride * MOD(hash_fmm_grid(idim), 2)
    end do

    icell = 1
    do idim=1,ndim
      nstride = 2**(idim-1)
      icell = icell + nstride * MOD(hash_fmm_cell(idim), 2)
    end do

    ! Check if we need to fetch neighbors & parent Taylor
    if (.not. all(hash_fmm_grid == prev_hash_fmm_grid)) then
      if (neighbors_cached) then
        do ind = 1, threetondim
          call unlock_cache(m, grid_nbors(ind))
        end do
      end if

      call get_threetondim_nbor_parent_cell(s, hash_fmm_grid, grid_nbors, ind_nbors, flush_cache=.false., fetch_cache=.true.)
      neighbors_cached = .true.
      prev_hash_fmm_grid = hash_fmm_grid
    end if

    ! Intermediate field & Near Field
    do ind = 1, threetondim
      igrid_nbor = grid_nbors(ind)
      if (igrid_nbor<=0) cycle
      jcell = ind_nbors(ind)
        cycle_flag = .false.
        do idim=1,ndim
          nstride = 2**(idim-1)
          cc_jcell_periodic(idim) = 2*m%grid(igrid_nbor)%ckey(idim) + MOD((jcell-1)/nstride, 2)
          if (cc_jcell_periodic(idim) < m%box_ckey_min(idim, ilev-1) .or. &
              cc_jcell_periodic(idim) >= m%box_ckey_max(idim, ilev-1)) then
            cycle_flag = .true.
        end if
      end do
      if (cycle_flag) cycle
      if (m%grid(igrid_nbor)%refined(jcell)) cycle
      do ifinecell=1, twotondim
        if (m%grid(ioct)%refined(ifinecell)) cycle
        D0 = D0_list(ifinecell, icell, ipcell, jcell, ind)
        m%phi(ifinecell, ioct) = m%phi(ifinecell, ioct) - D0 * m%rho(jcell,igrid_nbor) * vol
      end do
    end do
  end do

  deallocate(D0_list, cell_diff_list)
  call close_cache(mdl)
  end associate
end subroutine fmm_direct_coarsest
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_combined_direct(s, ilev, jlev)
  use amr_parameters, only: ndim, twotondim, threetondim
  use ramses_commons, only: ramses_t
  use nbors_utils
  use fmm_multipoles, only: pack_fetch_rho, unpack_fetch_rho
  use cache_commons
  use cache
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev

  integer :: ioct, idim, ind, icell, jcell, ifinecell, nstride
  real(kind=8) :: dx_loc, vol
  integer(kind=8), dimension(ndim) :: cc_jcell, cc_jcell_periodic, offset
  real(kind=8), dimension(ndim) :: diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        prev_hash_fmm_grid

  integer, dimension(1:threetondim) :: grid_nbors
  integer :: igrid_nbor
  type(msg_int4_small_realdp) :: dummy_rho
  logical :: cycle_flag, neighbors_cached

  real(kind=8) :: dist, D0
  integer(kind=8), dimension(ndim, threetondim) :: offset_list

  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
      [ &
        0, 1, 0, 1, 0, 1, 0, 1,  &
        0, 0, 1, 1, 0, 0, 1, 1,  &
        0, 0, 0, 0, 1, 1, 1, 1   &
      ], [twotondim, ndim] )

  ! =====================================================
  ! Runtime-allocated arrays depending on r%level_fmm_to_amr
  ! =====================================================
  real(kind=8), allocatable :: D0_list(:,:,:,:)
  real(kind=8), allocatable :: intermediate_diff_list(:,:,:,:,:), cell_diff_list(:,:,:,:)
  real(kind=8), allocatable :: far_diff_list(:,:,:)
  associate(r=>s%r, m=>s%m, mdl=>s%mdl)


  ! Open cache for multipoles
  call open_cache(mdl, m, pack_size=storage_size(dummy_rho)/32,&
                  pack=pack_fetch_rho, unpack=unpack_fetch_rho)

  hash_key(0) = ilev
  hash_fmm_grid(0) = ilev - 1
  prev_hash_fmm_grid(0) = ilev - 1
  hash_fmm_cell(0) = ilev
  prev_hash_fmm_grid(1:ndim) = -1 ! initialize

  dx_loc = r%boxlen / 2.0D0**ilev
  vol    = (dx_loc ** ndim) * twotondim

  neighbors_cached = .false.
  !ind  : index of source fmm_grid within the neighboring 3^n fmm_grid
  !jcell: index of source cell within the neighboring fmm_grid
  !igrid: index of oct within target fmm_grid
  !icell: index of target cell within target fmm_grid

  allocate(D0_list(twotondim, twotondim, twotondim, threetondim))

  allocate(intermediate_diff_list(ndim, twotondim, twotondim, twotondim, threetondim))
  allocate(cell_diff_list(threetondim, twotondim, twotondim, ndim))

  ! jcell to icell
  do ind = 1, threetondim
    do idim = 1, ndim
      offset_list(idim, ind) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids (-1, 0, 1)
    end do
    do jcell = 1, twotondim
      cc_jcell = 2 * offset_list(:,ind) + displacement_list(jcell,:) ! respect to grid left corner / fmm cell unit (-2, -1, 0, 1, 2, 3)
      offset   = 2 * (cc_jcell+ 0.5) ! offset by how many amr cells respect to grid left corner (-3, -1, 1, 3, 5, 7)
      do icell=1,twotondim
        cell_diff_list(ind, jcell, icell, :) = displacement_list(icell,:) - offset !(-7, -6, .., 4)
        do ifinecell=1, twotondim
          diff = (2 * displacement_list(icell,:) + displacement_list(ifinecell,:) - offset + 0.5) * dx_loc
          intermediate_diff_list(:, ifinecell, icell, jcell, ind) = diff
          dist = sqrt(sum(diff(:)**2))
          D0_list(ifinecell, icell, jcell, ind) = 1.0D0 / dist
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilev), m%tail(ilev)

    if (all(m%grid(ioct)%refined(1:twotondim))) cycle

    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)
    hash_fmm_grid(1:ndim) = m%grid(ioct)%ckey(1:ndim) / 2
    hash_fmm_cell(1:ndim) = m%grid(ioct)%ckey(1:ndim)

    icell = 1
    do idim=1,ndim
      nstride = 2**(idim-1)
      icell = icell + nstride * MOD(hash_key(idim), 2)
    end do

    ! Check if we need to fetch neighbors & parent Taylor
    if (.not. all(hash_fmm_grid == prev_hash_fmm_grid)) then
      if (neighbors_cached) then
        do ind = 1, threetondim
          call unlock_cache(m, grid_nbors(ind))
        end do
      end if

      call get_threetondim_nbor_grid(s, hash_fmm_grid, grid_nbors, flush_cache=.false., fetch_cache=.true.)
      neighbors_cached = .true.
      prev_hash_fmm_grid = hash_fmm_grid
    end if

    ! Intermediate field & Near Field
    do ind = 1, threetondim
      igrid_nbor = grid_nbors(ind)
      if (igrid_nbor<=0) cycle

      do jcell = 1, twotondim
        cycle_flag = .false.
        do idim=1,ndim
          nstride = 2**(idim-1)
          cc_jcell_periodic(idim) = 2*m%grid(igrid_nbor)%ckey(idim) + MOD((jcell-1)/nstride, 2)
          if (cc_jcell_periodic(idim) < m%box_ckey_min(idim, ilev) .or. &
              cc_jcell_periodic(idim) >= m%box_ckey_max(idim, ilev)) then
            cycle_flag = .true.
          end if
        end do
        if (cycle_flag) cycle
        if (m%grid(igrid_nbor)%refined(jcell)) cycle

        do ifinecell=1, twotondim
          if (m%grid(ioct)%refined(ifinecell)) cycle
          D0 = D0_list(ifinecell, icell, jcell, ind)
          m%phi(ifinecell, ioct) = m%phi(ifinecell, ioct) - D0 * m%rho(jcell,igrid_nbor) * vol
        end do
      end do
    end do
  end do

  deallocate(D0_list, intermediate_diff_list, cell_diff_list)
  call close_cache(mdl)
  end associate
end subroutine fmm_combined_direct
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fmm_amr_direct(s, ilev, jlev)
  use amr_parameters, only: ndim, twotondim, threetondim, nhilbert
  use ramses_commons, only: ramses_t
  use nbors_utils
  use fmm_multipoles, only: pack_fetch_rho, unpack_fetch_rho
  use cache_commons
  use cache
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev

  integer :: ioct, idim, ind, icell, jcell, nstride, nfine, nbox, jgrid, igrid, jfinecell
  integer :: i, j, k, iact, inear, nnear
  real(kind=8) :: phi, dx_loc, dxn
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_igrid, cc_jgrid, cc_fmm_cell, offset
  real(kind=8), dimension(ndim) :: diff, fine_diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_direct, prev_hash_fmm_cell, hash_fine
  integer(kind=8), dimension(ndim, threetondim) :: offset_list
  integer :: igrid_nbor, igrid_fine
  type(msg_int4_small_realdp) :: dummy_rho
  logical :: cycle_flag
  logical :: source_all_refined
  logical, dimension(twotondim) :: refined_target, refined_source
  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
    [ &
      0, 1, 0, 1, 0, 1, 0, 1,  &
      0, 0, 1, 1, 0, 0, 1, 1,  &
      0, 0, 0, 0, 1, 1, 1, 1   &
    ], [twotondim, ndim] )

  ! Arrays sized for all source cells: threetondim * (nfine/2)^ndim * twotondim
  real(kind=8), dimension(:,:,:), allocatable      :: mm_jcell_list
  real(kind=8), dimension(:,:,:,:), allocatable    :: mm_jfinecell_list
  real(kind=8), dimension(:,:,:,:,:), allocatable  :: inv_dist
  real(kind=8), dimension(:,:,:,:,:,:), allocatable:: nearest_inv_dist
  integer, dimension(:,:,:,:), allocatable         :: nearest_count
  integer, dimension(:,:,:,:,:), allocatable       :: nearest_idx
  integer, dimension(:), allocatable               :: unrefined_count, refined_count
  integer, dimension(:,:), allocatable             :: unrefined_idx, refined_idx

  associate(r=>s%r, m=>s%m, mdl=>s%mdl)


  ! Open cache for multipoles
  call open_cache(mdl, m, pack_size=storage_size(dummy_rho)/32,&
            pack=pack_fetch_rho, unpack=unpack_fetch_rho)

  hash_key(0) = ilev
  hash_fmm_grid(0) = ilev - r%level_fmm_to_amr
  hash_fmm_cell(0) = ilev - r%level_fmm_to_amr + 1
  hash_direct(0) = ilev
  hash_fine(0) = ilev + 1

  dx_loc = r%boxlen / 2.0D0**ilev
  dxn    = dx_loc**ndim
  nfine  = 2**r%level_fmm_to_amr
  nbox   = (nfine/2)**ndim
  ! Allocate arrays for all possible source cells
  allocate(mm_jcell_list(twotondim, nbox, threetondim))
  allocate(mm_jfinecell_list(twotondim, twotondim, nbox, threetondim))
  allocate(inv_dist(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(nearest_inv_dist(twotondim, twotondim, nbox, threetondim, twotondim, nbox))
  allocate(nearest_count(nbox, threetondim, twotondim, nbox))
  allocate(nearest_idx(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(unrefined_count(threetondim), refined_count(threetondim))
  allocate(unrefined_idx(nbox, threetondim), refined_idx(nbox, threetondim))

  prev_hash_fmm_cell = -huge(0_8)
  nearest_count = 0
  nearest_inv_dist = 0.0D0

  ! Get offset lists
  do ind = 1, threetondim
    do idim = 1, ndim
      offset_list(idim, ind) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids
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
            cc_jgrid(idim) = offset_list(idim, ind) * (nfine/2) + MOD((jgrid-1)/(nfine/2)**(idim-1), nfine/2)
          end do
          do jcell = 1, twotondim
            cc_jcell = 2 * cc_jgrid + displacement_list(jcell, :)
            if (all(cc_icell(1:ndim) == cc_jcell(1:ndim))) then
              inv_dist(jcell, jgrid, ind, icell, igrid) = 1 / dx_loc !! cap it to dx_loc instead of skipping
            else
              diff = (cc_icell - cc_jcell) * dx_loc
#if NDIM==3
              inv_dist(jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3))
#elif NDIM==2
              inv_dist(jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(diff(1)*diff(1) + diff(2)*diff(2))
#elif NDIM==1
              inv_dist(jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(diff(1)*diff(1))
#else
              inv_dist(jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(sum(diff(:)**2))
#endif
              if (is_direct_neighbor(cc_icell, cc_jcell)) then
                nearest_count(jgrid, ind, icell, igrid) = nearest_count(jgrid, ind, icell, igrid) + 1
                nearest_idx(nearest_count(jgrid, ind, icell, igrid), jgrid, ind, icell, igrid) = jcell
                do jfinecell = 1, twotondim
                  fine_diff = diff + (0.5 * (displacement_list(jfinecell, :) - 0.5)) * dx_loc
#if NDIM==3
                  nearest_inv_dist(jfinecell, jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(fine_diff(1)*fine_diff(1) + fine_diff(2)*fine_diff(2) + fine_diff(3)*fine_diff(3))
#elif NDIM==2
                  nearest_inv_dist(jfinecell, jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(fine_diff(1)*fine_diff(1) + fine_diff(2)*fine_diff(2))
#elif NDIM==1
                  nearest_inv_dist(jfinecell, jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(fine_diff(1)*fine_diff(1))
#else
                  nearest_inv_dist(jfinecell, jcell, jgrid, ind, icell, igrid) = 1.d0 / sqrt(sum(fine_diff(:)**2))
#endif
                end do
              end if
            end if
          end do
        end do
      end do
    end do
  end do

  ! Loop over octs at this level
  do ioct = m%head(ilev), m%tail(ilev)
    refined_target(:) = m%grid(ioct)%refined(1:twotondim)
    if (all(refined_target)) cycle

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
        unrefined_count(ind) = 0
        refined_count(ind) = 0

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
          cycle_flag = .false.
          do idim = 1, ndim
            if (r%periodic(idim)) then
              if (hash_direct(idim) < m%box_ckey_min(idim, ilev)) then
                hash_direct(idim) = m%box_ckey_max(idim, ilev) - 1
              end if
              if (hash_direct(idim) >= m%box_ckey_max(idim, ilev)) then
                hash_direct(idim) = m%box_ckey_min(idim, ilev)
              end if
            end if
            if (hash_direct(idim) < m%box_ckey_min(idim, ilev) .or. &
                hash_direct(idim) >= m%box_ckey_max(idim, ilev)) then
              cycle_flag = .true.
            end if
          end do
          if (cycle_flag) then
            mm_jcell_list(:, jgrid, ind) = 0.0d0
            cycle
          end if
          call get_grid(s, hash_direct, igrid_nbor, flush_cache = .false., fetch_cache = .true.)
          if (igrid_nbor .le. 0) then 
            mm_jcell_list(:, jgrid, ind) = 0.0d0
            cycle
          end if
          refined_source(:) = m%grid(igrid_nbor)%refined(1:twotondim)
          source_all_refined = all(refined_source)
          if (source_all_refined) then
            do jcell = 1, twotondim
              hash_fine(1:ndim) = 2 * hash_direct(1:ndim) + displacement_list(jcell, :)
              call get_grid(s, hash_fine, igrid_fine, flush_cache = .false., fetch_cache = .true.)
              if (igrid_fine .le. 0) then
                mm_jfinecell_list(:, jcell, jgrid, ind) = 0.0d0
              else
                mm_jfinecell_list(:, jcell, jgrid, ind) = m%rho(:,igrid_fine)*dxn/8
              end if
            end do
          end if
          mm_jcell_list(:, jgrid, ind) = m%rho(:,igrid_nbor)*dxn
          if (source_all_refined) then
            refined_count(ind) = refined_count(ind) + 1
            refined_idx(refined_count(ind), ind) = jgrid
          else
            unrefined_count(ind) = unrefined_count(ind) + 1
            unrefined_idx(unrefined_count(ind), ind) = jgrid
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
      if (refined_target(icell)) cycle
      phi = 0.0D0
      do ind = 1, threetondim
        do iact = 1, unrefined_count(ind)
          jgrid = unrefined_idx(iact, ind)
          do jcell = 1, twotondim
            phi = phi - mm_jcell_list(jcell, jgrid, ind) * inv_dist(jcell, jgrid, ind, icell, igrid)
          end do 
        end do 
        do iact = 1, refined_count(ind)
          jgrid = refined_idx(iact, ind)
          nnear = nearest_count(jgrid, ind, icell, igrid)
          do inear = 1, nnear
            jcell = nearest_idx(inear, jgrid, ind, icell, igrid)
            do jfinecell = 1, twotondim
              phi = phi - mm_jfinecell_list(jfinecell, jcell, jgrid, ind) * &
                          nearest_inv_dist(jfinecell, jcell, jgrid, ind, icell, igrid)
            end do
          end do
        end do
      end do
      m%phi(icell, ioct) = m%phi(icell, ioct) + phi
    end do
  end do ! end over all amr grids @ given ilev
  deallocate(mm_jcell_list, mm_jfinecell_list, inv_dist, nearest_inv_dist, nearest_count, nearest_idx, &
             unrefined_count, refined_count, unrefined_idx, refined_idx)
  call close_cache(mdl)
  end associate
end subroutine fmm_amr_direct
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fmm_amr_direct_taylor(s, ilev, jlev, use_merged)
  use amr_parameters, only: ndim, twotondim, threetondim, multipole_size
  use amr_commons, only: mesh_t
  use ramses_commons, only: ramses_t
  use nbors_utils
  use cache_commons
  use cache
  use fmm_taylor
  use fmm_multipoles, only: pack_fetch_multipole, unpack_fetch_multipole
  implicit none

  type(ramses_t) :: s
  integer :: ilev, jlev
  logical, optional, intent(in) :: use_merged

  integer :: ioct, idim, ind, icell, jcell, nstride, nfine, nbox, jgrid, igrid, iact, nact, ngrid, jcell_act
  integer :: i, j, k
  real(kind=8) :: phi, phi_out, dx_loc, dist, D0, D1, D2
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_igrid, cc_jgrid, cc_fmm_cell, offset
  real(kind=8), dimension(ndim) :: diff
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_fmm_grid, hash_fmm_cell, &
                                        hash_direct, prev_hash_fmm_cell
  integer(kind=8), dimension(ndim, threetondim) :: offset_list
  real(kind=8), dimension(1:multipole_size) :: multipole
  integer :: igrid_nbor
  type(msg_large_realdp) :: dummy_realdp
  logical :: cycle_flag
  logical, dimension(twotondim) :: refined_target
  integer, dimension(twotondim, ndim), parameter :: displacement_list = reshape( &
    [ &
      0, 1, 0, 1, 0, 1, 0, 1,  &
      0, 0, 1, 1, 0, 0, 1, 1,  &
      0, 0, 0, 0, 1, 1, 1, 1   &
    ], [twotondim, ndim] )

  ! Arrays sized for all source cells: threetondim * (nfine/2)^ndim * twotondim
  real(kind=8), dimension(:,:,:,:), allocatable    :: multipole_jcell_list
  real(kind=8), dimension(:,:,:,:,:), allocatable  :: D0_list, D1_list, D2_list
  real(kind=8), dimension(:,:,:,:,:,:), allocatable:: diff_list
  logical, dimension(:,:,:,:,:), allocatable       :: nearest_flags
  integer, dimension(:,:), allocatable             :: source_active_count
  integer, dimension(:,:,:), allocatable           :: source_active_idx
  integer, dimension(:), allocatable               :: active_grid_count
  integer, dimension(:,:), allocatable             :: active_grid_idx
  logical :: use_merged_local
  type(mesh_t), pointer :: m_source

  use_merged_local = .false.
  if (present(use_merged)) use_merged_local = use_merged
  if (use_merged_local .and. associated(s%m_fmm_merged)) then
    m_source => s%m_fmm_merged
  else
    m_source => s%m_fmm_list(jlev)
  end if

  associate(r=>s%r, m=>s%m, mdl=>s%mdl)


  ! Open cache for multipoles
  call open_cache(mdl, m_source, pack_size=storage_size(dummy_realdp)/32, pack=pack_fetch_multipole,unpack=unpack_fetch_multipole)

  hash_key(0) = ilev
  hash_fmm_grid(0) = ilev - r%level_fmm_to_amr
  hash_fmm_cell(0) = ilev - r%level_fmm_to_amr + 1
  hash_direct(0) = ilev 

  dx_loc = r%boxlen / 2.0D0**ilev
  nfine  = 2**r%level_fmm_to_amr
  nbox   = (nfine/2)**ndim

  ! Allocate arrays for all possible source cells
  allocate(multipole_jcell_list(multipole_size, twotondim, nbox, threetondim))
  allocate(D0_list(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(D1_list(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(D2_list(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(diff_list(ndim, twotondim, nbox, threetondim, twotondim, nbox))
  allocate(nearest_flags(twotondim, nbox, threetondim, twotondim, nbox))
  allocate(source_active_count(nbox, threetondim))
  allocate(source_active_idx(twotondim, nbox, threetondim))
  allocate(active_grid_count(threetondim))
  allocate(active_grid_idx(nbox, threetondim))

  prev_hash_fmm_cell = -huge(0_8)
  diff_list = 0.0D0
  source_active_count = 0
  active_grid_count = 0

  ! Get offset lists
  do ind = 1, threetondim
    do idim = 1, ndim
      offset_list(idim, ind) = MOD((ind-1)/3**(idim-1), 3) - 1 ! offset by how many fmm grids
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
            cc_jgrid(idim) = offset_list(idim, ind) * (nfine/2) + MOD((jgrid-1)/(nfine/2)**(idim-1), nfine/2)
          end do
          do jcell = 1, twotondim
            cc_jcell = 2 * cc_jgrid + displacement_list(jcell, :)
            if (all(cc_icell(1:ndim) == cc_jcell(1:ndim))) then
              ! Skip self-interaction to avoid dist=0
              nearest_flags(jcell, jgrid, ind, icell, igrid) = .true.
              diff_list(:, jcell, jgrid, ind, icell, igrid) = 0.0D0
              D0_list(jcell, jgrid, ind, icell, igrid) = 0.0D0
              D1_list(jcell, jgrid, ind, icell, igrid) = 0.0D0
              D2_list(jcell, jgrid, ind, icell, igrid) = 0.0D0
            else
              diff = (cc_icell - cc_jcell) * dx_loc
              diff_list(:, jcell, jgrid, ind, icell, igrid) = diff
#if NDIM==3
              dist = sqrt(diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3))
#elif NDIM==2
              dist = sqrt(diff(1)*diff(1) + diff(2)*diff(2))
#elif NDIM==1
              dist = sqrt(diff(1)*diff(1))
#else
              dist = sqrt(sum(diff(:)**2))
#endif
              D0_list(jcell, jgrid, ind, icell, igrid) = 1.0D0 / dist
              D1_list(jcell, jgrid, ind, icell, igrid) = -1.0D0 / dist**3
              D2_list(jcell, jgrid, ind, icell, igrid) = 3.0D0 / dist**5
              if (is_direct_neighbor(cc_icell, cc_jcell)) then
                nearest_flags(jcell, jgrid, ind, icell, igrid) = .true.
              else
                nearest_flags(jcell, jgrid, ind, icell, igrid) = .false.
              end if
            end if
          end do
        end do
      end do
    end do
  end do
  
  ! Loop over octs at this level
  do ioct = m%head(ilev), m%tail(ilev)
    refined_target(:) = m%grid(ioct)%refined(1:twotondim)
    if (all(refined_target)) cycle

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
        active_grid_count(ind) = 0
        source_active_count(:, ind) = 0
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
          cycle_flag = .false.
          do idim = 1, ndim
            if (r%periodic(idim)) then
              if (hash_direct(idim) < m%box_ckey_min(idim, ilev)) then
                hash_direct(idim) = m%box_ckey_max(idim, ilev) - 1
              end if
              if (hash_direct(idim) >= m%box_ckey_max(idim, ilev)) then
                hash_direct(idim) = m%box_ckey_min(idim, ilev)
              end if
            end if
            if (hash_direct(idim) < m%box_ckey_min(idim, ilev) .or. &
                hash_direct(idim) >= m%box_ckey_max(idim, ilev)) then
              cycle_flag = .true.
            end if
          end do
          if (cycle_flag) then
            cycle
          end if
          call get_grid(s, hash_direct, igrid_nbor, flush_cache = .false., fetch_cache = .true.)
          if (igrid_nbor .le. 0) then 
            cycle
          end if
          do jcell = 1, twotondim
#ifdef FMM
            multipole_jcell_list(:, jcell, jgrid, ind) = m_source%multipole(jcell, :, igrid_nbor)
            if (multipole_jcell_list(1, jcell, jgrid, ind) /= 0.0d0) then
              source_active_count(jgrid, ind) = source_active_count(jgrid, ind) + 1
              source_active_idx(source_active_count(jgrid, ind), jgrid, ind) = jcell
            end if
#endif
          end do
          if (source_active_count(jgrid, ind) > 0) then
            active_grid_count(ind) = active_grid_count(ind) + 1
            active_grid_idx(active_grid_count(ind), ind) = jgrid
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

    if (sum(active_grid_count) == 0) cycle

    ! Compute interactions for all cells in this AMR grid
    do icell = 1, twotondim
      if (refined_target(icell)) cycle
      phi = 0.0D0
      do ind = 1, threetondim
        ngrid = active_grid_count(ind)
        do iact = 1, ngrid
          jgrid = active_grid_idx(iact, ind)
          nact = source_active_count(jgrid, ind)
          do jcell = 1, nact
            jcell_act = source_active_idx(jcell, jgrid, ind)
            if (nearest_flags(jcell_act, jgrid, ind, icell, igrid)) cycle
            multipole = multipole_jcell_list(:, jcell_act, jgrid, ind)
            diff = diff_list(:, jcell_act, jgrid, ind, icell, igrid)
            D0 = D0_list(jcell_act, jgrid, ind, icell, igrid)
            D1 = D1_list(jcell_act, jgrid, ind, icell, igrid)
            D2 = D2_list(jcell_act, jgrid, ind, icell, igrid)
            call calc_phi_from_multipole(diff, D0, D1, D2, multipole, phi_out)
            phi = phi + phi_out
          end do
        end do 
      end do
      m%phi(icell, ioct) = m%phi(icell, ioct) + phi
    end do
  end do ! end over all amr grids @ given ilev
  deallocate(multipole_jcell_list, D0_list, D1_list, D2_list, nearest_flags, diff_list, source_active_count, &
             source_active_idx, active_grid_count, active_grid_idx)
  call close_cache(mdl)
  end associate
end subroutine fmm_amr_direct_taylor
!################################################################
!################################################################
!################################################################
!################################################################
logical function is_direct_neighbor(cc_icell, cc_jcell)
  use amr_parameters, only: ndim
  implicit none
  integer(kind=8), intent(in) :: cc_icell(ndim), cc_jcell(ndim)
  integer :: d, diff
  
  is_direct_neighbor = .true.
  do d = 1, ndim
     diff = abs(cc_icell(d) - cc_jcell(d))
     !if (ilev > 1) then
      !diff = min(diff, 2**ilev - diff)
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
subroutine init_flush_taylor(mesh,igrid,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  integer::igrid
  type(mesh_t)::mesh
  integer(kind=8),dimension(0:ndim)::hash_key

#ifdef FMM  
  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  mesh%multipole(:,:,igrid)=0.0
  mesh%taylor_coeff(:,:,igrid)=0.0
#endif
end subroutine init_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_taylor(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim,taylor_size
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_large_realdp)::msg
#ifdef FMM
  do ind=1,twotondim
    do ivar=1,taylor_size
      msg%realdp_fmm_taylor(ind, ivar)=mesh%taylor_coeff(ind, ivar, igrid)
    end do
  end do
#endif
  msg_array=transfer(msg,msg_array)
end subroutine pack_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_taylor(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim,taylor_size
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
  do ind=1,twotondim
    do ivar=1,taylor_size
      mesh%taylor_coeff(ind,ivar,igrid)=mesh%taylor_coeff(ind,ivar,igrid)+msg%realdp_fmm_taylor(ind,ivar)
    end do
  end do
#endif
end subroutine unpack_flush_taylor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_taylor(mesh,igrid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim,multipole_size
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_large_realdp)::msg
#ifdef FMM
  msg%realdp_fmm_multipole=mesh%multipole(:,:,igrid)
  msg%realdp_fmm_taylor=mesh%taylor_coeff(:,:,igrid)
#endif
  msg_array=transfer(msg,msg_array)
end subroutine pack_fetch_taylor
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_taylor(mesh,igrid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use amr_commons, only: mesh_t
  use cache_commons, only: msg_large_realdp
  integer::igrid
  type(mesh_t)::mesh
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  type(msg_large_realdp)::msg

  mesh%grid(igrid)%lev=hash_key(0)
  mesh%grid(igrid)%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)
#ifdef FMM
  mesh%multipole(:,:,igrid)=msg%realdp_fmm_multipole
  mesh%taylor_coeff(:,:,igrid)=msg%realdp_fmm_taylor
#endif
end subroutine unpack_fetch_taylor
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
!################################################################
subroutine dump_taylor(r, m, ilev)
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, mesh_t
  implicit none
  type(run_t) :: r
  type(mesh_t) :: m
  integer, intent(in) :: ilev

  integer :: ioct, icell, unit_debug
  character(len=256) :: filename

  ! Construct filename based on level
  write(filename, '(A,I0,A)') "./out_fmm/taylor_level", ilev, ".out"

  unit_debug = 999
  open(unit_debug, file=filename, status="replace")
#ifdef FMM
  do ioct = m%head(ilev), m%tail(ilev)
    write(unit_debug, '(3I6, 20E20.5)') m%grid(ioct)%ckey(1:ndim), m%taylor_coeff(:,:,ioct)
  end do
#endif
  close(unit_debug)
end subroutine dump_taylor
#endif
end module fmm_fine_commons

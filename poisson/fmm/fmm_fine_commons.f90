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
   call r_fmm_amr_direct(pst, pst%s%r%levelmin, 1)
   if(pst%s%r%verbose) print '(A,I2)','Direct Force Calculation done', pst%s%r%levelmin
    
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

  integer :: ioct, idim, pcell, icell, inbor, jcell, nstride
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, xx_pcell    ! box unit real coordinate
  real(kind=8) :: dx_loc
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor
  type(nbor), dimension(1:threetondim) :: grid_nbor
  integer, dimension(1:twotondim) :: ind_nbor_cells

  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp)::dummy_realdp
  real(kind=8), dimension(multipole_size) :: multipole
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor

  associate(r=>s%r, g=>s%g, m=>s%m)

  ! Open cache for multipoles
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,& 
            hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,& 
            pack=pack_fetch_taylor,unpack=unpack_fetch_taylor,& 
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel
  hash_nbor(0) = ilevel
  dx_loc = r%boxlen / 2.0D0**ilevel

  ! Loop over octs at this level
  do ioct = m%head_mg(ilevel), m%tail_mg(ilevel)
     hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim)

     call get_parent_cell(s, hash_key, m%mg_dict, gridp_parent, pcell, flush_cache=.false., fetch_cache=.true.)
     parent_taylor = gridp_parent%taylor_coeff(pcell, :)

     do idim = 1, ndim
        nstride = 2**(idim-1)
        xx_pcell(idim) = (hash_key(idim) + 0.5D0) *  (2 * dx_loc)
     end do

     ! Get neighboring grids. (3^n)
     call get_threetondim_nbor_grid(s, hash_key, m%mg_dict, grid_nbor, flush_cache=.false., fetch_cache=.true.)

     ! Loop over cells in current oct (target)
     do icell = 1, twotondim
        multipole(:) = 0.0D0
        temp_taylor(:) = 0.0D0

        ! get icell position in (0-2^ilevel)
        do idim = 1, ndim
          nstride = 2**(idim-1)
          cc_icell(idim) = 2*hash_key(idim) + MOD((icell-1)/nstride, 2)
          xx_icell(idim) = (cc_icell(idim) + 0.5D0) * dx_loc
        end do

        ! Add shifted taylor from parents
        call shift_taylor(parent_taylor, xx_icell - xx_pcell, temp_taylor)
        m%grid(ioct)%taylor_coeff(icell,:) = m%grid(ioct)%taylor_coeff(icell,:) + temp_taylor

        ! Loop over candidate neighbors
        do inbor = 1, threetondim
           gridp_nbor => grid_nbor(inbor)%p
           if (.not. associated(gridp_nbor)) cycle

           ! Loop over all cells in neighbor oct (sources)
           do jcell = 1, twotondim
              do idim = 1, ndim
                nstride = 2**(idim-1)
                cc_jcell(idim) = 2*hash_key(idim) + MOD((jcell-1)/nstride, 2)
                xx_jcell(idim) = (cc_jcell(idim) + 0.5D0) * dx_loc
              end do

              ! Skip near-field sources
              if (is_direct_neighbor(cc_icell, cc_jcell)) cycle

              ! Shift multipole from origin -> source center
              call shift_multipole(gridp_nbor%multipole(jcell,:), xx_jcell, multipole) !TODO: this is too much shifting

              ! Get taylor coeffs from local
              call calc_taylor(xx_icell, xx_jcell, multipole, temp_taylor)
              m%grid(ioct)%taylor_coeff(icell,:) = m%grid(ioct)%taylor_coeff(icell,:) + temp_taylor
           end do
        end do
     end do

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
  integer :: i, j, k, nfine
  real(kind=8) :: phi, dist2, fourpi
  integer(kind=8), dimension(ndim) :: cc_icell, cc_jcell, cc_fmm_cell ! cartesian coordinate
  real(kind=8), dimension(ndim) :: xx_icell, xx_jcell, xx_pcell       ! box unit real coordinate
  real(kind=8) :: dx_loc
  integer(kind=8), dimension(0:ndim) :: hash_key, hash_nbor
  integer, dimension(1:twotondim) :: ind_nbor_cells

  type(oct), pointer :: gridp_nbor, gridp_parent
  type(msg_large_realdp)::dummy_realdp
  real(kind=8), dimension(multipole_size) :: multipole
  real(kind=8), dimension(taylor_size) :: temp_taylor, parent_taylor

  associate(r=>s%r, g=>s%g, m=>s%m)
  fourpi = 4.D0*ACOS(-1.0D0)
  if(r%cosmo) fourpi = 1.5D0*g%omega_m*g%aexp

  ! Open cache for multipoles
  call open_cache(s,table=m%mg_dict,data_size=storage_size(m%grid(1))/32,& 
            hilbert=m%domain_mg, pack_size=storage_size(dummy_realdp)/32,& 
            pack=pack_fetch_taylor,unpack=unpack_fetch_taylor,& 
            init=init_flush_taylor, flush=pack_flush_taylor, combine=unpack_flush_taylor)

  hash_key(0) = ilevel - g%level_fmm_to_amr
  hash_nbor(0) = ilevel
  dx_loc = r%boxlen / 2.0D0**ilevel
  nfine = 2**g%level_fmm_to_amr

  ! Loop over octs at this level
  do ioct = m%head(ilevel), m%tail(ilevel)
    ! Get taylor from parent fmm grid
    hash_key(1:ndim) = m%grid(ioct)%ckey(1:ndim) / nfine
    pcell = 1

    ! Get fmm cell (not grid) number and its center position
    do idim = 1, ndim
      xx_pcell(idim) = (2 * m%grid(ioct)%ckey(idim) / nfine + 0.5d0) * (nfine / 2) * dx_loc
      if (btest(m%grid(ioct)%ckey(idim) / (2**(g%level_fmm_to_amr-1)), 0)) then
        pcell = pcell + 2**(idim-1)
      end if
    end do

    ! Get Taylor from Parents and calculate phi from above
    call get_grid(s, hash_key, m%mg_dict, gridp_parent, flush_cache=.false., fetch_cache=.true.)
    parent_taylor = gridp_parent%taylor_coeff(pcell, :)

    do icell = 1, twotondim
      phi = 0.0D0
      ! get icell position in (0-2^ilevel)
      do idim = 1, ndim
        nstride = 2**(idim-1)
        cc_icell(idim) = 2*m%grid(ioct)%ckey(idim) + MOD((icell-1)/nstride, 2)
        xx_icell(idim) = (cc_icell(idim) + 0.5D0) * dx_loc
      end do

      ! Now add taylor expansion terms to phi
      call calc_phi(parent_taylor, xx_icell - xx_pcell, phi)
      m%grid(ioct)%phi(icell) = fourpi * phi
    end do

    cc_fmm_cell = cc_icell / nfine
    ! Direct calculation for near neighbors
    ! Now get accumulate phi from direct calculation with neighbors (3^ndim * 2^(level_diff) * 2^ndim cells)  / executed 2^(3*ilevel) times. 
#if NDIM>2
    do k=-nfine/2,nfine
    hash_nbor(3) = nfine / 2 * cc_fmm_cell(3) + k
#endif
#if NDIM>1
    do j=-nfine/2,nfine
    hash_nbor(2) = nfine / 2 * cc_fmm_cell(2) / 2 + j
#endif
#if NDIM>0
    do i=-nfine/2,nfine
    hash_nbor(1) = nfine / 2 * cc_fmm_cell(1) / 2 + i
#endif
    ! Periodic boundary conditions
    do idim=1,ndim
      if(r%periodic(idim))then
        if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel))  hash_nbor(idim)=m%box_ckey_max(idim,ilevel)-1
        if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel)) hash_nbor(idim)=m%box_ckey_min(idim,ilevel)
      endif
    end do
    
    call get_grid(s,hash_nbor,m%grid_dict,gridp_nbor,flush_cache=.false., fetch_cache=.true.)

    ! Target cells
    do icell = 1, twotondim
      phi = 0.0D0
      do idim = 1, ndim
        nstride = 2**(idim-1)
        cc_icell(idim) = 2*m%grid(ioct)%ckey(idim) + MOD((icell-1)/nstride, 2)
        xx_icell(idim) = (cc_icell(idim) + 0.5D0) * dx_loc
      end do

      ! source cells
      do jcell = 1, twotondim
        dist2 = 0.0D0
        do idim = 1, ndim
          nstride = 2**(idim-1)
          cc_jcell(idim) = 2*hash_nbor(idim) + MOD((jcell-1)/nstride, 2)
          xx_jcell(idim) = (cc_jcell(idim) + 0.5D0) * dx_loc
          dist2 = dist2 + min(xx_icell(idim) - xx_jcell(idim), 1-(xx_icell(idim) - xx_jcell(idim))) ** 2
        end do

        if (all(hash_nbor(1:ndim) == m%grid(ioct)%ckey(1:ndim)) .and. icell == jcell) cycle
        print *, 'Dist: ', sqrt(dist2)
        phi = phi - (gridp_nbor%rho(jcell)-g%rho_tot) / sqrt(dist2)
      end do
    m%grid(ioct)%phi(icell) = m%grid(ioct)%phi(icell) + fourpi * phi
    end do
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
  call close_cache(s, m%mg_dict)
  end associate
end subroutine fmm_amr_direct
!################################################################
!################################################################
!################################################################
!################################################################
logical function is_direct_neighbor(cc_icell, cc_jcell)
  use amr_parameters, only: ndim
  implicit none
  integer(kind=8), intent(in) :: cc_icell(ndim), cc_jcell(ndim)
  integer :: d, n
  is_direct_neighbor = .true.
  do d = 1, ndim
     if (abs(cc_icell(d)-cc_jcell(d)) > 1) then
        is_direct_neighbor = .false.
        return
     end if
  end do
end function is_direct_neighbor
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
     do ind=1,twotondim
        msg%realdp_fmm_taylor(ind,ivar)=grid%taylor_coeff(ind,ivar)
     end do
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
     do ind=1,twotondim
        grid%taylor_coeff(ind,ivar)=grid%taylor_coeff(ind,ivar)+msg%realdp_fmm_taylor(ind,ivar)
     end do
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
#endif
end module fmm_fine_commons

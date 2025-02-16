module rt_params_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_read_rt_params(pst)
  use amr_parameters
  use hydro_parameters
  use rt_parameters
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
  implicit none
  type(pst_t)::pst

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  character(LEN=80)::infile
  logical::nml_ok

  !--------------------------------------------------
  ! RT namelist variables
  !--------------------------------------------------

  ! RT_PARAMS namelist
  logical::rt_advect=.false.           ! Advection of photons?                           !
  !logical::rt_smooth=.false.          ! Smooth the discrete RT update of op. splitting  !
  logical::rt_star=.false.             ! Activate radiation from star particles          !
  !logical::rt_sink=.false.            ! Activate radiation from sink particles          !
  real(dp)::rt_esc_frac=1d0            ! Photon escape fraction from stellar particles   !
  !character(LEN=10)::rt_flux_scheme='glf'                                               !
  !logical::rt_use_hll=.false.         ! Use hll flux (or the default glf)               !
  !logical::rt_is_outflow_bound=.false.! Make all boundaries=outflow for RT              !
  real(dp)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  real(dp)::rt_c_fraction=1d0          ! Lightspeed fraction for RT                      !
  !logical::rt_vsla=.false.            ! Are we using level variable light speed?        !
  integer::rt_nsubcycle=1              ! Maximum number of RT subcycles per hydro step   !
  logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  !logical::rt_isDiffuseUVsrc=.false.  ! UV emission from low-density cells              !
  !real(dp)::rt_UVsrc_nHmax=-1d0       ! Density threshold for UV emission               !
  !logical::upload_equilibrium_x=.true.! Enforce equilibrium xion when uploading         !
  !integer::heat_unresolved_HII=0      ! Subgrid model heating unresolved HII regions    !
  !integer::iHIIheat=6                 ! Var index for HII heating                       !
  !logical::cosmic_rays=.false.        ! Include cosmic ray ionisation                   !

  !character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                           !
  character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  !character(LEN=128)::uv_file=''      ! File containing stellar energy distributions    !

  ! RT source regions parameters----------------------------------------------------------
  integer                           ::rt_nsource=0
  character(LEN=10),dimension(1:MAXREGION)::rt_source_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_src_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_src_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_source=2.0
  integer, dimension(1:MAXREGION)   ::rt_src_group=1  
  integer, dimension(1:MAXREGION)   ::rt_src_trace_group=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_source=0.                      ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_source=0.                      !    Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_source=0.                      !    Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_source=0.                      !    Photon flux

  ! RT_GROUPS namelist---------------------------------------------------------------------
  integer::sedprops_update=-1           ! Update sedprops from stellar populations        !
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  ! logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count

  ! Group props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev)
  ! Indices nrtgrp, nion stand for photon group vs species (e.g. 1=H, 2=He).
  real(dp),dimension(nrtgrp,nion)::group_csn=0, group_cse=0      !    Cross sections (cm2)
  real(dp),dimension(nrtgrp)::group_egy=0                        !  Avg photon energy (ev)
  real(dp),dimension(nrtgrp)::group_L0=13.60                     ! Wavelength lower limits
  real(dp),dimension(nrtgrp)::group_L1=0                         ! Wavelength upper limits
  integer,dimension(nion)::spec2group=0                 ! Ion -> group # in recombinations
  real(dp),dimension(nrtgrp)::kappaAbs=0                         ! Dust absorption opacity
  real(dp),dimension(nrtgrp)::kappaSc=0                          ! Dust scattering opacity

  integer:: i

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/rt_params/rt_advect, rt_otsa, rt_c_fraction, rt_nsubcycle
  namelist/rt_sources/rt_star, rt_esc_frac                               &
       & ,rt_nsource, rt_source_type                                     &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center              &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z              &
       & ,rt_exp_source, rt_src_group, rt_src_trace_group                &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source
  namelist/rt_groups/group_csn, group_cse, group_egy, spec2group         &
       & ,group_L0, group_L1, kappaAbs, kappaSc, sed_dir, sedprops_update

  associate(s=>pst%s)

  !-------------------------------------------------
  ! Read the namelist file
  !-------------------------------------------------
  CALL getarg(1,infile)
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'File '//TRIM(infile)//' does not exist'
     call mdl_abort(s%mdl)
  end if
  open(1,file=infile)
  rewind(1)
  read(1,NML=rt_params,END=113)
113 continue
  rewind(1)
  read(1,NML=rt_groups,END=114)
114 continue
  rewind(1)
  read(1,NML=rt_sources,END=115)
115 continue
  close(1)

  ! Fill in all run parameters in corresponding structure
  s%r%rt_otsa=rt_otsa
  s%r%rt_advect=rt_advect
  s%r%rt_c_fraction=rt_c_fraction
  s%r%rt_nsubcycle=rt_nsubcycle
  s%r%rt_courant_factor=rt_courant_factor

  s%r%rt_star=rt_star
  s%r%rt_esc_frac=rt_esc_frac
  s%r%sed_dir=sed_dir
  s%r%rt_nsource=rt_nsource
  s%r%rt_source_type=rt_source_type
  s%r%rt_src_x_center=rt_src_x_center
  s%r%rt_src_y_center=rt_src_y_center
  s%r%rt_src_z_center=rt_src_z_center
  s%r%rt_src_length_x=rt_src_length_x
  s%r%rt_src_length_y=rt_src_length_y
  s%r%rt_src_length_z=rt_src_length_z
  s%r%rt_exp_source=rt_exp_source
  s%r%rt_src_group=rt_src_group
  s%r%rt_n_source=rt_n_source
  s%r%rt_u_source=rt_u_source
  s%r%rt_v_source=rt_v_source
  s%r%rt_w_source=rt_w_source

  s%r%sed_dir=sed_dir
  s%r%sedprops_update=sedprops_update
  s%r%group_csn=group_csn
  s%r%group_cse=group_cse
  s%r%group_egy=group_egy
  s%r%spec2group=spec2group
  s%r%group_L0=group_L0
  s%r%group_L1=group_L1
  s%r%kappaAbs=kappaAbs
  s%r%kappaSc=kappaSc

  if(s%r%isH2) then
    do i=1,nrtgrp
      if((s%r%group_L0(i) .ge. 11.2) .and. (s%r%group_L1(i) .le. 13.6)          &
        .and. (s%r%group_L0(i) .le. 13.6) .and. (s%r%group_L1(i) .ge. 11.2))then
          s%r%ssh2(i) = 4d2 ! H2 self-shielding factor
          s%r%isLW(i) = 1d0 ! Index for LW groups
      endif
    enddo
  endif

  end associate

end subroutine m_read_rt_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module rt_params_module

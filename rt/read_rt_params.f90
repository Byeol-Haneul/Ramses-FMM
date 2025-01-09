module rt_params_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_read_rt_params(pst)
  use amr_parameters
  use hydro_parameters
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
  use rt_parameters, only: nrtgroups, nrtvar
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
  !logical::rt_smooth=.false.           ! Smooth the discrete RT update of op. splitting  !
  real(dp)::units_Np=1.0               ! [#/cm^3]
  real(dp)::smallNp=1d-50              ! Floor value for photon number densities         !
  !real(dp)::rt_Tconst=-1               ! If pos. use this value for all T-depend. rates  !
  !logical::rt_isTconst=.false.         ! Const rates activated?                          !
  !logical::rt_star=.false.             ! Activate radiation from star particles          !
  !logical::rt_sink=.false.             ! Activate radiation from sink particles          !
  !real(dp)::rt_esc_frac=1d0            ! Escape fraction of light from stellar particles !
  !logical::rt_is_init_xion=.false.     ! Initialize ionization from T profile?           !
  !character(LEN=10)::rt_flux_scheme='glf'                                                !
  !logical::rt_use_hll=.false.          ! Use hll flux (or the default glf)               !
  !logical::rt_is_outflow_bound=.false. ! Make all boundaries=outflow for RT              !
  real(dp)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  real(dp)::rt_err_grad_n(nrtgroups)=-1. ! Photon number density gradient for refinement  !
  real(dp)::rt_floor_n(nrtgroups)=1d-10 ! Photon number density floor for refinement      !
  !real(dp)::rt_err_grad_xHI=-1.0       ! Ionization state gradient for refinement        !
  !real(dp)::rt_err_grad_xHII=-1.0      ! Ionization state gradient for refinement        !
  !real(dp)::rt_refine_aexp=-1.0        ! Start a for RT gradient refinement              !
  !real(dp)::rt_floor_xHI=1d-10         ! Ionization state floor for refinement           !
  !real(dp)::rt_floor_xHII=1d-10        ! Ionization state floor for refinement           !
  real(dp)::rt_c_fraction=1d0          ! Lightspeed fraction for RT        !
  !logical::rt_vsla=.false.            ! Are we using level variable light speed?        !
  integer::rt_nsubcycle=1              ! Maximum number of RT-steps during one hydro/    !
                                        ! gravity/etc timestep                            !
  !logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  !logical::rt_isDiffuseUVsrc=.false.  ! UV emission from low-density cells              !
  !real(dp)::rt_UVsrc_nHmax=-1d0       ! Density threshold for UV emission               !
  !logical::upload_equilibrium_x=.true. ! Enforce equilibrium xion when uploading         !
  !integer::heat_unresolved_HII=0      ! Subgrid model heating unresolved HII regions    !
  !integer::iHIIheat=6                 ! Var index for HII heating                       !
  !logical::cosmic_rays=.false.         ! Include cosmic ray ionisation                   !

  !character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                            !
  !character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  !character(LEN=128)::uv_file=''       ! File containing stellar energy distributions    !

  ! Initial condition RT regions parameters-----------------------------------------------
  integer                           ::rt_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::rt_region_type='square'
  real(dp),dimension(1:MAXREGION)   ::rt_reg_x_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_y_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_z_center=0.
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_reg_length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::rt_exp_region=2.0
  integer,dimension(1:MAXREGION)    ::rt_reg_group=1
  real(dp),dimension(1:MAXREGION)   ::rt_n_region=0.                    ! Photon density
  real(dp),dimension(1:MAXREGION)   ::rt_u_region=0.                    !    Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_v_region=0.                    !    Photon flux
  real(dp),dimension(1:MAXREGION)   ::rt_w_region=0.                    !    Photon flux

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
  ! integer::sedprops_update=-1                     ! Update sedprops from star populations
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  ! logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count

  ! Group props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev)
  ! Indexes nrtgroups, nions stand for photon group vs species (e.g. 1=H, 2=He).
  real(dp),dimension(nrtgroups,nions)::group_csn=0, group_cse=0  !    Cross sections (cm2)
  real(dp),dimension(nrtgroups)::group_egy=0                     !  Avg photon energy (ev)
  real(dp),dimension(nrtgroups)::groupL0=13.60                   ! Wavelength lower limits
  real(dp),dimension(nrtgroups)::groupL1=0                       ! Wavelength upper limits
  integer,dimension(nions)::spec2group=0                ! Ion -> group # in recombinations
  real(dp),dimension(nrtgroups)::kappaAbs=0                      ! Dust absorption opacity
  real(dp),dimension(nrtgroups)::kappaSc=0                       ! Dust scattering opacity

  ! NEQ_CHEM namelist---------------------------------------------------------------------
  logical::is_init_xion=.false.                    ! Initialize ionization from T profile?
  logical::isHe=.true.                             !      He ionization fractions tracked?
  logical::isH2=.false.                            !                           H2 tracked?
  real(dp)::X
  real(dp)::Y

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/rt_params/rt_advect, rt_c_fraction, rt_nsubcycle              &
       & ,rt_err_grad_n, rt_floor_n                                      &
       ! RT regions (for initialization)                                 &
       & ,units_np, smallnp, rt_nregion, rt_region_type                  &
       & ,rt_reg_x_center, rt_reg_y_center, rt_reg_z_center              &
       & ,rt_reg_length_x, rt_reg_length_y, rt_reg_length_z              &
       & ,rt_exp_region, rt_reg_group                                    &
       & ,rt_n_region, rt_u_region, rt_v_region, rt_w_region             &
       ! RT source regions (for every timestep)                          &
       & ,rt_nsource, rt_source_type                                     &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center              &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z              &
       & ,rt_exp_source, rt_src_group,   rt_src_trace_group              &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source             
  
  namelist/rt_groups/group_csn, group_cse, group_egy, spec2group         &
       & ,groupL0, groupL1, kappaAbs, kappaSc

  namelist/neq_chem/isHe, isH2, X, Y, is_init_xion

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
  read(1,NML=neq_chem,END=115)
115 continue
  close(1)

  ! Fill in all run parameters in corresponding structure

  ! rt_params
  s%r%rt_advect=rt_advect
  s%r%rt_c_fraction=rt_c_fraction
  s%r%rt_nsubcycle=rt_nsubcycle
  s%r%rt_err_grad_n=rt_err_grad_n
  s%r%rt_floor_n=rt_floor_n
  s%r%rt_courant_factor=rt_courant_factor
  s%r%units_np=units_np
  s%r%smallnp=smallnp
  s%r%rt_nregion=rt_nregion
  s%r%rt_region_type=rt_region_type
  s%r%rt_reg_x_center=rt_reg_x_center
  s%r%rt_reg_y_center=rt_reg_y_center
  s%r%rt_reg_z_center=rt_reg_z_center
  s%r%rt_reg_length_x=rt_reg_length_x
  s%r%rt_reg_length_y=rt_reg_length_y
  s%r%rt_reg_length_z=rt_reg_length_z
  s%r%rt_exp_region=rt_exp_region
  s%r%rt_reg_group=rt_reg_group
  s%r%rt_n_region=rt_n_region
  s%r%rt_u_region=rt_u_region
  s%r%rt_v_region=rt_v_region
  s%r%rt_w_region=rt_w_region
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

  ! rt_groups
  s%r%group_csn=group_csn
  s%r%group_cse=group_cse
  s%r%group_egy=group_egy
  s%r%spec2group=spec2group
  s%r%groupL0=groupL0
  s%r%groupL1=groupL1
  s%r%kappaAbs=kappaAbs
  s%r%kappaSc=kappaSc

  ! neq_chem
  s%r%is_init_xion=is_init_xion
  s%r%isHe=isH2
  s%r%isH2=isH2
  s%r%X=X
  s%r%Y=Y

  end associate

end subroutine m_read_rt_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module rt_params_module

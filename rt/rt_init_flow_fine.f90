!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_rt_init_flow_fine(pst,ilevel)
  use amr_commons, only: ndim
  use amr_parameters, only: nvector, twotondim
  use ramses_commons, only: pst_t
  use rt_input_condinit_module, only: r_rt_input_condinit
  use rt_parameters, only: nrtgrp
  implicit none
  type(pst_t)::pst
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to input a given initial 
  ! condition for RT variables into the exisiting AMR structure from an
  ! analytical model given as an external routine.
  !--------------------------------------------------------------------
  character(len=80)::filename
  logical::ok_file1,ok_file2,ok_file

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m)

  if(s%m%noct_tot(ilevel)==0)return
  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering init_rt_fine for level ',I2)

  ! Use internal-defined or user-defined functions.
  ! We always call the condinit routine, even in cosmological simulations, 
  ! just to initialisethe RT variables (to small values)
  if(s%r%verbose)write(*,*)'Initialising RT variables'
  call r_rt_input_condinit(pst,ilevel,1)

  filename=TRIM(s%r%initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  filename=TRIM(s%r%initfile(ilevel))//'/ic_deltab'
  INQUIRE(file=filename,exist=ok_file2)
  ok_file=ok_file1.or.ok_file2
  if(ok_file)then
     ! No initialization necessary for photons
  else if (s%r%filetype=='gadget')then
     ! No initialization necessary for photons
  endif

  end associate

end subroutine m_rt_init_flow_fine
!###############################################
!###############################################
!###############################################
!###############################################
module update_rt_c_module
  type :: in_broadcast_rt_c
     real(kind=8)::rt_c
  end type in_broadcast_rt_c
contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine update_rt_var(pst, ilevel)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  integer::ilevel
  logical:: cont

  ! Continue if ilevel=levelmin, for regular coarse level updates
  ! Also continue if rt_star=.true. and rt_advect=.false. and nstar_tot>0, 
  ! and if so, activate RT on all cpus, and do stellar p
  cont=.false.
  if(ilevel .eq. pst%s%r%levelmin) cont=.true.

  if(.not. pst%s%r%rt_advect            &
     .and. pst%s%r%rt_star              &
     .and. pst%s%star%npart_tot .gt. 0) then
        ! Need to activate rt_advect on all cpus
        cont=.true.
  endif

  if(cont) call r_update_rt_var(pst, pst%s%g%nstep_coarse, 1)

end subroutine update_rt_var
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_update_rt_var(pst, nstep_coarse, input_size)
  use mdl_module
  use coolrates_module, only: update_rt_c, update_coolrates_tables
  use neq_cooling_module, only: updateRTGroups_CoolConstants, update_metal_cooling
  use ramses_commons, only: pst_t
  use output_rt_module, only: output_photon_emission_stats
  use SED_module, only: update_SED_group_props
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::nstep_coarse, rID

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_UPDATE_RT_VAR,pst%iUpper+1,input_size,0,nstep_coarse)
     call r_update_rt_var(pst%pLower, nstep_coarse,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     ! Update reduced speed of light
     if(r%cosmo)call update_rt_c(r, g, pst%s%tables)

     ! Update Compton heating
     if(r%cosmo)call update_coolrates_tables(r, pst%s%tables, dble(g%aexp))

     ! Check if radiadion advection should be turned on
     if(.not. pst%s%r%rt_advect            &
        .and. pst%s%r%rt_star              &
        .and. pst%s%star%npart_tot .gt. 0) then
        if(pst%s%g%myid==1) then
          write(*,*) '*****************************************'
          write(*,*) 'Stellar RT turned on at a=',pst%s%g%aexp
          write(*,*) '*****************************************'
        endif
        pst%s%r%rt_advect=.true.
        ! Update cross sections based on star properties
        if(pst%s%r%sedprops_update .gt. 0) then
           call update_SED_group_props(r, g, s%SED, s%star)
        endif
     else if(r%star.and.r%rt_advect .and. r%rt_star .and. r%sedprops_update .gt. 0  &
        .and. mod(nstep_coarse,r%sedprops_update)==0)  then
           ! Update cross sections based on evolving star properties
           call update_SED_group_props(r, g, s%SED, s%star)
     endif

     ! Update radiation heating and cooling constants
     if(r%cosmo.or.(r%star.and.r%rt))then
        call updateRTGroups_CoolConstants(r, s%tables)
     endif

     ! Update UV background constants for metal cooling
     if(r%cosmo)call update_metal_cooling(r, s%tables, dble(g%aexp))

     if(r%rt_advect .and. r%rt_emission_stats .and. r%rt_star) &
        call output_photon_emission_stats(r,g)

     if(g%myid==1) write(*,*)'Time dependent RT quantities updated'
  endif

  end associate

end subroutine r_update_rt_var
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module update_rt_c_module

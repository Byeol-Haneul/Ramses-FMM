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
recursive subroutine r_update_rt_var(pst)
  use mdl_module
  use coolrates_module, only: update_rt_c, update_coolrates_tables
  use neq_cooling_module, only: updateRTGroups_CoolConstants, update_metal_cooling
  use ramses_commons, only: pst_t
  use SED_module, only: update_SED_group_props
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::rID

  associate(s=>pst%s,r=>pst%s%r,m=>pst%s%m,g=>pst%s%g)

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_UPDATE_RT_VAR,pst%iUpper+1)
     call r_update_rt_var(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     ! Update reduced speed of light
     if(r%cosmo)call update_rt_c(r, g, pst%s%tables)

     ! Update Compton heating
     if(r%cosmo)call update_coolrates_tables(r, pst%s%tables, dble(g%aexp))

     ! Update cross sections based on evolving star properties
     if(r%star.and.r%rt .and. r%rt_star .and. r%sedprops_update .gt. 0  &
         .and. mod(g%nstep_coarse,r%sedprops_update)==0)  then

                        call update_SED_group_props(r, g, s%SED, s%star)
     endif

     ! Update radiation heating and cooling constants
     if(r%cosmo.or.(r%star.and.r%rt))then
        call updateRTGroups_CoolConstants(r, s%tables)
     endif

     ! Update UV background constants for metal cooling
     if(r%cosmo)call update_metal_cooling(r, s%tables, dble(g%aexp))

     if(g%myid==1) write(*,*)'Time dependent RT quantities updated'
  endif

  end associate

end subroutine r_update_rt_var
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module update_rt_c_module

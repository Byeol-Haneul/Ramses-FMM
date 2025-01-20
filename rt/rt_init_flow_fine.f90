!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_rt_init_flow_fine(pst,ilevel)
  use amr_commons, only: ndim
  use amr_parameters, only: nvector, twotondim
  use ramses_commons, only: pst_t
  use rt_input_condinit_module, only: r_rt_input_condinit
  use rt_parameters, only: nrtgroups
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
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_update_rt_c(pst)
  use constants, only: clight
  use amr_commons, only: run_t, global_t, dp
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !type(run_t)::r
  !type(global_t)::g
  !-------------------------------------------------------------------------
  ! Update the speed of light for radiative transfer, in code units.
  ! This cannot be just a constant, since scale_v changes with time in
  ! cosmological simulations.
  !-------------------------------------------------------------------------
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  type(in_broadcast_rt_c)::in_broadcast

  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)

  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  g%rt_c_cgs = clight * r%rt_c_fraction
  g%rt_c = g%rt_c_cgs/scale_v

  ! Broadcast updated rt_c to all CPUs
  in_broadcast%rt_c=g%rt_c
  call r_broadcast_rt_c(pst,in_broadcast,storage_size(in_broadcast)/32)

  end associate

end subroutine m_update_rt_c
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_broadcast_rt_c(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(in_broadcast_rt_c)::input

  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BROADCAST_RT_C,pst%iUpper+1,input_size,0,input)
     call r_broadcast_rt_c(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%g%rt_c  = input%rt_c
     pst%s%g%rt_c2 = input%rt_c**2
  endif

end subroutine r_broadcast_rt_c
!##############################################################
!##############################################################
!##############################################################
!##############################################################
recursive subroutine r_update_rt_var(pst)
  use mdl_module
  use coolrates_module, only: update_rt_c, update_coolrates_tables
  use rt_cooling_module, only: updateRTGroups_CoolConstants, update_metal_cooling
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_UPDATE_RT_VAR,pst%iUpper+1)
     call r_update_rt_var(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     ! Update reduced speed of light
     call update_rt_c(pst%s%r, pst%s%g, pst%s%tables)
     ! Update Compton heating
     call update_coolrates_tables(pst%s%r, pst%s%tables, dble(pst%s%g%aexp))
     ! Update radiation heating and cooling constants
     call updateRTGroups_CoolConstants(pst%s%r, pst%s%tables)
     ! Update UV background constants for metal cooling
     call update_metal_cooling(pst%s%r, pst%s%tables, dble(pst%s%g%aexp))
  endif

end subroutine r_update_rt_var
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module update_rt_c_module

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_rt_fine(pst,ilevel)
  use ramses_commons, only: pst_t
  use input_rt_condinit_module, only: r_input_rt_condinit
  use rt_parameters, only: nrtgroups
  implicit none
  type(pst_t)::pst
  integer::ilevel,igrp,idim
  !--------------------------------------------------------------------
  ! This routine is the master procedure to input a given initial 
  ! condition for RT variables into the exisiting AMR structure from an
  ! analytical model given as an external routine.
  !--------------------------------------------------------------------
  character(len=80)::filename
  logical::ok_file1,ok_file2,ok_file

  associate(s=>pst%s)
  
  if(s%m%noct_tot(ilevel)==0)return
  if(s%r%verbose)write(*,111)ilevel
111 format(' Entering init_rt_fine for level ',I2)
  ! First initialize everything to zero or small values
  do igrp = 1, nrtgroups
     !rtuold(:,1+(igrp-1)*(ndim+1)) = smallnp
     !do idim = 1, ndim
     !   rtuold(ind,1+idim+(igrp-1)*(ndim+1)) = 0.0
     !end do
  end do
  
  filename=TRIM(s%r%initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  filename=TRIM(s%r%initfile(ilevel))//'/ic_deltab'
  INQUIRE(file=filename,exist=ok_file2)
  ok_file=ok_file1.or.ok_file2
  if(ok_file)then
     ! No initialization necessary for photons
  else if (s%r%filetype=='gadget')then
     ! No initialization necessary for photons
  else
     ! Use internal-defined or user-defined functions
     if(s%r%verbose)write(*,*)'Computing RT initial conditions from analytical model'
     call r_input_rt_condinit(pst,ilevel,1)
  endif

  end associate

end subroutine m_init_rt_fine
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
!*************************************************************************
subroutine m_update_rt_c(pst)

! Update the speed of light for radiative transfer, in code units.
! This cannot be just a constant, since scale_v changes with time in
! cosmological simulations.
!-------------------------------------------------------------------------
  use rt_parameters,only: rt_c, rt_c2, rt_c_cgs
  use amr_parameters, only: clight
  use amr_commons, only: run_t, global_t, dp
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !type(run_t)::r
  !type(global_t)::g
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  type(in_broadcast_rt_c)::in_broadcast
  integer::i
!-------------------------------------------------------------------------
  associate(r=>pst%s%r,g=>pst%s%g,m=>pst%s%m,p=>pst%s%p,mdl=>pst%s%mdl)
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rt_c_cgs = clight * r%rt_c_fraction
  rt_c=rt_c_cgs/scale_v

  ! Broadcast updated rt_c to all CPUs
  in_broadcast%rt_c=rt_c
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
  use rt_parameters,only: rt_c, rt_c2, rt_c_cgs
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
     rt_c   =input%rt_c
     rt_c2  =rt_c**2
  endif

end subroutine r_broadcast_rt_c
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end module update_rt_c_module

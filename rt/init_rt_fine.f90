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


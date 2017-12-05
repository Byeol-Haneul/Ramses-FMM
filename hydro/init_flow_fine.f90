!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_flow_fine(r,g,m,p,mdl,ilevel)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine is the master procedure to input a given initial 
  ! condition into the exisiting AMR structure, either from a file
  ! of from an analytical model given as an external routine.
  !--------------------------------------------------------------------
  character(len=80)::filename
  logical::ok_file1,ok_file2,ok_file

  if(m%noct_tot(ilevel)==0)return
  if(r%verbose)write(*,111)ilevel
111 format(' Entering init_flow_fine for level ',I2)
  
  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  if(r%hydro)then
     filename=TRIM(r%initfile(ilevel))//'/ic_d'
     INQUIRE(file=filename,exist=ok_file1)
     filename=TRIM(r%initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
     ok_file=ok_file1.or.ok_file2
     if(ok_file)then
        ! Read external grafic files 
        if(r%verbose)write(*,*)'Reading initial conditions from grafic file'
        call r_input_hydro_grafic(r,g,m,p,mdl,g%ncpu,1,0,ilevel)
     else
        ! Use internal-defined or user-defined functions
        if(r%verbose)write(*,*)'Computing initial conditions from analytical model'
        call r_input_hydro_condinit(r,g,m,p,mdl,g%ncpu,1,0,ilevel)
     endif
  endif

end subroutine m_init_flow_fine
!###############################################
!###############################################
!###############################################
!###############################################


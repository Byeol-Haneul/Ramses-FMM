!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part(r,g,m,p,mdl)
  use amr_commons, only: run_t,global_t,mesh_t
  use pm_commons, only: part_t
  use mdl_commons, only: mdl_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(mdl_t)::mdl
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from many different initial conditions file formats.
  !--------------------------------------------------------------------
  select case(r%filetype)
  case ('grafic')
    call m_input_part_grafic(r,g,m,p,mdl)
  case('ascii')
    call m_input_part_ascii(r,g,m,p,mdl)
  case('restart')
    call m_input_part_restart(r,g,m,p,mdl)
  case('gadget')
    write(*,*)'Gadget format not supported yet'
    call mdl_abort
  case DEFAULT
    write(*,*) 'Unsupported format file ' // r%filetype
    call mdl_abort
  end select
end subroutine m_input_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

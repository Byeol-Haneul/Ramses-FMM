!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_adaptive(pst)
  use ramses_commons, only: pst_t
  use flag_utils, only: m_flag_fine
  use refine_utils, only: m_refine_fine
  use upload_module, only: m_upload_fine
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  use input_part_zoom_module, only: m_input_part_zoom
  use input_hydro_grafic_module, only: r_input_refmap_grafic
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  integer::istep,ilevel

  write(*,*)'Building initial adaptive grid'

  do istep=pst%s%r%levelmin,pst%s%r%nlevelmax+1

     call m_refine_fine(pst,pst%s%r%levelmin)

     do ilevel=pst%s%r%nlevelmax,pst%s%r%levelmin,-1
        if(pst%s%r%hydro)then
           call m_init_flow_fine(pst,ilevel)
           call m_upload_fine(pst,ilevel)
        endif
        if(pst%s%r%filetype=='grafic_zoom'.and.pst%s%r%ivar_refine==0)then
           call r_input_refmap_grafic(pst,ilevel,1)
        endif
     end do

#ifdef GRAV
     if(pst%s%r%filetype.NE.'grafic_zoom')then
        call m_rho_fine(pst,pst%s%r%levelmin)
     endif
#endif

     do ilevel=pst%s%r%nlevelmax,pst%s%r%levelmin,-1
        call m_flag_fine(pst,ilevel,2)
     end do

  end do

  if(pst%s%r%filetype=='grafic_zoom'.and.pst%s%r%pic)then
     call m_input_part_zoom(pst)
  endif

  call write_screen(pst%s%r,pst%s%m)

end subroutine m_init_refine_adaptive
!###############################################
!###############################################
!###############################################
!###############################################

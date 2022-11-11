!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_adaptive(pst)
  use ramses_commons, only: pst_t
  use flag_utils, only:m_flag_fine
  use refine_utils, only: m_refine_fine
  use upload_module, only: m_upload_fine
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
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
     end do

#ifdef GRAV
     call m_rho_fine(pst,pst%s%r%levelmin)
#endif
     
     do ilevel=pst%s%r%nlevelmax,pst%s%r%levelmin,-1
        call m_flag_fine(pst,ilevel,2)
     end do

  end do

  call write_screen(pst%s%r,pst%s%m)

end subroutine m_init_refine_adaptive
!###############################################
!###############################################
!###############################################
!###############################################

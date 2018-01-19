!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_adaptive(s)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  !--------------------------------------------------------------------
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  integer::istep,ilevel
  
  write(*,*)'Building initial adaptive grid'

  do istep=s%r%levelmin,s%r%nlevelmax+1

     call m_refine_fine(s,s%r%levelmin)

     do ilevel=s%r%nlevelmax,s%r%levelmin,-1
        if(s%r%hydro)then
           call m_init_flow_fine(s,ilevel)
           call m_upload_fine(s,ilevel)
        endif
     end do

     call m_rho_fine(s,s%r%levelmin)

     do ilevel=s%r%nlevelmax,s%r%levelmin,-1
        call m_flag_fine(s,ilevel,2)
     end do

  end do

  do ilevel=s%r%levelmin,s%r%nlevelmax
     call write_screen(s%m,ilevel)
  end do

end subroutine m_init_refine_adaptive
!###############################################
!###############################################
!###############################################
!###############################################

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_refine_adaptive(r,g,m,p,mdl)
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
  ! This routine is the master procedure to set the base grid
  ! and initialize all cell-based variables within it.
  !--------------------------------------------------------------------
  integer::istep,ilevel
  
  write(*,*)'Building initial adaptive grid'

  do istep=r%levelmin,r%nlevelmax+1

     call m_refine_fine(r,g,m,p,mdl,r%levelmin)

     do ilevel=r%nlevelmax,r%levelmin,-1
        if(r%hydro)then
           call m_init_flow_fine(r,g,m,p,mdl,ilevel)
           call m_upload_fine(r,g,m,p,mdl,ilevel)
        endif
     end do

     call m_rho_fine(r,g,m,p,mdl,r%levelmin)

     do ilevel=r%nlevelmax,r%levelmin,-1
        call m_flag_fine(r,g,m,p,mdl,ilevel,2)
     end do

  end do

  do ilevel=r%levelmin,r%nlevelmax
     call write_screen(m,ilevel)
  end do

end subroutine m_init_refine_adaptive
!###############################################
!###############################################
!###############################################
!###############################################

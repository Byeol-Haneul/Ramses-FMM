module cleanup_fmm_module
#ifdef GRAV
contains
recursive subroutine r_cleanup_fmm(pst, ilevel)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::rID, ilevel
    
  if(pst%nLower>0) then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEANUP_FMM,pst%iUpper+1,1,0,ilevel)
     call r_cleanup_fmm(pst%pLower, ilevel)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call m_cleanup_fmm(pst%s%m_fmm_list(ilevel))
  endif
end subroutine r_cleanup_fmm

subroutine m_cleanup_fmm(m)
  use amr_commons, only: mesh_t
  use hash
  implicit none
  type(mesh_t)::m

  ! Reset the MG hash table
  call reset_entire_hash(m%grid_dict,.false.)
  
  ! Restore AMR grid array into its original state
  m%head=1
  m%tail=0
  m%noct=0
  m%ifree=1
  m%noct_tot=0
  m%noct_used=0

end subroutine m_cleanup_fmm
#endif
end module cleanup_fmm_module

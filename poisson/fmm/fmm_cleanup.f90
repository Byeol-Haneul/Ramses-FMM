module cleanup_fmm_module
#ifdef FMM
contains
recursive subroutine r_cleanup_fmm(pst)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer::rID
  
  if(pst%s%r%verbose) write(*,*)'Entering cleanup_fmm'
  
  if(pst%nLower>0) then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEANUP_FMM,pst%iUpper+1)
     call r_cleanup_fmm(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call m_cleanup_fmm(pst%s%m)
  endif
end subroutine r_cleanup_fmm

subroutine m_cleanup_fmm(m)
  use amr_commons, only: mesh_t
  use hash
  implicit none
  type(mesh_t)::m

  integer :: ilev

   ! Deallocate processor boundary array
   deallocate(m%head_mg,m%tail_mg,m%noct_mg)
   do ilev=1,size(m%domain_mg)
      call m%domain_mg(ilev)%destroy
   end do
   deallocate(m%domain_mg)

  ! Reset the MG hash table
  call reset_entire_hash(m%mg_dict,.false.)
  
  ! Restore AMR grid array into its original state
  m%ifree=m%ifree_mg
  m%noct_used=m%ifree-1

end subroutine m_cleanup_fmm
#endif
end module cleanup_fmm_module

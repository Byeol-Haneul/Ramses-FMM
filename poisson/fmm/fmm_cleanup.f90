module cleanup_fmm_module
#ifdef GRAV
contains
recursive subroutine r_cleanup_fmm(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
#if defined(_CUDA) && defined(WITHOUTMPI)
  use gpu_runner, only: gpu_clean_fmm, gpu_clean_fmm_merged
#endif
  implicit none
  type(pst_t)::pst
  integer::ilevel
  integer,VALUE::input_size
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLEANUP_FMM,pst%iUpper+1,input_size,0,ilevel)
     call r_cleanup_fmm(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#if defined(_CUDA) && defined(WITHOUTMPI)
     if (pst%s%m%data_on_device) then
        if (ilevel > 0) then
           call gpu_clean_fmm(pst%s, ilevel)
        else
           call gpu_clean_fmm_merged(pst%s)
        end if
        return
     end if
#endif
     if (ilevel > 0) then
        call m_cleanup_fmm(pst%s%m_fmm_list(ilevel))
     else
        if (associated(pst%s%m_fmm_merged)) call m_cleanup_fmm(pst%s%m_fmm_merged)
     end if
  endif

end subroutine r_cleanup_fmm

recursive subroutine r_cleanup_fmm_merge(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst

  call r_cleanup_fmm(pst, 0, 1)
end subroutine r_cleanup_fmm_merge

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

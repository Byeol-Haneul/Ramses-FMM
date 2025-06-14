module gpu_manager
  use cudafor
  use nvtx
  use gpu_utils
  use gpu_runner
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_grid_device(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_GRID_DEVICE,pst%iUpper+1)
     call r_set_grid_device(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call nvtxStartRange("Copy entire grid from host to device", color=5)!red
     grid_device = pst%s%m%grid
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()
  endif

end subroutine r_set_grid_device
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module gpu_manager

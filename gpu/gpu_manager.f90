module gpu_manager
  use cudafor
  use nvtx
  use gpu_utils
  use gpu_runner
  use amr_parameters, only: threetondim
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_set_grid_device(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use refine_device, only: insert_hash_kernel, update_nbor_array
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID
  integer::head_idx, num_octs, num_threads, num_blocks, ind

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_GRID_DEVICE,pst%iUpper+1)
     call r_set_grid_device(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else

     ! Copy grid from host to device
     call nvtxStartRange("Copy entire mesh from host to device", color=5)!red
     mesh_device%head = pst%s%m%head
     mesh_device%tail = pst%s%m%tail
     mesh_device%noct = pst%s%m%noct
     mesh_device%ifree = pst%s%m%ifree

     mesh_device%head_cache = 1
     mesh_device%tail_cache = 0
     mesh_device%noct_cache = 0
     mesh_device%ifree_cache = 1

     mesh_device%grid = pst%s%m%grid
     mesh_device%flag1 = pst%s%m%flag1
     mesh_device%uold = pst%s%m%uold
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     ! Insert entire grid in the device hash table
     call nvtxStartRange("Insert grid in hash table", color=5)!red
     head_idx = 1
     num_octs = pst%s%m%ifree - 1
     num_threads = 128
     num_blocks = (num_octs + num_threads - 1) / num_threads
     call insert_hash_kernel<<<num_blocks, num_threads>>>(mesh_device, head_idx, num_octs)
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     ! Compute nbor grids array
     call nvtxStartRange("Compute nbor array", color=5)!red
     head_idx = 1
     num_octs = pst%s%m%ifree - 1
     num_threads = 128
     num_blocks = (num_octs + num_threads - 1) / num_threads
     do ind = 1, threetondim
        call update_nbor_array<<<num_blocks, num_threads>>>(mesh_device, head_idx, num_octs, ind)
     end do
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

  endif

end subroutine r_set_grid_device
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module gpu_manager

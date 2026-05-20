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
  use refine_device, only: insert_hash_kernel, update_nbor_array
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID
  integer::head_idx, num_octs, num_subgrids, num_threads, num_blocks, ind
  
  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SET_GRID_DEVICE,pst%iUpper+1)
     call r_set_grid_device(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else

     ! Copy grid from host to device
     call nvtxStartRange("Copy entire mesh from host to device", color=5)!red
     grid = pst%s%m%grid
     flag1 = pst%s%m%flag1
#ifdef HYDRO
     uold = pst%s%m%uold
#endif
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     ! Insert entire grid in the device hash table
     head_idx = 1
     num_octs = pst%s%m%ifree - 1
     num_threads = 128
     num_blocks = (num_octs + num_threads - 1) / num_threads

     call nvtxStartRange("Insert grid in hash table", color=5)!red
     call insert_hash_kernel<<<num_blocks, num_threads>>>(grid, hash_key, hash_val, pst%s%m%hash_size, &
          & ckey_max, key_off, head_idx, num_octs)
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     ! Compute nbor grids array for coarse level only
     head_idx = 1
     num_subgrids = pst%s%m%noct(pst%s%r%levelmin) / nsubgridtondim
     num_threads = 128
     num_blocks = (num_subgrids + num_threads - 1) / num_threads
     call nvtxStartRange("Compute nbor array", color=5)!red
     do ind = 1, subgridsize
        call update_nbor_array<<<num_blocks, num_threads>>>(nbor, grid, hash_key, hash_val, pst%s%m%hash_size, &
             & ckey_max, key_off, box_ckey_min, box_ckey_max, periodic, head_idx, num_subgrids, ind)
     end do
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     pst%s%m%data_on_device=.true.

  endif

end subroutine r_set_grid_device
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_transfer_grid_host(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_TRANSFER_GRID_HOST,pst%iUpper+1)
     call r_transfer_grid_host(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else

     ! Copy grid from device to host
     call nvtxStartRange("Copy entire mesh from device to host", color=5)!red
     pst%s%m%grid = grid
     pst%s%m%flag1 = flag1
#ifdef HYDRO
     pst%s%m%uold = uold
#endif
#ifdef GRAV
     pst%s%m%f = f
     pst%s%m%phi = phi
#endif
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

  endif

end subroutine r_transfer_grid_host
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module gpu_manager

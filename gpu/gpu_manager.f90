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

#ifdef _CUDA
     ! Allocate and copy particle arrays (DM only in phase 1; gpu_part_prompt.md §3, §9).
     ! pst%s%p is value-typed (not a pointer) — gate on npart_max.
     if (pst%s%p%npart_max > 0) then
        call nvtxStartRange("Copy particles from host to device", color=5)!red
        if (allocated(xp))      deallocate(xp)
        if (allocated(vp))      deallocate(vp)
        if (allocated(fp))      deallocate(fp)
        if (allocated(mp))      deallocate(mp)
        if (allocated(levelp))  deallocate(levelp)
        if (allocated(sortp))   deallocate(sortp)
        if (allocated(workp))   deallocate(workp)
        if (allocated(idp))     deallocate(idp)
        if (allocated(hkey_part))   deallocate(hkey_part)
        if (allocated(bucket_part)) deallocate(bucket_part)
        if (allocated(cell_part_count)) deallocate(cell_part_count)
        if (allocated(cell_part_head))  deallocate(cell_part_head)
        if (allocated(cell_part_idx))   deallocate(cell_part_idx)

        allocate(xp(1:pst%s%p%npart_max, 1:ndim))
        allocate(vp(1:pst%s%p%npart_max, 1:ndim))
        allocate(fp(1:pst%s%p%npart_max, 1:ndim))
        allocate(mp(1:pst%s%p%npart_max))
        allocate(levelp(1:pst%s%p%npart_max))
        allocate(sortp(1:pst%s%p%npart_max))
        allocate(workp(1:pst%s%p%npart_max))
        allocate(idp(1:pst%s%p%npart_max))
        allocate(hkey_part(1:pst%s%p%npart_max))
        allocate(bucket_part(1:pst%s%p%npart_max))
        allocate(cell_part_count(1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(cell_part_head (1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(cell_part_idx  (1:pst%s%p%npart_max))

        ! Host -> device (only the fields populated by phase-1 host code)
        if (allocated(pst%s%p%xp))     xp     = pst%s%p%xp
        if (allocated(pst%s%p%vp))     vp     = pst%s%p%vp
        if (allocated(pst%s%p%fp))     fp     = pst%s%p%fp
        if (allocated(pst%s%p%mp))     mp     = pst%s%p%mp
        if (allocated(pst%s%p%levelp)) levelp = pst%s%p%levelp
        if (allocated(pst%s%p%sortp))  sortp  = pst%s%p%sortp
        if (allocated(pst%s%p%workp))  workp  = pst%s%p%workp
        if (allocated(pst%s%p%idp))    idp    = pst%s%p%idp
        call GPU_Error_Check(__FILE__, __LINE__)
        call nvtxEndRange()
     endif
#endif

     ! Insert entire grid in the device hash table
     call nvtxStartRange("Insert grid in hash table", color=5)!red
     head_idx = 1
     num_octs = pst%s%m%ifree - 1
     num_threads = 128
     num_blocks = (num_octs + num_threads - 1) / num_threads
     call insert_hash_kernel<<<num_blocks, num_threads>>>(grid, hash_key, hash_val, pst%s%m%hash_size, &
          & ckey_max, key_off, head_idx, num_octs)
     call GPU_Error_Check(__FILE__, __LINE__)
     call nvtxEndRange()

     ! Compute nbor grids array for coarse level only
     call nvtxStartRange("Compute nbor array", color=5)!red
     head_idx = 1
     num_subgrids = pst%s%m%noct(pst%s%r%levelmin) / nsubgridtondim
     num_threads = 128
     num_blocks = (num_subgrids + num_threads - 1) / num_threads
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

#ifdef _CUDA
     call gpu_to_host_part(pst)
#endif

  endif

end subroutine r_transfer_grid_host
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef _CUDA
!> Copy device particle arrays back to host. Required so PART_DUMP / I/O paths
!> that read pst%s%p%xp etc. see the post-kernel state. See §0.5 / §10.
!> No-op when no particle arrays were allocated on device (npart_max==0).
subroutine gpu_to_host_part(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst

  if (pst%s%p%npart_max <= 0) return
  if (.not. allocated(xp)) return

  call nvtxStartRange("Copy particles from device to host", color=5)!red
  if (allocated(pst%s%p%xp))     pst%s%p%xp     = xp
  if (allocated(pst%s%p%vp))     pst%s%p%vp     = vp
  if (allocated(pst%s%p%fp))     pst%s%p%fp     = fp
  if (allocated(pst%s%p%mp))     pst%s%p%mp     = mp
  if (allocated(pst%s%p%levelp)) pst%s%p%levelp = levelp
  if (allocated(pst%s%p%sortp))  pst%s%p%sortp  = sortp
  if (allocated(pst%s%p%workp))  pst%s%p%workp  = workp
  if (allocated(pst%s%p%idp))    pst%s%p%idp    = idp
  call GPU_Error_Check(__FILE__, __LINE__)
  call nvtxEndRange()

end subroutine gpu_to_host_part
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module gpu_manager

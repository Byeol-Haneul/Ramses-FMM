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
     ! Optional fields (zp/tp/tm/jp/idm/idt/size_p/charge) mirror the host
     ! `allocated(p%X)` gate used by split_part (pm/rho_fine.f90:1413-1494):
     ! allocate device storage iff the host counterpart is allocated.
     if (pst%s%p%npart_max > 0) then
        call nvtxStartRange("Copy particles from host to device", color=5)!red
        ! Mandatory arrays (always allocated when npart_max > 0)
        if (allocated(xp))      deallocate(xp)
        if (allocated(vp))      deallocate(vp)
        if (allocated(mp))      deallocate(mp)
        if (allocated(levelp))  deallocate(levelp)
        if (allocated(sortp))   deallocate(sortp)
        if (allocated(workp))   deallocate(workp)
        if (allocated(idp))     deallocate(idp)
        ! Optional DM fields (allocated only when host has them)
        if (allocated(fp))      deallocate(fp)
        if (allocated(jp))      deallocate(jp)
        if (allocated(zp))      deallocate(zp)
        if (allocated(tp))      deallocate(tp)
        if (allocated(tm))      deallocate(tm)
        if (allocated(size_p))  deallocate(size_p)
        if (allocated(charge))  deallocate(charge)
        if (allocated(idm))     deallocate(idm)
        if (allocated(idt))     deallocate(idt)
        ! Per-particle / per-cell scratch
        if (allocated(hkey_part))   deallocate(hkey_part)
        if (allocated(bucket_part)) deallocate(bucket_part)
        if (allocated(cell_part_count)) deallocate(cell_part_count)
        if (allocated(cell_part_head))  deallocate(cell_part_head)
        if (allocated(cell_part_idx))   deallocate(cell_part_idx)
        ! Out-of-place gather/scatter scratch for split/sort kernels (§4.1, §13).
        ! See the typed-scratch table in gpu/gpu_runner.cuf next to the decl.
        if (allocated(xp_swap))  deallocate(xp_swap)
        if (allocated(mp_swap))  deallocate(mp_swap)
        if (allocated(idp_swap)) deallocate(idp_swap)

        allocate(xp(1:pst%s%p%npart_max, 1:ndim))
        allocate(vp(1:pst%s%p%npart_max, 1:ndim))
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
        ! Gather/scatter scratch — no host->device copy (kernel-overwritten).
        allocate(xp_swap (1:pst%s%p%npart_max, 1:ndim))
        allocate(mp_swap (1:pst%s%p%npart_max))
        allocate(idp_swap(1:pst%s%p%npart_max))

        ! Mandatory host -> device copies
        xp     = pst%s%p%xp
        vp     = pst%s%p%vp
        mp     = pst%s%p%mp
        levelp = pst%s%p%levelp
        sortp  = pst%s%p%sortp
        workp  = pst%s%p%workp
        idp    = pst%s%p%idp

        ! Optional fields: allocate-and-copy only if host has them.
        if (allocated(pst%s%p%fp)) then
           allocate(fp(1:pst%s%p%npart_max, 1:ndim))
           fp = pst%s%p%fp
        endif
        if (allocated(pst%s%p%jp)) then
           allocate(jp(1:pst%s%p%npart_max, 1:ndim))
           jp = pst%s%p%jp
        endif
        if (allocated(pst%s%p%zp)) then
           allocate(zp(1:pst%s%p%npart_max))
           zp = pst%s%p%zp
        endif
        if (allocated(pst%s%p%tp)) then
           allocate(tp(1:pst%s%p%npart_max))
           tp = pst%s%p%tp
        endif
        if (allocated(pst%s%p%tm)) then
           allocate(tm(1:pst%s%p%npart_max))
           tm = pst%s%p%tm
        endif
        if (allocated(pst%s%p%size)) then
           allocate(size_p(1:pst%s%p%npart_max))
           size_p = pst%s%p%size
        endif
        if (allocated(pst%s%p%charge)) then
           allocate(charge(1:pst%s%p%npart_max))
           charge = pst%s%p%charge
        endif
        if (allocated(pst%s%p%idm)) then
           allocate(idm(1:pst%s%p%npart_max))
           idm = pst%s%p%idm
        endif
        if (allocated(pst%s%p%idt)) then
           allocate(idt(1:pst%s%p%npart_max))
           idt = pst%s%p%idt
        endif
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
  ! Mandatory mirrors
  pst%s%p%xp     = xp
  pst%s%p%vp     = vp
  pst%s%p%mp     = mp
  pst%s%p%levelp = levelp
  pst%s%p%sortp  = sortp
  pst%s%p%workp  = workp
  pst%s%p%idp    = idp
  ! Optional mirrors — copy iff both sides are allocated.
  if (allocated(fp)     .and. allocated(pst%s%p%fp))     pst%s%p%fp     = fp
  if (allocated(jp)     .and. allocated(pst%s%p%jp))     pst%s%p%jp     = jp
  if (allocated(zp)     .and. allocated(pst%s%p%zp))     pst%s%p%zp     = zp
  if (allocated(tp)     .and. allocated(pst%s%p%tp))     pst%s%p%tp     = tp
  if (allocated(tm)     .and. allocated(pst%s%p%tm))     pst%s%p%tm     = tm
  if (allocated(size_p) .and. allocated(pst%s%p%size))   pst%s%p%size   = size_p
  if (allocated(charge) .and. allocated(pst%s%p%charge)) pst%s%p%charge = charge
  if (allocated(idm)    .and. allocated(pst%s%p%idm))    pst%s%p%idm    = idm
  if (allocated(idt)    .and. allocated(pst%s%p%idt))    pst%s%p%idt    = idt
  call GPU_Error_Check(__FILE__, __LINE__)
  call nvtxEndRange()

end subroutine gpu_to_host_part
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module gpu_manager

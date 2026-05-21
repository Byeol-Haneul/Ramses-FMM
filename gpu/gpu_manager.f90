module gpu_manager
  use cudafor
  use nvtx
  use gpu_utils
  use gpu_runner
  use gpu_part_state
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
     ! Copy particle arrays host→device (size = host array capacity, not npart_max).
     if (pst%s%r%part .and. allocated(pst%s%p%xp)) then
        call nvtxStartRange("Copy particles from host to device", color=5)!red
        ! Mandatory arrays (always allocated when host xp is allocated)
        if (allocated(xp))      deallocate(xp)
        if (allocated(vp))      deallocate(vp)
        if (allocated(mp))      deallocate(mp)
        if (allocated(levelp))  deallocate(levelp)
        if (allocated(sortp))   deallocate(sortp)
        if (allocated(workp))   deallocate(workp)
        if (allocated(idp))     deallocate(idp)
        ! Optional DM fields (allocated only when host has them)
        if (allocated(fp))      deallocate(fp)
        if (allocated(phip))    deallocate(phip)
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
        if (allocated(src_icell_part))  deallocate(src_icell_part)
        if (allocated(src_igrid_part))  deallocate(src_igrid_part)
        if (allocated(dest_oct_per_cell))   deallocate(dest_oct_per_cell)
        if (allocated(dest_icell_per_cell)) deallocate(dest_icell_per_cell)
        if (allocated(weight_part_rho))  deallocate(weight_part_rho)
        if (allocated(weight_part_nref)) deallocate(weight_part_nref)
        if (allocated(multipole_q_dev))  deallocate(multipole_q_dev)
        ! Gather/scatter scratch for split/sort.
        if (allocated(xp_swap))  deallocate(xp_swap)
        if (allocated(mp_swap))  deallocate(mp_swap)
        if (allocated(idp_swap)) deallocate(idp_swap)

        ! Allocate device mirrors to host capacity. All rank-2 host particle
        ! arrays share the same leading dim (r%npartmax); size(p%xp,1) is the
        ! canonical anchor for the per-particle scratch arrays.
        allocate(xp(1:size(pst%s%p%xp, 1), 1:ndim))
        allocate(vp(1:size(pst%s%p%vp, 1), 1:ndim))
        allocate(mp(1:size(pst%s%p%mp)))
        allocate(levelp(1:size(pst%s%p%levelp)))
        allocate(sortp(1:size(pst%s%p%sortp)))
        allocate(workp(1:size(pst%s%p%workp)))
        allocate(idp(1:size(pst%s%p%idp)))
        allocate(hkey_part  (1:size(pst%s%p%xp, 1)))
        allocate(bucket_part(1:size(pst%s%p%xp, 1)))
        allocate(cell_part_count(1:twotondim, 1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(cell_part_head (1:twotondim, 1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(cell_part_idx  (1:size(pst%s%p%xp, 1)))
        allocate(src_icell_part (1:size(pst%s%p%xp, 1)))
        allocate(src_igrid_part (1:size(pst%s%p%xp, 1)))
        allocate(dest_oct_per_cell   (1:threetondim, 1:twotondim, 1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(dest_icell_per_cell (1:threetondim, 1:twotondim, 1:pst%s%m%ngridmax+pst%s%m%ncachemax))
        allocate(weight_part_rho (1:size(pst%s%p%xp, 1)))
        allocate(weight_part_nref(1:size(pst%s%p%xp, 1)))
        allocate(multipole_q_dev (1:ndim+1))
        ! Gather/scatter scratch — kernel-overwritten, no host->device copy.
        allocate(xp_swap (1:size(pst%s%p%xp, 1), 1:ndim))
        allocate(mp_swap (1:size(pst%s%p%mp)))
        allocate(idp_swap(1:size(pst%s%p%idp)))

        ! Mandatory host -> device copies (shapes match by construction above)
        xp     = pst%s%p%xp
        vp     = pst%s%p%vp
        mp     = pst%s%p%mp
        levelp = pst%s%p%levelp
        sortp  = pst%s%p%sortp
        workp  = pst%s%p%workp
        idp    = pst%s%p%idp


        allocate(fp(1:size(pst%s%p%xp, 1), 1:ndim))
        if (allocated(pst%s%p%fp)) then
           fp(1:size(pst%s%p%fp, 1), 1:ndim) = pst%s%p%fp
        else
           fp = 0d0
        endif
        if (allocated(pst%s%p%phip)) then
           allocate(phip(1:size(pst%s%p%phip)))
           phip = pst%s%p%phip
        endif
        if (allocated(pst%s%p%jp)) then
           allocate(jp(1:size(pst%s%p%jp, 1), 1:ndim))
           jp = pst%s%p%jp
        endif
        if (allocated(pst%s%p%zp)) then
           allocate(zp(1:size(pst%s%p%zp)))
           zp = pst%s%p%zp
        endif
        if (allocated(pst%s%p%tp)) then
           allocate(tp(1:size(pst%s%p%tp)))
           tp = pst%s%p%tp
        endif
        if (allocated(pst%s%p%tm)) then
           allocate(tm(1:size(pst%s%p%tm)))
           tm = pst%s%p%tm
        endif
        if (allocated(pst%s%p%size)) then
           allocate(size_p(1:size(pst%s%p%size)))
           size_p = pst%s%p%size
        endif
        if (allocated(pst%s%p%charge)) then
           allocate(charge(1:size(pst%s%p%charge)))
           charge = pst%s%p%charge
        endif
        if (allocated(pst%s%p%idm)) then
           allocate(idm(1:size(pst%s%p%idm)))
           idm = pst%s%p%idm
        endif
        if (allocated(pst%s%p%idt)) then
           allocate(idt(1:size(pst%s%p%idt)))
           idt = pst%s%p%idt
        endif
        call GPU_Error_Check(__FILE__, __LINE__)
        call nvtxEndRange()
     endif
#endif

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
!> Copy device particle arrays back to host for host-side I/O readers.
!> Call only at explicit host-consumer boundaries.
subroutine gpu_to_host_part(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst

  if (.not. allocated(xp)) return

  call nvtxStartRange("Copy particles from device to host", color=5)!red
  ! Mandatory mirrors (shapes match by construction in r_set_grid_device)
  pst%s%p%xp     = xp
  pst%s%p%vp     = vp
  pst%s%p%mp     = mp
  pst%s%p%levelp = levelp
  pst%s%p%sortp  = sortp
  pst%s%p%workp  = workp
  pst%s%p%idp    = idp
  ! Optional mirrors — copy iff both sides are allocated.
  if (allocated(fp)     .and. allocated(pst%s%p%fp))     pst%s%p%fp     = fp
  if (allocated(phip)   .and. allocated(pst%s%p%phip))   pst%s%p%phip   = phip
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

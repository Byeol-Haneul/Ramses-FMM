# develop sync with rteyssie/mini-ramses (2026-05-26)

Merged `rteyssie/develop` (`724cc5f0`) into `ermoseley/mini-ramses-dev` branch `develop`.

- **Merge commit:** `20a33aaf`
- **Previous develop tip:** `28575111`
- **Upstream range merged:** `ef3a9101` … `724cc5f0` (5 commits)

`develop` is now aligned with upstream `develop` (no remaining commits on `rteyssie/develop` that are not in this branch). The fork still carries local-only history (zero-weight atomicAdd guards, prior merge commits, etc.).

---

## Upstream commits absorbed

| Commit | Summary |
|--------|---------|
| `ef3a9101` | NVTX markers in `adaptive_loop.f90`; NVHPC `-cudalib=nvtx`; Makefile `CUDA_ARCH` / `-Minstrument` hooks |
| `080f20fa` | Merge PR #55 (burlen NVTX instrumentation) |
| `83511388` | Major PIC GPU refactor (`gpu_part.cuf`, `gpu_part_state.cuf`, `gpu_scan.cuf`, `init_part.f90`) |
| `10400bd3` | Compilation fixes in `gpu_manager.f90`, `gpu_part.cuf`, `gpu_part_state.cuf` |
| `724cc5f0` | `CUDA_PROFILE=1` build flags moved under Makefile conditional |

---

## Files changed on develop (net vs `28575111`)

| File | Δ lines | What changed |
|------|---------|--------------|
| `gpu/gpu_part.cuf` | −625 / +565 | Largest change — see below |
| `gpu/gpu_scan.cuf` | +145 | New `real(dp)` warp/block scan kernels (`warp_scan_dp`, `block_scan_dp`, …) |
| `gpu/gpu_runner.cuf` | −44 / +58 | NVTX range cleanup / instrumentation alignment |
| `gpu/gpu_part_state.cuf` | +26 | Device state trimmed; `src_part` layout documented |
| `gpu/gpu_manager.f90` | ±32 | Optional tracer/sink device arrays commented out (`!!$`) |
| `bin/Makefile` | +10 | `CUDA_ARCH`, `CUDA_PROFILE=0`; profile flags only when `CUDA_PROFILE=1` |
| `amr/adaptive_loop.f90` | +15 | Per coarse-step NVTX range (`step_N`) when `_CUDA` |
| `pm/init_part.f90` | +11 | Init aligned with simplified GPU particle state |

---

## GPU PIC / CIC deposit (`gpu/gpu_part.cuf`)

Upstream refactor (Romain / burlen):

- **Removed** `bucket_part_boundary`, old multipole device buffer, and **`reset_cell_part_count_kernel`**.
- **`gpu_cic_part` rewritten:** uses level-local grid slice (`head_grids:tail_grids`), zeros `cell_part_count` via device assignment (`cell_part_count(...) = 0`), then build → flatten → scan → place → `cic_part_warp_kernel`.
- **Multipole:** `multipole_q_kernel` now writes to a local `q_d(1:ndim+1)` device array; host accumulates into `sim%g%multipole%q`.
- **`part_device` imports simplified:** drops unused optional particle fields from module `use`; scan capacity helper `ensure_scan_capacity_part` retained.
- **NVTX:** finer-grained `nvtxStartRange` / `nvtxEndRange` around CIC pipeline stages.

### Kept from mini-ramses-dev (conflict resolution)

1. **Zero-weight atomicAdd guards** in `cic_part_warp_kernel` (from `01e6a34d`):
   - Fast path: skip when `my_rho /= 0d0` / `my_nref /= 0d0`.
   - Slow path: outer `w > 0d0`; inner `my_nref /= 0d0` for `mass_cut_refine`.
2. **`real(..., kind=dp)` casts** on atomic payloads (upstream NPRE=4 fix) — both applied together.
3. **Extended `src_part` bit-layout comments** in `gpu_part.cuf` and `gpu_part_state.cuf` (upstream had a one-line note).
4. **Dropped** local-only `reset_cell_part_count_kernel` — obsolete after upstream switched to direct array zeroing in `gpu_cic_part`.

---

## Build / profiling

- **`bin/Makefile`:** `CUDA_ARCH = sm_80,sm_90` (overrideable); `-mp` removed from default NVHPC flags.
- **`CUDA_PROFILE=1`:** adds `-gpu=lineinfo,ptxinfo -Minstrument -traceback -gopt` for automated NVTX instrumentation.
- **NVTX library:** linked via `-cudalib=nvtx` on NVHPC builds.

---

## Conflict details (manual resolution)

| File | Conflict | Resolution |
|------|----------|------------|
| `gpu/gpu_part.cuf` | Upstream removed `reset_cell_part_count_kernel` + shortened header comment vs local kernel + long comment | Keep long comment; remove dead reset kernel |
| `gpu/gpu_part_state.cuf` | Short vs long `src_part` documentation | Keep long bit-layout comment |

No conflict markers remain; `cic_part_warp_kernel` atomics auto-merged cleanly (guards + `real(...,kind=dp)`).

---

## Suggested follow-up

```bash
git push origin develop
```

Rebuild GPU binary after pull (`make -f Makefile.a100` or your site Makefile). Re-run DMO smoke / bit-exact rho compare if you rely on the CIC deposit path.

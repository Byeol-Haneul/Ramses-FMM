# Halo Gas Analysis

This folder contains a gas-focused, AMR-agnostic isolated-halo analysis pipeline.

## 1) Per-cell IC comparison

`gas_ic_compare.py` compares each output AMR leaf cell against the analytic halo initial-condition equations in `patch/init/halo/condinit.f90`, evaluated at that cell center.

It reports gas drift in:
- density (`rho`)
- thermal pressure (`P`)
- radial velocity (`v_r`, expected near 0)
- azimuthal velocity mismatch (`v_phi - v_phi,IC`)

Outputs:
- `gas_ic_summary.csv`: one row per output snapshot
- `gas_ic_level_stats.csv`: per-AMR-level residual summary
- `gas_ic_radial_profile_noutXXXXX.csv`: radial profile at a chosen output
- `gas_ic_timeseries.png`
- `gas_ic_radial_noutXXXXX.png`

Example:

```bash
python analyze/halo/gas_ic_compare.py \
  --run-dir /path/to/halo_run \
  --namelist namelist/amr_tests/isolated_halo.nml \
  --profile-nout last
```

If your run directory already has `output_XXXXX/namelist.txt`, `--namelist` is optional.

## 2) nstepmax=1 smoke namelist

`prepare_nstep1_namelist.py` creates a cluster smoke-test namelist to validate the pipeline quickly.

Example:

```bash
python analyze/halo/prepare_nstep1_namelist.py \
  --input namelist/amr_tests/isolated_halo.nml \
  --output tmp/isolated_halo_nstep1.nml
```

Then run your cluster executable with the generated namelist and analyze the output:

```bash
python analyze/halo/gas_ic_compare.py --run-dir /path/to/smoke_run
```


# Spheres Convergence (MG vs FMM)

This folder contains scripts to run uniform-spheres convergence studies and compare multigrid (MG) vs FMM.

## Files

- `run_spheres.sh`: Slurm-capable runner for generating namelists, running solver(s), and launching analysis.
- `generate_spheres_namelists.py`: creates level-specific namelists from `namelist/spheres.nml`.
- `spheres_convergence.py`: computes `log10(h)` vs `log10(mean L2 phi error)` and plots convergence curves.

## Typical usage

Run both solvers for levels 6..10 and compare:

```bash
analyze/spheres/run_spheres.sh \
  --solver both \
  --level-start 6 \
  --level-end 10 \
  --nprocs 96 \
  --run-root /path/to/runs/spheres_compare
```

Run only FMM:

```bash
analyze/spheres/run_spheres.sh --solver fmm --level-start 6 --level-end 10
```

Run only MG:

```bash
analyze/spheres/run_spheres.sh --solver mg --level-start 6 --level-end 10
```

Analyze existing runs directly:

```bash
python3 analyze/spheres/spheres_convergence.py \
  --case FMM=/path/to/runs/spheres_compare/fmm \
  --case MG=/path/to/runs/spheres_compare/mg \
  --level-min 6 --level-max 10 --nout 2
```

## Notes

- `nout=2` is the default analysis target (post-solve state).
- Generated per-level runs are stored as:
  - `<run-root>/fmm/level_XX/`
  - `<run-root>/mg/level_XX/`
- Comparison outputs:
  - `<run-root>/spheres_convergence.csv`
  - `<run-root>/spheres_convergence.png`

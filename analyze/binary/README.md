# Binary Analysis Toolkit

This folder contains analysis and movie-making tools for close-binary MG/FMM tests.

## Scripts
- `binary_mg_fmm_compare.py`: MG vs FMM direct comparison for one run root.
- `binary_vs_analytic.py`: compare MG/FMM trajectories against analytic circular orbit.
- `binary_amr_depth_compare.py`: aggregate diagnostics across `levelmax` cases.
- `binary_lmax_compilation.py`: multi-level diagnostic figure.
- `binary_bg6_grid_movie.py`: MG-vs-FMM AMR-grid movie builder for BG6 runs.
  - Supports fixed colorbars with `--global-limits` or manual `--vmin/--vmax`.
  - Overlays analytic particle locations as hollow circles.
- `binary_bg6_six_panel_movie.py`: synchronized 6-panel movie (L6/L7/L8 × MG/FMM) with one fixed colorbar.
- `run_all_bg6.sh`: one-command pipeline for BG6 L6/L7/L8 analyses + movies.

## Notebooks
- `binary_fixed_dt_diagnostics.ipynb`
- `compare_movies.ipynb` is intentionally kept at `analyze/compare_movies.ipynb`.

## Output Layout
By default, generated artifacts go to:
- `analyze/binary/results/bg6_runs/` (per-level symlinked run roots + per-level analyses)
- `analyze/binary/results/analysis_fixed_dt_bg6/` (cross-level summaries/figures)
- `analyze/binary/results/movies_bg6/` (GIF/MP4 + frames)
- `analyze/binary/results/logs/` (pipeline logs)

## Quick Start
From repository root:

```bash
analyze/binary/run_all_bg6.sh
```

Use custom paths:

```bash
analyze/binary/run_all_bg6.sh \
  --binary-root /abs/path/to/binary_test \
  --output-root /abs/path/to/analyze/binary/results
```

Run just movies:

```bash
analyze/binary/run_all_bg6.sh --skip-analyses
```

Skip the combined 6-panel movie:

```bash
analyze/binary/run_all_bg6.sh --skip-six-panel
```

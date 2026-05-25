# Analysis Layout

Analysis tools are organized by test/problem under `analyze/`.

## Problem Folders
- `analyze/binary/`: close-binary MG/FMM analysis and movie tooling.
- `analyze/halo/`: isolated-halo comparison and setup utilities.
- `analyze/spheres/`: spheres test orchestration and convergence analysis.

## Shared/Legacy Files
- `analyze/compare_movies.ipynb`: movie comparison notebook shared across workflows.
- Other top-level scripts/notebooks are legacy/general utilities and can be migrated into per-test folders over time.

The namelist block `&CLUMP_PARAMS` is used to specify parameters controlling the clump and halo finder in the code. Clump finding is performed before each snapshot is outputted to disk. To activate clump finding, set `clump_finder=.true.` in the `&RUN_PARAMS` namelist.

| Variable name   | Fortran type | Default value       | Description                                                                                    |
|:----------------|:-------------|:--------------------|:-----------------------------------------------------------------------------------------------|
| `clump_info`    | `logical`    | `.false.`           | Turn clump finder diagnostics on or off.                                                       |
| `output_clump`  | `logical`    | `.false.`           | Output clump information to disk.                                                              |
| `output_peak`   | `logical`    | `.false.`           | Output grid-based density and peak IDs fields.                                                 |
| `relevance_threshold`  | `real` | 2                  | Set the relevance threshold to remove noisy peaks. |
| `density_threshold`  | `real` | -1                   | Set the density threshold in code units above which peaks are detected. |
| `saddle_threshold`  | `real` | -1                    | Set the saddle point density threshold in code units above which peaks are merged into haloes. |
| `mass_threshold`  | `real` | 0                       | Set the mass threshold in code units abiove which clumps and haloes are discarded. |




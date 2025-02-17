The namelist block `&RT_PARAMS` is used to specify parameters controlling the radiative transfer M1 solver.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `rt_advect`         | `logical`   | `.false.`  | Turn on and off radiation advection. It is off by default and turned on automatically once radiation sources appear for the first time in the simulation. |
| `rt_otsa`           | `logical`   | `.false.`  | Turn on and off the on the spot approximation for ionizing radiation. |
| `rt_smooth`         | `logical`   | `.false.`  | Turn on and off rt smoothing. This means we add the radiation fluxes progressively during the thermo-chemistry step. This is usually faster but not strictly conservative anymore. |
| `rt_c_fraction`     | `real`      | 1.0        | Fraction of the true speed of light used to define the adopted reduced speed of light. |
| `rt_nsubcycle`      | `integer`   | 1          | Maximum number of radiation substeps allowed per hydro steps. This is usually faster but not strictly conservative anymore in case of AMR. |



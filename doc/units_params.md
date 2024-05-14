The namelist block `&UNITS_PARAMS` is used to specify the unit system to convert code units into cgs units and vice-versa. This is the default unit system used in the code. In most cases, units are irrelevant unless you need to use the cooling routines or any of the small-scale physics subgrid models.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `units_density`     | `real`  | 1 | This is the scaling factor used to convert densities from code units to cgs units. Units are therefore [g/cm^3]. |
| `units_time`        | `real`  | 1 | This is the scaling factor used to convert times from code units to cgs units. Units are therefore [s]. When `poisson=.true.`, units of time must be set so that `G=1`. If you are not sure, please use predefined unit systems at compilation time using for example `UNITS=COSMO` or `UNITS=MERGER`. Units are defined in the file `amr/units.f90`. |
| `units_length`      | `real`  | 1 | This is the scaling factor used to convert length scales from code units to cgs units. Units are therefore [cm]. |




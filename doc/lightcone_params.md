The namelist block `&LIGHTCONE_PARAMS` is used to specify parameters controlling the light cone output in cosmological simulations.

| Variable name         | Fortran type | Default value | Description                                                                                    |
|:----------------------|:-------------|:--------------|:-----------------------------------------------------------------------------------------------|
| `lightcone`           | `logical`    | `.false.`     | Activate light cone output for dark matter particles if set to .true.                         |
| `cone_opening_angle_y`| `real`       | 0.0           | Field of view half opening angle in the y-direction (degrees).                                |
| `cone_opening_angle_z`| `real`       | 0.0           | Field of view half opening angle in the z-direction (degrees).                                |
| `cone_theta`          | `real`       | 0.0           | Polar angle (theta) defining the line of sight axis (degrees).                                |
| `cone_phi`            | `real`       | 0.0           | Azimuthal angle (phi) defining the line of sight axis (degrees).                              |
| `cone_z_max`          | `real`       | 0.0           | Maximum redshift to consider for the light cone output.                                       |
| `cone_z_min`          | `real`       | 0.0           | Minimum redshift to consider for the light cone output.                                       |
| `cone_observer`       | `real` array | 0.5,0.5,0.5   | Position of the observer in the simulation box (x,y,z coordinates in box units).              |

**Note:** If both `lightcone=.true.` and `merger_tree=.true.` (from `&CLUMP_PARAMS`) are set, then merger tree particles are also outputted to their own light cone files using the same geometric configuration.

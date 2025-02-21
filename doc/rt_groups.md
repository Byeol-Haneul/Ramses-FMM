The namelist block `&RT_GROUPS` is used to specify parameters controlling the radiative transfer multigroups' properrties.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `group_L0`    | `real array`   | `(\13.60,24.59,54.42\)`  | Lower energy boundaries, in eV, of each photon group. Used for example when calculating SED model emission from stellar particles. |
| `group_L1`    | `real array`   | `(\24.59,54.42,0.0\)`    | Upper energy boundaries, in eV, of each photon group. A value of 0.0 is used to represent infinity. |
| `spec2group`  | `integer array` | `(\1,2,3\)`  | Determines, for each recombining species (HII, HeII, HeIII) which photon group the recombination photons are injected into. Note that recombination emission must be activated with `rt_otsa=.false.` in namelist block `&RT_PARAMS`. |
| `sed_dir`     | `character(len=128)`   | `" "`  | Directory containing spectral energy distribution (SED) model for the stellar emission. This can also be set by the environment variable `RAMSES_SED_DIR`. !
| `sedprops_update`  | `integer` | `-1`  | Frequency (in units of coarse timestep) of photon group updates according to the SED model. The default value of `-1` means that the update is never done. |
| `group_egy` | `real array`   | `(\18.85,35.079,65.666\)` | Average photon energies (eV) for each group. These can either be set manually or left to RAMSES to derive from SED models (with `sedprops_update>0`). |
| `group_csn(1,:)` | `real array`   | `(\3.0d-18,0.0,0.0\)`| 2D array (matrix) representing number-weighted average ionisation cross sections (cm^2) between each group (first index) and species (second index). These can either be set manually or left to RAMSES to derive from SED models (with `sedprops_update>0`). |
| `group_csn(2,:)` |                | `(\5.7d-19,4.5d-18,0.0\),`|  |
| `group_csn(3,:)` |                | `(\7.9d-20,1.2d-18,1.1d-18\),`|  |
| `group_cse(1,;)` | `real array`   | `(\2.8d-18,0.0,0.0\)`| 2D array (matrix) representing energy-weighted average ionisation cross sections (cm^2) between each group (first index) and species (second index). These can either be set manually or left to RAMSES to derive from SED models (with `sedprops_update>0`). |
| `group_cse(2,;)` |                | `(\5.0d-19,4.1d-18,0.0\)`|  |
| `group_cse(3,:)` |                | `(\7.4d-20,1.1d-18,1.0d-18\)`|  |
| `kappaAbs`       | `real array`   | `(\0,0,0\)`  | Dust absorption coefficient (opacity with Planck mean) for each group. The opacity scales with the local metallicity. If `is_kIR_T=.true.`, the IR opacity also scales with the local gas temperature. |
| `kappaSc`        | `real array`   | `(\0,0,0\)`  | Dust scattering coefficient (opacity with Rosseland mean) for each group. In general, this is only for the IR photon group (which is usually the first group). The IT opacity scales with the local metallicity. If `is_kIR_T=.true.`, it also scales with the local gas temperature. |





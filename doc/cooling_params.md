The namelist block `&COOLING_PARAMS` is used to specify parameters controlling gas cooling. 

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `cooling`           | `logical`  | `.false.`  | Turn galaxy formation cooling processes on or off. |
| `cooling_ism`       | `logical`  | `.false.`  | Turn ISM cooling processes on or off. |
| `metal`             | `logical`  | `.false.`  | Turn metal cooling on or off. This also adds metallicity as an additional passive scalar. You must set `NVAR>5` to account for this additional hydro variable. |
| `isothermal`        | `logical`  | `.false.`  | Turn on or off a density-dependent equation of state (EoS). This forces the gas temperature to follow a polytropic EoS whose parameters can be set using the variables defined| `eos_type`          | `int`      | 1  | Type of EoS used in case `isothermal=.true.`. `eos_type=1` correspnds to a strict isothermal EoS. `eos_type=2` corresponds to a polytropic EoS. `eos_type=3` corresponds to the sum of an isothermal and a polytropic EoS, `eos_type=4` corresponds to the max of an isothermal and a polytropic EoS. |
| `eos_nH`            | `real`      | 1d50  | Characteristic density (in units of [H/cc]) used in the polytropic EoS. |
| `eos_T2`            | `real`      | 10    | Temperature (T/mu in units of [K]) of the isothermal EoS or temperature of the polytropic EoS at the characteristic density. |
| `eos_index`         | `real`      | 1     | Power law density dependance of the polytropic EoS. `eos_index=1` is equivalent to the isothermal case `eos_type=1`. |
| `haardt_madau`      | `logical`  | .false. | Turn on or off the default UV background model in the code. |
| `self_shielding`    | `logical`  | .false. | Turn on or off a simple self-shielding model to reduce the UV radiation using the following attenuation factor: `exp(-nH/0.01)`. |
| `J21`               | `real`      | 0      | Nornalization of the UV background Spectral Energy Distribution at 13.6 eV in units of `10^21 erg`. |
| `a_spec`            | `real`      | -1     | Power law of the UV background Spectral Energy Distribution as `J(nu)=J21*(nu/13.6)^(a_spec)`. |
| `z_reion`           | `real`      | 8.5    | Reionization redshift that marks the starting eoph of the adopted UV background. |
| `z_ave`             | `real`      | 0      | Average metallicity (in solar units) used in the cooling function in case `metal=.false.`. If `metal=.true.`, it is used as the default value in the initial conditions for the corresonding passive scalar. Note that 1 solar unit corresponds to the passive scalar `Z=0.02`. |
| `T2max`            | `real`      | 1d50    | Maximum temperature (T/mu in units of [K]) allowed in the cooling routine. If for some reason `T2>T2max`, then the code sets `T2=T2max` in the corresponding cell. |







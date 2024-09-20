This namelist block, called `&OUTPUT_PARAMS`, is used to set up the output strategy of data to disk.
 
| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `nfile=1     ` | `integer` | Number of files per family (amr, hydro, grav, part, star...) per snapshot. These files are stored in directories called `output_00001`, `output_00002`, etc. Default value is `nfile=1`, even if one uses a large number of processors. `nfile=-1` will set `nfile=ncpu`, which is the maximum allowed value. |
| `tend=10     ` | `real`    | End time of the simulation, in time code units. |
| `delta_tout=0` | `real`    | Frequency of outputs in time code units. |
| `aend=0`       | `real`    | End time of the simulation, in expansion factor (for cosmology runs only). |
| `delta_aout=0` | `real`    | Frequency of outputs in expansion factor (for cosmology runs only). |
| `noutput=1`    | `integer` | Number of specified output times.  At least one output time should be given, corresponding to the end of the simulation. Note that this output strategy should no be used in conjuction with `tend` or `aend`.|
| `tout=0.0,0.0,0.0,` | `real array` | Value of specified output times. |
| `aout=1.1,1.1,1.1,` | `real array` | Value of specified output expansion factor (for cosmology runs only). `aout=1.0` means "present epoch" or "zero redshift". |
| `bkp_time_hrs=2` | `real` | Wallclock time in hours between two consecutive backup files. Backup files are used to perform a checkpoint/restart of a simulation that ended prematurely.|
| `bkp_modulo=0`   | `integer` | Backup files are written using a numbering scheme that cycles modulo the prescribed number. Default value of 0 means never. A value of 1 corresponds to only 1 backup file at all time, with the risk of loosing it if something goes wrong during the write. Preferred values are 2 or 3. | 
| `run_time_hrs=0` | `real` | Esimated wallclock time of the current job in hours. This is used to dump a last backup file just before the end of the simulation. Default value of 0 means this is not used.|
| `bkp_last_min=10` | `real` | Dump the last backup file just before the end of the simulation using the prescribed time in minutes. |
| `foutput=N` | `integer` | Frequency of outputs in units of main time steps. If N=0 then no output files are generated. |

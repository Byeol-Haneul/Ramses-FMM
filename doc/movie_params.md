This namelist block, called `&MOVIE_PARAMS`, is used to configure the generation of movie frames during simulation runs. The movie system outputs a series of 2D projection images that can later be combined into movies or animations. Each frame is written to a separate directory with the format `movie1/`, `movie2/`, etc., containing projection maps of various physical quantities.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `movie=.false.` | `logical` | Enable or disable movie frame generation. If `false`, no movie frames are created. |
| `tstartmov=0.0` | `real` | Start time for movie generation in code time units. Frames before this time are skipped. |
| `astartmov=0.0` | `real` | Start time for movie generation in expansion factor (cosmological runs only). Frames before this expansion are skipped. |
| `tendmov=0.0` | `real` | End time for movie generation in code time units. If set to 0, movie generation continues until the end of the simulation. |
| `aendmov=0.0` | `real` | End time for movie generation in expansion factor (for cosmological runs only). If set to 0, movie generation continues until the end of the simulation. |
| `imovout=0` | `integer` | Total number of movie frames to generate. If set to 0, frames are generated based on time intervals. |
| `imov=1` | `integer` | Initial movie frame counter. Usually starts at 1. |
| `nw_frame=512` | `integer` | Width of each movie frame in pixels. |
| `nh_frame=512` | `integer` | Height of each movie frame in pixels. |
| `levelmax_frame=0` | `integer` | Maximum AMR level to use for movie frame generation. If 0, uses the current maximum level. |
| `ivar_frame=1` | `integer` | Hydro variable index used together with `varmin_frame`/`varmax_frame` to filter cells in/out of the projection (e.g. only project cells whose density falls in a band). |
| `xcentre_frame=0.0,0.0,...` | `real array` | X-values for a cubic polynomial in `aexp` determining camera trajectory. 4 values per projection axis. |
| `ycentre_frame=0.0,0.0,...` | `real array` | The same as `xcentre_frame`, but for the Y-direction. |
| `zcentre_frame=0.0,0.0,...` | `real array` | The same as `xcentre_frame`, but for the Z-direction. |
| `deltax_frame=0.0,0.0,...` | `real array` | Extent of the movie frame in the horizontal direction for each projection axis. Two values per projection axis, first specifies a comoving width, second specifies a physical width. `deltax_frame` always corresponds to the first index of the `*.map` file, and so which direction this corresponds to is projection axis dependent. |
| `deltay_frame=0.0,0.0,...` | `real array` | The same as `deltax_frame`, but for the second index of the `*.map' file. |
| `deltaz_frame=0.0,0.0,...` | `real array` | The same as `deltax_frame`, but for the Z-direction. Corresponds to the thickness of the frame.|
| `proj_axis='z'` | `character` | Projection axis for the movie frames. Can be 'x', 'y', 'z', or combinations like 'xy' for multiple projections (which will create a movie directory for each axis). Maximum of 5 projection axes. |
| `zoom_only=.false.` | `logical` | Legacy global flag; superseded by `zoom_only_frame`. Kept for backwards compatibility. |
| `zoom_only_frame=.false.,...` | `logical array` | Per-projection (length 5). When `.true.` and `mass_cut_refine>0` (set in `&REFINE_PARAMS`), the dark matter map only includes particles lighter than `mass_cut_refine` — the standard "show high-resolution particles only" mode for zoom simulations. |
| `movie_vars=0,0,0,...` | `integer array` | Array specifying which variables to include in movie frames. 1 means include, 0 means exclude. Indexed by variable number; usually set indirectly via `movie_vars_txt`. |
| `movie_vars_txt='','','',...` | `character array` | Text labels for variables to include in movie frames. Valid options: `'dens'`, `'vx'`, `'vy'`, `'vz'`, `'temp'`, `'dm'` (dark matter mass), `'stars'` (stellar mass), `'var6'`, `'var7'`, ... (passive scalars `6..NVAR`), `'Fp1'`, `'Fp2'`, ... (radiation photon groups when compiled with RT). |

### Camera parameters

These are length-5 arrays (one entry per projection in `proj_axis`). Default values give no rotation, no perspective, and unit smoothing — set them only if you need a rotating or zooming camera.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `theta_camera=0.0,...` | `real array` | Initial camera angle in the screen-X / screen-Y plane (degrees). |
| `phi_camera=0.0,...` | `real array` | Initial camera angle in the screen-Y / line-of-sight plane (degrees). |
| `dtheta_camera=0.0,...` | `real array` | Total change in `theta_camera` (degrees) over the rotation interval. |
| `dphi_camera=0.0,...` | `real array` | Total change in `phi_camera` (degrees) over the rotation interval. |
| `tstart_theta_camera=0.0,...` | `real array` | Time/expansion at which the θ rotation starts (code units / `aexp`). |
| `tstart_phi_camera=0.0,...` | `real array` | Time/expansion at which the φ rotation starts. |
| `tend_theta_camera=0.0,...` | `real array` | Time/expansion at which the θ rotation ends. If `≤0`, defaults to `tendmov`/`aendmov`. |
| `tend_phi_camera=0.0,...` | `real array` | As `tend_theta_camera`, for φ. |
| `focal_camera=0.0,...` | `real array` | Focal length used for the perspective projection (code units). |
| `dist_camera=0.0,...` | `real array` | Distance from the camera to the frame centre (code units). If `≤0`, defaults to `boxlen`. |
| `ddist_camera=0.0,...` | `real array` | Linear drift in `dist_camera` over the rotation interval (currently unused — reserved). |
| `perspective_camera=.false.,...` | `logical array` | If `.true.` on a projection, apply a perspective transform using `focal_camera`/`dist_camera`. Otherwise the projection is orthographic (default). |

### Rendering options

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `smooth_frame=1.0,...` | `real array` | Per-projection multiplier on the cell half-size used by the shader. Higher values make cells "bleed" further into neighbouring pixels — useful for sparsely-populated frames. |
| `shader_frame='square',...` | `character array` | Cell shape for the hydro projection. `'square'` (default) projects each cell as an axis-aligned rectangle weighted by pixel/cell intersection volume. `'sphere'` projects only the disc inscribed in the cell. `'cube'` falls back to `'square'` (the full rotated-cube shader from RAMSES is not ported). |
| `method_frame='mean_mass',...` | `character array` | Reduction operator for hydro variables. `'mean_mass'`, `'mean_dens'`, `'mean_vol'` produce mass/density/volume-weighted means (the column-weight map is computed once and used as the divisor). `'sum'` accumulates the raw value with no normalisation. `'min'` and `'max'` take the per-pixel extremum across all CPU subdomains. Particle (DM/stars) maps are always summed regardless of `method_frame`. |
| `varmin_frame=-1d60,...` | `real array` | Lower bound on `ivar_frame` for cells to enter the projection. |
| `varmax_frame=1d60,...` | `real array` | Upper bound on `ivar_frame` for cells to enter the projection. |

## Technical Details

### Frame Size Calculation
The physical extent of each movie frame is calculated as:
```
frame_size = comoving_width + physical_width/aexp
```

Where:
- **comoving_width**: Remains constant regardless of expansion
- **physical_width**: Scales with the expansion factor `aexp`
- **aexp**: Current expansion factor (1.0 for non-cosmological runs)

### Parameter Indexing
For multiple projection axes (e.g., `proj_axis='xy'`):
- **First projection axis** (index 1): uses `deltax_frame(1:2)`, `deltay_frame(1:2)`, `deltaz_frame(1:2)`
- **Second projection axis** (index 2): uses `deltax_frame(3:4)`, `deltay_frame(3:4)`, `deltaz_frame(3:4)`
- And so on...

Per-projection arrays (`theta_camera`, `shader_frame`, `method_frame`, `varmin_frame`, …) are indexed directly: `theta_camera(1)` for the first projection, `theta_camera(2)` for the second, etc.

### Dark matter and star projections

Particle maps are written even when `hydro=.false.` (e.g. cosmological DMO runs). They iterate over all locally-resident particles of the relevant type — `pst%s%p` for `'dm'`, `pst%s%star` for `'stars'` — and bin their mass into the per-pixel map. Reductions across MPI ranks are always summed for particle maps. The DM map honours `zoom_only_frame` + `mass_cut_refine`; the stars map has no mass filter.

### Method semantics

- `'mean_mass'`, `'mean_dens'`: numerator is `Σ dvol·ρ·uvar`, denominator is `Σ dvol·ρ` (column mass).
- `'mean_vol'`: numerator is `Σ dvol·uvar`, denominator is `Σ dvol` (intersected volume).
- `'sum'`: map is `Σ uvar` with no normalisation.
- `'min'` / `'max'`: map is the per-pixel min/max of `uvar`. Empty pixels (no contribution) are reset to 0 by the master before writing.

## Output Format

Movie frames are written as binary files in the following format:
- Each frame is stored in a separate directory (e.g., `movie1/`, `movie2/`)
- Files include projection maps of density (`'dens'`), velocity components (`'vx'`, `'vy'`, `'vz'`), temperature (`'temp'`), dark matter mass (`'dm'`), star mass (`'stars'`), and other selected variables
- Each file contains: time/expansion factor, frame dimensions, and the 2D projection data
- Files can be processed with external tools (e.g. `amr2vid.py`) to create animations

## Example Usage (from coeur.nml test problem)

```fortran
&MOVIE_PARAMS
movie=.true.
tendmov=5.62
imovout=1600
nw_frame=512
nh_frame=512
levelmax_frame=16
xcentre_frame=0.5,0.,0.,0.
ycentre_frame=0.5,0.,0.,0.
zcentre_frame=0.5,0.,0.,0.
deltax_frame=0.01,0.
deltay_frame=0.01,0.
deltaz_frame=0.01,0.
proj_axis='z'
movie_vars_txt='dens','vx'
/
```

This configuration generates 1600 movie frames until time 5.62, with 512x512 pixel resolution, centered at (0.5, 0.5, 0.5), showing density and x-velocity projections along the z-axis.

## Example: Cosmological DMO movie

```fortran
&MOVIE_PARAMS
movie=.true.
aendmov=1.0
imovout=50
nw_frame=1024
nh_frame=1024
xcentre_frame=0.5,0.,0.,0.
ycentre_frame=0.5,0.,0.,0.
zcentre_frame=0.5,0.,0.,0.
deltax_frame=1.0,0.
deltay_frame=1.0,0.
deltaz_frame=1.0,0.
proj_axis='z'
movie_vars_txt='dm'
/
```

This writes a 1024×1024 dark-matter column-mass projection along the z-axis at 50 expansion-factor steps between `astartmov=0` and `aendmov=1`. For a zoom-in simulation, add `zoom_only_frame=.true.` and set `mass_cut_refine` in `&REFINE_PARAMS` to exclude the coarse outer particles.

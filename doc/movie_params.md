This namelist block, called `&MOVIE_PARAMS`, is used to configure the generation of movie frames during simulation runs. The movie system outputs a series of 2D projection images that can later be combined into movies or animations. Each frame is written to a separate directory with the format `movie1/`, `movie2/`, etc., containing projection maps of various physical quantities.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `movie=.false.` | `logical` | Enable or disable movie frame generation. If `false`, no movie frames are created. |
| `tendmov=0.0` | `real` | End time for movie generation in code time units. If set to 0, movie generation continues until the end of the simulation. |
| `aendmov=0.0` | `real` | End time for movie generation in expansion factor (for cosmological runs only). If set to 0, movie generation continues until the end of the simulation. |
| `imovout=0` | `integer` | Total number of movie frames to generate. If set to 0, frames are generated based on time intervals. |
| `imov=1` | `integer` | Initial movie frame counter. Usually starts at 1. |
| `nw_frame=512` | `integer` | Width of each movie frame in pixels. |
| `nh_frame=512` | `integer` | Height of each movie frame in pixels. |
| `levelmax_frame=0` | `integer` | Maximum AMR level to use for movie frame generation. If 0, uses the current maximum level. |
| `ivar_frame=1` | `integer` | Variable index for frame generation (legacy parameter). |
| `xcentre_frame=0.0,0.0,...` | `real array` | X-values for a cubic polynomial in `aexp` determining camera trajectory. 4 values per projection axis. |
| `ycentre_frame=0.0,0.0,...` | `real array` | The same as `xcentre_frame`, but for the Y-direction. |
| `zcentre_frame=0.0,0.0,...` | `real array` | The same as `xcentre_frame`, but for the Z-direction. |
| `deltax_frame=0.0,0.0,...` | `real array` | Extent of the movie frame in the horizontal direction for each projection axis. Two values per projection axis, first specifies a comoving width, second specifies a physical width. `deltax_frame` always corresponds to the first index of the `*.map` file, and so which direction this corresponds to is projection axis dependent. |
| `deltay_frame=0.0,0.0,...` | `real array` | The same as `deltax_frame`, but for the second index of the `*.map' file. |
| `deltaz_frame=0.0,0.0,...` | `real array` | The same as `deltax_frame`, but for the Z-direction. Corresponds to the thickness of the frame.|
| `proj_axis='z'` | `character` | Projection axis for the movie frames. Can be 'x', 'y', 'z', or combinations like 'xy' for multiple projections (which will create a movie directory for each axis). Maximum of 5 projection axes. |
| `zoom_only=.false.` | `logical` | If true, only generate zoomed-in frames. Not yet implemented.|
| `movie_vars=0,0,0,...` | `integer array` | Array specifying which variables to include in movie frames. 1 means include, 0 means exclude. |
| `movie_vars_txt='','','',...` | `character array` | Text labels for variables to include in movie frames. Valid options include 'dens', 'vx', 'vy', 'vz', 'temp', 'dm', 'stars', 'var6', 'var7', etc. |

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

## Output Format

Movie frames are written as binary files in the following format:
- Each frame is stored in a separate directory (e.g., `movie1/`, `movie2/`)
- Files include projection maps of density (`'dens'`), velocity components (`'vx'`, `'vy'`, `'vz'`), temperature (`'temp'`), and other selected variables
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
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
| `xcentre_frame=0.0,0.0,...` | `real array` | X-coordinate of the center of the movie frame for each projection axis. Array of up to 20 values. |
| `ycentre_frame=0.0,0.0,...` | `real array` | Y-coordinate of the center of the movie frame for each projection axis. Array of up to 20 values. |
| `zcentre_frame=0.0,0.0,...` | `real array` | Z-coordinate of the center of the movie frame for each projection axis. Array of up to 20 values. |
| `deltax_frame=0.0,0.0,...` | `real array` | Physical size of the movie frame in the X-direction for each projection axis. Array of up to 10 values. |
| `deltay_frame=0.0,0.0,...` | `real array` | Physical size of the movie frame in the Y-direction for each projection axis. Array of up to 10 values. |
| `deltaz_frame=0.0,0.0,...` | `real array` | Physical size of the movie frame in the Z-direction for each projection axis. Array of up to 10 values. |
| `proj_axis='z'` | `character` | Projection axis for the movie frames. Can be 'x', 'y', 'z', or combinations like 'xy' for multiple projections. |
| `zoom_only=.false.` | `logical` | If true, only generate zoomed-in frames around specified centers. |
| `movie_vars=0,0,0,...` | `integer array` | Array specifying which variables to include in movie frames. 1 means include, 0 means exclude. |
| `movie_vars_txt='','','',...` | `character array` | Text labels for variables to include in movie frames. Valid options include 'dens', 'vx', 'vy', 'vz', 'temp', 'dm', 'stars', 'var6', 'var7', etc. |

## Output Format

Movie frames are written as binary files in the following format:
- Each frame is stored in a separate directory (e.g., `movie1/`, `movie2/`)
- Files include projection maps of density, velocity components, temperature, and other selected variables
- Each file contains: time/expansion factor, frame dimensions, and the 2D projection data
- Files can be processed with external tools to create animations

## Example Usage

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

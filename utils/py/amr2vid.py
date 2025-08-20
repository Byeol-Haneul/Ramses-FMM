"""
AMR2VID - Generate movies from RAMSES snapshots using amr2img.py

This script loops through a sequence of RAMSES output snapshots and generates
individual frames using amr2img.py, then combines them into a movie.
Supports parallel processing with MPI for faster generation of large sequences.
"""

import os
import sys
import argparse
import subprocess
import glob
import numpy as np
from pathlib import Path

# Try to import MPI, but don't fail if not available
try:
    from mpi4py import MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False
    print("Warning: mpi4py not available. Running in serial mode.")

def find_output_directories(path, prefix="output_"):
    """Find all output directories matching the pattern."""
    if path is None:
        path = "./"
    
    # Look for directories matching output_XXXXX pattern (supports 5-digit numbers)
    pattern = os.path.join(path, f"{prefix}*")
    dirs = glob.glob(pattern)
    
    # Filter to only include directories
    dirs = [d for d in dirs if os.path.isdir(d)]
    
    # Extract numbers and sort
    output_numbers = []
    for d in dirs:
        try:
            # Extract number from directory name (e.g., "output_00001" -> 1, "output_00301" -> 301)
            # Handle both 4-digit and 5-digit formats
            dirname = os.path.basename(d)
            if dirname.startswith(prefix):
                num_str = dirname[len(prefix):]  # Remove "output_" prefix
                num = int(num_str)
                output_numbers.append((num, d))
        except (ValueError, IndexError):
            continue
    
    # Sort by output number
    output_numbers.sort(key=lambda x: x[0])
    
    return output_numbers

def generate_frame(output_num, output_dir, args, frame_dir):
    """Generate a single frame using amr2img.py."""
    
    # Ensure frame directory exists
    frame_dir = Path(frame_dir)
    frame_dir.mkdir(exist_ok=True)
    
    # Build the amr2img command
    cmd = ["python", "amr2img.py", str(output_num), "--no-display"]
    
    # Add all the arguments from amr2img
    if args.path:
        cmd.extend(["--path", args.path])
    if args.log:
        cmd.extend(["--log"])
    if args.prefix:
        cmd.extend(["--prefix", args.prefix])
    if args.col:
        cmd.extend(["--col", args.col])
    if args.min:
        cmd.extend(["--min", str(args.min)])
    if args.max:
        cmd.extend(["--max", str(args.max)])
    if args.var:
        cmd.extend(["--var", str(args.var)])
    if args.xcen:
        cmd.extend(["--xcen", str(args.xcen)])
    if args.ycen:
        cmd.extend(["--ycen", str(args.ycen)])
    if args.zcen:
        cmd.extend(["--zcen", str(args.zcen)])
    if args.rad:
        cmd.extend(["--rad", str(args.rad)])
    if args.clump:
        cmd.extend(["--clump"])
    if args.sink:
        cmd.extend(["--sink"])
    if args.dir:
        cmd.extend(["--dir", args.dir])
    if args.grid:
        cmd.extend(["--grid"])
    
    # Set output filename for this frame
    frame_filename = f"frame_{output_num:05d}.png"
    frame_path = frame_dir / frame_filename
    cmd.extend(["--out", str(frame_path)])
    
    # Run amr2img
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
        if result.returncode == 0:
            print(f"Generated frame {output_num:05d}")
            return str(frame_path)
        else:
            print(f"Error generating frame {output_num:05d}: {result.stderr}")
            return None
    except Exception as e:
        print(f"Exception generating frame {output_num:05d}: {e}")
        return None

def detect_frame_pattern(frame_dir):
    """Detect the frame numbering pattern from existing frames."""
    frame_pattern = str(Path(frame_dir) / "frame_*.png")
    frames = sorted(glob.glob(frame_pattern))
    
    if not frames:
        return None, None
    
    # Extract the numbering pattern from the first few frames
    sample_frames = frames[:min(5, len(frames))]
    patterns = []
    
    for frame in sample_frames:
        basename = os.path.basename(frame)
        if basename.startswith("frame_"):
            # Extract the number part
            number_part = basename[6:-4]  # Remove "frame_" and ".png"
            patterns.append(len(number_part))
    
    if not patterns:
        return None, None
    
    # Use the most common pattern length
    pattern_length = max(set(patterns), key=patterns.count)
    
    # Determine the format string
    if pattern_length == 1:
        format_str = "frame_%d.png"
    elif pattern_length == 2:
        format_str = "frame_%02d.png"
    elif pattern_length == 3:
        format_str = "frame_%03d.png"
    elif pattern_length == 4:
        format_str = "frame_%04d.png"
    elif pattern_length == 5:
        format_str = "frame_%05d.png"
    else:
        format_str = f"frame_%0{pattern_length}d.png"
    
    return format_str, len(frames)

def create_movie(frame_dir, output_movie, fps=30, quality="high"):
    """Create a movie from the generated frames using ffmpeg."""
    
    # Check if ffmpeg is available
    try:
        subprocess.run(["ffmpeg", "-version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: ffmpeg not found. Please install ffmpeg to create movies.")
        print("Frames have been generated in:", frame_dir)
        return False
    
    # Detect frame pattern automatically
    format_str, num_frames = detect_frame_pattern(frame_dir)
    if format_str is None:
        print("No frames found to create movie.")
        return False
    
    print(f"Detected frame pattern: {format_str}")
    print(f"Number of frames: {num_frames}")
    
    # Try different encoding options based on what's available
    encoding_options = [
        # Option 1: H.264 with libx264
        {
            "name": "H.264 (libx264)",
            "cmd": [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(Path(frame_dir) / format_str),
                "-c:v", "libx264",
                "-crf", "23" if quality == "high" else "28",
                "-pix_fmt", "yuv420p",
                output_movie
            ]
        },
        # Option 2: H.264 with h264_videotoolbox (macOS)
        {
            "name": "H.264 (videotoolbox)",
            "cmd": [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(Path(frame_dir) / format_str),
                "-c:v", "h264_videotoolbox",
                "-b:v", "10M" if quality == "high" else "5M",
                "-pix_fmt", "yuv420p",
                output_movie
            ]
        },
        # Option 3: MPEG-4 (more widely supported)
        {
            "name": "MPEG-4",
            "cmd": [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(Path(frame_dir) / format_str),
                "-c:v", "mpeg4",
                "-q:v", "2" if quality == "high" else "5",
                "-pix_fmt", "yuv420p",
                output_movie
            ]
        }
    ]
    
    print(f"Creating movie: {output_movie}")
    print(f"Using {num_frames} frames at {fps} fps")
    
    # Try each encoding option until one works
    for option in encoding_options:
        print(f"Trying {option['name']}...")
        try:
            result = subprocess.run(option['cmd'], capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                print(f"Movie created successfully using {option['name']}: {output_movie}")
                return True
            else:
                print(f"Failed with {option['name']}: {result.stderr}")
                continue
        except subprocess.TimeoutExpired:
            print(f"Timeout with {option['name']}")
            continue
        except Exception as e:
            print(f"Exception with {option['name']}: {e}")
            continue
    
    print("All encoding options failed. Frames are available for manual processing.")
    print("You can try creating the movie manually with:")
    print(f"ffmpeg -framerate {fps} -i {frame_dir}/frame_%05d.png -c:v mpeg4 -q:v 2 {output_movie}")
    return False

def main():
    parser = argparse.ArgumentParser(
        description="Generate movies from RAMSES snapshots using amr2img.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate movie from outputs 1-100
  python amr2vid.py 1 100
  
  # Generate movie with custom settings
  python amr2vid.py 1 100 --fps 60 --quality high --log --var 1
  
  # Run in parallel with MPI (4 processes)
  mpirun -np 4 python amr2vid.py 1 100 --parallel
  
  # Create movie from existing frames (skip generation)
  python amr2vid.py 1 100 --ffmpeg-only --fps 30 --quality high
  
  # Create movie from frames in custom directory
  python amr2vid.py 1 100 --ffmpeg-only --frame-dir my_frames --output my_movie.mp4
        """
    )
    
    # Range arguments
    parser.add_argument("start", type=int, help="starting output number")
    parser.add_argument("end", type=int, help="ending output number")
    
    # Movie-specific arguments
    parser.add_argument("--fps", type=int, default=30, help="frames per second (default: 30)")
    parser.add_argument("--quality", choices=["low", "medium", "high"], default="medium", 
                       help="movie quality (default: medium)")
    parser.add_argument("--output", help="output movie filename (default: movie.mp4)")
    parser.add_argument("--frame-dir", help="directory to store frames (default: frames/)")
    parser.add_argument("--parallel", action="store_true", help="use MPI parallel processing")
    parser.add_argument("--keep-frames", action="store_true", help="keep individual frames after movie creation")
    parser.add_argument("--ffmpeg-only", action="store_true", help="skip frame generation, create movie from existing frames")
    
    # All amr2img arguments
    parser.add_argument("--path", help="specify a path")
    parser.add_argument("--log", help="plot log variable", action="store_true")
    parser.add_argument("--prefix", help="specify a file prefix")
    parser.add_argument("--col", help="choose the color map")
    parser.add_argument("--min", help="specify a minimum variable value for colorbar")
    parser.add_argument("--max", help="specify a maximum variable value for colorbar")
    parser.add_argument("--var", help="specify a variable number")
    parser.add_argument("--xcen", help="specify the image center x-coordinate")
    parser.add_argument("--ycen", help="specify the image center y-coordinate")
    parser.add_argument("--zcen", help="specify the image center z-coordinate")
    parser.add_argument("--rad", help="specify the image radius")
    parser.add_argument("--clump", help="specify if clumps are overplotted")
    parser.add_argument("--sink", help="specify if sinks are overplotted")
    parser.add_argument("--dir", help="specify the projection axis")
    parser.add_argument("--grid", help="overlay the AMR grid", action="store_true")
    
    args = parser.parse_args()
    
    # Check MPI availability
    if args.parallel and not MPI_AVAILABLE:
        print("Error: --parallel requested but mpi4py not available.")
        print("Install mpi4py with: pip install mpi4py")
        sys.exit(1)
    
    # Initialize MPI if using parallel mode
    if args.parallel:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        is_root = (rank == 0)
    else:
        comm = None
        rank = 0
        size = 1
        is_root = True
    
    # Set default values
    if args.output is None:
        args.output = "movie.mp4"
    if args.frame_dir is None:
        args.frame_dir = "frames"
    
    # Create frame directory (all processes need to do this)
    frame_dir = Path(args.frame_dir).resolve()  # Use absolute path
    frame_dir.mkdir(exist_ok=True)
    if is_root:
        print(f"Frame directory: {frame_dir}")
    
    # If ffmpeg-only mode, skip frame generation and just create movie
    if args.ffmpeg_only:
        if is_root:
            print("ffmpeg-only mode: skipping frame generation")
            # Check if frames exist
            frame_pattern = str(frame_dir / "frame_*.png")
            existing_frames = sorted(glob.glob(frame_pattern))
            if not existing_frames:
                print(f"No frames found in {frame_dir}")
                print("Use pattern: frame_00001.png, frame_00002.png, etc.")
                sys.exit(1)
            
            print(f"Found {len(existing_frames)} existing frames")
            # Create movie from existing frames
            success = create_movie(frame_dir, args.output, args.fps, args.quality)
            if success:
                print(f"Movie created successfully: {args.output}")
            else:
                print("Failed to create movie")
            sys.exit(0)
    
    # Find output directories
    if is_root:
        output_dirs = find_output_directories(args.path)
        # Filter to requested range
        output_dirs = [(num, dir_path) for num, dir_path in output_dirs 
                      if args.start <= num <= args.end]
        
        if not output_dirs:
            print(f"No output directories found in range {args.start}-{args.end}")
            sys.exit(1)
        
        print(f"Found {len(output_dirs)} output directories")
        for num, dir_path in output_dirs:
            print(f"  {num:5d}: {dir_path}")
    
    # Broadcast output directories to all processes
    if args.parallel:
        if is_root:
            output_nums = [num for num, _ in output_dirs]
        else:
            output_nums = None
        output_nums = comm.bcast(output_nums, root=0)
        output_dirs = [(num, None) for num in output_nums]  # Only need numbers for parallel processing
    else:
        output_nums = [num for num, _ in output_dirs]
    
    # Distribute work among processes
    if args.parallel:
        # Simple round-robin distribution
        my_outputs = [num for i, num in enumerate(output_nums) if i % size == rank]
        print(f"Process {rank}: processing {len(my_outputs)} outputs")
    else:
        my_outputs = output_nums
    
    # Generate frames
    generated_frames = []
    for output_num in my_outputs:
        frame_path = generate_frame(output_num, None, args, frame_dir)
        if frame_path:
            generated_frames.append(frame_path)
    
    # Gather all generated frames
    if args.parallel:
        all_frames = comm.gather(generated_frames, root=0)
        if is_root:
            generated_frames = [frame for sublist in all_frames for frame in sublist]
            generated_frames.sort()  # Sort by frame number
    else:
        generated_frames.sort()
    
    # Create movie (only on root process)
    if is_root and generated_frames:
        print(f"Generated {len(generated_frames)} frames")
        
        # Create movie
        success = create_movie(frame_dir, args.output, args.fps, args.quality)
        
        # Clean up frames if requested
        if success and not args.keep_frames:
            print("Cleaning up frame files...")
            for frame_file in generated_frames:
                try:
                    os.remove(frame_file)
                except OSError:
                    pass
            try:
                os.rmdir(frame_dir)
            except OSError:
                pass
            print("Frame cleanup complete")
        elif args.keep_frames:
            print(f"Frames kept in: {frame_dir}")
    
    if args.parallel:
        MPI.Finalize()

if __name__ == "__main__":
    main()

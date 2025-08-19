#!/usr/bin/env python3
"""
Uniform initial condition generator for mini-ramses in GRAFIC format.

This script writes unformatted Fortran records compatible with the mini-ramses
"grafic" file reader. It supports both 3D and 2D outputs.

Usage examples:

  - 3D (n = 2^6 = 64):
      python3 uniform.py 6 --size 1.0 --ndim 3 \
        --outdir /Users/moseley/ramses-development/mini-ramses-main/ics/ic_uniform/ic_uniform_6

  - 2D (writes (n, n, 1) arrays):
      python3 uniform.py 6 --size 1.0 --ndim 2 \
        --outdir /Users/moseley/ramses-development/mini-ramses-main/ics/ic_uniform/ic_uniform_6_d2

Options:
  --rho, --p0, --u0, --v0, --w0 control uniform hydro fields.
  --bz writes uniform magnetic Bz boundary files (left/right) when nonzero.
  --seed controls random particle ID shuffle reproducibly.

Notes:
  - GRAFIC headers encode [n1, n2, n3, dx, ...], with dx = size / n1.
  - For 2D, n3 = 1 and slices are still written as z-slice records.
  - Particle IDs are a random permutation 1..(n1*n2*n3), reshaped to (n1,n2,n3)
    and written as int64.
"""
import argparse
import os
from pathlib import Path

import numpy as np

import grafic


def write_array(filename: str, array: np.ndarray, box_size_cu: float, *, as_int64: bool = False) -> None:
    """Write a 3D numpy array to a GRAFIC file.

    Args:
        filename: Destination filename (created/overwritten).
        array: 3D array shaped (n1, n2, n3). For 2D output, use n3=1.
        box_size_cu: Physical box size in code units used to set dx in header.
        as_int64: If True, write data as int64; otherwise float32.
    """
    g = grafic.Grafic()
    g.set_data(np.asarray(array))
    g.make_header(box_size_cu)
    if as_int64:
        g.write_int64(filename)
    else:
        g.write_float(filename)


def main():
    parser = argparse.ArgumentParser(description="Generate uniform GRAFIC ICs (hydro + optional Bz + particle IDs)")
    parser.add_argument("lvl", type=int, help="Refinement level (grid size is 2^lvl)")
    parser.add_argument("--size", type=float, default=1.0, help="Box size in code units (default: 1.0)")
    parser.add_argument("--ndim", type=int, default=3, help="Output dimensionality: 2 or 3 (default: 3)")
    parser.add_argument("--rho", type=float, default=1.0, help="Uniform density")
    parser.add_argument("--p0", type=float, default=1.0, help="Uniform pressure")
    parser.add_argument("--u0", type=float, default=0.0, help="Uniform vx")
    parser.add_argument("--v0", type=float, default=0.0, help="Uniform vy")
    parser.add_argument("--w0", type=float, default=0.0, help="Uniform vz")
    parser.add_argument("--bz", type=float, default=0.0, help="Uniform Bz (writes left/right boundary fields)")
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Output directory for ICs. Default: ./ic_uniform/ic_uniform_<lvl>",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for particle IDs")
    parser.add_argument("--particles", action="store_true", help="Also write particle ICs (IDs and velocities)")

    args = parser.parse_args()

    lvl = int(args.lvl)
    L = float(args.size)
    n = 2 ** lvl
    if args.ndim == 3:
        n1, n2, n3 = n, n, n
    elif args.ndim == 2:
        n1, n2, n3 = n, n, 1
    else:
        raise SystemExit("--ndim must be 2 or 3")

    # Uniform hydro fields
    d = np.full((n1, n2, n3), float(args.rho), dtype=np.float32)
    p = np.full((n1, n2, n3), float(args.p0), dtype=np.float32)
    u = np.full((n1, n2, n3), float(args.u0), dtype=np.float32)
    v = np.full((n1, n2, n3), float(args.v0), dtype=np.float32)
    w = np.full((n1, n2, n3), float(args.w0), dtype=np.float32)

    # Optional particles
    if args.particles:
        rng = np.random.default_rng(args.seed)
        ids = np.arange(1, (n1*n2*n3) + 1, dtype=np.int64)
        rng.shuffle(ids)
        ids = ids.reshape((n1, n2, n3))
        velcx = np.zeros((n1, n2, n3), dtype=np.float32)
        velcy = np.zeros((n1, n2, n3), dtype=np.float32)
        velcz = np.zeros((n1, n2, n3), dtype=np.float32)

    # Determine output dir
    tag = f"ic_uniform_{lvl}_{args.ndim}d"
    if args.outdir is None:
        outdir = Path("ic_uniform") / tag
    else:
        outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    # Write hydro primitives
    write_array("ic_d", d, L)
    write_array("ic_u", u, L)
    write_array("ic_v", v, L)
    write_array("ic_w", w, L)
    write_array("ic_p", p, L)

    # Optional magnetic fields (left/right boundaries)
    Bx = 0.0
    By = 0.0
    Bz = float(args.bz)
    if Bx != 0.0 or By != 0.0 or Bz != 0.0:
        write_array("ic_bxleft",  np.full((n1, n2, n3), Bx, dtype=np.float32), L)
        write_array("ic_bxright", np.full((n1, n2, n3), Bx, dtype=np.float32), L)
        write_array("ic_byleft",  np.full((n1, n2, n3), By, dtype=np.float32), L)
        write_array("ic_byright", np.full((n1, n2, n3), By, dtype=np.float32), L)
        write_array("ic_bzleft",  np.full((n1, n2, n3), Bz, dtype=np.float32), L)
        write_array("ic_bzright", np.full((n1, n2, n3), Bz, dtype=np.float32), L)

    # Particles
    if args.particles:
        write_array("ic_particle_ids", ids, L, as_int64=True)
        write_array("ic_velcx", velcx, L)
        write_array("ic_velcy", velcy, L)
        write_array("ic_velcz", velcz, L)

    print(f"Uniform ICs written to: {outdir}")


if __name__ == "__main__":
    main()



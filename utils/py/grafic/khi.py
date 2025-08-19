#!/usr/bin/env python3
"""
Kelvin–Helmholtz (KHI) initial condition generator for mini-ramses in GRAFIC format.

This script writes unformatted Fortran records compatible with the mini-ramses
"grafic" file reader. It supports both 3D and 2D outputs.

Usage examples:

  - 3D (n = 2^6 = 64):
      python3 khi.py 6 --size 1.0 --ndim 3 \
        --outdir /Users/moseley/ramses-development/mini-ramses-main/ics/ic_khi/ic_khi_6

  - 2D (writes (n, n, 1) arrays):
      python3 khi.py 6 --size 1.0 --ndim 2 \
        --outdir /Users/moseley/ramses-development/mini-ramses-main/ics/ic_khi/ic_khi_6_d2

Important details:
  - Output files: ic_d, ic_u, ic_v, ic_w, ic_p (and optional magnetic boundaries
    if you extend this script).
  - Header encodes [n1, n2, n3, dx, ...] with dx = size / n1. For 2D, n3 = 1
    and data are still written as a sequence of z-slices.
  - Pressure and other fields must be written as full 3D arrays to ensure each
    slice record has the correct length; do not write (1, n, 1) views.
"""
import argparse
import os
from pathlib import Path

import numpy as np

import grafic


def write_array(filename: str, array: np.ndarray, box_size_cu: float) -> None:
    """Write a 3D numpy array to a GRAFIC file.

    Args:
        filename: Destination filename (created/overwritten).
        array: 3D array shaped (n1, n2, n3). For 2D output, use n3=1.
        box_size_cu: Physical box size in code units used to set dx in header.
    """
    g = grafic.Grafic()
    g.set_data(np.asarray(array))
    g.make_header(box_size_cu)
    g.write_float(filename)


def generate_kh_fields(
    lvl_max: int,
    box_size_cu: float,
    u0: float = 1.0,
    shear_thickness: float = 0.02,
    density_outer: float = 1.0,
    density_inner: float = 1.0,
    pressure0: float = 1.0,
    pressure_inner: float = 1.0,
    perturb_eps: float = 0.01,
    perturb_sigma: float = 0.02,
    random_seed: int | None = 42,
    ndim: int = 3,
):
    """Generate smoothed KHI hydro primitives.

    Produces density, velocity components (u, v, w) and pressure arrays with a
    tanh shear in y and a small vy sinusoidal perturbation localized near two
    interfaces.

    Args:
        lvl_max: Refinement level; linear resolution n = 2**lvl_max.
        box_size_cu: Domain size in code units.
        u0: Velocity jump magnitude between layers.
        shear_thickness: Half-thickness of the shear transition (fraction of box).
        density_outer: Density in outer slabs (top/bottom).
        density_inner: Density in inner slab (middle region).
        pressure0: Pressure in outer slabs.
        pressure_inner: Pressure in inner slab.
        perturb_eps: Amplitude of vy perturbation.
        perturb_sigma: Gaussian width of perturbation (fraction of box).
        random_seed: RNG seed for any randomized content (not currently used).
        ndim: 2 or 3; 2D returns arrays with shape (n, n, 1).

    Returns:
        Tuple (d, u, v, w, p), each shaped (n, n, n3) with n3 = 1 if ndim=2.
    """
    """
    Generate 3D Kelvin–Helmholtz initial conditions with smoothed shear layers.

    Velocity is along x, varying in y with two tanh interfaces at y=0.25 and y=0.75.
    Adds a small vy sinusoidal perturbation localized near the shear layers.qq
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    n = 2 ** int(lvl_max)
    L = float(box_size_cu)
    if ndim == 3:
        n1, n2, n3 = n, n, n
    elif ndim == 2:
        n1, n2, n3 = n, n, 1
    else:
        raise ValueError("ndim must be 2 or 3")

    # Coordinates in code units [0, L)
    x = (np.arange(n1, dtype=np.float32) + 0.5) * (L / n1)
    y = (np.arange(n2, dtype=np.float32) + 0.5) * (L / n2)
    # z not used explicitly; data arrays carry n3

    # Normalize to [0,1) for profile definitions
    yn = y / L
    X = x / L

    # Smoothed shear profile: outer (top+bottom) = U1, inner (middle slab) = U2
    U1 = -0.5 * u0
    U2 = +0.5 * u0
    a = float(shear_thickness)
    tanh_low = np.tanh((yn - 0.25) / a)
    tanh_high = np.tanh((yn - 0.75) / a)
    vx_1d = U1 + 0.5 * (U2 - U1) * (tanh_low - tanh_high)  # shape (n2,)

    # Densities: choose uniform by default; allow contrast if desired
    rho_outer = float(density_outer)
    rho_inner = float(density_inner)
    rho_1d = rho_outer + 0.5 * (rho_inner - rho_outer) * (tanh_low - tanh_high)

    # Pressures: choose uniform by default; allow contrast if desired
    p0 = float(pressure0)
    p_inner = float(pressure_inner)
    p_1d = p0 + 0.5 * (p_inner - p0) * (tanh_low - tanh_high)

    # Small vy perturbation to seed KH, localized near interfaces
    sig = float(perturb_sigma)
    gauss = np.exp(-((yn - 0.25) ** 2) / (2.0 * sig * sig)) + np.exp(
        -((yn - 0.75) ** 2) / (2.0 * sig * sig)
    )  # shape (n,)
    sinus = np.sin(2.0 * np.pi * X, dtype=np.float32)  # shape (n,)

    # Allocate 3D fields
    d = np.empty((n1, n2, n3), dtype=np.float32)
    p = np.empty((n1, n2, n3), dtype=np.float32)
    u = np.empty((n1, n2, n3), dtype=np.float32)
    v = np.empty((n1, n2, n3), dtype=np.float32)
    w = np.zeros((n1, n2, n3), dtype=np.float32)

    # Broadcast 1D profiles
    # Axis order is (x, y, z) to match grafic writer expectations (Fortran order handled inside)
    vx = vx_1d[None, :, None]
    rho = rho_1d[None, :, None]
    d[:, :, :] = rho
    u[:, :, :] = vx
    p[:, :, :] = p_1d[None, :, None]

    # vy perturbation: eps * sin(2π x/L) * gaussian(y)
    vy = perturb_eps * sinus[:, None, None] * gauss[None, :, None]
    v[:, :, :] = vy.astype(np.float32)

    return d, u, v, w, p


def main():
    parser = argparse.ArgumentParser(description="Generate Kelvin–Helmholtz ICs (GRAFIC format)")
    parser.add_argument("lvl", type=int, help="Refinement level (grid size is 2^lvl)")
    parser.add_argument("--size", type=float, default=1.0, help="Box size in code units (default: 1.0)")
    parser.add_argument("--ndim", type=int, default=3, help="Output dimensionality: 2 or 3 (default: 3)")
    parser.add_argument("--u0", type=float, default=1.0, help="Velocity jump magnitude (center vs outer)")
    parser.add_argument("--thickness", type=float, default=0.02, help="Shear layer half-thickness (fraction of box)")
    parser.add_argument("--rho_outer", type=float, default=1.0, help="Outer slab density")
    parser.add_argument("--rho_inner", type=float, default=1.0, help="Inner slab density")
    parser.add_argument("--p0", type=float, default=2.5, help="Uniform pressure")
    parser.add_argument("--eps", type=float, default=0.01, help="Perturbation amplitude for vy")
    parser.add_argument("--sigma", type=float, default=0.02, help="Perturbation Gaussian width (fraction of box)")
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Output directory for ICs. Default: ./ic_khi/ic_khi_<lvl>",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed (not critical; for reproducibility)")
    parser.add_argument("--pressure_inner", type=float, default=1.0, help="Inner slab pressure")

    args = parser.parse_args()

    lvl = int(args.lvl)
    size = float(args.size)

    d, u, v, w, p = generate_kh_fields(
        lvl_max=lvl,
        box_size_cu=size,
        u0=args.u0,
        shear_thickness=args.thickness,
        density_outer=args.rho_outer,
        density_inner=args.rho_inner,
        pressure0=args.p0,
        perturb_eps=args.eps,
        perturb_sigma=args.sigma,
        random_seed=args.seed,
        ndim=args.ndim,
    )

    # Determine output dir
    tag = f"ic_khi_{lvl}_{args.ndim}d"
    if args.outdir is None:
        outdir = Path("ic_khi") / tag
    else:
        outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    # Write primitive hydro fields expected by mini-ramses (grafic input)
    write_array("ic_d", d, size)
    write_array("ic_u", u, size)
    write_array("ic_v", v, size)
    write_array("ic_w", w, size)
    write_array("ic_p", p, size)

    print(f"Kelvin–Helmholtz ICs written to: {outdir}")


if __name__ == "__main__":
    main()



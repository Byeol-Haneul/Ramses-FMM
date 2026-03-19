#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import argparse
import time
from mpi4py import MPI


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Deposit particles onto a 256^3 PM grid using CIC, then compute the "
            "direct gravitational force from PM cells onto a 256x256 slice."
        )
    )
    p.add_argument("--ic-file", type=str, required=True,
                   help="Path to ic_part file with columns: x y z vx vy vz m")
    p.add_argument("--box-cut", type=float, default=300.0,
                   help="Half-box size. Domain is [-box_cut, +box_cut]^3")
    p.add_argument("--nmesh", type=int, default=256,
                   help="PM mesh size per dimension (default: 256)")
    p.add_argument("--z-cell", type=int, default=128,
                   help="Slice index in z (default: 128)")
    p.add_argument("--G", type=float, default=1.0,
                   help="Gravitational constant in your chosen units")
    p.add_argument("--eps", type=float, default=None,
                   help="Softening length. Default = 0.5 * cell size")
    p.add_argument("--out-dir", type=str, default="./txt",
                   help="Output directory")
    p.add_argument("--prefix", type=str, default="pm_direct_force_slice",
                   help="Output filename prefix")
    p.add_argument("--save-density", action="store_true",
                   help="Also save the deposited PM density/mass grid")
    return p.parse_args()


def split_rows(n_rows, size, rank):
    counts = np.full(size, n_rows // size, dtype=int)
    counts[: n_rows % size] += 1
    starts = np.zeros(size, dtype=int)
    starts[1:] = np.cumsum(counts[:-1])
    return starts[rank], starts[rank] + counts[rank], counts, starts


def cic_deposit_3d(x, y, z, m, box_cut, nmesh):
    """
    CIC deposit onto a uniform 3D grid.
    Returns cell mass array of shape (nmesh, nmesh, nmesh).
    """
    dx = 2.0 * box_cut / nmesh
    grid = np.zeros((nmesh, nmesh, nmesh), dtype=np.float64)

    # Convert physical coordinates to cell coordinates
    # Cell centers are at (i + 0.5) * dx - box_cut
    # For CIC, use fractional position relative to cell-center indexing
    gx = (x + box_cut) / dx - 0.5
    gy = (y + box_cut) / dx - 0.5
    gz = (z + box_cut) / dx - 0.5

    i0 = np.floor(gx).astype(np.int64)
    j0 = np.floor(gy).astype(np.int64)
    k0 = np.floor(gz).astype(np.int64)

    tx = gx - i0
    ty = gy - j0
    tz = gz - k0

    for ox in (0, 1):
        wx = (1.0 - tx) if ox == 0 else tx
        ix = i0 + ox

        valid_x = (ix >= 0) & (ix < nmesh)
        if not np.any(valid_x):
            continue

        for oy in (0, 1):
            wy = (1.0 - ty) if oy == 0 else ty
            jy = j0 + oy

            valid_xy = valid_x & (jy >= 0) & (jy < nmesh)
            if not np.any(valid_xy):
                continue

            for oz in (0, 1):
                wz = (1.0 - tz) if oz == 0 else tz
                kz = k0 + oz

                valid = valid_xy & (kz >= 0) & (kz < nmesh)
                if not np.any(valid):
                    continue

                w = wx[valid] * wy[valid] * wz[valid]
                np.add.at(
                    grid,
                    (ix[valid], jy[valid], kz[valid]),
                    m[valid] * w
                )

    return grid


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = parse_args()

    nmesh = args.nmesh
    z_cell = args.z_cell
    box_cut = args.box_cut
    G = args.G
    dx = 2.0 * box_cut / nmesh
    eps = args.eps if args.eps is not None else 0.5 * dx

    if not (0 <= z_cell < nmesh):
        if rank == 0:
            raise ValueError(f"z_cell must satisfy 0 <= z_cell < {nmesh}")
        return

    if rank == 0:
        t0 = time.time()
        print("=" * 72)
        print("PM-CIC deposition + direct force from PM cells onto z-slice")
        print("=" * 72)
        print(f"ic_file   = {args.ic_file}")
        print(f"box_cut   = {box_cut}")
        print(f"nmesh     = {nmesh}")
        print(f"z_cell    = {z_cell}")
        print(f"mpi_size  = {size}")
        print(f"G         = {G}")
        print(f"dx_cell   = {dx}")
        print(f"eps       = {eps}")
        print("=" * 72)

    # --------------------------------------------------
    # Load particles on root
    # --------------------------------------------------
    if rank == 0:
        data = np.loadtxt(args.ic_file)
        x, y, z, vx, vy, vz, m = data.T

        mask = (
            (x >= -box_cut) & (x <= box_cut) &
            (y >= -box_cut) & (y <= box_cut) &
            (z >= -box_cut) & (z <= box_cut)
        )

        x = np.ascontiguousarray(x[mask], dtype=np.float64)
        y = np.ascontiguousarray(y[mask], dtype=np.float64)
        z = np.ascontiguousarray(z[mask], dtype=np.float64)
        m = np.ascontiguousarray(m[mask], dtype=np.float64)

        npart = len(m)

        print(f"Particles kept inside box: {npart}")

        # --------------------------------------------------
        # CIC deposit on root
        # --------------------------------------------------
        t_dep0 = time.time()
        mass_grid = cic_deposit_3d(x, y, z, m, box_cut, nmesh)
        t_dep1 = time.time()

        print(f"CIC deposition done in {t_dep1 - t_dep0:.2f} s")

        # Build compact source list from nonzero cells
        src_idx = np.nonzero(mass_grid > 0.0)
        m_cells = mass_grid[src_idx]

        ix, iy, iz = src_idx
        x_cells = -box_cut + (ix + 0.5) * dx
        y_cells = -box_cut + (iy + 0.5) * dx
        z_cells = -box_cut + (iz + 0.5) * dx

        nsrc = len(m_cells)

        print(f"Nonzero PM source cells: {nsrc}")

        x_eval = -box_cut + (np.arange(nmesh) + 0.5) * dx
        y_eval = -box_cut + (np.arange(nmesh) + 0.5) * dx
        z_eval = -box_cut + (z_cell + 0.5) * dx

        print(f"Slice z_phys = {z_eval}")
    else:
        x_cells = None
        y_cells = None
        z_cells = None
        m_cells = None
        x_eval = None
        y_eval = None
        z_eval = None
        nsrc = None
        mass_grid = None

    # --------------------------------------------------
    # Broadcast compact source list + eval coords
    # --------------------------------------------------
    nsrc = comm.bcast(nsrc, root=0)
    z_eval = comm.bcast(z_eval, root=0)

    if rank != 0:
        x_cells = np.empty(nsrc, dtype=np.float64)
        y_cells = np.empty(nsrc, dtype=np.float64)
        z_cells = np.empty(nsrc, dtype=np.float64)
        m_cells = np.empty(nsrc, dtype=np.float64)
        x_eval = np.empty(nmesh, dtype=np.float64)
        y_eval = np.empty(nmesh, dtype=np.float64)

    comm.Bcast(x_cells, root=0)
    comm.Bcast(y_cells, root=0)
    comm.Bcast(z_cells, root=0)
    comm.Bcast(m_cells, root=0)
    comm.Bcast(x_eval, root=0)
    comm.Bcast(y_eval, root=0)

    # --------------------------------------------------
    # Row decomposition over x
    # --------------------------------------------------
    i0, i1, counts, starts = split_rows(nmesh, size, rank)
    local_nx = i1 - i0

    if rank == 0:
        print("Row distribution:")
        for r in range(size):
            s = starts[r]
            e = starts[r] + counts[r] - 1
            if counts[r] > 0:
                print(f"  rank {r}: rows {s} to {e}")
            else:
                print(f"  rank {r}: no rows")

    # --------------------------------------------------
    # Local force computation
    # --------------------------------------------------
    local_fx = np.zeros((local_nx, nmesh), dtype=np.float64)
    local_fy = np.zeros((local_nx, nmesh), dtype=np.float64)
    local_fz = np.zeros((local_nx, nmesh), dtype=np.float64)

    # These depend only on y and z for the whole local block
    dy = y_eval[:, None] - y_cells[None, :]
    dz = z_eval - z_cells
    dy2 = dy * dy
    dz2 = dz * dz
    eps2 = eps * eps

    for local_i, global_i in enumerate(range(i0, i1)):
        xi = x_eval[global_i]
        dxv = xi - x_cells
        dx2 = dxv * dxv

        r2 = dx2[None, :] + dy2 + dz2[None, :] + eps2
        inv_r3 = 1.0 / (r2 * np.sqrt(r2))

        # g = -G * sum_j m_j * (r - r_j) / |r-r_j|^3
        local_fx[local_i, :] = -G * np.sum(m_cells[None, :] * dxv[None, :] * inv_r3, axis=1)
        local_fy[local_i, :] = -G * np.sum(m_cells[None, :] * dy * inv_r3, axis=1)
        local_fz[local_i, :] = -G * np.sum(m_cells[None, :] * dz[None, :] * inv_r3, axis=1)

        if local_nx > 0 and (local_i % max(1, local_nx // 10) == 0 or local_i == local_nx - 1):
            print(f"[rank {rank}] done {local_i + 1}/{local_nx} local rows", flush=True)

    # --------------------------------------------------
    # Gather on root
    # --------------------------------------------------
    recvcounts = counts * nmesh
    displs = starts * nmesh

    if rank == 0:
        fx = np.empty((nmesh, nmesh), dtype=np.float64)
        fy = np.empty((nmesh, nmesh), dtype=np.float64)
        fz = np.empty((nmesh, nmesh), dtype=np.float64)

        recvbuf_fx = [fx.ravel(), recvcounts, displs, MPI.DOUBLE]
        recvbuf_fy = [fy.ravel(), recvcounts, displs, MPI.DOUBLE]
        recvbuf_fz = [fz.ravel(), recvcounts, displs, MPI.DOUBLE]
    else:
        fx = fy = fz = None
        recvbuf_fx = recvbuf_fy = recvbuf_fz = None

    comm.Gatherv(local_fx.ravel(), recvbuf_fx, root=0)
    comm.Gatherv(local_fy.ravel(), recvbuf_fy, root=0)
    comm.Gatherv(local_fz.ravel(), recvbuf_fz, root=0)

    # --------------------------------------------------
    # Save on root
    # --------------------------------------------------
    if rank == 0:
        os.makedirs(args.out_dir, exist_ok=True)

        X, Y = np.meshgrid(x_eval, y_eval, indexing="ij")
        Fmag = np.sqrt(fx**2 + fy**2 + fz**2)

        out_force = os.path.join(
            args.out_dir,
            f"{args.prefix}_nmesh{nmesh}_zcell{z_cell}.txt"
        )

        np.savetxt(
            out_force,
            np.column_stack([
                X.ravel(),
                Y.ravel(),
                fx.ravel(),
                fy.ravel(),
                fz.ravel(),
                Fmag.ravel()
            ]),
            header="x  y  Fx  Fy  Fz  |F|",
            fmt="%.10e"
        )

        if args.save_density:
            out_rho = os.path.join(
                args.out_dir,
                f"{args.prefix}_massgrid_nmesh{nmesh}.npy"
            )
            np.save(out_rho, mass_grid)

        t1 = time.time()
        print("=" * 72)
        print(f"Saved force slice to {out_force}")
        if args.save_density:
            print(f"Saved PM mass grid to {out_rho}")
        print(f"Elapsed time: {t1 - t0:.2f} s")
        print("=" * 72)


if __name__ == "__main__":
    main()

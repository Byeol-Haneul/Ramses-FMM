#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import argparse
import time
from mpi4py import MPI


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute direct particle-sum potential and force on a 2D slice using mpi4py."
    )
    p.add_argument("--ic-file", type=str, required=True,
                   help="Path to ic_part file with columns: x y z vx vy vz m")
    p.add_argument("--box-cut", type=float, default=300.0,
                   help="Half-box size. Domain is [-box_cut, +box_cut]")
    p.add_argument("--lvl", type=int, default=9,
                   help="Level used to define grid size: n_points = 2**(lvl-1)")
    p.add_argument("--z-cell", type=int, default=None,
                   help="Slice index in z. Default: center slice")
    p.add_argument("--G", type=float, default=1.0,
                   help="Gravitational constant in your chosen units")
    p.add_argument("--out-dir", type=str, default="./txt",
                   help="Output directory")
    p.add_argument("--prefix", type=str, default="direct_force_phi_slice",
                   help="Output filename prefix")
    return p.parse_args()


def split_rows(n_rows, size, rank):
    counts = np.full(size, n_rows // size, dtype=int)
    counts[: n_rows % size] += 1
    starts = np.zeros(size, dtype=int)
    starts[1:] = np.cumsum(counts[:-1])
    return starts[rank], starts[rank] + counts[rank], counts, starts


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = parse_args()

    if rank == 0:
        t0 = time.time()
        print("=" * 60)
        print("MPI direct particle-sum slice")
        print("=" * 60)
        print(f"ic_file   = {args.ic_file}")
        print(f"box_cut   = {args.box_cut}")
        print(f"lvl       = {args.lvl}")
        print(f"mpi_size  = {size}")
        print(f"G         = {args.G}")
        print("=" * 60)

    ic_file = args.ic_file
    box_cut = args.box_cut
    G = args.G
    lvl = args.lvl
    n_points = 2 ** (lvl - 1)
    n_cells = n_points
    z_cell = args.z_cell if args.z_cell is not None else n_cells // 2

    # -------------------------------
    # Load particle data on root
    # -------------------------------
    if rank == 0:
        data = np.loadtxt(ic_file)
        x, y, z, vx, vy, vz, m = data.T

        mask = (
            (x >= -box_cut) & (x <= box_cut) &
            (y >= -box_cut) & (y <= box_cut) &
            (z >= -box_cut) & (z <= box_cut)
        )

        x_particles = np.ascontiguousarray(x[mask], dtype=np.float64)
        y_particles = np.ascontiguousarray(y[mask], dtype=np.float64)
        z_particles = np.ascontiguousarray(z[mask], dtype=np.float64)
        m_particles = np.ascontiguousarray(m[mask], dtype=np.float64)

        npart = len(m_particles)
        dx_cell = 2.0 * box_cut / n_cells
        x_eval_phys = -box_cut + (np.arange(n_points) + 0.5) * dx_cell
        y_eval_phys = -box_cut + (np.arange(n_points) + 0.5) * dx_cell
        z_eval_phys = -box_cut + (z_cell + 0.5) * dx_cell

        print(f"Particles kept inside box: {npart}")
        print(f"n_points  = {n_points}")
        print(f"z_cell    = {z_cell}")
        print(f"dx_cell   = {dx_cell}")
        print(f"z_phys    = {z_eval_phys}")
    else:
        x_particles = None
        y_particles = None
        z_particles = None
        m_particles = None
        x_eval_phys = None
        y_eval_phys = None
        z_eval_phys = None
        npart = None

    # -------------------------------
    # Broadcast metadata and arrays
    # -------------------------------
    npart = comm.bcast(npart, root=0)
    z_eval_phys = comm.bcast(z_eval_phys, root=0)

    if rank != 0:
        x_particles = np.empty(npart, dtype=np.float64)
        y_particles = np.empty(npart, dtype=np.float64)
        z_particles = np.empty(npart, dtype=np.float64)
        m_particles = np.empty(npart, dtype=np.float64)
        x_eval_phys = np.empty(n_points, dtype=np.float64)
        y_eval_phys = np.empty(n_points, dtype=np.float64)

    comm.Bcast(x_particles, root=0)
    comm.Bcast(y_particles, root=0)
    comm.Bcast(z_particles, root=0)
    comm.Bcast(m_particles, root=0)
    comm.Bcast(x_eval_phys, root=0)
    comm.Bcast(y_eval_phys, root=0)

    # -------------------------------
    # Row decomposition
    # -------------------------------
    i0, i1, counts, starts = split_rows(n_points, size, rank)
    local_nx = i1 - i0

    if rank == 0:
        print("Row distribution:")
        for r in range(size):
            s = starts[r]
            e = starts[r] + counts[r] - 1
            print(f"  rank {r}: rows {s} to {e}" if counts[r] > 0 else f"  rank {r}: no rows")

    # -------------------------------
    # Local computation
    # -------------------------------
    local_phi = np.zeros((local_nx, n_points), dtype=np.float64)
    local_fx = np.zeros((local_nx, n_points), dtype=np.float64)
    local_fy = np.zeros((local_nx, n_points), dtype=np.float64)
    local_fz = np.zeros((local_nx, n_points), dtype=np.float64)

    dz = z_eval_phys - z_particles
    dz2 = dz ** 2
    dy = y_eval_phys[:, None] - y_particles[None, :]
    dy2 = dy ** 2

    for local_i, global_i in enumerate(range(i0, i1)):
        xi = x_eval_phys[global_i]

        dx = xi - x_particles
        dx2 = dx ** 2

        r2 = dx2[None, :] + dy2 + dz2[None, :]
        r = np.sqrt(r2)

        # avoid singularity
        r[r == 0.0] = 1e-12

        inv_r = 1.0 / r
        inv_r3 = inv_r / r2
        inv_r3[r2 == 0.0] = 0.0

        local_phi[local_i, :] = -G * np.sum(m_particles[None, :] * inv_r, axis=1)

        local_fx[local_i, :] = -G * np.sum(
            m_particles[None, :] * dx[None, :] * inv_r3, axis=1
        )
        local_fy[local_i, :] = -G * np.sum(
            m_particles[None, :] * dy * inv_r3, axis=1
        )
        local_fz[local_i, :] = -G * np.sum(
            m_particles[None, :] * dz[None, :] * inv_r3, axis=1
        )

        if local_nx > 0 and (local_i % max(1, local_nx // 10) == 0 or local_i == local_nx - 1):
            print(f"[rank {rank}] done {local_i + 1}/{local_nx} local rows", flush=True)

    # -------------------------------
    # Gather on root
    # -------------------------------
    recvcounts = counts * n_points
    displs = starts * n_points

    send_fx = local_fx.ravel()
    send_fy = local_fy.ravel()
    send_fz = local_fz.ravel()
    send_phi = local_phi.ravel()

    if rank == 0:
        fx = np.empty((n_points, n_points), dtype=np.float64)
        fy = np.empty((n_points, n_points), dtype=np.float64)
        fz = np.empty((n_points, n_points), dtype=np.float64)
        phi = np.empty((n_points, n_points), dtype=np.float64)

        recv_fx = [fx.ravel(), recvcounts, displs, MPI.DOUBLE]
        recv_fy = [fy.ravel(), recvcounts, displs, MPI.DOUBLE]
        recv_fz = [fz.ravel(), recvcounts, displs, MPI.DOUBLE]
        recv_phi = [phi.ravel(), recvcounts, displs, MPI.DOUBLE]
    else:
        fx = fy = fz = phi = None
        recv_fx = recv_fy = recv_fz = recv_phi = None

    comm.Gatherv(send_fx, recv_fx, root=0)
    comm.Gatherv(send_fy, recv_fy, root=0)
    comm.Gatherv(send_fz, recv_fz, root=0)
    comm.Gatherv(send_phi, recv_phi, root=0)

    # -------------------------------
    # Save on root
    # -------------------------------
    if rank == 0:
        os.makedirs(args.out_dir, exist_ok=True)

        out_file = os.path.join(
            args.out_dir,
            f"{args.prefix}_lvl{lvl}_zcell{z_cell}.txt"
        )

        X, Y = np.meshgrid(x_eval_phys, y_eval_phys, indexing="ij")
        Z = np.full_like(X, z_eval_phys)

        np.savetxt(
            out_file,
            np.column_stack([
                X.ravel(),
                Y.ravel(),
                Z.ravel(),
                fx.ravel(),
                fy.ravel(),
                fz.ravel(),
                phi.ravel()
            ]),
            header="x  y  z  fx  fy  fz  phi",
            fmt="%.10e"
        )

        t1 = time.time()
        print("=" * 60)
        print(f"Saved direct slice to {out_file}")
        print(f"Elapsed time: {t1 - t0:.2f} s")
        print("=" * 60)


if __name__ == "__main__":
    main()
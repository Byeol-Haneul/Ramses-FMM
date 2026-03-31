from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class PhiRegionStats:
    count: int
    mean: float
    median: float
    p95: float
    max: float


def _compute_fine_raster_indices(u, v, dx, box):
    xmin, xmax, ymin, ymax = box
    fine_dx = float(np.min(dx))
    nx = int(round((xmax - xmin) / fine_dx))
    ny = int(round((ymax - ymin) / fine_dx))

    id_grid = -np.ones((ny, nx), dtype=np.int32)
    lev_grid = -np.ones((ny, nx), dtype=np.int16)

    return fine_dx, nx, ny, id_grid, lev_grid


def audit_phi_amr_slice(
    s,
    *,
    box,
    sphere_centers,
    sphere_radii,
    surface_factor=1.0,
):
    """Classify a slice into AMR-interface / edge / interior / sphere-surface regions."""
    u = np.asarray(s["u"])
    v = np.asarray(s["v"])
    dx = np.asarray(s["dx"])
    level = np.asarray(s.get("level_plot", s["level"]))
    phi_sim = np.asarray(s["phi_sim"])
    phi_ana = np.asarray(s["phi_ana"])
    relerr = np.asarray(s["phi_relerr"])
    abserr = np.abs(relerr)

    fine_dx, nx, ny, id_grid, lev_grid = _compute_fine_raster_indices(u, v, dx, box)
    xmin, xmax, ymin, ymax = box

    for i, (uu, vv, dd, lev) in enumerate(zip(u, v, dx, level)):
        ix0 = max(int(round((uu - dd / 2.0 - xmin) / fine_dx)), 0)
        ix1 = min(int(round((uu + dd / 2.0 - xmin) / fine_dx)), nx)
        iy0 = max(int(round((vv - dd / 2.0 - ymin) / fine_dx)), 0)
        iy1 = min(int(round((vv + dd / 2.0 - ymin) / fine_dx)), ny)
        id_grid[iy0:iy1, ix0:ix1] = i
        lev_grid[iy0:iy1, ix0:ix1] = lev

    interface_grid = np.zeros_like(lev_grid, dtype=bool)
    jump_x = np.not_equal(lev_grid[:, 1:], lev_grid[:, :-1])
    jump_x &= (lev_grid[:, 1:] >= 0) & (lev_grid[:, :-1] >= 0)
    interface_grid[:, 1:] |= jump_x
    interface_grid[:, :-1] |= jump_x

    jump_y = np.not_equal(lev_grid[1:, :], lev_grid[:-1, :])
    jump_y &= (lev_grid[1:, :] >= 0) & (lev_grid[:-1, :] >= 0)
    interface_grid[1:, :] |= jump_y
    interface_grid[:-1, :] |= jump_y

    interface_ids = np.unique(id_grid[interface_grid & (id_grid >= 0)])
    interface_cell = np.zeros(len(u), dtype=bool)
    interface_cell[interface_ids] = True

    edge_grid = np.zeros_like(lev_grid, dtype=bool)
    edge_grid[0, :] = True
    edge_grid[-1, :] = True
    edge_grid[:, 0] = True
    edge_grid[:, -1] = True

    edge_ids = np.unique(id_grid[edge_grid & (id_grid >= 0)])
    edge_cell = np.zeros(len(u), dtype=bool)
    edge_cell[edge_ids] = True

    x3d = np.asarray(s["x3d"])
    y3d = np.asarray(s["y3d"])
    z3d = np.asarray(s["z3d"])
    near_surface = np.zeros(len(u), dtype=bool)
    for center, radius in zip(sphere_centers, sphere_radii):
        rr = np.sqrt((x3d - center[0]) ** 2 + (y3d - center[1]) ** 2 + (z3d - center[2]) ** 2)
        near_surface |= np.abs(rr - radius) <= surface_factor * dx

    interior_cell = ~(interface_cell | edge_cell)

    return {
        "u": u,
        "v": v,
        "dx": dx,
        "level": level,
        "phi_sim": phi_sim,
        "phi_ana": phi_ana,
        "relerr": relerr,
        "abserr": abserr,
        "fine_dx": fine_dx,
        "box": box,
        "id_grid": id_grid,
        "lev_grid": lev_grid,
        "interface_grid": interface_grid,
        "interface_cell": interface_cell,
        "edge_cell": edge_cell,
        "interior_cell": interior_cell,
        "near_surface": near_surface,
    }


def _region_stats(values):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return PhiRegionStats(0, np.nan, np.nan, np.nan, np.nan)
    return PhiRegionStats(
        count=int(finite.size),
        mean=float(np.mean(finite)),
        median=float(np.median(finite)),
        p95=float(np.quantile(finite, 0.95)),
        max=float(np.max(finite)),
    )


def summarize_phi_amr_audit(audit):
    """Return a dict of region summaries based on |phi relative error|."""
    abserr = np.asarray(audit["abserr"])
    interface = np.asarray(audit["interface_cell"])
    edge = np.asarray(audit["edge_cell"])
    interior = np.asarray(audit["interior_cell"])
    surface = np.asarray(audit["near_surface"])

    summaries = {
        "all": _region_stats(abserr),
        "interface": _region_stats(abserr[interface]),
        "interior": _region_stats(abserr[interior]),
        "edge": _region_stats(abserr[edge]),
        "surface": _region_stats(abserr[surface]),
        "interface_not_surface": _region_stats(abserr[interface & ~surface]),
        "interior_not_surface": _region_stats(abserr[interior & ~surface]),
        "interior_surface": _region_stats(abserr[interior & surface]),
    }
    return summaries


def print_phi_amr_audit_summary(audit, *, top_n=10):
    summaries = summarize_phi_amr_audit(audit)
    print(f"slice cells: {len(audit['u'])}, fine dx: {audit['fine_dx']:.6f}")
    for key, stat in summaries.items():
        if stat.count == 0:
            print(f"{key:>22s}: none")
            continue
        print(
            f"{key:>22s}: n={stat.count:5d} "
            f"mean={stat.mean:.3e} median={stat.median:.3e} "
            f"p95={stat.p95:.3e} max={stat.max:.3e}"
        )

    mask = np.isfinite(audit["abserr"])
    idx = np.argsort(np.where(mask, audit["abserr"], -1.0))[-top_n:][::-1]
    print("top phi-error cells:")
    for i in idx:
        print(
            f"  x={audit['u'][i]:.6f} y={audit['v'][i]:.6f} "
            f"dx={audit['dx'][i]:.6f} L={int(audit['level'][i])} "
            f"rel={audit['relerr'][i]:+.3e} "
            f"interface={bool(audit['interface_cell'][i])} "
            f"edge={bool(audit['edge_cell'][i])} "
            f"surface={bool(audit['near_surface'][i])}"
        )


def plot_phi_amr_audit(audit, *, figsize=(15, 5)):
    box = audit["box"]
    size = np.maximum(4, 7000.0 * audit["dx"] / (box[1] - box[0]))

    abs_rel = np.abs(audit["relerr"])
    finite_abs = abs_rel[np.isfinite(abs_rel)]
    if finite_abs.size:
        err_lim = float(np.quantile(finite_abs, 0.995))
        if err_lim == 0.0:
            err_lim = 1.0e-16
    else:
        err_lim = 1.0e-16

    interface_abs = abs_rel[audit["interface_cell"]]
    interface_abs = interface_abs[np.isfinite(interface_abs)]
    if interface_abs.size:
        interface_lim = float(np.quantile(interface_abs, 0.995))
        if interface_lim == 0.0:
            interface_lim = 1.0e-16
    else:
        interface_lim = 1.0e-16

    fig, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)

    sc0 = axes[0].scatter(
        audit["u"], audit["v"], c=audit["level"], s=size,
        marker="s", cmap="viridis", linewidths=0
    )
    axes[0].set_title("AMR level")

    sc1 = axes[1].scatter(
        audit["u"], audit["v"], c=audit["relerr"], s=size,
        marker="s", cmap="RdBu_r", linewidths=0,
        vmin=-err_lim, vmax=err_lim
    )
    axes[1].set_title("phi relative error")

    sc2 = axes[2].scatter(
        audit["u"][audit["interface_cell"]],
        audit["v"][audit["interface_cell"]],
        c=audit["abserr"][audit["interface_cell"]],
        s=np.maximum(4, 7000.0 * audit["dx"][audit["interface_cell"]] / (box[1] - box[0])),
        marker="s", cmap="magma", linewidths=0,
        vmin=0.0, vmax=interface_lim
    )
    axes[2].set_title("|phi relerr| on interface cells")

    for ax in axes:
        ax.set_xlim(box[0], box[1])
        ax.set_ylim(box[2], box[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(sc0, ax=axes[0], fraction=0.046)
    fig.colorbar(sc1, ax=axes[1], fraction=0.046)
    fig.colorbar(sc2, ax=axes[2], fraction=0.046)
    return fig, axes


def plot_phi_amr_region_histogram(audit, *, figsize=(7, 4)):
    fig, ax = plt.subplots(figsize=figsize)
    abserr = np.asarray(audit["abserr"])
    regions = [
        ("interior", audit["interior_cell"], "tab:blue"),
        ("interface", audit["interface_cell"], "tab:orange"),
        ("edge", audit["edge_cell"], "tab:green"),
        ("surface", audit["near_surface"], "tab:red"),
    ]

    for label, mask, color in regions:
        vals = abserr[np.asarray(mask)]
        vals = vals[np.isfinite(vals) & (vals > 0.0)]
        if vals.size == 0:
            continue
        bins = np.logspace(np.log10(np.min(vals)), np.log10(np.max(vals)), 40)
        ax.hist(vals, bins=bins, histtype="step", lw=1.5, label=f"{label} (n={vals.size})", color=color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("|phi relative error|")
    ax.set_ylabel("count")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    return fig, ax

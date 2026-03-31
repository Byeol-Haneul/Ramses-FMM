from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "utils/py/miniramses.py").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root from the current path.")


try:
    HERE = Path(__file__).resolve().parent
except NameError:
    HERE = Path.cwd()

REPO_ROOT = find_repo_root(HERE)
if str(REPO_ROOT / "utils/py") not in sys.path:
    sys.path.append(str(REPO_ROOT / "utils/py"))

import miniramses as ram

from amr_phi_interface_audit import audit_phi_amr_slice, summarize_phi_amr_audit


PLOT_BOX = [1.0, 3.0, 1.0, 3.0]
SLICE_AXIS = "z"
SLICE_CENTER = 2.0
SPHERE_CENTERS = (
    np.array([1.7, 2.0, 2.0]),
    np.array([2.2, 2.0, 2.0]),
)
SPHERE_RADII = (0.1, 0.2)


@dataclass
class CompareSummary:
    ncell: int
    time_mg: float
    time_fmm: float
    time_delta: float
    mean_delta: float
    median_delta: float
    std_delta: float
    max_abs_delta: float
    mean_abs_rel: float
    median_abs_rel: float
    p95_abs_rel: float
    max_abs_rel: float


def slice_cells(c, axis=SLICE_AXIS, center=SLICE_CENTER):
    axis_to_index = {"x": 0, "y": 1, "z": 2}
    axis_idx = axis_to_index[axis]
    mask = np.abs(c.x[axis_idx] - center) <= c.dx / 2.0

    data = {
        "x3d": np.asarray(c.x[0][mask], dtype=np.float64),
        "y3d": np.asarray(c.x[1][mask], dtype=np.float64),
        "z3d": np.asarray(c.x[2][mask], dtype=np.float64),
        "dx": np.asarray(c.dx[mask], dtype=np.float64),
        "level": np.asarray(c.level[mask], dtype=np.int32),
        "level_plot": np.asarray(c.level[mask] + 1, dtype=np.int32),
        "phi": np.asarray(c.g[0][mask], dtype=np.float64),
    }

    if axis == "x":
        data["u"] = data["y3d"]
        data["v"] = data["z3d"]
    elif axis == "y":
        data["u"] = data["x3d"]
        data["v"] = data["z3d"]
    else:
        data["u"] = data["x3d"]
        data["v"] = data["y3d"]
    return data


def load_output(run_dir: Path, nout: int):
    c = ram.rd_cell(str(nout), path=str(run_dir), prefix="grav")
    info = ram.rd_info(str(nout), path=str(run_dir))
    return slice_cells(c), info


def structured_keys(s):
    n = len(s["phi"])
    key = np.empty(
        n,
        dtype=[
            ("x", np.int64),
            ("y", np.int64),
            ("z", np.int64),
            ("dx", np.int64),
            ("level", np.int32),
        ],
    )
    scale = 10**12
    key["x"] = np.rint(s["x3d"] * scale).astype(np.int64)
    key["y"] = np.rint(s["y3d"] * scale).astype(np.int64)
    key["z"] = np.rint(s["z3d"] * scale).astype(np.int64)
    key["dx"] = np.rint(s["dx"] * scale).astype(np.int64)
    key["level"] = s["level"]
    return key


def align_slices(mg, fmm):
    key_mg = structured_keys(mg)
    key_fmm = structured_keys(fmm)

    order_mg = np.argsort(key_mg, order=("level", "dx", "z", "y", "x"))
    order_fmm = np.argsort(key_fmm, order=("level", "dx", "z", "y", "x"))

    if len(order_mg) != len(order_fmm):
        raise ValueError(f"Cell-count mismatch: MG={len(order_mg)} FMM={len(order_fmm)}")
    if not np.array_equal(key_mg[order_mg], key_fmm[order_fmm]):
        raise ValueError("MG/FMM slices do not share the same cell keys.")

    mg_aligned = {k: np.asarray(v)[order_mg] for k, v in mg.items()}
    fmm_aligned = {k: np.asarray(v)[order_fmm] for k, v in fmm.items()}
    return mg_aligned, fmm_aligned


def build_compare_slice(mg, fmm):
    delta = fmm["phi"] - mg["phi"]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(np.abs(mg["phi"]) > 0.0, delta / mg["phi"], np.nan)
    compare = {
        "x3d": mg["x3d"],
        "y3d": mg["y3d"],
        "z3d": mg["z3d"],
        "u": mg["u"],
        "v": mg["v"],
        "dx": mg["dx"],
        "level": mg["level"],
        "level_plot": mg["level_plot"],
        "phi_sim": fmm["phi"],
        "phi_ana": mg["phi"],
        "phi_relerr": rel,
        "phi_delta": delta,
    }
    return compare


def summarize_compare(info_mg, info_fmm, compare):
    delta = compare["phi_delta"]
    rel = compare["phi_relerr"]
    finite_rel = np.abs(rel[np.isfinite(rel)])
    return CompareSummary(
        ncell=len(delta),
        time_mg=float(info_mg.time),
        time_fmm=float(info_fmm.time),
        time_delta=float(info_fmm.time - info_mg.time),
        mean_delta=float(np.mean(delta)),
        median_delta=float(np.median(delta)),
        std_delta=float(np.std(delta)),
        max_abs_delta=float(np.max(np.abs(delta))),
        mean_abs_rel=float(np.mean(finite_rel)),
        median_abs_rel=float(np.median(finite_rel)),
        p95_abs_rel=float(np.quantile(finite_rel, 0.95)),
        max_abs_rel=float(np.max(finite_rel)),
    )


def print_summary(summary, audit):
    region = summarize_phi_amr_audit(audit)
    print(f"cells compared: {summary.ncell}")
    print(
        "times: "
        f"mg={summary.time_mg:.16e} "
        f"fmm={summary.time_fmm:.16e} "
        f"delta={summary.time_delta:+.3e}"
    )
    print(
        "raw phi delta (FMM-MG): "
        f"mean={summary.mean_delta:+.3e} "
        f"median={summary.median_delta:+.3e} "
        f"std={summary.std_delta:.3e} "
        f"maxabs={summary.max_abs_delta:.3e}"
    )
    print(
        "relative phi delta vs MG: "
        f"meanabs={summary.mean_abs_rel:.3e} "
        f"medianabs={summary.median_abs_rel:.3e} "
        f"p95abs={summary.p95_abs_rel:.3e} "
        f"maxabs={summary.max_abs_rel:.3e}"
    )
    for label in ("interface", "interior", "edge", "surface", "interface_not_surface", "interior_not_surface"):
        stat = region[label]
        if stat.count == 0:
            print(f"{label}: none")
        else:
            print(
                f"{label}: n={stat.count} "
                f"meanabs={stat.mean:.3e} "
                f"medianabs={stat.median:.3e} "
                f"p95abs={stat.p95:.3e} "
                f"maxabs={stat.max:.3e}"
            )


def plot_compare_maps(compare, *, figsize=(15, 5)):
    size = np.maximum(4, 7000.0 * compare["dx"] / (PLOT_BOX[1] - PLOT_BOX[0]))

    delta = compare["phi_delta"]
    rel = compare["phi_relerr"]

    finite_delta = np.abs(delta[np.isfinite(delta)])
    delta_lim = float(np.quantile(finite_delta, 0.995)) if finite_delta.size else 1.0
    if delta_lim == 0.0:
        delta_lim = 1.0e-16

    finite_rel = np.abs(rel[np.isfinite(rel)])
    rel_lim = float(np.quantile(finite_rel, 0.995)) if finite_rel.size else 1.0
    if rel_lim == 0.0:
        rel_lim = 1.0e-16

    fig, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)

    sc0 = axes[0].scatter(compare["u"], compare["v"], c=compare["level_plot"], s=size, marker="s", cmap="viridis", linewidths=0)
    axes[0].set_title("AMR level")

    sc1 = axes[1].scatter(compare["u"], compare["v"], c=delta, s=size, marker="s", cmap="RdBu_r", linewidths=0, vmin=-delta_lim, vmax=delta_lim)
    axes[1].set_title("phi delta (FMM-MG)")

    sc2 = axes[2].scatter(compare["u"], compare["v"], c=rel, s=size, marker="s", cmap="RdBu_r", linewidths=0, vmin=-rel_lim, vmax=rel_lim)
    axes[2].set_title("relative phi delta vs MG")

    for ax in axes:
        ax.set_xlim(PLOT_BOX[0], PLOT_BOX[1])
        ax.set_ylim(PLOT_BOX[2], PLOT_BOX[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(sc0, ax=axes[0], fraction=0.046)
    fig.colorbar(sc1, ax=axes[1], fraction=0.046)
    fig.colorbar(sc2, ax=axes[2], fraction=0.046)
    return fig, axes


def plot_compare_histogram(compare, audit, *, figsize=(7, 4)):
    rel = compare["phi_relerr"]
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for label, mask, color in (
        ("interior", audit["interior_cell"], "tab:blue"),
        ("interface", audit["interface_cell"], "tab:orange"),
        ("edge", audit["edge_cell"], "tab:green"),
        ("surface", audit["near_surface"], "tab:red"),
    ):
        vals = np.abs(rel[np.asarray(mask)])
        vals = vals[np.isfinite(vals) & (vals > 0.0)]
        if vals.size == 0:
            continue
        bins = np.logspace(np.log10(np.min(vals)), np.log10(np.max(vals)), 40)
        ax.hist(vals, bins=bins, histtype="step", lw=1.5, label=f"{label} (n={vals.size})", color=color)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("|(phi_fmm-phi_mg)/phi_mg|")
    ax.set_ylabel("count")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    return fig, ax


def save_plots(compare, audit, save_dir: Path):
    save_dir.mkdir(parents=True, exist_ok=True)
    fig, _ = plot_compare_maps(compare)
    fig.savefig(save_dir / "mg_fmm_phi_compare_maps.png", dpi=180)
    plt.close(fig)

    fig, _ = plot_compare_histogram(compare, audit)
    fig.savefig(save_dir / "mg_fmm_phi_compare_hist.png", dpi=180)
    plt.close(fig)


def main():
    run_root = REPO_ROOT / "tmp" / "mg_vs_fmm"
    mg_dir = run_root / "mg"
    fmm_dir = run_root / "fmm"
    save_dir = run_root / "analysis_phi_compare"

    mg, info_mg = load_output(mg_dir, 2)
    fmm, info_fmm = load_output(fmm_dir, 2)
    mg, fmm = align_slices(mg, fmm)
    compare = build_compare_slice(mg, fmm)
    audit = audit_phi_amr_slice(
        compare,
        box=PLOT_BOX,
        sphere_centers=SPHERE_CENTERS,
        sphere_radii=SPHERE_RADII,
        surface_factor=1.0,
    )
    summary = summarize_compare(info_mg, info_fmm, compare)

    save_dir.mkdir(parents=True, exist_ok=True)
    with (save_dir / "summary.txt").open("w") as fh:
        old_stdout = sys.stdout
        try:
            sys.stdout = fh
            print_summary(summary, audit)
        finally:
            sys.stdout = old_stdout

    print_summary(summary, audit)
    save_plots(compare, audit, save_dir)


if __name__ == "__main__":
    main()

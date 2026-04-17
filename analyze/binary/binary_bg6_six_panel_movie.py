#!/usr/bin/env python3
"""
Build one synchronized 6-panel BG6 movie:
  rows = levels (e.g., lmax 6, 7, 8)
  cols = solvers (MG, FMM)

The colorbar is fixed for the full movie.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from binary_bg6_grid_movie import (
    REPO_ROOT,
    common_outputs,
    get_box_from_particles,
    infer_value_limits,
    read_frame,
    stitch_gif,
)


def level_num(level: str) -> int:
    token = level.split("_")[-1]
    if token.startswith("l"):
        token = token[1:]
    return int(token)


def build_cases(binary_root: Path, levels: list[str], box_margin: float):
    cases = []
    for level in levels:
        run_mg = binary_root / "mg" / f"equal_mass_binary_amr_{level}"
        run_fmm = binary_root / "fmm" / f"equal_mass_binary_amr_{level}"
        if not run_mg.exists() or not run_fmm.exists():
            raise FileNotFoundError(f"Missing run directories for level {level}: {run_mg}, {run_fmm}")

        outs = common_outputs(run_mg, run_fmm)
        box = get_box_from_particles([run_mg, run_fmm], outputs=[outs[0], outs[-1]], margin=box_margin)
        cases.append(
            {
                "level": level,
                "run_mg": run_mg,
                "run_fmm": run_fmm,
                "outputs": outs,
                "box": box,
            }
        )
    return cases


def common_nouts_across_cases(cases: list[dict]) -> list[int]:
    sets = [set(case["outputs"]) for case in cases]
    shared = sorted(set.intersection(*sets))
    if not shared:
        raise RuntimeError("No common output indices shared across all requested levels.")
    return shared


def infer_global_limits(cases: list[dict], nouts: list[int], z_center: float) -> tuple[float, float]:
    mins = []
    maxs = []
    for case in cases:
        vmin, vmax = infer_value_limits(
            [case["run_mg"], case["run_fmm"]],
            nouts,
            box=case["box"],
            z_center=z_center,
        )
        mins.append(vmin)
        maxs.append(vmax)
    return float(min(mins)), float(max(maxs))


def render_six_panel_frame(
    cases: list[dict],
    nout: int,
    npix: int,
    z_center: float,
    vmin: float,
    vmax: float,
    out_png: Path,
    cmap: str = "magma",
):
    nrows = len(cases)
    fig, axes = plt.subplots(nrows, 2, figsize=(12.2, max(4.2, 3.35 * nrows)), constrained_layout=True)
    if nrows == 1:
        axes = np.array([axes])

    im = None

    for i, case in enumerate(cases):
        level = case["level"]
        lnum = level_num(level)

        mg = read_frame(case["run_mg"], nout, box=case["box"], npix=npix, z_center=z_center, with_grid=True)
        fmm = read_frame(case["run_fmm"], nout, box=case["box"], npix=npix, z_center=z_center, with_grid=True)
        for j, (solver, frame_data) in enumerate((("MG", mg), ("FMM", fmm))):
            ax = axes[i, j]
            box = case["box"]
            im = ax.imshow(
                frame_data["image"],
                origin="lower",
                extent=[box[0], box[1], box[2], box[3]],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )
            if frame_data["segments"]:
                lc = LineCollection(
                    frame_data["segments"],
                    colors=(0, 0, 0, 0.72),
                    linewidths=frame_data["widths"],
                )
                ax.add_collection(lc)

            ax.scatter(
                frame_data["particle_x"],
                frame_data["particle_y"],
                s=28,
                c=["#00D9FF", "#FF4D4D"],
                edgecolors="white",
                linewidths=0.55,
                zorder=6,
            )
            ax.scatter(
                frame_data["analytic_x"],
                frame_data["analytic_y"],
                s=120,
                facecolors="none",
                edgecolors=["#00D9FF", "#FF4D4D"],
                linewidths=1.5,
                zorder=7,
            )
            ax.set_aspect("equal")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title(f"L{lnum} {solver} | out={nout:05d} | t={frame_data['time']:.6f}", fontsize=10.5)

    assert im is not None
    fig.suptitle(
        f"Binary BG6 6-Panel | out={nout:05d} | z={z_center:.3f} | circles=analytic",
        fontsize=12,
    )
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92, pad=0.012)
    cbar.set_label(r"$\Phi$")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Build synchronized 6-panel BG6 MG/FMM movie.")
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["bg6_l6", "bg6_l7", "bg6_l8"],
        help="Level suffixes; one row per level.",
    )
    parser.add_argument("--npix", type=int, default=560, help="Per-panel raster resolution.")
    parser.add_argument("--fps", type=int, default=12, help="GIF frame rate.")
    parser.add_argument("--z-center", type=float, default=1.0, help="Z slice center.")
    parser.add_argument("--box-margin", type=float, default=0.08, help="Margin around particle orbit for panel boxes.")
    parser.add_argument(
        "--binary-root",
        type=Path,
        default=REPO_ROOT / "binary_test",
        help="Input root containing mg/ and fmm/ run folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "analyze" / "binary" / "results" / "movies_bg6" / "six_panel",
        help="Output directory for 6-panel movie and frames.",
    )
    parser.add_argument("--vmin", type=float, default=None, help="Manual fixed colorbar minimum.")
    parser.add_argument("--vmax", type=float, default=None, help="Manual fixed colorbar maximum.")
    return parser.parse_args()


def main():
    args = parse_args()
    binary_root = args.binary_root.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if (args.vmin is None) != (args.vmax is None):
        raise ValueError("Provide both --vmin and --vmax together, or neither.")
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        raise ValueError(f"Expected --vmin < --vmax, got {args.vmin} >= {args.vmax}")

    levels = list(args.levels)
    cases = build_cases(binary_root, levels, box_margin=args.box_margin)
    nouts = common_nouts_across_cases(cases)

    if args.vmin is not None and args.vmax is not None:
        vmin, vmax = float(args.vmin), float(args.vmax)
        print(f"Using manual fixed limits: vmin={vmin:.6e}, vmax={vmax:.6e}")
    else:
        vmin, vmax = infer_global_limits(cases, nouts, z_center=args.z_center)
        print(f"Using global fixed limits: vmin={vmin:.6e}, vmax={vmax:.6e}")

    level_tag = "_".join(levels)
    frame_dir = outdir / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    for stale in frame_dir.glob("frame_*.png"):
        stale.unlink()

    print(f"Rendering {len(nouts)} synchronized frames for levels={levels}")
    frame_paths = []
    for i, nout in enumerate(nouts, start=1):
        out_png = frame_dir / f"frame_{i:04d}.png"
        render_six_panel_frame(
            cases=cases,
            nout=nout,
            npix=args.npix,
            z_center=args.z_center,
            vmin=vmin,
            vmax=vmax,
            out_png=out_png,
        )
        frame_paths.append(out_png)
        if i % 10 == 0 or i == len(nouts):
            print(f"  rendered {i}/{len(nouts)} frames")

    out_gif = outdir / f"binary_{level_tag}_mg_vs_fmm_grid_6panel.gif"
    stitch_gif(frame_paths, out_gif=out_gif, fps=args.fps)

    print("Done.")
    print(f"  GIF:    {out_gif}")
    print(f"  Frames: {frame_dir}")


if __name__ == "__main__":
    main()

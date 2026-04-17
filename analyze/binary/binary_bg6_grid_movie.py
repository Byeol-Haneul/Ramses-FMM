#!/usr/bin/env python3
"""
Build side-by-side MG vs FMM movies for BG6 binary runs with AMR grid overlays.

Outputs, for each level:
  - PNG frame sequence
  - GIF animation stitched from frames
"""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image


plt.rcParams["font.family"] = "serif"
plt.rcParams["figure.dpi"] = 150
plt.rcParams.update(
    {
        "font.size": 12,
        "legend.fontsize": 11,
    }
)


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "utils/py/miniramses.py").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root containing utils/py/miniramses.py")


REPO_ROOT = find_repo_root(Path.cwd())
UTILS_PY = REPO_ROOT / "utils/py"
if str(UTILS_PY) not in sys.path:
    sys.path.append(str(UTILS_PY))

import miniramses as ram  # noqa: E402

ORBIT_REF_CACHE: dict[str, dict[str, float | np.ndarray]] = {}


def quiet_call(fn, *args, **kwargs):
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        return fn(*args, **kwargs)


def list_outputs(run_dir: Path) -> list[int]:
    outs = []
    for p in run_dir.glob("output_*"):
        if not p.is_dir():
            continue
        try:
            outs.append(int(p.name.split("_")[-1]))
        except ValueError:
            continue
    outs.sort()
    if not outs:
        raise FileNotFoundError(f"No output_* directories in {run_dir}")
    return outs


def common_outputs(run_mg: Path, run_fmm: Path) -> list[int]:
    mg = set(list_outputs(run_mg))
    fmm = set(list_outputs(run_fmm))
    shared = sorted(mg & fmm)
    if not shared:
        raise RuntimeError(f"No common outputs between {run_mg} and {run_fmm}")
    return shared


def get_box_from_particles(run_dirs: Iterable[Path], outputs: Iterable[int], margin: float) -> list[float]:
    outputs = list(outputs)
    if not outputs:
        raise ValueError("outputs must contain at least one snapshot index.")

    mins = []
    maxs = []
    for run_dir in run_dirs:
        c0 = quiet_call(ram.rd_cell, outputs[0], path=str(run_dir), prefix="grav")
        xmin = float(np.min(c0.x[0] - c0.dx / 2))
        xmax = float(np.max(c0.x[0] + c0.dx / 2))
        ymin = float(np.min(c0.x[1] - c0.dx / 2))
        ymax = float(np.max(c0.x[1] + c0.dx / 2))

        for nout in outputs:
            p = quiet_call(ram.rd_part, nout, path=str(run_dir), silent=True)
            px = np.asarray(p.pos[0], dtype=np.float64)
            py = np.asarray(p.pos[1], dtype=np.float64)
            px, py = map_xy_to_plot(px, py, xmin, xmax, ymin, ymax)

            mins.append(np.array([np.min(px), np.min(py)]))
            maxs.append(np.array([np.max(px), np.max(py)]))
    mins = np.min(np.stack(mins), axis=0)
    maxs = np.max(np.stack(maxs), axis=0)
    return [
        float(mins[0] - margin),
        float(maxs[0] + margin),
        float(mins[1] - margin),
        float(maxs[1] + margin),
    ]


def amr_slice(c, z_center: float):
    mask = np.abs(c.x[2] - z_center) <= c.dx / 2
    return c.x[0][mask], c.x[1][mask], c.dx[mask], c.g[0][mask], c.level[mask]


def clip_to_box(x, y, dx, val, lev, box):
    xmin, xmax, ymin, ymax = box
    in_box = (x + dx / 2 > xmin) & (x - dx / 2 < xmax) & (y + dx / 2 > ymin) & (y - dx / 2 < ymax)
    return x[in_box], y[in_box], dx[in_box], val[in_box], lev[in_box]


def rasterize_amr(x, y, dx, val, box, npix: int):
    xmin, xmax, ymin, ymax = box
    sx = npix / (xmax - xmin)
    sy = npix / (ymax - ymin)
    img = np.full((npix, npix), np.nan, dtype=np.float32)

    order = np.argsort(-dx)  # coarse -> fine overwrite
    for i in order:
        x0 = x[i] - dx[i] / 2
        x1 = x[i] + dx[i] / 2
        y0 = y[i] - dx[i] / 2
        y1 = y[i] + dx[i] / 2
        if x1 <= xmin or x0 >= xmax or y1 <= ymin or y0 >= ymax:
            continue

        x0c = max(x0, xmin)
        x1c = min(x1, xmax)
        y0c = max(y0, ymin)
        y1c = min(y1, ymax)

        ix0 = max(int(math.floor((x0c - xmin) * sx)), 0)
        ix1 = min(int(math.ceil((x1c - xmin) * sx)), npix)
        iy0 = max(int(math.floor((y0c - ymin) * sy)), 0)
        iy1 = min(int(math.ceil((y1c - ymin) * sy)), npix)
        if ix1 > ix0 and iy1 > iy0:
            img[iy0:iy1, ix0:ix1] = val[i]
    return img


def build_grid_segments(x, y, dx, lev, box, npix: int):
    xmin, xmax, ymin, ymax = box
    pixel_size = max((xmax - xmin) / npix, (ymax - ymin) / npix)
    order = np.argsort(-dx)

    edges = set()
    segments = []
    widths = []

    for i in order:
        if dx[i] <= pixel_size * 1.5:
            continue

        x0 = x[i] - dx[i] / 2
        x1 = x[i] + dx[i] / 2
        y0 = y[i] - dx[i] / 2
        y1 = y[i] + dx[i] / 2
        if x1 <= xmin or x0 >= xmax or y1 <= ymin or y0 >= ymax:
            continue

        x0c = max(x0, xmin)
        x1c = min(x1, xmax)
        y0c = max(y0, ymin)
        y1c = min(y1, ymax)

        lw = 0.8 * 2 ** (-0.4 * (int(lev[i]) - 1))
        lw = max(lw, 0.05)

        cell_edges = [
            ((x0c, y0c), (x1c, y0c)),
            ((x1c, y0c), (x1c, y1c)),
            ((x1c, y1c), (x0c, y1c)),
            ((x0c, y1c), (x0c, y0c)),
        ]

        for edge in cell_edges:
            key = tuple(sorted(edge))
            if key not in edges:
                edges.add(key)
                segments.append(edge)
                widths.append(lw)

    return segments, widths


def ordered_particle_state(run_dir: Path, nout: int):
    p = quiet_call(ram.rd_part, nout, path=str(run_dir), silent=True)
    order = np.argsort(p.birth_id.astype(np.int64))
    pos = np.asarray(p.pos[:, order], dtype=np.float64)
    vel = np.asarray(p.vel[:, order], dtype=np.float64)
    mass = np.asarray(p.mass[order], dtype=np.float64)
    return pos, vel, mass


def ordered_particles(run_dir: Path, nout: int):
    pos, _, _ = ordered_particle_state(run_dir, nout)
    return pos[0], pos[1]


def map_xy_to_plot(
    x_raw: np.ndarray,
    y_raw: np.ndarray,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
):
    x = np.asarray(x_raw, dtype=np.float64)
    y = np.asarray(y_raw, dtype=np.float64)

    # Some runs store particles in normalized [0,1], while AMR cells are in shifted/scaled coordinates.
    if np.nanmax(x) <= 1.000001 and np.nanmin(x) >= -1e-6:
        x = xmin + x * (xmax - xmin)
    if np.nanmax(y) <= 1.000001 and np.nanmin(y) >= -1e-6:
        y = ymin + y * (ymax - ymin)
    return x, y


def get_orbit_reference(run_dir: Path):
    key = str(run_dir.resolve())
    cached = ORBIT_REF_CACHE.get(key)
    if cached is not None:
        return cached

    n0 = list_outputs(run_dir)[0]
    pos0, vel0, mass = ordered_particle_state(run_dir, n0)
    if mass.size != 2:
        raise ValueError(f"Expected exactly 2 particles in {run_dir}, found {mass.size}")

    m1 = float(mass[0])
    m2 = float(mass[1])
    mtot = m1 + m2
    com0 = np.sum(pos0[:2] * mass[None, :], axis=1) / mtot
    rel_pos = pos0[:2, 1] - pos0[:2, 0]
    rel_vel = vel0[:2, 1] - vel0[:2, 0]
    sep0 = float(np.linalg.norm(rel_pos))
    phi0 = float(np.arctan2(rel_pos[1], rel_pos[0]))

    cross = float(rel_pos[0] * rel_vel[1] - rel_pos[1] * rel_vel[0])
    omega = 0.0 if sep0 <= 0.0 else cross / (sep0 * sep0)
    # Fallback to circular estimate when initial relative velocity is tiny.
    if abs(omega) < 1e-14 and sep0 > 0.0:
        sign = 1.0 if cross >= 0.0 else -1.0
        omega = sign * float(np.sqrt(max(mtot, 0.0) / (sep0**3)))

    ref = {
        "m1": m1,
        "m2": m2,
        "com0": com0,
        "sep0": sep0,
        "phi0": phi0,
        "omega": omega,
    }
    ORBIT_REF_CACHE[key] = ref
    return ref


def analytic_positions_raw(run_dir: Path, time: float):
    ref = get_orbit_reference(run_dir)
    m1 = float(ref["m1"])
    m2 = float(ref["m2"])
    mtot = m1 + m2
    sep0 = float(ref["sep0"])
    phi = float(ref["phi0"]) + float(ref["omega"]) * float(time)
    com0 = np.asarray(ref["com0"], dtype=np.float64)

    xrel = sep0 * np.cos(phi)
    yrel = sep0 * np.sin(phi)
    x1 = com0[0] - (m2 / mtot) * xrel
    y1 = com0[1] - (m2 / mtot) * yrel
    x2 = com0[0] + (m1 / mtot) * xrel
    y2 = com0[1] + (m1 / mtot) * yrel
    return np.array([x1, x2], dtype=np.float64), np.array([y1, y2], dtype=np.float64)


def read_frame(run_dir: Path, nout: int, box, npix: int, z_center: float, with_grid: bool):
    c = quiet_call(ram.rd_cell, nout, path=str(run_dir), prefix="grav")
    info = quiet_call(ram.rd_info, nout, path=str(run_dir))

    xmin = float(np.min(c.x[0] - c.dx / 2))
    xmax = float(np.max(c.x[0] + c.dx / 2))
    ymin = float(np.min(c.x[1] - c.dx / 2))
    ymax = float(np.max(c.x[1] + c.dx / 2))

    x, y, dx, val, lev = amr_slice(c, z_center=z_center)
    x, y, dx, val, lev = clip_to_box(x, y, dx, val, lev, box=box)
    img = rasterize_amr(x, y, dx, val, box=box, npix=npix)

    segments, widths = [], []
    if with_grid:
        segments, widths = build_grid_segments(x, y, dx, lev, box=box, npix=npix)

    px, py = ordered_particles(run_dir, nout)
    px, py = map_xy_to_plot(px, py, xmin, xmax, ymin, ymax)

    ax_raw, ay_raw = analytic_positions_raw(run_dir, float(info.time))
    ax, ay = map_xy_to_plot(ax_raw, ay_raw, xmin, xmax, ymin, ymax)

    return {
        "image": img,
        "segments": segments,
        "widths": widths,
        "time": float(info.time),
        "particle_x": px,
        "particle_y": py,
        "analytic_x": ax,
        "analytic_y": ay,
    }


def infer_value_limits(run_dirs: list[Path], outputs: list[int], box, z_center: float, qlo=1.0, qhi=99.7):
    vals = []
    sample_ids = np.unique(np.rint(np.linspace(0, len(outputs) - 1, min(9, len(outputs)))).astype(int))
    sampled_outputs = [outputs[i] for i in sample_ids]
    for run_dir in run_dirs:
        for nout in sampled_outputs:
            c = quiet_call(ram.rd_cell, nout, path=str(run_dir), prefix="grav")
            x, y, dx, val, lev = amr_slice(c, z_center=z_center)
            _, _, _, v, _ = clip_to_box(x, y, dx, val, lev, box=box)
            if v.size:
                vals.append(v)
    if not vals:
        return -1.0, 1.0
    allv = np.concatenate(vals)
    vmin = float(np.percentile(allv, qlo))
    vmax = float(np.percentile(allv, qhi))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(allv))
        vmax = float(np.nanmax(allv))
    if vmin == vmax:
        vmax = vmin + 1e-12
    return vmin, vmax


def render_pair_frame(
    mg_frame,
    fmm_frame,
    level_label: str,
    nout: int,
    box,
    z_center: float,
    vmin: float,
    vmax: float,
    out_png: Path,
    cmap: str = "magma",
):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.8), constrained_layout=True)
    panel_data = [("MG", mg_frame, axes[0]), ("FMM", fmm_frame, axes[1])]

    for name, data, ax in panel_data:
        im = ax.imshow(
            data["image"],
            origin="lower",
            extent=[box[0], box[1], box[2], box[3]],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )

        if data["segments"]:
            lc = LineCollection(data["segments"], colors=(0, 0, 0, 0.72), linewidths=data["widths"])
            ax.add_collection(lc)

        ax.scatter(
            data["particle_x"],
            data["particle_y"],
            s=36,
            c=["#00D9FF", "#FF4D4D"],
            edgecolors="white",
            linewidths=0.6,
            zorder=6,
        )
        ax.scatter(
            data["analytic_x"],
            data["analytic_y"],
            s=190,
            facecolors="none",
            edgecolors=["#00D9FF", "#FF4D4D"],
            linewidths=2.0,
            zorder=7,
        )
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"{name} | output_{nout:05d} | t={data['time']:.6f}")

    dt = abs(mg_frame["time"] - fmm_frame["time"])
    fig.suptitle(
        f"Binary BG6 {level_label} | grav slice z={z_center:.3f} | circles=analytic | |Δt|={dt:.2e}",
        y=1.01,
        fontsize=12,
    )
    cbar = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02)
    cbar.set_label(r"$\Phi$")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)


def stitch_gif(frame_paths: list[Path], out_gif: Path, fps: int):
    images = [Image.open(p).convert("P", palette=Image.Palette.ADAPTIVE) for p in frame_paths]
    duration_ms = int(1000 / max(fps, 1))
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(
        out_gif,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
        disposal=2,
    )
    for img in images:
        img.close()


def make_level_movie(
    level: str,
    run_mg: Path,
    run_fmm: Path,
    output_root: Path,
    npix: int,
    fps: int,
    z_center: float,
    box_margin: float,
    fixed_limits: tuple[float, float] | None = None,
):
    outputs = common_outputs(run_mg, run_fmm)
    box = get_box_from_particles([run_mg, run_fmm], outputs=[outputs[0], outputs[-1]], margin=box_margin)
    if fixed_limits is None:
        vmin, vmax = infer_value_limits([run_mg, run_fmm], outputs, box=box, z_center=z_center)
    else:
        vmin, vmax = fixed_limits

    level_out = output_root / level
    frame_dir = level_out / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{level}] outputs={len(outputs)} box={box} vmin={vmin:.6e} vmax={vmax:.6e}")
    frame_paths = []
    for i, nout in enumerate(outputs, start=1):
        mg_frame = read_frame(run_mg, nout, box=box, npix=npix, z_center=z_center, with_grid=True)
        fmm_frame = read_frame(run_fmm, nout, box=box, npix=npix, z_center=z_center, with_grid=True)
        out_png = frame_dir / f"frame_{i:04d}.png"
        render_pair_frame(
            mg_frame=mg_frame,
            fmm_frame=fmm_frame,
            level_label=level.upper(),
            nout=nout,
            box=box,
            z_center=z_center,
            vmin=vmin,
            vmax=vmax,
            out_png=out_png,
        )
        frame_paths.append(out_png)
        if i % 10 == 0 or i == len(outputs):
            print(f"  rendered {i}/{len(outputs)} frames")

    out_gif = level_out / f"binary_{level}_mg_vs_fmm_grid.gif"
    stitch_gif(frame_paths, out_gif=out_gif, fps=fps)
    print(f"[{level}] wrote GIF: {out_gif}")
    return out_gif, frame_dir


def parse_args():
    parser = argparse.ArgumentParser(description="Build BG6 MG-vs-FMM AMR-grid comparison movies.")
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["bg6_l6", "bg6_l7", "bg6_l8"],
        help="Run suffix levels to render.",
    )
    parser.add_argument("--npix", type=int, default=720, help="Per-panel raster resolution.")
    parser.add_argument("--fps", type=int, default=12, help="GIF frame rate.")
    parser.add_argument("--z-center", type=float, default=1.0, help="Z slice center.")
    parser.add_argument("--box-margin", type=float, default=0.08, help="Margin around particle orbit for box.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "analyze" / "binary" / "results" / "movies_bg6",
        help="Output root directory for movies and frames.",
    )
    parser.add_argument(
        "--binary-root",
        type=Path,
        default=REPO_ROOT / "binary_test",
        help="Input root containing mg/ and fmm/ run folders.",
    )
    parser.add_argument(
        "--global-limits",
        action="store_true",
        help="Use one shared colorbar range across all rendered levels.",
    )
    parser.add_argument("--vmin", type=float, default=None, help="Manual fixed colorbar minimum.")
    parser.add_argument("--vmax", type=float, default=None, help="Manual fixed colorbar maximum.")
    return parser.parse_args()


def main():
    args = parse_args()
    binary_root = args.binary_root.resolve()
    output_root = args.output_dir.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"Repo root: {REPO_ROOT}")
    print(f"Binary root: {binary_root}")
    print(f"Output root: {output_root}")

    if (args.vmin is None) != (args.vmax is None):
        raise ValueError("Provide both --vmin and --vmax together, or neither.")
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        raise ValueError(f"Expected --vmin < --vmax, got {args.vmin} >= {args.vmax}")

    cases = []
    for level in args.levels:
        run_mg = binary_root / "mg" / f"equal_mass_binary_amr_{level}"
        run_fmm = binary_root / "fmm" / f"equal_mass_binary_amr_{level}"
        if not run_mg.exists() or not run_fmm.exists():
            raise FileNotFoundError(f"Missing run directories for level {level}: {run_mg}, {run_fmm}")
        outputs = common_outputs(run_mg, run_fmm)
        box = get_box_from_particles([run_mg, run_fmm], outputs=[outputs[0], outputs[-1]], margin=args.box_margin)
        cases.append(
            {
                "level": level,
                "run_mg": run_mg,
                "run_fmm": run_fmm,
                "outputs": outputs,
                "box": box,
            }
        )

    fixed_limits = None
    if args.vmin is not None and args.vmax is not None:
        fixed_limits = (args.vmin, args.vmax)
        print(f"Using manual fixed limits: vmin={fixed_limits[0]:.6e}, vmax={fixed_limits[1]:.6e}")
    elif args.global_limits:
        mins = []
        maxs = []
        for case in cases:
            vmin, vmax = infer_value_limits(
                [case["run_mg"], case["run_fmm"]],
                case["outputs"],
                box=case["box"],
                z_center=args.z_center,
            )
            mins.append(vmin)
            maxs.append(vmax)
        fixed_limits = (float(min(mins)), float(max(maxs)))
        print(f"Using global fixed limits: vmin={fixed_limits[0]:.6e}, vmax={fixed_limits[1]:.6e}")

    outputs = []
    for case in cases:
        level = case["level"]
        gif_path, frame_dir = make_level_movie(
            level=level,
            run_mg=case["run_mg"],
            run_fmm=case["run_fmm"],
            output_root=output_root,
            npix=args.npix,
            fps=args.fps,
            z_center=args.z_center,
            box_margin=args.box_margin,
            fixed_limits=fixed_limits,
        )
        outputs.append((level, gif_path, frame_dir))

    print("\nDone.")
    for level, gif_path, frame_dir in outputs:
        print(f"  {level}:")
        print(f"    GIF:    {gif_path}")
        print(f"    Frames: {frame_dir}")


if __name__ == "__main__":
    main()

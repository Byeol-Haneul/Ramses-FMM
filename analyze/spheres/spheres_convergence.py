#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
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
    raise FileNotFoundError("Could not locate repository root from current path.")


try:
    HERE = Path(__file__).resolve().parent
except NameError:
    HERE = Path.cwd()

REPO_ROOT = find_repo_root(HERE)
if str(REPO_ROOT / "utils/py") not in sys.path:
    sys.path.append(str(REPO_ROOT / "utils/py"))

import miniramses as ram


G = 1.0
RHO = 1.0
SPHERE1_CENTER = np.array([1.7, 2.0, 2.0], dtype=np.float64)
SPHERE2_CENTER = np.array([2.2, 2.0, 2.0], dtype=np.float64)
SPHERE1_RADIUS = 0.1
SPHERE2_RADIUS = 0.2


@dataclass
class ConvergenceRow:
    case: str
    level: int
    ncell: int
    sim_time: float
    h: float
    mean_l2_abs: float
    mean_l2_rel: float
    rho_phi_dv_sum: float
    rho_phi_dv_abs_sum: float
    rho_phi_dv_rms: float


def _uniform_sphere_mass(radius: float, density: float = RHO) -> float:
    return (4.0 / 3.0) * np.pi * radius**3 * density


def _sphere_offsets(x: np.ndarray, y: np.ndarray, z: np.ndarray, center: np.ndarray):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    r = np.maximum(r, 1.0e-14)
    return dx, dy, dz, r


def phi_uniform_sphere(r: np.ndarray, radius: float, mass: float) -> np.ndarray:
    phi_out = -G * mass / r
    phi_in = -1.5 * G * mass / radius + 0.5 * G * mass * r**2 / radius**3
    return np.where(r >= radius, phi_out, phi_in)


def analytic_double_sphere_phi(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    m1 = _uniform_sphere_mass(SPHERE1_RADIUS)
    m2 = _uniform_sphere_mass(SPHERE2_RADIUS)

    _, _, _, r1 = _sphere_offsets(x, y, z, SPHERE1_CENTER)
    _, _, _, r2 = _sphere_offsets(x, y, z, SPHERE2_CENTER)

    return phi_uniform_sphere(r1, SPHERE1_RADIUS, m1) + phi_uniform_sphere(r2, SPHERE2_RADIUS, m2)


def analytic_double_sphere_rho(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    _, _, _, r1 = _sphere_offsets(x, y, z, SPHERE1_CENTER)
    _, _, _, r2 = _sphere_offsets(x, y, z, SPHERE2_CENTER)
    rho1 = np.where(r1 <= SPHERE1_RADIUS, RHO, 0.0)
    rho2 = np.where(r2 <= SPHERE2_RADIUS, RHO, 0.0)
    return rho1 + rho2


def load_grav(run_dir: Path, nout: int, prefix: str):
    c = ram.rd_cell(str(nout), path=str(run_dir), prefix=prefix)
    info = ram.rd_info(str(nout), path=str(run_dir))
    return c, info


def compute_row(case: str, level: int, run_dir: Path, nout: int, prefix: str) -> ConvergenceRow:
    c, info = load_grav(run_dir, nout=nout, prefix=prefix)

    x = np.asarray(c.x[0], dtype=np.float64)
    y = np.asarray(c.x[1], dtype=np.float64)
    z = np.asarray(c.x[2], dtype=np.float64)
    phi_sim = np.asarray(c.g[0], dtype=np.float64)
    phi_ana = analytic_double_sphere_phi(x, y, z)
    rho_ana = analytic_double_sphere_rho(x, y, z)
    dx = np.asarray(c.dx, dtype=np.float64)
    dvol = dx ** int(c.ndim)

    err = phi_sim - phi_ana
    mean_l2_abs = float(np.sqrt(np.mean(err**2)))

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(np.abs(phi_ana) > 0.0, err / phi_ana, np.nan)
    finite_rel = rel[np.isfinite(rel)]
    mean_l2_rel = float(np.sqrt(np.mean(finite_rel**2))) if finite_rel.size else float("nan")

    weighted = rho_ana * err * dvol
    rho_phi_dv_sum = float(np.sum(weighted))
    rho_phi_dv_abs_sum = float(np.sum(np.abs(weighted)))
    rho_phi_dv_rms = float(np.sqrt(np.mean(weighted**2)))

    h = float(np.min(dx))

    return ConvergenceRow(
        case=case,
        level=level,
        ncell=int(c.ncell),
        sim_time=float(info.time),
        h=h,
        mean_l2_abs=mean_l2_abs,
        mean_l2_rel=mean_l2_rel,
        rho_phi_dv_sum=rho_phi_dv_sum,
        rho_phi_dv_abs_sum=rho_phi_dv_abs_sum,
        rho_phi_dv_rms=rho_phi_dv_rms,
    )


def parse_levels(levels_arg: list[int] | None, level_min: int, level_max: int) -> list[int]:
    if levels_arg:
        levels = sorted(set(levels_arg))
    else:
        levels = list(range(level_min, level_max + 1))
    if not levels:
        raise ValueError("No levels requested.")
    return levels


def parse_cases(case_args: list[str], runs_root: Path | None, label: str) -> list[tuple[str, Path]]:
    cases: list[tuple[str, Path]] = []

    for item in case_args:
        if "=" not in item:
            raise ValueError(f"Invalid --case '{item}'. Use NAME=PATH.")
        name, path = item.split("=", 1)
        name = name.strip()
        path = path.strip()
        if not name or not path:
            raise ValueError(f"Invalid --case '{item}'. Use NAME=PATH.")
        cases.append((name, Path(path).resolve()))

    if not cases and runs_root is not None:
        cases.append((label, runs_root.resolve()))

    if not cases:
        raise ValueError("No run cases provided. Use --case NAME=PATH or --runs-root.")

    return cases


def write_csv(rows: list[ConvergenceRow], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "case",
                "level",
                "ncell",
                "time",
                "h",
                "mean_l2_abs",
                "mean_l2_rel",
                "rho_phi_dv_sum",
                "rho_phi_dv_abs_sum",
                "rho_phi_dv_rms",
                "log10_h",
                "log10_rho_phi_dv_abs_sum",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.case,
                    row.level,
                    row.ncell,
                    f"{row.sim_time:.16e}",
                    f"{row.h:.16e}",
                    f"{row.mean_l2_abs:.16e}",
                    f"{row.mean_l2_rel:.16e}",
                    f"{row.rho_phi_dv_sum:.16e}",
                    f"{row.rho_phi_dv_abs_sum:.16e}",
                    f"{row.rho_phi_dv_rms:.16e}",
                    f"{math.log10(row.h):.16e}",
                    f"{math.log10(row.rho_phi_dv_abs_sum):.16e}",
                ]
            )


def plot_convergence(rows: list[ConvergenceRow], plot_path: Path) -> dict[str, float | None]:
    cases = sorted({row.case for row in rows})
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    slopes: dict[str, float | None] = {}

    for idx, case in enumerate(cases):
        subset = sorted((r for r in rows if r.case == case), key=lambda r: r.level)
        logh = np.log10(np.array([r.h for r in subset], dtype=np.float64))
        loge = np.log10(np.array([r.rho_phi_dv_abs_sum for r in subset], dtype=np.float64))
        color = cmap(idx % 10)

        slope = None
        if len(subset) >= 2:
            slope, intercept = np.polyfit(logh, loge, 1)
            fit_x = np.linspace(np.min(logh), np.max(logh), 100)
            fit_y = slope * fit_x + intercept
            ax.plot(fit_x, fit_y, "--", color=color, lw=1.2, alpha=0.9)

        label = case if slope is None else f"{case} (slope={slope:.3f})"
        ax.plot(logh, loge, "o-", color=color, lw=1.6, label=label)

        for row, xval, yval in zip(subset, logh, loge):
            ax.annotate(f"L{row.level}", (xval, yval), textcoords="offset points", xytext=(4, 4), fontsize=8, color=color)

        slopes[case] = None if slope is None else float(slope)

    ax.set_xlabel("log10(h)")
    ax.set_ylabel("log10(sum |rho_ana*(phi-phi_ana)*dV|)")
    ax.set_title("Uniform Spheres Convergence (rho-weighted): MG vs FMM")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    return slopes


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute convergence using rho_ana*(phi-phi_ana)*dV for one or more uniform-spheres runs."
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="Case in NAME=PATH form. Can be passed multiple times.",
    )
    parser.add_argument(
        "--runs-root",
        type=Path,
        default=None,
        help="Single-case root directory containing level_XX runs (legacy mode).",
    )
    parser.add_argument(
        "--label",
        type=str,
        default="FMM",
        help="Label used with --runs-root single-case mode.",
    )
    parser.add_argument(
        "--levels",
        type=int,
        nargs="+",
        default=None,
        help="Explicit levels to analyze. If omitted, uses --level-min..--level-max.",
    )
    parser.add_argument("--level-min", type=int, default=6, help="Min level for auto-range mode.")
    parser.add_argument("--level-max", type=int, default=10, help="Max level for auto-range mode.")
    parser.add_argument("--nout", type=int, default=2, help="Output number to analyze (default: 2).")
    parser.add_argument("--prefix", type=str, default="grav", help="RAMSES field prefix (default: grav).")
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("analyze/spheres/spheres_convergence.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("analyze/spheres/spheres_convergence.png"),
        help="Output plot path.",
    )
    args = parser.parse_args()

    levels = parse_levels(args.levels, args.level_min, args.level_max)
    cases = parse_cases(args.case, args.runs_root, args.label)

    rows: list[ConvergenceRow] = []
    for case_name, case_root in cases:
        for level in levels:
            run_dir = case_root / f"level_{level:02d}"
            out_dir = run_dir / f"output_{args.nout:05d}"
            if not out_dir.exists():
                raise FileNotFoundError(
                    f"Missing {out_dir}. Ensure run '{case_name}' completed and --nout is correct."
                )
            rows.append(compute_row(case_name, level, run_dir, nout=args.nout, prefix=args.prefix))

    rows.sort(key=lambda r: (r.case, r.level))
    if any((row.h <= 0.0 or row.rho_phi_dv_abs_sum <= 0.0) for row in rows):
        raise ValueError("Found non-positive h or rho-weighted absolute-sum error; cannot compute log10 values.")

    csv_path = (REPO_ROOT / args.csv).resolve() if not args.csv.is_absolute() else args.csv
    plot_path = (REPO_ROOT / args.plot).resolve() if not args.plot.is_absolute() else args.plot

    write_csv(rows, csv_path)
    slopes = plot_convergence(rows, plot_path)

    print("Convergence rows:")
    for row in rows:
        print(
            f"  case={row.case:>4s} level={row.level:2d} "
            f"h={row.h:.6e} mean_l2_abs={row.mean_l2_abs:.6e} "
            f"mean_l2_rel={row.mean_l2_rel:.6e} "
            f"rho_phi_dv_sum={row.rho_phi_dv_sum:.6e} "
            f"rho_phi_dv_abs_sum={row.rho_phi_dv_abs_sum:.6e} "
            f"log10(h)={math.log10(row.h):.6f} "
            f"log10(err)={math.log10(row.rho_phi_dv_abs_sum):.6f}"
        )

    for case, slope in slopes.items():
        if slope is None:
            print(f"Case {case}: slope unavailable (need >= 2 levels)")
        else:
            print(f"Case {case}: best-fit slope (log10(err) vs log10(h)) = {slope:.6f}")

    print(f"Wrote CSV : {csv_path}")
    print(f"Wrote plot: {plot_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

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

G = 1.0

COLORS = {"MG": "#4C5B7A", "FMM": "#C06C2B"}


def list_outputs(run_dir: Path) -> list[int]:
    return sorted(int(path.name.split("_")[-1]) for path in run_dir.glob("output_*") if path.is_dir())


def ordered_particle_state(part) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(part.birth_id.astype(np.int64))
    pos = np.asarray(part.pos[:, order], dtype=np.float64)
    vel = np.asarray(part.vel[:, order], dtype=np.float64)
    mass = np.asarray(part.mass[order], dtype=np.float64)
    return pos, vel, mass


def load_solver_history(run_root: Path, solver: str) -> pd.DataFrame:
    solver_dir = run_root / solver
    rows = []
    for nout in list_outputs(solver_dir):
        info = ram.rd_info(nout, path=str(solver_dir))
        part = ram.rd_part(nout, path=str(solver_dir), silent=True)
        pos, vel, mass = ordered_particle_state(part)
        rel_pos = pos[:, 1] - pos[:, 0]
        rel_vel = vel[:, 1] - vel[:, 0]
        com = np.sum(pos * mass[None, :], axis=1) / np.sum(mass)
        kinetic = float(0.5 * np.sum(mass * np.sum(vel**2, axis=0)))
        potential = float(-G * mass[0] * mass[1] / np.linalg.norm(rel_pos))
        angmom = np.sum(np.cross(pos.T, (mass[None, :] * vel).T), axis=0)
        rows.append(
            {
                "time": float(info.time),
                "x1": float(pos[0, 0]),
                "y1": float(pos[1, 0]),
                "z1": float(pos[2, 0]),
                "x2": float(pos[0, 1]),
                "y2": float(pos[1, 1]),
                "z2": float(pos[2, 1]),
                "vx1": float(vel[0, 0]),
                "vy1": float(vel[1, 0]),
                "vz1": float(vel[2, 0]),
                "vx2": float(vel[0, 1]),
                "vy2": float(vel[1, 1]),
                "vz2": float(vel[2, 1]),
                "sep": float(np.linalg.norm(rel_pos)),
                "phase": float(np.arctan2(rel_pos[1], rel_pos[0])),
                "com_x": float(com[0]),
                "com_y": float(com[1]),
                "com_z": float(com[2]),
                "m1": float(mass[0]),
                "m2": float(mass[1]),
                "kinetic": kinetic,
                "potential": potential,
                "energy": kinetic + potential,
                "lx": float(angmom[0]),
                "ly": float(angmom[1]),
                "lz": float(angmom[2]),
            }
        )
    df = pd.DataFrame(rows).sort_values("time").reset_index(drop=True)
    df["solver"] = solver.upper()
    df["lmag"] = np.sqrt(df["lx"] ** 2 + df["ly"] ** 2 + df["lz"] ** 2)
    return df


def analytic_reference(df: pd.DataFrame) -> dict[str, np.ndarray | float]:
    t = df["time"].to_numpy()
    m1 = float(df["m1"].iloc[0])
    m2 = float(df["m2"].iloc[0])
    mtot = m1 + m2
    r0 = np.array([df["x2"].iloc[0] - df["x1"].iloc[0], df["y2"].iloc[0] - df["y1"].iloc[0]])
    v0 = np.array([df["vx2"].iloc[0] - df["vx1"].iloc[0], df["vy2"].iloc[0] - df["vy1"].iloc[0]])
    sep0 = float(np.linalg.norm(r0))
    omega = float(np.sqrt(G * mtot / sep0**3))
    sign = float(np.sign(r0[0] * v0[1] - r0[1] * v0[0]))
    phi0 = float(np.arctan2(r0[1], r0[0]))
    phi = phi0 + sign * omega * t
    xrel = sep0 * np.cos(phi)
    yrel = sep0 * np.sin(phi)
    com_x0 = float(df["com_x"].iloc[0])
    com_y0 = float(df["com_y"].iloc[0])
    phase_ref = np.unwrap(phi)
    return {
        "sep0": sep0,
        "period": 2.0 * np.pi / omega,
        "phase_ref": phase_ref,
        "x1": com_x0 - 0.5 * xrel,
        "y1": com_y0 - 0.5 * yrel,
        "x2": com_x0 + 0.5 * xrel,
        "y2": com_y0 + 0.5 * yrel,
    }


def add_error_columns(df: pd.DataFrame, ref: dict[str, np.ndarray | float]) -> None:
    phase = np.unwrap(df["phase"].to_numpy())
    df["phase_err"] = phase - ref["phase_ref"]
    df["sep_rel_err"] = (df["sep"] - ref["sep0"]) / ref["sep0"]
    body1_err = np.sqrt((df["x1"].to_numpy() - ref["x1"]) ** 2 + (df["y1"].to_numpy() - ref["y1"]) ** 2)
    body2_err = np.sqrt((df["x2"].to_numpy() - ref["x2"]) ** 2 + (df["y2"].to_numpy() - ref["y2"]) ** 2)
    df["body_err"] = 0.5 * (body1_err + body2_err)

    com0 = df.loc[0, ["com_x", "com_y", "com_z"]].to_numpy(dtype=float)
    df["com_r"] = np.sqrt(
        (df["com_x"] - com0[0]) ** 2
        + (df["com_y"] - com0[1]) ** 2
        + (df["com_z"] - com0[2]) ** 2
    )

    e0 = float(df["energy"].iloc[0])
    l0 = float(df["lmag"].iloc[0])
    df["energy_rel_drift"] = (df["energy"] - e0) / abs(e0)
    df["lmag_rel_drift"] = (df["lmag"] - l0) / abs(l0)


def build_summary(*histories: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for df in histories:
        rows.append(
            {
                "solver": df["solver"].iloc[0],
                "rms_body_err": float(np.sqrt(np.mean(df["body_err"] ** 2))),
                "max_body_err": float(np.max(df["body_err"])),
                "rms_sep_rel_err": float(np.sqrt(np.mean(df["sep_rel_err"] ** 2))),
                "max_sep_rel_err": float(np.max(np.abs(df["sep_rel_err"]))),
                "rms_phase_err": float(np.sqrt(np.mean(df["phase_err"] ** 2))),
                "max_phase_err": float(np.max(np.abs(df["phase_err"]))),
                "final_phase_err": float(df["phase_err"].iloc[-1]),
                "max_com_drift": float(np.max(df["com_r"])),
                "max_abs_energy_rel_drift": float(np.max(np.abs(df["energy_rel_drift"]))),
                "max_abs_lmag_rel_drift": float(np.max(np.abs(df["lmag_rel_drift"]))),
            }
        )
    return pd.DataFrame(rows)


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "font.family": "serif",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titlepad": 10.0,
        }
    )


def style_axis(ax) -> None:
    ax.grid(alpha=0.22, linewidth=0.8)
    ax.tick_params(direction="out", length=4)


def positive_for_log(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    positive = arr[arr > 0.0]
    floor = 1e-16 if positive.size == 0 else max(np.min(positive) * 0.5, 1e-16)
    return np.maximum(arr, floor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    run_root = args.run_root.resolve()
    outdir = (args.output_dir or (run_root / "analysis")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mg = load_solver_history(run_root, "mg")
    fmm = load_solver_history(run_root, "fmm")
    ref = analytic_reference(mg)
    add_error_columns(mg, ref)
    add_error_columns(fmm, ref)
    summary = build_summary(mg, fmm)
    summary.to_csv(outdir / "binary_vs_analytic_summary.csv", index=False)

    apply_plot_style()
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(ref["x1"], ref["y1"], color="black", lw=2, label="Analytic body 1")
    ax.plot(ref["x2"], ref["y2"], color="black", lw=2, alpha=0.45, label="Analytic body 2")
    ax.plot(mg["x1"], mg["y1"], color=COLORS["MG"], label="MG body 1")
    ax.plot(fmm["x1"], fmm["y1"], color=COLORS["FMM"], label="FMM body 1")
    ax.set_title("Orbital Track vs Analytic")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axis("equal")
    style_axis(ax)
    ax.legend(loc="best", fontsize=9)

    ax = axes[0, 1]
    ax.plot(mg["time"], mg["sep_rel_err"], color=COLORS["MG"], lw=2, label="MG")
    ax.plot(fmm["time"], fmm["sep_rel_err"], color=COLORS["FMM"], lw=2, label="FMM")
    ax.axhline(0.0, color="black", lw=1, alpha=0.45)
    ax.set_title("Relative Separation Error")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$(r-r_0)/r_0$")
    style_axis(ax)
    ax.legend(loc="best")

    ax = axes[0, 2]
    ax.plot(mg["time"], mg["phase_err"], color=COLORS["MG"], lw=2, label="MG")
    ax.plot(fmm["time"], fmm["phase_err"], color=COLORS["FMM"], lw=2, label="FMM")
    ax.axhline(0.0, color="black", lw=1, alpha=0.45)
    ax.set_title("Phase Error vs Analytic")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\Delta \phi$ [rad]")
    style_axis(ax)
    ax.legend(loc="best")

    ax = axes[1, 0]
    ax.plot(mg["time"], mg["body_err"], color=COLORS["MG"], lw=2, label="MG")
    ax.plot(fmm["time"], fmm["body_err"], color=COLORS["FMM"], lw=2, label="FMM")
    ax.set_title("Mean Body Position Error")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\frac{|{\bf x}_1-{\bf x}_{1,\rm ana}| + |{\bf x}_2-{\bf x}_{2,\rm ana}|}{2}$")
    style_axis(ax)
    ax.legend(loc="best")

    ax = axes[1, 1]
    ax.plot(mg["time"], positive_for_log(mg["com_r"]), color=COLORS["MG"], lw=2, label="MG")
    ax.plot(fmm["time"], positive_for_log(fmm["com_r"]), color=COLORS["FMM"], lw=2, label="FMM")
    ax.set_yscale("log")
    ax.set_title("Center-of-Mass Drift")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$|{\bf x}_{\rm COM}-{\bf x}_{{\rm COM},0}|$ (log)")
    style_axis(ax)
    ax.legend(loc="best")

    ax = axes[1, 2]
    ax.plot(mg["time"], mg["energy_rel_drift"], color=COLORS["MG"], lw=2, label=r"MG $\Delta E/|E_0|$")
    ax.plot(fmm["time"], fmm["energy_rel_drift"], color=COLORS["FMM"], lw=2, label=r"FMM $\Delta E/|E_0|$")
    ax.plot(
        mg["time"],
        mg["lmag_rel_drift"],
        color=COLORS["MG"],
        lw=2,
        alpha=0.55,
        linestyle=":",
        label=r"MG $\Delta L/|L_0|$",
    )
    ax.plot(
        fmm["time"],
        fmm["lmag_rel_drift"],
        color=COLORS["FMM"],
        lw=2,
        alpha=0.55,
        linestyle=":",
        label=r"FMM $\Delta L/|L_0|$",
    )
    ax.axhline(0.0, color="black", lw=1, alpha=0.45)
    ax.set_title("Energy and Angular Momentum Drift")
    ax.set_xlabel("Time")
    ax.set_ylabel("relative drift")
    style_axis(ax)
    ax.legend(loc="best", fontsize=9)

    fig.suptitle(f"Equal-Mass Binary vs Analytic: {run_root.name}", fontsize=14)
    fig.savefig(outdir / "binary_vs_analytic.png", bbox_inches="tight")

    print(summary.to_string(index=False))
    print(outdir / "binary_vs_analytic.png")
    print(outdir / "binary_vs_analytic_summary.csv")


if __name__ == "__main__":
    main()

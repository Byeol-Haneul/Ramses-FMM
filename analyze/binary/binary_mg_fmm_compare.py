# %% [markdown]
# Equal-mass binary MG/FMM comparison.
#
# Open this file in Jupyter as a percent-format notebook, or run it as a
# regular Python script. It reads the latest run from `amr_test/runs/latest`
# by default and compares the orbital dynamics measured from the particle
# snapshots.

# %%
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

try:
    from IPython.display import display
except ImportError:
    def display(obj):
        print(obj)

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

plt.rcParams["figure.dpi"] = 160
plt.rcParams["font.family"] = "serif"
plt.rcParams.update({"font.size": 12, "legend.fontsize": 10})

G = 1.0
RUN_ROOT = REPO_ROOT / "amr_test" / "runs" / "latest"
SOLVERS = ("mg", "fmm")
SAVE_DIR = RUN_ROOT / "analysis"
SAVE_DIR.mkdir(parents=True, exist_ok=True)


# %%
def list_outputs(run_dir: Path) -> list[int]:
    outputs = []
    for path in sorted(run_dir.glob("output_*")):
        if path.is_dir():
            outputs.append(int(path.name.split("_")[-1]))
    if not outputs:
        raise FileNotFoundError(f"No output_* directories found in {run_dir}")
    return outputs


def ordered_particle_state(part) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(part.birth_id.astype(np.int64))
    pos = np.asarray(part.pos[:, order], dtype=np.float64)
    vel = np.asarray(part.vel[:, order], dtype=np.float64)
    mass = np.asarray(part.mass[order], dtype=np.float64)
    ids = np.asarray(part.birth_id[order], dtype=np.int64)
    return pos, vel, mass, ids


def snapshot_record(run_dir: Path, nout: int) -> dict[str, float]:
    info = ram.rd_info(nout, path=str(run_dir))
    part = ram.rd_part(nout, path=str(run_dir), silent=True)
    if part.npart != 2:
        raise ValueError(f"Expected exactly 2 particles in {run_dir}/output_{nout:05d}, found {part.npart}")

    pos, vel, mass, ids = ordered_particle_state(part)
    m_tot = float(np.sum(mass))
    com_pos = np.sum(pos * mass[None, :], axis=1) / m_tot
    com_vel = np.sum(vel * mass[None, :], axis=1) / m_tot

    rel_pos = pos[:, 1] - pos[:, 0]
    rel_vel = vel[:, 1] - vel[:, 0]
    sep = float(np.linalg.norm(rel_pos))

    kinetic = float(0.5 * np.sum(mass * np.sum(vel**2, axis=0)))
    potential = float(-G * mass[0] * mass[1] / sep)
    angmom = np.sum(np.cross(pos.T, (mass[None, :] * vel).T), axis=0)

    return {
        "nout": int(nout),
        "time": float(info.time),
        "boxlen": float(info.boxlen),
        "id1": int(ids[0]),
        "id2": int(ids[1]),
        "m1": float(mass[0]),
        "m2": float(mass[1]),
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
        "sep": sep,
        "phase": float(np.arctan2(rel_pos[1], rel_pos[0])),
        "com_x": float(com_pos[0]),
        "com_y": float(com_pos[1]),
        "com_z": float(com_pos[2]),
        "com_vx": float(com_vel[0]),
        "com_vy": float(com_vel[1]),
        "com_vz": float(com_vel[2]),
        "kinetic": kinetic,
        "potential": potential,
        "energy": kinetic + potential,
        "lx": float(angmom[0]),
        "ly": float(angmom[1]),
        "lz": float(angmom[2]),
    }


def load_solver_history(run_dir: Path, solver: str) -> pd.DataFrame:
    solver_dir = run_dir / solver
    records = [snapshot_record(solver_dir, nout) for nout in list_outputs(solver_dir)]
    df = pd.DataFrame.from_records(records).sort_values("time").reset_index(drop=True)
    df["solver"] = solver.upper()
    df["phase_unwrapped"] = np.unwrap(df["phase"].to_numpy())
    df["com_v"] = np.sqrt(df["com_vx"] ** 2 + df["com_vy"] ** 2 + df["com_vz"] ** 2)
    df["lmag"] = np.sqrt(df["lx"] ** 2 + df["ly"] ** 2 + df["lz"] ** 2)

    com0 = df.loc[0, ["com_x", "com_y", "com_z"]].to_numpy(dtype=float)
    df["com_dx"] = df["com_x"] - com0[0]
    df["com_dy"] = df["com_y"] - com0[1]
    df["com_dz"] = df["com_z"] - com0[2]
    df["com_r"] = np.sqrt(df["com_dx"] ** 2 + df["com_dy"] ** 2 + df["com_dz"] ** 2)
    df["x1_relcom0"] = df["x1"] - com0[0]
    df["y1_relcom0"] = df["y1"] - com0[1]
    df["x2_relcom0"] = df["x2"] - com0[0]
    df["y2_relcom0"] = df["y2"] - com0[1]

    e0 = float(df["energy"].iloc[0])
    l0 = float(df["lmag"].iloc[0])
    sep0 = float(df["sep"].iloc[0])
    mtot = float(df["m1"].iloc[0] + df["m2"].iloc[0])
    period = 2.0 * np.pi * np.sqrt(sep0**3 / (G * mtot))

    df["energy_rel_drift"] = (df["energy"] - e0) / abs(e0)
    df["lmag_rel_drift"] = (df["lmag"] - l0) / abs(l0)
    df["sep_rel_dev"] = (df["sep"] - sep0) / sep0
    df["expected_phase"] = df["time"] * (2.0 * np.pi / period)
    df["phase_residual"] = df["phase_unwrapped"] - df["expected_phase"]
    df.attrs["period"] = period
    df.attrs["sep0"] = sep0
    df.attrs["mtot"] = mtot
    return df


def build_summary(*histories: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for df in histories:
        rows.append(
            {
                "solver": df["solver"].iloc[0],
                "n_outputs": int(len(df)),
                "time_max": float(df["time"].max()),
                "orbits_completed": float(df["phase_unwrapped"].iloc[-1] / (2.0 * np.pi)),
                "sep0": float(df["sep"].iloc[0]),
                "sep_mean": float(df["sep"].mean()),
                "max_abs_sep_rel_dev": float(np.max(np.abs(df["sep_rel_dev"]))),
                "max_com_drift": float(np.max(df["com_r"])),
                "final_abs_energy_rel_drift": float(abs(df["energy_rel_drift"].iloc[-1])),
                "max_abs_energy_rel_drift": float(np.max(np.abs(df["energy_rel_drift"]))),
                "final_abs_lmag_rel_drift": float(abs(df["lmag_rel_drift"].iloc[-1])),
                "max_abs_lmag_rel_drift": float(np.max(np.abs(df["lmag_rel_drift"]))),
                "max_abs_phase_residual": float(np.max(np.abs(df["phase_residual"]))),
                "period_from_ic": float(df.attrs["period"]),
            }
        )
    return pd.DataFrame(rows)


def compare_histories(mg: pd.DataFrame, fmm: pd.DataFrame) -> pd.DataFrame:
    mg_small = mg[["time", "sep", "com_r", "energy_rel_drift", "lmag_rel_drift", "phase_unwrapped"]].sort_values("time")
    fmm_small = fmm[["time", "sep", "com_r", "energy_rel_drift", "lmag_rel_drift", "phase_unwrapped"]].sort_values("time")
    tol = 0.5 * min(np.median(np.diff(mg_small["time"])), np.median(np.diff(fmm_small["time"])))
    merged = pd.merge_asof(
        mg_small,
        fmm_small,
        on="time",
        direction="nearest",
        tolerance=tol,
        suffixes=("_mg", "_fmm"),
    ).dropna()
    merged["sep_diff"] = merged["sep_fmm"] - merged["sep_mg"]
    merged["phase_diff"] = merged["phase_unwrapped_fmm"] - merged["phase_unwrapped_mg"]
    merged["com_r_diff"] = merged["com_r_fmm"] - merged["com_r_mg"]
    merged["energy_drift_diff"] = merged["energy_rel_drift_fmm"] - merged["energy_rel_drift_mg"]
    merged["lmag_drift_diff"] = merged["lmag_rel_drift_fmm"] - merged["lmag_rel_drift_mg"]
    return merged


# %%
mg = load_solver_history(RUN_ROOT, "mg")
fmm = load_solver_history(RUN_ROOT, "fmm")
summary = build_summary(mg, fmm)
comparison = compare_histories(mg, fmm)

print(f"Using run root: {RUN_ROOT}")
display(summary)
display(comparison.head())


# %%
fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)

for df in (mg, fmm):
    label = df["solver"].iloc[0]
    axes[0, 0].plot(df["time"], df["sep"], label=label)
    axes[0, 1].plot(df["time"], df["sep_rel_dev"], label=label)
    axes[0, 2].plot(df["time"], df["phase_unwrapped"] / (2.0 * np.pi), label=label)
    axes[1, 0].plot(df["time"], df["com_r"], label=label)
    axes[1, 1].plot(df["time"], df["energy_rel_drift"], label=label)
    axes[1, 2].plot(df["time"], df["lmag_rel_drift"], label=label)

axes[0, 0].set_title("Separation")
axes[0, 1].set_title("Relative Separation Deviation")
axes[0, 2].set_title("Completed Orbits")
axes[1, 0].set_title("COM Drift")
axes[1, 1].set_title("Relative Energy Drift")
axes[1, 2].set_title("Relative Angular Momentum Drift")

for ax in axes.flat:
    ax.set_xlabel("Time")
    ax.grid(alpha=0.3)

axes[0, 0].set_ylabel("r")
axes[0, 1].set_ylabel("(r-r0)/r0")
axes[0, 2].set_ylabel("phase / 2pi")
axes[1, 0].set_ylabel("|r_com|")
axes[1, 1].set_ylabel("dE / |E0|")
axes[1, 2].set_ylabel("dL / |L0|")
axes[0, 0].legend(loc="best")

fig.savefig(SAVE_DIR / "binary_mg_fmm_dynamics.png", bbox_inches="tight")
plt.show()


# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

axes[0].plot(mg["x1_relcom0"], mg["y1_relcom0"], label="MG body 1")
axes[0].plot(mg["x2_relcom0"], mg["y2_relcom0"], label="MG body 2")
axes[0].plot(fmm["x1_relcom0"], fmm["y1_relcom0"], "--", label="FMM body 1")
axes[0].plot(fmm["x2_relcom0"], fmm["y2_relcom0"], "--", label="FMM body 2")
axes[0].set_title("Orbital Tracks")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
axes[0].axis("equal")
axes[0].grid(alpha=0.3)
axes[0].legend(loc="best")

axes[1].plot(comparison["time"], comparison["sep_diff"])
axes[1].set_title("FMM - MG Separation")
axes[1].set_xlabel("Time")
axes[1].set_ylabel("dr")
axes[1].grid(alpha=0.3)

axes[2].plot(comparison["time"], comparison["phase_diff"])
axes[2].set_title("FMM - MG Phase")
axes[2].set_xlabel("Time")
axes[2].set_ylabel("dphase")
axes[2].grid(alpha=0.3)

fig.savefig(SAVE_DIR / "binary_mg_fmm_solver_deltas.png", bbox_inches="tight")
plt.show()


# %%
summary.to_csv(SAVE_DIR / "binary_solver_summary.csv", index=False)
comparison.to_csv(SAVE_DIR / "binary_solver_deltas.csv", index=False)

print("Saved analysis products to", SAVE_DIR)


# %%
if __name__ == "__main__":
    print(summary.to_string(index=False))

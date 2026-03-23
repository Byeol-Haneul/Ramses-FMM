#!/usr/bin/env python3
"""Aggregate close-binary MG/FMM comparisons across AMR depth."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def load_case(run_root: Path, levelmax: int) -> dict[str, float]:
    summary = pd.read_csv(run_root / "analysis" / "binary_solver_summary.csv")
    deltas = pd.read_csv(run_root / "analysis" / "binary_solver_deltas.csv")

    mg = summary.loc[summary["solver"] == "MG"].iloc[0]
    fmm = summary.loc[summary["solver"] == "FMM"].iloc[0]
    max_phase_diff = float(deltas["phase_diff"].abs().max())
    final_phase_diff = float(abs(deltas["phase_diff"].iloc[-1]))
    fmm_orbits = float(abs(fmm["orbits_completed"]))

    return {
        "levelmax": levelmax,
        "mg_sep_dev": float(mg["max_abs_sep_rel_dev"]),
        "fmm_sep_dev": float(fmm["max_abs_sep_rel_dev"]),
        "fmm_energy_drift": float(fmm["final_abs_energy_rel_drift"]),
        "fmm_lmag_drift": float(fmm["final_abs_lmag_rel_drift"]),
        "max_phase_diff": max_phase_diff,
        "final_phase_diff": final_phase_diff,
        "fmm_orbits": fmm_orbits,
        "phase_diff_per_fmm_orbit": max_phase_diff / fmm_orbits if fmm_orbits else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", nargs=2, metavar=("LEVELMAX", "RUN_ROOT"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for levelmax_str, run_root_str in args.run:
        rows.append(load_case(Path(run_root_str), int(levelmax_str)))

    df = pd.DataFrame(rows).sort_values("levelmax").reset_index(drop=True)
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "binary_amr_depth_summary.csv"
    png_path = outdir / "binary_amr_depth_summary.png"
    df.to_csv(csv_path, index=False)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(df["levelmax"], df["phase_diff_per_fmm_orbit"], marker="o", color="#c44900")
    axes[0].set_xlabel("levelmax")
    axes[0].set_ylabel("max phase diff / FMM orbit [rad]")
    axes[0].set_title("Phase Mismatch")

    axes[1].plot(df["levelmax"], df["fmm_sep_dev"], marker="o", color="#2a6f97")
    axes[1].plot(df["levelmax"], df["mg_sep_dev"], marker="s", color="#6c757d")
    axes[1].set_xlabel("levelmax")
    axes[1].set_ylabel("max |sep rel dev|")
    axes[1].set_title("Separation Drift")
    axes[1].legend(["FMM", "MG"])

    axes[2].plot(df["levelmax"], df["fmm_lmag_drift"], marker="o", color="#588157")
    axes[2].set_xlabel("levelmax")
    axes[2].set_ylabel("final |L drift|")
    axes[2].set_title("Angular Momentum Drift")

    for ax in axes:
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(png_path, dpi=180)

    print(df.to_string(index=False))
    print(f"Saved {csv_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()

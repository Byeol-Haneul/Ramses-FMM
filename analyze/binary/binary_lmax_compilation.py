#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from binary_vs_analytic import (
    add_error_columns,
    analytic_reference,
    apply_plot_style,
    load_solver_history,
    positive_for_log,
    style_axis,
)


LEVEL_COLORS = {5: "#4C5B7A", 6: "#2A9D8F", 7: "#C06C2B", 8: "#9C6644"}
SOLVER_LINESTYLES = {"MG": "-", "FMM": "--"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", nargs=2, metavar=("LMAX", "RUN_ROOT"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    apply_plot_style()
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 8.4), constrained_layout=True)

    for lmax_str, run_root_str in args.run:
        lmax = int(lmax_str)
        run_root = Path(run_root_str).resolve()
        mg = load_solver_history(run_root, "mg")
        fmm = load_solver_history(run_root, "fmm")
        ref = analytic_reference(mg)
        add_error_columns(mg, ref)
        add_error_columns(fmm, ref)

        for df in (mg, fmm):
            solver = df["solver"].iloc[0]
            color = LEVEL_COLORS.get(lmax, "#7B7B7B")
            linestyle = SOLVER_LINESTYLES[solver]
            axes[0, 0].plot(df["time"], df["sep_rel_err"], color=color, linestyle=linestyle, lw=2.2)
            axes[0, 1].plot(df["time"], df["phase_err"], color=color, linestyle=linestyle, lw=2.2)
            axes[0, 2].plot(df["time"], df["body_err"], color=color, linestyle=linestyle, lw=2.2)
            axes[1, 0].plot(df["time"], positive_for_log(df["com_r"]), color=color, linestyle=linestyle, lw=2.2)
            axes[1, 1].plot(df["time"], df["energy_rel_drift"], color=color, linestyle=linestyle, lw=2.2)
            axes[1, 2].plot(df["time"], df["lmag_rel_drift"], color=color, linestyle=linestyle, lw=2.2)

    axes[0, 0].set_title("Relative Separation Error")
    axes[0, 0].set_xlabel("Time")
    axes[0, 0].set_ylabel(r"$(r-r_0)/r_0$")

    axes[0, 1].set_title("Phase Error")
    axes[0, 1].set_xlabel("Time")
    axes[0, 1].set_ylabel(r"$\Delta \phi$ [rad]")

    axes[0, 2].set_title("Mean Body Position Error")
    axes[0, 2].set_xlabel("Time")
    axes[0, 2].set_ylabel(r"$\frac{|{\bf x}_1-{\bf x}_{1,\rm ana}| + |{\bf x}_2-{\bf x}_{2,\rm ana}|}{2}$")

    axes[1, 0].set_title("Center-of-Mass Drift")
    axes[1, 0].set_xlabel("Time")
    axes[1, 0].set_ylabel(r"$|{\bf x}_{\rm COM}-{\bf x}_{{\rm COM},0}|$ (log)")
    axes[1, 0].set_yscale("log")

    axes[1, 1].set_title("Relative Energy Drift")
    axes[1, 1].set_xlabel("Time")
    axes[1, 1].set_ylabel(r"$\Delta E/|E_0|$")

    axes[1, 2].set_title("Relative Angular Momentum Drift")
    axes[1, 2].set_xlabel("Time")
    axes[1, 2].set_ylabel(r"$\Delta L/|L_0|$")

    for ax in axes.flat:
        ax.axhline(0.0, color="black", lw=0.9, alpha=0.18)
        style_axis(ax)

    plotted_levels = sorted({int(lmax_str) for lmax_str, _ in args.run})
    color_handles = [
        plt.Line2D([0], [0], color=LEVEL_COLORS.get(level, "#7B7B7B"), lw=2.6, label=fr"$l_{{\max}}={level}$")
        for level in plotted_levels
    ]
    solver_handles = [
        plt.Line2D([0], [0], color="black", lw=2.6, linestyle=SOLVER_LINESTYLES[solver], label=solver)
        for solver in ("MG", "FMM")
    ]

    leg1 = axes[0, 0].legend(handles=color_handles, loc="upper right", title="Resolution")
    axes[0, 0].add_artist(leg1)
    axes[0, 1].legend(handles=solver_handles, loc="upper right", title="Solver")

    fig.suptitle(r"Equal-Mass Binary Diagnostics, $l_{\min}=5$", fontsize=15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    print(args.output)


if __name__ == "__main__":
    main()

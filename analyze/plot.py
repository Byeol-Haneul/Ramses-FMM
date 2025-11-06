# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "serif"
plt.rcParams["figure.dpi"] = 150
plt.rcParams.update({"font.size": 10})


def load_ic(ic_name):
    names = ["x", "y", "z", "fx", "fy", "fz", "phi"]
    fmm_file = f"../out_fmm/{ic_name}_fmm.out"
    mg_file = f"../out_mg/{ic_name}_mg.out"
    fmm = pd.read_csv(fmm_file, header=None, sep=r"\s+", names=names)
    mg = pd.read_csv(mg_file, header=None, sep=r"\s+", names=names)
    for df in [fmm, mg]:
        df["x"] -= df["x"].min()
        df["y"] -= df["y"].min()
        df["z"] -= df["z"].min()
    return fmm, mg


def plot_2d(fmm, mg, ic_name, z_slice):
    fmm_slice = fmm[fmm["z"] == z_slice]
    mg_slice = mg[mg["z"] == z_slice]

    fmm_grid = fmm_slice.pivot(index="y", columns="x", values="phi")
    mg_grid = mg_slice.pivot(index="y", columns="x", values="phi")

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    im0 = axs[0].imshow(mg_grid.values, origin="lower", cmap="viridis")
    axs[0].set_title(f"MG  (z={z_slice})")
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(fmm_grid.values, origin="lower", cmap="viridis")
    axs[1].set_title(f"FMM  (z={z_slice})")
    plt.colorbar(im1, ax=axs[1])

    diff = mg_grid.values - fmm_grid.values
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = np.where(mg_grid.values != 0, diff / mg_grid.values, 0)
    im2 = axs[2].imshow(rel_diff, origin="lower", cmap="coolwarm")
    axs[2].set_title("Relative Difference")
    plt.colorbar(im2, ax=axs[2])

    plt.tight_layout()

    # Save figure
    save_dir = "./fig"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(f"{save_dir}/{ic_name}_2d_{z_slice}.png")
    plt.show()


def plot_profile(fmm, mg, z_slice):
    fmm_slice = fmm[fmm["z"] == z_slice]
    mg_slice = mg[mg["z"] == z_slice]

    # --- Convert to 2D grids ---
    fmm_grid = fmm_slice.pivot(index="y", columns="x", values="phi").values
    mg_grid = mg_slice.pivot(index="y", columns="x", values="phi").values

    ny, nx = fmm_grid.shape
    dx = 1.0 / 128.0  # cell size
    L = 1.0
    xc_phys = L / 2.0
    yc_phys = L / 2.0

    # --- Physical coordinates of cell centers ---
    y_idx, x_idx = np.indices(fmm_grid.shape)
    x_phys = (x_idx + 0.5) * dx
    y_phys = (y_idx + 0.5) * dx

    # --- Physical radial distance from particle ---
    r = np.sqrt((x_phys - xc_phys) ** 2 + (y_phys - yc_phys) ** 2)

    # --- Flatten arrays ---
    r_flat = r.flatten()
    fmm_flat = fmm_grid.flatten()
    mg_flat = mg_grid.flatten()

    # --- Mask out the singularity at the center ---
    mask = r_flat > 0
    r_fit = r_flat[mask]
    fmm_fit = fmm_flat[mask]
    mg_fit = mg_flat[mask]

    # --- Sort by radius ---
    sort_idx = np.argsort(r_fit)
    r_fit = r_fit[sort_idx]
    fmm_fit = fmm_fit[sort_idx]
    mg_fit = mg_fit[sort_idx]

    # --- Model: pure 1/r ---
    def phi_model(r, A):
        return A / r

    # --- Fit both ---
    popt_fmm, _ = curve_fit(phi_model, r_fit, fmm_fit)
    popt_mg, _ = curve_fit(phi_model, r_fit, mg_fit)

    A_fmm = popt_fmm[0]
    A_mg = popt_mg[0]

    print(f"FMM fit: A = {A_fmm:.6e}")
    print(f"MG  fit: A = {A_mg:.6e}")

    # --- Prepare smooth r range ---
    r_plot = np.logspace(np.log10(dx), np.log10(r_fit.max()), 300)
    phi_fmm_fit = phi_model(r_plot, A_fmm)
    phi_mg_fit = phi_model(r_plot, A_mg)

    # --- Residuals vs 1/r ---
    fmm_resid = 100 * (fmm_fit + 1 / r_fit) / (-1 / r_fit)
    mg_resid = 100 * (mg_fit + 1 / r_fit) / (-1 / r_fit)

    # --- Plot main + residuals ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7, 9), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    # --- Main log–log plot ---
    ax1.loglog(r_fit, np.abs(fmm_fit), "o", markersize=3, alpha=0.5, label="|FMM data|")
    ax1.loglog(r_fit, np.abs(mg_fit), "o", markersize=3, alpha=0.5, label="|MG data|")

    ref_r = np.array([r_plot.min(), r_plot.max()])
    ref_phi = 1 / ref_r
    ax1.loglog(ref_r, ref_phi, "k--", lw=1, label="∝ 1/r ref")

    ax1.set_ylabel("|φ|")
    ax1.set_title("Radial potential profile and residuals (log–log)")
    ax1.legend()
    ax1.grid(True, which="both", ls="--", alpha=0.3)

    # --- Residual plot ---
    ax2.plot(r_fit, fmm_resid, "r-", alpha=0.6, label="FMM residual (%)")
    ax2.plot(r_fit, mg_resid, "b-", alpha=0.6, label="MG residual (%)")
    ax2.axhline(0, color="k", lw=1)
    ax2.set_xscale("log")
    ax2.set_xlabel("r")
    ax2.set_ylabel("Residual [%]")
    ax2.legend()
    ax2.grid(True, which="both", ls="--", alpha=0.3)
    ax2.set_ylim(-10, )

    print(max(fmm_resid))

    plt.tight_layout()
    os.makedirs("./fig", exist_ok=True)
    plt.savefig("./fig/single_profile.png", dpi=200)
    plt.show()


def main():
    if len(sys.argv) < 3:
        print("Usage: python plot.py ic_name z_slice")
        return

    ic = sys.argv[1]
    z_slice = int(sys.argv[2])

    fmm, mg = load_ic(ic)
    plot_2d(fmm, mg, ic, z_slice)
    if ic == "single" and z_slice == 63:
        plot_profile(fmm, mg, z_slice)


if __name__ == "__main__":
    main()

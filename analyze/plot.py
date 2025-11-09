# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "serif"
plt.rcParams["figure.dpi"] = 150
plt.rcParams.update({"font.size": 10})


def load_data(filename1, filename2, out_dir="../out"):
    names = ["x", "y", "z", "fx", "fy", "fz", "phi"]

    file1 = os.path.join(out_dir, filename1)
    file2 = os.path.join(out_dir, filename2)

    df1 = pd.read_csv(file1, header=None, sep=r"\s+", names=names)
    df2 = pd.read_csv(file2, header=None, sep=r"\s+", names=names)

    for df in [df1, df2]:
        df["x"] -= df["x"].min()
        df["y"] -= df["y"].min()
        df["z"] -= df["z"].min()

    return df1, df2, filename1, filename2


def plot_2d(df1, df2, label1, label2, z_slice):
    df1_slice = df1[df1["z"] == z_slice]
    df2_slice = df2[df2["z"] == z_slice]

    grid1 = df1_slice.pivot(index="y", columns="x", values="phi")
    grid2 = df2_slice.pivot(index="y", columns="x", values="phi")

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    im0 = axs[0].imshow(grid1.values, origin="lower", cmap="viridis")
    axs[0].set_title(f"{label1}  (z={z_slice})")
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(grid2.values, origin="lower", cmap="viridis")
    axs[1].set_title(f"{label2}  (z={z_slice})")
    plt.colorbar(im1, ax=axs[1])

    diff = grid2.values - grid1.values
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = np.where(grid1.values != 0, diff / grid1.values, 0)
    im2 = axs[2].imshow(rel_diff, origin="lower", cmap="coolwarm")
    axs[2].set_title("Relative Difference")
    plt.colorbar(im2, ax=axs[2])

    plt.tight_layout()
    os.makedirs("./fig", exist_ok=True)
    plt.savefig(f"./fig/compare_2d_{z_slice}_{label1}_vs_{label2}.png", dpi=200)
    plt.show()


def plot_profile(df1, df2, label1, label2, z_slice):
    df1_slice = df1[df1["z"] == z_slice]
    df2_slice = df2[df2["z"] == z_slice]

    grid1 = df1_slice.pivot(index="y", columns="x", values="phi").values
    grid2 = df2_slice.pivot(index="y", columns="x", values="phi").values

    dx = 1.0 / 128.0
    L = 1.0
    xc_phys = L / 2.0
    yc_phys = L / 2.0

    y_idx, x_idx = np.indices(grid1.shape)
    x_phys = (x_idx + 0.5) * dx
    y_phys = (y_idx + 0.5) * dx
    r = np.sqrt((x_phys - xc_phys) ** 2 + (y_phys - yc_phys) ** 2)

    r_flat = r.flatten()
    vals1 = grid1.flatten()
    vals2 = grid2.flatten()

    mask = r_flat > 0
    r_fit = r_flat[mask]
    vals1 = vals1[mask]
    vals2 = vals2[mask]

    sort_idx = np.argsort(r_fit)
    r_fit = r_fit[sort_idx]
    vals1 = vals1[sort_idx]
    vals2 = vals2[sort_idx]

    # Residuals relative to the analytical 1/r profile
    resid1 = 100 * (vals1 + 1 / r_fit) / (-1 / r_fit)
    resid2 = 100 * (vals2 + 1 / r_fit) / (-1 / r_fit)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7, 9), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )

    # Plot raw profiles
    ax1.loglog(r_fit, np.abs(vals1), "o", markersize=3, alpha=0.5, label=f"|{label1}|")
    ax1.loglog(r_fit, np.abs(vals2), "o", markersize=3, alpha=0.5, label=f"|{label2}|")

    # Reference 1/r line
    ref_r = np.array([r_fit.min(), r_fit.max()])
    ref_phi = 1 / ref_r
    ax1.loglog(ref_r, ref_phi, "k--", lw=1, label="∝ 1/r ref")

    ax1.set_ylabel("|φ|")
    ax1.set_title("Radial potential profile and residuals (log–log)")
    ax1.legend()
    ax1.grid(True, which="both", ls="--", alpha=0.3)

    # Residuals
    ax2.plot(r_fit, resid1, "r-", alpha=0.6, label=f"{label1} residual (%)")
    ax2.plot(r_fit, resid2, "b-", alpha=0.6, label=f"{label2} residual (%)")
    ax2.axhline(0, color="k", lw=1)
    ax2.set_xscale("log")
    ax2.set_xlabel("r")
    ax2.set_ylabel("Residual [%]")
    ax2.legend()
    ax2.grid(True, which="both", ls="--", alpha=0.3)
    ax2.set_ylim(-10, 10)

    plt.tight_layout()
    os.makedirs("./fig", exist_ok=True)
    plt.savefig(f"./fig/profile_{label1}_vs_{label2}.png", dpi=200)
    plt.show()


def main():
    if len(sys.argv) < 4:
        print("Usage: python plot.py <file1> <file2> <z_slice>")
        return

    file1, file2, z_slice = sys.argv[1], sys.argv[2], int(sys.argv[3])

    df1, df2, label1, label2 = load_data(file1, file2)
    plot_2d(df1, df2, label1, label2, z_slice)

    # optional: profile for specific z
    if z_slice == 63:
        plot_profile(df1, df2, label1, label2, z_slice)


if __name__ == "__main__":
    main()

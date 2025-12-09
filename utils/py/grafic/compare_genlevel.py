#!/usr/bin/env python3
"""
Compare multi-level ICs generated with turb_genlevel.py.

Generates:
  1) A 2x2 panel of mid-plane slices of the z-momentum (rho * w)
     for levels 8, 7, 6, 5 (generated at gen_level=8).
  2) An overlaid isotropic 3D power spectrum of the same z-momentum
     component for each level.

Defaults assume ICs live under:
  ic_turb/ic_turb_lvl{lvl}_gen8_3d
relative to the repository root.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import grafic


LEVELS = [8, 7, 6, 5]
GEN_LEVEL = 8
IC_DIR_TMPL = "ic_turb/ic_turb_lvl{lvl}_gen{gen_level}_3d"


def generate_k_grid(n1: int, n2: int, n3: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return kx, ky, kz, kmag arrays for an FFT grid with unit box length."""
    kx = np.fft.fftfreq(n1, d=1.0 / n1).reshape(n1, 1, 1)
    ky = np.fft.fftfreq(n2, d=1.0 / n2).reshape(1, n2, 1)
    kz = np.fft.fftfreq(n3, d=1.0 / n3).reshape(1, 1, n3)
    kmag = np.sqrt(kx * kx + ky * ky + kz * kz)
    return kx, ky, kz, kmag


def read_field(path: Path, name: str) -> np.ndarray:
    g = grafic.Grafic()
    g.read(str(path / name))
    return np.asarray(g.data, dtype=np.float32)


def load_level(path: Path) -> dict[str, np.ndarray]:
    rho = read_field(path, "ic_d")
    u = read_field(path, "ic_u")
    v = read_field(path, "ic_v")
    w = read_field(path, "ic_w")
    return {
        "rho": rho,
        "u": u,
        "v": v,
        "w": w,
        "mz": rho * w,
    }


def power_spectrum_isotropic(field: np.ndarray, kmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute isotropic 3D power spectrum of a scalar field."""
    n1, n2, n3 = field.shape
    kx, ky, kz, kmag = generate_k_grid(n1, n2, n3)
    fft = np.fft.fftn(field)
    power = (np.abs(fft) ** 2) / (n1 * n2 * n3) ** 2  # volume-normalized

    kint = np.rint(kmag).astype(np.int64)
    mask = kint <= kmax

    weights = power[mask].ravel()
    bins = kint[mask].ravel()
    ps_sum = np.bincount(bins, weights=weights, minlength=kmax + 1)
    ps_count = np.bincount(bins, minlength=kmax + 1)
    ps = ps_sum / np.maximum(ps_count, 1)
    k = np.arange(kmax + 1, dtype=np.float64)
    return k, ps


def plot_slices(data_by_level: dict[int, dict[str, np.ndarray]], outfile: Path) -> None:
    levels_sorted = sorted(data_by_level.keys(), reverse=True)
    slices = []
    for lvl in levels_sorted:
        arr = data_by_level[lvl]["mz"]
        k_mid = arr.shape[2] // 2
        slices.append(arr[:, :, k_mid])
    vmax = max(np.abs(s).max() for s in slices)
    vmin = -vmax

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    for ax, lvl, sli in zip(axes.ravel(), levels_sorted, slices):
        im = ax.imshow(sli.T, origin="lower", cmap="RdBu_r", vmin=vmin, vmax=vmax)
        ax.set_title(f"Level {lvl} (n={sli.shape[0]})")
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.9)
    cbar.set_label("rho * w")
    fig.suptitle("Mid-plane z-momentum slices (rho * w)", fontsize=14)
    fig.savefig(outfile, dpi=150)
    plt.close(fig)


def plot_power_spectra(ps_by_level: dict[int, tuple[np.ndarray, np.ndarray]], outfile: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    for lvl, (k, ps) in sorted(ps_by_level.items(), reverse=True):
        mask = k > 0  # drop k=0 DC component
        ax.loglog(k[mask], ps[mask], label=f"lvl {lvl} (n={2**lvl})")
    ax.set_xlabel("k (grid units)")
    ax.set_ylabel("P(k) of rho*w")
    ax.set_title("Isotropic power spectra (rho*w)")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Compare multi-level ICs generated with turb_genlevel.py")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[4],
        help="Repository root containing ic_turb/ (default: auto-detected)",
    )
    parser.add_argument(
        "--gen_level",
        type=int,
        default=GEN_LEVEL,
        help="Generation level used to create the ICs (default: 8)",
    )
    args = parser.parse_args()

    root = args.root
    gen_level = args.gen_level
    data_by_level: dict[int, dict[str, np.ndarray]] = {}
    ps_by_level: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    # Load fields
    for lvl in LEVELS:
        ic_dir = root / IC_DIR_TMPL.format(lvl=lvl, gen_level=gen_level)
        if not ic_dir.exists():
            raise SystemExit(f"Missing IC directory: {ic_dir}")
        data = load_level(ic_dir)
        data_by_level[lvl] = data

    # Determine common kmax (Nyquist of coarsest grid)
    min_n = min(data["mz"].shape[0] for data in data_by_level.values())
    kmax_common = min_n // 2

    for lvl, data in data_by_level.items():
        k, ps = power_spectrum_isotropic(data["mz"], kmax=kmax_common)
        ps_by_level[lvl] = (k, ps)

    out_dir = Path(__file__).resolve().parent
    plot_slices(data_by_level, out_dir / "compare_genlevel_slices.png")
    plot_power_spectra(ps_by_level, out_dir / "compare_genlevel_power.png")
    print("Saved plots:")
    print(out_dir / "compare_genlevel_slices.png")
    print(out_dir / "compare_genlevel_power.png")


if __name__ == "__main__":
    main()


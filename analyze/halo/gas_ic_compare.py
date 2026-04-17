#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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
    raise FileNotFoundError("Could not locate repository root from current path.")


try:
    HERE = Path(__file__).resolve().parent
except NameError:
    HERE = Path.cwd()

REPO_ROOT = find_repo_root(HERE)
UTILS_PY = REPO_ROOT / "utils/py"
if str(UTILS_PY) not in sys.path:
    sys.path.append(str(UTILS_PY))

import miniramses as ram


M_H = 1.6605390e-24


@dataclass
class HaloParams:
    halo_center: np.ndarray
    v_200: float = 150.0
    concentration: float = 10.0
    baryon_fraction: float = 0.15
    spin_lambda: float = 0.04
    halo_nmin: float = 1e-6
    halo_tmin: float = 1e6
    halo_eps: float = 0.01
    x_h: float = 0.76


@dataclass
class LeafHydroCells:
    x: np.ndarray
    dx: np.ndarray
    level: np.ndarray
    u: np.ndarray
    ndim: int
    nvar: int
    ncell: int


@dataclass
class HydroFieldMap:
    rho: int
    vx: int
    vy: int
    vz: int
    pressure: int


class HaloICModel:
    def __init__(
        self,
        params: HaloParams,
        boxlen: float,
        unit_d: float,
        unit_l: float,
        unit_t: float,
        pressure_grid_size: int = 20000,
    ):
        if pressure_grid_size < 512:
            raise ValueError("pressure_grid_size should be >= 512.")
        self.params = params
        self.boxlen = float(boxlen)
        self.unit_d = float(unit_d)
        self.unit_l = float(unit_l)
        self.unit_t = float(unit_t)
        self.center = np.array(
            [
                self.params.halo_center[0] + 0.5 * self.boxlen,
                self.params.halo_center[1] + 0.5 * self.boxlen,
                self.params.halo_center[2] + 0.5 * self.boxlen,
            ],
            dtype=np.float64,
        )

        self.scale_v = self.unit_l / self.unit_t
        self.scale_t2 = M_H / 1.3806490e-16 * self.scale_v**2
        self.scale_nh = self.params.x_h / M_H * self.unit_d

        self.hsmall = 0.7
        self.pi = math.acos(-1.0)
        self.hub = self.hsmall * 100.0 / 1.0e3
        self.rhocrit = 3.0 * self.hub**2 / (8.0 * self.pi)

        self.v200 = float(self.params.v_200)
        self.r200 = self.v200 / self.hsmall
        self.m200_msun = (self.r200 * self.hsmall / 1.63e-2) ** 3 / self.hsmall
        self.c = float(self.params.concentration)
        self.fb = float(self.params.baryon_fraction)
        self.rs = self.r200 / self.c

        self.eps_dimless = float(self.params.halo_eps) / self.rs
        self.eps_phys = float(self.params.halo_eps)
        self.rmax_dimless = 2.0 * self.r200 / self.rs
        self.rmax_phys = self.rmax_dimless * self.rs

        self.rhos = (
            self.rhocrit
            * 200.0
            / 3.0
            * self.c**3
            / (math.log(1.0 + self.c) - self.c / (1.0 + self.c))
        )
        self.m200_internal = self.m200_msun / 2.3262e5
        self.j_max = self.params.spin_lambda * self.v200 * self.r200
        self.dmin = self.params.halo_nmin / self.scale_nh
        self.pmin = self.dmin * self.params.halo_tmin / self.scale_t2

        self._build_pressure_table(pressure_grid_size)

    def _mass_enclosed(self, rr: np.ndarray) -> np.ndarray:
        return 4.0 * self.pi * self.rhos * self.rs**3 * (np.log1p(rr) - rr / (1.0 + rr))

    def _rho_nfw(self, rr: np.ndarray) -> np.ndarray:
        return self.fb * self.rhos / rr / (1.0 + rr) ** 2

    def _build_pressure_table(self, ngrid: int) -> None:
        lnr = np.linspace(np.log(self.eps_dimless), np.log(self.rmax_dimless), ngrid, dtype=np.float64)
        rr = np.exp(lnr)
        rho = self._rho_nfw(rr)
        mass = self._mass_enclosed(rr)
        integrand = rho * mass / (rr * self.rs)
        dln = np.diff(lnr)
        seg = 0.5 * (integrand[:-1] + integrand[1:]) * dln
        p = np.zeros_like(rr)
        p[:-1] = np.cumsum(seg[::-1])[::-1]
        self._lnr_grid = lnr
        self._p_grid = p

    def evaluate(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        center_override: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        center = self.center if center_override is None else np.asarray(center_override, dtype=np.float64)
        xx = np.asarray(x, dtype=np.float64) - center[0]
        yy = np.asarray(y, dtype=np.float64) - center[1]
        zz = np.asarray(z, dtype=np.float64) - center[2]

        rc = np.sqrt(xx**2 + yy**2)
        rr_dimless = np.sqrt(xx**2 + yy**2 + zz**2) / self.rs
        rr_dimless = np.maximum(rr_dimless, self.eps_dimless)
        rc = np.maximum(rc, self.eps_phys)

        rho = self._rho_nfw(rr_dimless)
        rho = np.maximum(rho, self.dmin)

        mass = self._mass_enclosed(rr_dimless)
        vphi = self.j_max * mass / self.m200_internal / (rr_dimless * self.rs)
        vx = -vphi * yy / rc
        vy = +vphi * xx / rc
        vz = np.zeros_like(vx)

        lnr = np.log(rr_dimless)
        pressure = np.interp(lnr, self._lnr_grid, self._p_grid, left=self._p_grid[0], right=0.0)
        pressure = np.maximum(pressure, self.pmin)

        rr_phys = rr_dimless * self.rs
        return {
            "rho": rho,
            "pressure": pressure,
            "vx": vx,
            "vy": vy,
            "vz": vz,
            "vphi": vphi,
            "r_phys": rr_phys,
            "rc_phys": rc,
            "xx": xx,
            "yy": yy,
            "zz": zz,
        }


def parse_fortran_float(text: str) -> float:
    return float(text.strip().replace("D", "E").replace("d", "e"))


def parse_fortran_float_list(text: str) -> list[float]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    return [parse_fortran_float(p) for p in parts]


def extract_namelist_blocks(text: str, block_name: str) -> list[str]:
    pattern = re.compile(rf"(?is)&\s*{re.escape(block_name)}\b(.*?)/")
    return [m.group(1) for m in pattern.finditer(text)]


def parse_halo_params_from_namelist(path: Path) -> HaloParams:
    params = HaloParams(halo_center=np.zeros(3, dtype=np.float64))
    if (not path.exists()) or (not path.is_file()):
        return params

    content = path.read_text()
    halo_blocks = extract_namelist_blocks(content, "HALO_PARAMS")
    run_blocks = extract_namelist_blocks(content, "RUN_PARAMS")

    halo_center_indexed = params.halo_center.copy()
    for block in halo_blocks:
        for raw_line in block.splitlines():
            line = raw_line.split("!")[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip().lower()
            value = value.strip().rstrip(",")
            if key.startswith("halo_center("):
                m = re.search(r"halo_center\((\d+)\)", key)
                if m:
                    idx = int(m.group(1)) - 1
                    if 0 <= idx < 3:
                        halo_center_indexed[idx] = parse_fortran_float(value)
            elif key == "halo_center":
                vals = parse_fortran_float_list(value)
                for i in range(min(3, len(vals))):
                    halo_center_indexed[i] = vals[i]
            elif key == "v_200":
                params.v_200 = parse_fortran_float(value)
            elif key == "concentration":
                params.concentration = parse_fortran_float(value)
            elif key == "baryon_fraction":
                params.baryon_fraction = parse_fortran_float(value)
            elif key == "lambda":
                params.spin_lambda = parse_fortran_float(value)
            elif key == "halo_eps":
                params.halo_eps = parse_fortran_float(value)
            elif key == "halo_nmin":
                params.halo_nmin = parse_fortran_float(value)
            elif key == "halo_tmin":
                params.halo_tmin = parse_fortran_float(value)

    params.halo_center = halo_center_indexed
    for block in run_blocks:
        for raw_line in block.splitlines():
            line = raw_line.split("!")[0].strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip().lower()
            value = value.strip().rstrip(",")
            if key == "x_h":
                params.x_h = parse_fortran_float(value)
    return params


def read_hydro_var_names(output_dir: Path, hydro_prefix: str) -> list[str]:
    header_path = output_dir / f"{hydro_prefix}_header.txt"
    if not header_path.exists():
        return []
    names: list[str] = []
    for line in header_path.read_text().splitlines():
        m = re.search(r"variable\s*#\s*\d+\s*:\s*(.+)$", line.strip(), re.IGNORECASE)
        if m:
            names.append(m.group(1).strip().lower())
    return names


def build_field_map(nvar: int, names: list[str]) -> HydroFieldMap:
    rho = 0
    pressure = min(4, nvar - 1)
    vx, vy, vz = 1, 2, 3

    if names:
        for i, name in enumerate(names[:nvar]):
            if "density" in name:
                rho = i
            if "velocity_x" in name or name == "vx":
                vx = i
            if "velocity_y" in name or name == "vy":
                vy = i
            if "velocity_z" in name or name == "vz":
                vz = i
            if "pressure" in name and "total" not in name:
                pressure = i

    if nvar < 5:
        raise ValueError(f"Hydro file has nvar={nvar}, expected at least 5.")
    return HydroFieldMap(rho=rho, vx=vx, vy=vy, vz=vz, pressure=pressure)


def discover_outputs(run_root: Path, hydro_prefix: str = "hydro") -> list[int]:
    out = []
    for p in run_root.glob("output_*"):
        if not p.is_dir():
            continue
        m = re.fullmatch(r"output_(\d+)", p.name)
        if not m:
            continue
        if (p / f"{hydro_prefix}.00001").exists():
            out.append(int(m.group(1)))
    return sorted(out)


def resolve_run_root(path: Path, hydro_prefix: str) -> tuple[Path, list[int] | None]:
    p = path.resolve()
    if re.fullmatch(r"output_(\d+)", p.name) and (p / f"{hydro_prefix}.00001").exists():
        return p.parent, [int(p.name.split("_")[1])]
    return p, None


def load_leaf_hydro_cells(nout: int, run_root: Path, hydro_prefix: str = "hydro") -> LeafHydroCells:
    amr = ram.rd_amr(nout, path=str(run_root))
    hydro = ram.rd_hydro(nout, path=str(run_root), prefix=hydro_prefix)

    ndim = amr[0].ndim
    if ndim != 3:
        raise ValueError(f"This halo analysis expects 3D data; got ndim={ndim}.")
    nvar = hydro[0].nvar
    boxlen = float(amr[0].boxlen)
    nlevelmax = len(amr)

    offsets = np.zeros((ndim, 2**ndim), dtype=np.float64)
    offsets[0, :] = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5]
    offsets[1, :] = [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5]
    offsets[2, :] = [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5]

    x_chunks: list[np.ndarray] = []
    u_chunks: list[np.ndarray] = []
    dx_chunks: list[np.ndarray] = []
    level_chunks: list[np.ndarray] = []

    for ilev in range(nlevelmax):
        dx = 0.5 * boxlen / (2**ilev)
        for ind in range(2**ndim):
            leaf_mask = ~amr[ilev].refined[ind]
            nc = int(np.count_nonzero(leaf_mask))
            if nc == 0:
                continue

            xg = amr[ilev].xg[:, leaf_mask]
            xc = (2.0 * xg + 1.0 + offsets[:, ind][:, None]) * dx
            uc = hydro[ilev].u[:, ind, leaf_mask]

            x_chunks.append(xc)
            u_chunks.append(uc)
            dx_chunks.append(np.full(nc, dx, dtype=np.float64))
            level_chunks.append(np.full(nc, ilev, dtype=np.int16))

    x = np.concatenate(x_chunks, axis=1) if x_chunks else np.empty((ndim, 0))
    u = np.concatenate(u_chunks, axis=1) if u_chunks else np.empty((nvar, 0))
    dx = np.concatenate(dx_chunks) if dx_chunks else np.empty((0,), dtype=np.float64)
    level = np.concatenate(level_chunks) if level_chunks else np.empty((0,), dtype=np.int16)
    return LeafHydroCells(x=x, dx=dx, level=level, u=u, ndim=ndim, nvar=nvar, ncell=int(dx.size))


def align_positions(
    x: np.ndarray,
    center: np.ndarray,
    boxlen: float,
    enabled: bool,
    shift_trials: Iterable[float],
) -> tuple[np.ndarray, np.ndarray]:
    xout = np.array(x, copy=True, dtype=np.float64)
    applied = np.zeros(x.shape[0], dtype=np.float64)
    trials = [0.0] if not enabled else [float(t) * boxlen for t in shift_trials]

    for d in range(x.shape[0]):
        med = float(np.nanmedian(xout[d]))
        best_shift = trials[0]
        best_err = abs((med + best_shift) - center[d])
        for s in trials[1:]:
            err = abs((med + s) - center[d])
            if err < best_err:
                best_err = err
                best_shift = s
        xout[d] += best_shift
        applied[d] = best_shift

    return xout, applied


def weighted_rms(values: np.ndarray, weights: np.ndarray) -> float:
    wsum = float(np.sum(weights))
    if wsum <= 0.0:
        return float("nan")
    return float(np.sqrt(np.sum(weights * values**2) / wsum))


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    wsum = float(np.sum(weights))
    if wsum <= 0.0:
        return float("nan")
    return float(np.sum(weights * values) / wsum)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    if cdf[-1] <= 0:
        return float("nan")
    cdf /= cdf[-1]
    return float(np.interp(q, cdf, v))


def compute_radial_profile(
    r: np.ndarray,
    mass: np.ndarray,
    abs_rel_rho: np.ndarray,
    abs_rel_p: np.ndarray,
    abs_vr: np.ndarray,
    abs_dvphi: np.ndarray,
    bins: int,
    rmax: float,
) -> pd.DataFrame:
    valid = (r > 0.0) & np.isfinite(r) & np.isfinite(mass) & (mass > 0.0) & (r <= rmax)
    if not np.any(valid):
        return pd.DataFrame(
            columns=[
                "r_lo",
                "r_hi",
                "r_mid",
                "ncell",
                "mass",
                "mean_abs_rel_rho",
                "mean_abs_rel_pressure",
                "mean_abs_vr",
                "mean_abs_dvphi",
            ]
        )

    rv = r[valid]
    mv = mass[valid]
    q = max(np.min(rv), 1e-12)
    edges = np.geomspace(q, rmax, bins + 1)
    rows = []
    for i in range(bins):
        lo = edges[i]
        hi = edges[i + 1]
        mask = (rv >= lo) & (rv < hi)
        if not np.any(mask):
            continue
        mm = mv[mask]
        wsum = np.sum(mm)
        rows.append(
            {
                "r_lo": lo,
                "r_hi": hi,
                "r_mid": math.sqrt(lo * hi),
                "ncell": int(np.count_nonzero(mask)),
                "mass": float(wsum),
                "mean_abs_rel_rho": float(np.sum(abs_rel_rho[valid][mask] * mm) / wsum),
                "mean_abs_rel_pressure": float(np.sum(abs_rel_p[valid][mask] * mm) / wsum),
                "mean_abs_vr": float(np.sum(abs_vr[valid][mask] * mm) / wsum),
                "mean_abs_dvphi": float(np.sum(abs_dvphi[valid][mask] * mm) / wsum),
            }
        )
    return pd.DataFrame(rows)


def analyze_snapshot(
    run_root: Path,
    nout: int,
    model: HaloICModel,
    hydro_prefix: str,
    align: bool,
    shift_trials: Iterable[float],
    radial_bins: int,
    radial_rmax: float | None,
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    output_dir = run_root / f"output_{nout:05d}"
    info = ram.rd_info(nout, path=str(run_root))
    cells = load_leaf_hydro_cells(nout=nout, run_root=run_root, hydro_prefix=hydro_prefix)
    names = read_hydro_var_names(output_dir, hydro_prefix)
    fmap = build_field_map(cells.nvar, names)

    _, shifts = align_positions(
        cells.x,
        center=model.center,
        boxlen=float(info.boxlen),
        enabled=align,
        shift_trials=shift_trials,
    )
    # Shift the analytic frame instead of shifting output cells.
    # If x_aligned = x_raw + shift, then using raw x is equivalent to
    # evaluating IC at center' = center - shift.
    analytic_center = model.center - shifts

    rho = cells.u[fmap.rho].astype(np.float64)
    vx = cells.u[fmap.vx].astype(np.float64)
    vy = cells.u[fmap.vy].astype(np.float64)
    vz = cells.u[fmap.vz].astype(np.float64)
    p = cells.u[fmap.pressure].astype(np.float64)

    ref = model.evaluate(cells.x[0], cells.x[1], cells.x[2], center_override=analytic_center)
    vol = cells.dx**3
    mass = np.maximum(rho, 0.0) * vol

    r = ref["r_phys"]
    rc = np.maximum(ref["rc_phys"], model.eps_phys)
    vr = (ref["xx"] * vx + ref["yy"] * vy + ref["zz"] * vz) / np.maximum(r, model.eps_phys)
    vphi = (-ref["yy"] * vx + ref["xx"] * vy) / rc

    rel_rho = (rho - ref["rho"]) / np.maximum(ref["rho"], model.dmin)
    rel_p = (p - ref["pressure"]) / np.maximum(ref["pressure"], model.pmin)
    dvphi = vphi - ref["vphi"]

    abs_rel_rho = np.abs(rel_rho)
    abs_rel_p = np.abs(rel_p)
    abs_vr = np.abs(vr)
    abs_dvphi = np.abs(dvphi)

    summary = {
        "nout": int(nout),
        "time": float(info.time),
        "ncell": int(cells.ncell),
        "boxlen": float(info.boxlen),
        "shift_x": float(shifts[0]),
        "shift_y": float(shifts[1]),
        "shift_z": float(shifts[2]),
        "analytic_center_x": float(analytic_center[0]),
        "analytic_center_y": float(analytic_center[1]),
        "analytic_center_z": float(analytic_center[2]),
        "gas_mass": float(np.sum(mass)),
        "mean_vr_mass": weighted_mean(vr, mass),
        "rms_abs_vr_mass": weighted_rms(vr, mass),
        "p99_abs_vr_mass": weighted_quantile(abs_vr, mass, 0.99),
        "mean_dvphi_mass": weighted_mean(dvphi, mass),
        "rms_abs_dvphi_mass": weighted_rms(dvphi, mass),
        "p99_abs_dvphi_mass": weighted_quantile(abs_dvphi, mass, 0.99),
        "rms_rel_rho_mass": weighted_rms(rel_rho, mass),
        "p99_abs_rel_rho_mass": weighted_quantile(abs_rel_rho, mass, 0.99),
        "rms_rel_pressure_mass": weighted_rms(rel_p, mass),
        "p99_abs_rel_pressure_mass": weighted_quantile(abs_rel_p, mass, 0.99),
        "rms_rel_rho_vol": weighted_rms(rel_rho, vol),
        "rms_rel_pressure_vol": weighted_rms(rel_p, vol),
        "rms_abs_vr_vol": weighted_rms(vr, vol),
        "rms_abs_dvphi_vol": weighted_rms(dvphi, vol),
    }

    level_rows = []
    for lev in np.unique(cells.level):
        m = cells.level == lev
        msum = float(np.sum(mass[m]))
        if msum <= 0.0:
            continue
        level_rows.append(
            {
                "nout": nout,
                "level": int(lev),
                "ncell": int(np.count_nonzero(m)),
                "mass_fraction": float(msum / np.sum(mass)),
                "rms_abs_vr_mass": weighted_rms(vr[m], mass[m]),
                "rms_abs_dvphi_mass": weighted_rms(dvphi[m], mass[m]),
                "rms_rel_rho_mass": weighted_rms(rel_rho[m], mass[m]),
                "rms_rel_pressure_mass": weighted_rms(rel_p[m], mass[m]),
            }
        )
    level_df = pd.DataFrame(level_rows)

    profile_rmax = radial_rmax if radial_rmax is not None else model.rmax_phys
    radial_df = compute_radial_profile(
        r=r,
        mass=mass,
        abs_rel_rho=abs_rel_rho,
        abs_rel_p=abs_rel_p,
        abs_vr=abs_vr,
        abs_dvphi=abs_dvphi,
        bins=radial_bins,
        rmax=profile_rmax,
    )
    return summary, radial_df, level_df


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "font.size": 11,
            "font.family": "serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def plot_timeseries(summary: pd.DataFrame, out_png: Path) -> None:
    if summary.empty:
        return
    apply_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    x = summary["time"].to_numpy()

    ax = axes[0]
    ax.plot(x, np.abs(summary["rms_rel_rho_mass"]), lw=2, label=r"$\rho$")
    ax.plot(x, np.abs(summary["rms_rel_pressure_mass"]), lw=2, label=r"$P$")
    ax.set_title("Mass-Weighted RMS Relative Error")
    ax.set_xlabel("Time")
    ax.set_ylabel("RMS")
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1]
    ax.plot(x, np.abs(summary["rms_abs_vr_mass"]), lw=2, label=r"$|v_r|$")
    ax.plot(x, np.abs(summary["rms_abs_dvphi_mass"]), lw=2, label=r"$|v_\phi-v_{\phi,0}|$")
    ax.set_title("Mass-Weighted Velocity Drift")
    ax.set_xlabel("Time")
    ax.set_ylabel("Velocity")
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[2]
    ax.plot(x, summary["mean_vr_mass"], lw=2, label=r"$\langle v_r \rangle_M$")
    ax.plot(x, summary["mean_dvphi_mass"], lw=2, label=r"$\langle \Delta v_\phi \rangle_M$")
    ax.axhline(0.0, color="black", lw=1, alpha=0.4)
    ax.set_title("Net Signed Drift")
    ax.set_xlabel("Time")
    ax.set_ylabel("Velocity")
    ax.grid(alpha=0.25)
    ax.legend()

    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def plot_radial_profile(radial_df: pd.DataFrame, out_png: Path, nout: int) -> None:
    if radial_df.empty:
        return
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    r = radial_df["r_mid"].to_numpy()

    ax = axes[0]
    ax.plot(r, radial_df["mean_abs_rel_rho"], lw=2, label=r"$|\Delta \rho/\rho_0|$")
    ax.plot(r, radial_df["mean_abs_rel_pressure"], lw=2, label=r"$|\Delta P/P_0|$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Radius")
    ax.set_ylabel("Mass-weighted mean abs relative error")
    ax.set_title(f"Gas Scalar Drift (nout={nout:05d})")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1]
    ax.plot(r, radial_df["mean_abs_vr"], lw=2, label=r"$|v_r|$")
    ax.plot(r, radial_df["mean_abs_dvphi"], lw=2, label=r"$|v_\phi-v_{\phi,0}|$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Radius")
    ax.set_ylabel("Mass-weighted mean abs velocity drift")
    ax.set_title(f"Gas Motion Drift (nout={nout:05d})")
    ax.grid(alpha=0.25)
    ax.legend()

    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def pick_profile_nout(nouts: list[int], spec: str) -> int:
    if not nouts:
        raise ValueError("No outputs available.")
    if spec == "last":
        return nouts[-1]
    if spec == "first":
        return nouts[0]
    return int(spec)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare RAMSES halo gas state against patch/init/halo/condinit.f90 ground truth."
    )
    parser.add_argument("--run-dir", type=Path, required=True, help="Directory containing output_XXXXX folders.")
    parser.add_argument("--nout", type=int, nargs="*", default=None, help="Specific output numbers to analyze.")
    parser.add_argument("--namelist", type=Path, default=None, help="Namelist path for HALO/RUN params.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Where analysis files are written.")
    parser.add_argument("--hydro-prefix", type=str, default="hydro")
    parser.add_argument("--x-h", type=float, default=None, help="Override hydrogen mass fraction X_H.")
    parser.add_argument("--pressure-grid-size", type=int, default=20000)
    parser.add_argument("--radial-bins", type=int, default=80)
    parser.add_argument(
        "--radial-rmax",
        type=float,
        default=None,
        help="Max radius for radial profile (default: 2*r200 from halo model).",
    )
    parser.add_argument("--profile-nout", type=str, default="last", help="'first', 'last', or integer nout.")
    parser.add_argument("--align", dest="align", action="store_true")
    parser.add_argument("--no-align", dest="align", action="store_false")
    parser.set_defaults(align=True)
    parser.add_argument(
        "--shift-trials",
        type=float,
        nargs="+",
        default=[0.0, -0.5, 0.5, -1.0, 1.0],
        help="Multiples of boxlen tested for coordinate alignment.",
    )
    args = parser.parse_args()

    run_root, single_nout = resolve_run_root(args.run_dir, args.hydro_prefix)
    discovered = discover_outputs(run_root, args.hydro_prefix)
    if not discovered and single_nout is None:
        raise FileNotFoundError(f"No output_XXXXX directories with '{args.hydro_prefix}.00001' under {run_root}")

    if args.nout:
        nouts = sorted(args.nout)
    elif single_nout is not None:
        nouts = single_nout
    else:
        nouts = discovered
    if not nouts:
        raise ValueError("No outputs selected for analysis.")

    first_output_dir = run_root / f"output_{nouts[0]:05d}"
    info0 = ram.rd_info(nouts[0], path=str(run_root))

    namelist_path: Path | None = args.namelist
    if namelist_path is None:
        candidate = first_output_dir / "namelist.txt"
        namelist_path = candidate if candidate.exists() else None

    halo_params = (
        parse_halo_params_from_namelist(namelist_path)
        if namelist_path is not None
        else HaloParams(np.zeros(3))
    )
    if args.x_h is not None:
        halo_params.x_h = float(args.x_h)

    model = HaloICModel(
        params=halo_params,
        boxlen=float(info0.boxlen),
        unit_d=float(info0.unit_d),
        unit_l=float(info0.unit_l),
        unit_t=float(info0.unit_t),
        pressure_grid_size=args.pressure_grid_size,
    )

    outdir = (args.output_dir or (run_root / "analysis/halo")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, float]] = []
    level_dfs: list[pd.DataFrame] = []
    profile_target = pick_profile_nout(nouts, args.profile_nout)
    radial_profile = pd.DataFrame()

    print(f"[halo] run_root={run_root}")
    if namelist_path:
        print(f"[halo] namelist={namelist_path}")
    print(f"[halo] outputs={nouts}")
    print(f"[halo] model center={model.center}, rs={model.rs:.6g}, r200={model.r200:.6g}, rmax={model.rmax_phys:.6g}")

    for nout in nouts:
        print(f"[halo] analyzing output_{nout:05d}")
        summary, radial_df, level_df = analyze_snapshot(
            run_root=run_root,
            nout=nout,
            model=model,
            hydro_prefix=args.hydro_prefix,
            align=args.align,
            shift_trials=args.shift_trials,
            radial_bins=args.radial_bins,
            radial_rmax=args.radial_rmax,
        )
        summary_rows.append(summary)
        if not level_df.empty:
            level_dfs.append(level_df)
        if nout == profile_target:
            radial_profile = radial_df.copy()

    summary_df = pd.DataFrame(summary_rows).sort_values("nout").reset_index(drop=True)
    summary_path = outdir / "gas_ic_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[halo] wrote {summary_path}")

    if level_dfs:
        level_all = pd.concat(level_dfs, ignore_index=True)
        level_path = outdir / "gas_ic_level_stats.csv"
        level_all.to_csv(level_path, index=False)
        print(f"[halo] wrote {level_path}")

    if not radial_profile.empty:
        radial_path = outdir / f"gas_ic_radial_profile_nout{profile_target:05d}.csv"
        radial_profile.to_csv(radial_path, index=False)
        print(f"[halo] wrote {radial_path}")

    ts_png = outdir / "gas_ic_timeseries.png"
    plot_timeseries(summary_df, ts_png)
    print(f"[halo] wrote {ts_png}")

    if not radial_profile.empty:
        radial_png = outdir / f"gas_ic_radial_nout{profile_target:05d}.png"
        plot_radial_profile(radial_profile, radial_png, profile_target)
        print(f"[halo] wrote {radial_png}")

    print("[halo] done")


if __name__ == "__main__":
    main()

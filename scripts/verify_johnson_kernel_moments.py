#!/usr/bin/env python3
"""
Cross-check MC trinomial kernel moments vs mini-ramses-dev-2 Johnson SU/SB fitters (mirror move_fine.f90).

Requires: numpy
Optional: scipy (compare SU solve to johnsonsu if installed)

Usage:
  python3 scripts/verify_johnson_kernel_moments.py
"""

from __future__ import annotations

import math

import numpy as np


def mc_kernel_skewness_paper(pr: float, pl: float) -> float:
    cm = pr - pl
    cp = pr + pl
    var_y = cp - cm * cm
    if var_y <= 0:
        return 0.0
    mu3 = cm * (1.0 - 3.0 * cp + 2.0 * cm * cm)
    return mu3 / var_y**1.5


def mc_kernel_kurtosis_excess(pr: float, pl: float) -> float:
    cm = pr - pl
    cp = pr + pl
    var_y = cp - cm * cm
    if var_y <= 0:
        return 0.0
    m1, m2, m3, m4 = cm, cp, cm, cp
    mu4 = m4 - 4 * m3 * m1 + 6 * m2 * m1 * m1 - 3 * m1**4
    return mu4 / (var_y * var_y) - 3.0


def johnson_su_moments_std(eta: float, delta: float):
    """Returns (skew, excess_kurt, ok). Mirrors Fortran johnson_su_moments_std."""
    if delta <= 1e-15:
        return 0.0, 0.0, False
    omega = math.exp(1.0 / (delta * delta))
    if omega > 1e200 or math.isnan(omega):
        return 0.0, 0.0, False
    if omega <= 1.0:
        return 0.0, 0.0, False
    lam2 = 2.0 / ((omega - 1.0) * (omega * math.cosh(2.0 * eta) + 1.0))
    if lam2 <= 0 or math.isnan(lam2):
        return 0.0, 0.0, False
    lam = math.sqrt(lam2)
    var1 = 0.5 * lam2 * (omega - 1.0) * (omega * math.cosh(2.0 * eta) + 1.0)
    if abs(var1 - 1.0) > 1e-6:
        return 0.0, 0.0, False
    skew = (
        -lam**3
        * math.sqrt(omega)
        * (omega - 1.0) ** 2
        * (omega * (omega + 2.0) * math.sinh(3.0 * eta) + 3.0 * math.sinh(eta))
        / 4.0
    )
    K1 = omega**2 * (omega**4 + 2 * omega**3 + 3 * omega**2 - 3) * math.cosh(4.0 * eta)
    K2 = 4.0 * omega**2 * (omega + 2.0) * math.cosh(3.0 * eta)
    K3 = 3.0 * (2.0 * omega + 1.0)
    k_ex = lam**4 * (omega - 1.0) ** 2 * (K1 + K2 + K3) / 8.0 - 3.0
    return skew, k_ex, True


def main() -> None:
    rng = np.random.default_rng(42)
    worst = 0.0
    for _ in range(500):
        pr = float(rng.random())
        pl = float(rng.random())
        if pr + pl > 1.0:
            continue
        var_y = (pr + pl) - (pr - pl) ** 2
        if var_y < 0.08:
            continue
        g1 = mc_kernel_skewness_paper(pr, pl)
        kex = mc_kernel_kurtosis_excess(pr, pl)
        # Monte Carlo trinomial vs formulas
        n = 500_000
        y = rng.choice([-1.0, 0.0, 1.0], size=n, p=[pl, 1 - pr - pl, pr])
        ys = (y - (pr - pl)) / math.sqrt((pr + pl) - (pr - pl) ** 2)
        g1_mc = float(np.mean(ys**3))
        kex_mc = float(np.mean(ys**4) - 3.0)
        worst = max(worst, abs(g1 - g1_mc), abs(kex - kex_mc))
    print(f"Trinomial vs MC sample max error (skew/kurt): {worst:.3e}")

    # Spot-check SU moments routine
    sk, ku, ok = johnson_su_moments_std(0.3, 1.1)
    print(f"SU moments sample eta=0.3 delta=1.1 -> skew={sk:.6f} k_ex={ku:.6f} ok={ok}")

    try:
        from scipy import stats

        # Compare one SU solve by scanning eta,delta to hit target (coarse)
        g1_t, kex_t = 0.2, -0.5
        best = 1e300
        for eta in np.linspace(-2, 2, 41):
            for delta in np.linspace(0.3, 2.0, 35):
                sk, ku, ok_m = johnson_su_moments_std(float(eta), float(delta))
                if not ok_m:
                    continue
                err = abs(sk - g1_t) + abs(ku - kex_t)
                if err < best:
                    best = err
        print(f"Best coarse SU match error for target (0.2,-0.5): {best:.4f}")
    except ImportError:
        print("scipy not installed; skipped optional comparison.")


if __name__ == "__main__":
    main()

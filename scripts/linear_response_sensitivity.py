"""How can the rotation curve set the pattern-speed sensitivity? Linear epicyclic response to a rotating m = 2 bar.

A star on a near-circular orbit at radius R, forced by an m = 2 bar rotating at Omega_b, responds with an amplitude
proportional to 1 / D, D = kappa^2 - 4 (Omega - Omega_b)^2 (D = 0 at the inner and outer Lindblad resonances). The
fractional change of that response per fractional change of the pattern speed is
    s(R) = | d ln(1/D) / d ln Omega_b | = | 8 Omega_b (Omega - Omega_b) / D |.
It depends on the rotation curve only (Omega and kappa of the axisymmetrised potential), not on the stars or the
integration. Large s: the predicted kinematics there change strongly with Omega_b; it diverges at a resonance (where
linear theory fails and orbits are trapped, but the response is strongest). Times the bar forcing Q_T(R), it indicates
where the pattern-speed information comes from.

Usage: python linear_response_sensitivity.py
Output: results/linear_response_sensitivity.csv and a printed table
"""
import numpy as np
import pandas as pd

from common import RESULTS_DIR
from resonances_vs_pattern_speed import curves, R
from compare_bar_strength_and_resonances import forces, R as RQ
import steady_state_likelihood as ssl

CASES = [("t7 N-body family", "nbody", 26.0), ("Sormani+22", "sormani22", 32.5), ("Hunter+24", "mw", 32.5)]
RS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]


def main():
    rows = []
    for lab, gal, ob in CASES:
        om, kap = curves(gal)
        D = kap ** 2 - 4 * (om - ob) ** 2
        s = np.abs(8 * ob * (om - ob) / D)
        _, qt = forces(ssl.bar_aligned_potential(1.0, 1.0, gal))
        sv = np.interp(RS, R, s)
        qv = np.interp(RS, RQ, qt)
        band = (R >= 2) & (R <= 6)
        print(f"{lab:18s} Omega_b {ob:4.1f}: s(R) at R = " + " ".join(f"{r:g}:{v:5.1f}" for r, v in zip(RS, sv)))
        print(f"{'':18s}             Q_T * s:          " + " ".join(f"{r:g}:{v:5.2f}" for r, v in zip(RS, sv * qv)))
        print(f"{'':18s}             median s over 2-6 kpc {np.median(s[band]):.1f}; fraction of 2-6 kpc with s > 5: "
              f"{np.mean(s[band] > 5):.2f}; D = 0 (Lindblad resonance) at R = "
              + ", ".join(f"{R[i]:.2f}" for i in np.where(np.diff(np.sign(D)) != 0)[0]))
        for r, v, q in zip(RS, sv, qv):
            rows.append(dict(family=lab, Omega_b=ob, R=r, s=v, Q_T=q, Q_T_times_s=v * q))
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "linear_response_sensitivity.csv", index=False)


if __name__ == "__main__":
    main()

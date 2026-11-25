"""Do resonance diagnostics separate a galaxy that constrains Omega_b (N-body t = 2, 10%) from one that does not
(Sormani+22, 38%)? (user, 2026-11-25: Auriga galaxies with Sormani+22-like ILR and ILR/CR still constrain 8-22%.)

For the 12 N-body snapshots (own finite-difference Omega_b; trial-family potential dm_spherical + Nbody_baryon_axi from
the portable drive) and the five Milky Way models (Omega_b 30.53): ILR; (Omega - kappa/2)/Omega_b at R = 2, 3, 4, 5 kpc;
mean |(Omega - kappa/2)/Omega_b - 1| over 2-6 kpc; median and fraction > 5 of the linear sensitivity
s = |8 Omega_b (Omega - Omega_b) / (kappa^2 - 4 (Omega - Omega_b)^2)| over 2-6 kpc; d ln v_c / d ln R over 2-6 kpc;
and the sqrt(2N)-tempered Omega_b width. Spearman of each diagnostic against the width over the 12 snapshots.
Result (2026-11-25): t = 2 and Sormani+22 are alike in every diagnostic (gap 0.51 / 0.56, median s 0.91 / 1.09) yet
constrain 10% / 38%; none correlates with width within the snapshots (|rho| <= 0.26, p >= 0.42).
Usage: python frequency_gap_nbody_vs_milky_way.py   (run from a directory holding plot_milky_way_rotation_curves.py,
       plot_rotation_curves_mocks_vs_milky_way.py and plot_resonances_vs_constraint_simulations_and_milky_way.py)
Output: results/frequency_gap_nbody_vs_milky_way.csv
"""
import csv
import os
import sys

import numpy as np
import agama
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_milky_way_rotation_curves as mw
from plot_rotation_curves_mocks_vs_milky_way import NB
from plot_resonances_vs_constraint_simulations_and_milky_way import NBODY, USB, REPO

R = mw.R
BAND = (R >= 2) & (R <= 6)


def stats(vc, ob):
    om, kap = mw.freqs(vc)
    g = (om - kap / 2) / ob
    s = np.abs(8 * ob * (om - ob) / (kap ** 2 - 4 * (om - ob) ** 2))
    il = [r for r in mw.crossings(om - kap / 2, ob) if r < (mw.crossings(om, ob) or [99])[0]]
    return dict(ILR=il[-1] if il else np.nan, g2=np.interp(2, R, g), g3=np.interp(3, R, g), g4=np.interp(4, R, g),
                g5=np.interp(5, R, g), gap=np.mean(np.abs(g[BAND] - 1)), s_med=np.median(s[BAND]),
                s_gt5=np.mean(s[BAND] > 5), slope=np.polyfit(np.log(R[BAND]), np.log(vc[BAND]), 1)[0])


def main():
    rows = []
    for t in sorted(NB):
        d = f"{USB}/nbody_potential" if t == 7 else f"{USB}/nbody_t{t:g}_potential"
        pot = agama.Potential(agama.Potential(f"{d}/dm_spherical.ini"), agama.Potential(f"{d}/Nbody_baryon_axi.ini"))
        r = stats(mw.vc_phi_avg(pot, R), NB[t][0])
        r.update(galaxy=f"N-body t={t:g}", Omega_b=NB[t][0], width_pct=NBODY[t][1])
        rows.append(r)
    for lab, faxi, _, _, _ in mw.MODELS:
        r = stats(mw.vc_axi(agama.Potential(faxi), R), mw.OMB)
        r.update(galaxy=lab, Omega_b=mw.OMB, width_pct={"Sormani": 38.1, "Hunter": 43.3}.get(lab.split()[0], np.nan))
        rows.append(r)
    nb = [r for r in rows if r["galaxy"].startswith("N-body")]
    for k in ("ILR", "gap", "s_med", "s_gt5", "slope"):
        rho, p = spearmanr([r[k] for r in nb], [r["width_pct"] for r in nb])
        print(f"Spearman over the 12 N-body, {k} vs width: {rho:+.2f} (p {p:.2f})")
    keys = ["galaxy", "Omega_b", "ILR", "g2", "g3", "g4", "g5", "gap", "s_med", "s_gt5", "slope", "width_pct"]
    with open(os.path.join(REPO, "results", "frequency_gap_nbody_vs_milky_way.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: r[k] for k in keys} for r in rows])


if __name__ == "__main__":
    main()

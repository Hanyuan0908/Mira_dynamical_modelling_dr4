"""Are the cluster-built Auriga trial-family "axisymmetric" files axisymmetric? (user, 2026-11-26)

The 13 Auriga potential folders were copied from the cluster (/data/hz420/.../mock/au<N>_potential, built there
2026-11-04/05) to the portable drive. The N-body families built on the cluster on 2026-10-29/30 had "_axi" files
identical to the barred ones (PROGRESS 2026-11-25). Same test here, per galaxy:
  - m=2/m=0 of the midplane radial force at R = 1..5 kpc of Au<N>_baryon_axi.ini (0 to round-off if axisymmetric),
  - the ratio m=2(baryon_axi)/m=2(baryon_full) at R = 3, 4, 5 kpc (1 = the "axi" file carries the whole bar),
  - the m=2 phase difference between the two (0 = the same bar orientation),
  - the trial potential EXACTLY as the pipeline builds it (steady_state_likelihood.bar_aligned_potential: the barred
    file turned by 27 deg onto x, the "axi" file not turned, merged into one mmax-8 CylSpline), at S = 1 and
    bar_frac = 0 (1e-4), 0.5, 1, 1.5, 2: the m=2/m=0 of F_R at R = 4 kpc relative to bar_frac = 1, and how far the
    bar is turned from bar_frac = 1. With correct files: 0, 0.5, 1, 1.5, 2 and no turning.
    (A first version omitted the 27 deg turn and wrongly found the bar frozen; withdrawn, PROGRESS 2026-11-26.)
The laptop-built Au18 family (mock/au18_potential in the repo, 2026-10-27) is included as the reference.
Usage: python check_auriga_family_axisymmetry.py  ->  results/auriga_family_axisymmetry_check.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import agama

from common import ROOT, RESULTS_DIR
import galaxies
import steady_state_likelihood as ssl

agama.setUnits(mass=1, length=1, velocity=1)
USB = "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/mock"
HALOS = [5, 9, 10, 12, 14, 17, 18, 20, 22, 24, 26, 27, 28]
PH = np.linspace(0, 2 * np.pi, 72, endpoint=False)


def m2(pot, R):
    """Complex m=2 coefficient and m=0 mean of the midplane radial force on a ring."""
    X = np.column_stack([R * np.cos(PH), R * np.sin(PH), np.zeros_like(PH)])
    F = pot.force(X)
    FR = F[:, 0] * np.cos(PH) + F[:, 1] * np.sin(PH)
    return 2 * np.mean(FR * np.exp(-2j * PH)), np.mean(FR)


def check(label, d, prefix):
    full, axi = agama.Potential(f"{d}/{prefix}_baryon_full.ini"), agama.Potential(f"{d}/{prefix}_baryon_axi.ini")
    fixed = agama.Potential(agama.Potential(f"{d}/dm.ini"), agama.Potential(f"{d}/bh.ini"))
    r = dict(galaxy=label, folder=d)
    r["axi_m2_over_m0_max_R1to5"] = max(abs(m2(axi, R)[0]) / abs(m2(axi, R)[1]) for R in (1, 2, 3, 4, 5))
    for R in (3, 4, 5):
        ca, cf = m2(axi, R)[0], m2(full, R)[0]
        r[f"m2_axi_over_full_R{R}"] = abs(ca) / abs(cf)
        r[f"phase_diff_deg_R{R}"] = np.degrees(0.5 * np.angle(ca / cf)) if abs(ca) > 1e-12 * abs(cf) else np.nan
    g = f"au{prefix[2:]}"
    galaxies.galaxy(g)
    galaxies.GALAXIES[g]["pot_dir"] = Path(d)
    ssl._CACHE.pop(g, None)
    res = {}
    for f in (1.0, 1e-4, 0.5, 1.5, 2.0):
        c, c0 = m2(ssl.bar_aligned_potential(f, 1.0, g), 4.0)
        res[f] = (abs(c) / abs(c0), c)
    for f in (1e-4, 0.5, 1.0, 1.5, 2.0):
        k = "0" if f == 1e-4 else f"{f:g}"
        r[f"pipeline_bar_force_rel_R4_barfrac{k}"] = res[f][0] / res[1.0][0]
        r[f"pipeline_bar_turned_deg_R4_barfrac{k}"] = np.degrees(0.5 * np.angle(res[f][1] / res[1.0][1])) if res[f][0] > 1e-6 else np.nan
    galaxies.GALAXIES[g]["pot_dir"] = galaxies.MOCK_DIR / f"{g}_potential"
    ssl._CACHE.pop(g, None)
    return r


def main():
    rows = [check("Au18 (laptop-built, 2026-10-27; reference)", ROOT / "mock" / "au18_potential", "Au18")]
    for h in HALOS:
        rows.append(check(f"Au{h} (cluster-built)", f"{USB}/au{h}_potential", f"Au{h}"))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "auriga_family_axisymmetry_check.csv", index=False)
    print(f"{'galaxy':44s} axi m2/m0 | m2(axi)/m2(full) R3 R4 R5 | phase | pipeline bar force at bar_frac 0/0.5/1.5/2 (x truth) | turned [deg]")
    for _, r in df.iterrows():
        print(f"{r.galaxy:44s} {r.axi_m2_over_m0_max_R1to5:9.2e} | {r.m2_axi_over_full_R3:6.3f} {r.m2_axi_over_full_R4:6.3f} "
              f"{r.m2_axi_over_full_R5:6.3f} | {r.phase_diff_deg_R4:6.2f} | "
              + " ".join(f"{r[f'pipeline_bar_force_rel_R4_barfrac{k}']:.3f}" for k in ("0", "0.5", "1.5", "2")) + " | "
              + " ".join(f"{r[f'pipeline_bar_turned_deg_R4_barfrac{k}']:+5.1f}" for k in ("0", "0.5", "1.5", "2")), flush=True)


if __name__ == "__main__":
    main()

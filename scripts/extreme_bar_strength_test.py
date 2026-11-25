"""Sanity test of the bar strength (user, 2026-11-25): does the pipeline respond to bar_frac at all, down to no bar?

Everything else fixed: (|Omega_b|, phi_b, S) at the run's best point (OSARG in Sormani+22) or at the truth (the
steady-state mock in Hunter+24: 32.5, 25, 1, truth bar_frac = 1). bar_frac is stepped through 0 (no bar: the
family's axisymmetric potential only), 0.25, 0.5, 1, 1.5, 2, 3, 4. At each value the stars are integrated (1 Gyr,
50 Myr) and scored exactly as the pipeline (data fitted at 27 deg with the 300-bootstrap covariance; library sfratio,
K = 10, rotated to 27 deg; scheme B) on the default grid and on geo20 |z| < 2. Recorded:
  - chi2 and its exact split into the four fields (d . C^-1 d per block, summed by field);
  - the <v_R> m = 2 amplitude sqrt(cos2^2 + sin2^2) of the library and of the data, median over the default-grid blocks
    with 1 < R < 4 kpc, |z| < 1 (the bar's kinematic signature);
  - the potential's bar force, m=2/m=0 of the midplane F_R at R = 3 kpc (should scale exactly with bar_frac).
Expectations if the code is right: at bar_frac 0 the library quadrupole vanishes and chi2 is far worse than near the
best; the library quadrupole and the bar force grow monotonically with bar_frac.
Usage: python extreme_bar_strength_test.py --case osarg_s22 | h24mock
Output: results/extreme_bar_strength_<case>.csv and a printed table
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from scan_omega_at_fixed_bar import CASES
from bimodality_spatial_origin import block_terms, FIELDS

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
FRACS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
FIXED = {"osarg_s22": dict(Omega=32.06, phi=20.46, S=0.949, note="the run's best evaluated point (bar_frac 1.495 there)"),
         "h24mock": dict(Omega=32.5, phi=25.0, S=1.0, note="the truth (bar_frac 1)")}
GRIDS = [((0.3, 3.0, 1.0, None), "default"), ((0.3, 2.0, 2.0, "geo20"), "geo20 |z|<2")]


def bar_force_m2(pot, R=3.0):
    ph = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    X = np.column_stack([R * np.cos(ph), R * np.sin(ph), np.zeros_like(ph)])
    F = pot.force(X)
    FR = F[:, 0] * np.cos(ph) + F[:, 1] * np.sin(ph)
    return 2 * abs(np.mean(FR * np.exp(-2j * ph))) / abs(np.mean(FR))


def vr_quad(coef, ring_edges, nl):
    """Median <v_R> m=2 amplitude over blocks with 1 < R < 4 kpc in the first |z| layer."""
    rc = 0.5 * (ring_edges[:-1] + ring_edges[1:])
    amps = [np.hypot(coef[i, 0, 1], coef[i, 0, 2]) for i in range(len(coef))
            if i % nl == 0 and 1 < rc[i // nl] < 4 and np.all(np.isfinite(coef[i, 0, 1:3]))]
    return float(np.median(amps))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, choices=list(FIXED))
    args = ap.parse_args()
    c, fx = CASES[args.case], FIXED[args.case]
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    ssl.ACCURACY = old.ACCURACY
    data = {}
    for gargs, gl in GRIDS:
        fs.set_grid(*gargs)
        cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
        data[gl] = (cd0.reshape(fs.NB, 6, fs.NK), cov0, ok, fs.RING_EDGES.copy(), fs.NL, fs.NB)
    fs.set_grid()
    print(f"[{args.case}] {c['label']}: Omega {fx['Omega']}, phi {fx['phi']}, S {fx['S']} fixed ({fx['note']}); "
          f"{len(w):,} stars, sigma_DM {sigma_dm:.4f}", flush=True)
    rows = []
    for f in FRACS:
        pot = ssl.bar_aligned_potential(max(f, 1e-9), fx["S"], c["galaxy"])
        _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(fx["phi"])), -abs(fx["Omega"]))
        traj = np.ascontiguousarray(traj[:, ::2])
        r = dict(bar_frac=f, bar_force_m2_R3=bar_force_m2(pot))
        for gargs, gl in GRIDS:
            fs.set_grid(*gargs)
            cd, cov, ok, edges, nl, nb = data[gl]
            cl, _ = library_fit_empirical(traj, fx["phi"], sigma_dm, dat, ok)
            cl = (cl.reshape(nb, -1) @ fs.rotation(BAR_ANGLE_DEG - fx["phi"]).T).reshape(nb, 6, fs.NK)
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            x2, nd = fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[B])
            bt = block_terms(cd, cl, cov, valid, fs.SCHEMES[B][1])
            r[f"chi2_{gl}"] = x2
            r[f"ndata_{gl}"] = nd
            for k, fld in enumerate(FIELDS.values()):
                r[f"chi2_{gl}_{fld}"] = float(np.nansum(bt[:, k]))
            if gl == "default":
                r["vR_m2_lib"] = vr_quad(cl, edges, nl)
                r["vR_m2_data"] = vr_quad(cd, edges, nl)
        fs.set_grid()
        rows.append(r)
        print(f"   bar_frac {f:4.2f}: bar force m2/m0 at 3 kpc {r['bar_force_m2_R3']:.4f} | <v_R> m=2 (1-4 kpc, |z|<1) library "
              f"{r['vR_m2_lib']:5.2f} vs data {r['vR_m2_data']:5.2f} km/s | chi2 default {r['chi2_default']:7.1f} (N {r['ndata_default']}), "
              f"geo20 {r['chi2_geo20 |z|<2']:6.1f} | default by field: "
              + ", ".join(f"{fld} {r[f'chi2_default_{fld}']:.0f}" for fld in FIELDS.values()), flush=True)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"extreme_bar_strength_{args.case}.csv", index=False)


if __name__ == "__main__":
    main()

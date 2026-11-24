"""The temperature T = sigma_2 (chi2 jitter between two near-identical models) for each Fourier grid.

A coarser grid has fewer data terms, so both the chi2 jitter and the old sqrt(2N) shrink; comparing grids at the
sqrt(2N) tempering would favour coarse grids for no reason. This measures each grid's own jitter, the user's
definition of T (2026-11-22): N_RAND runs at the best evaluated point of the run, all four parameters perturbed by
the same tiny Gaussians as chi2_jitter_small_perturbations.py (seed 2123), plus the point itself. Trajectories are
integrated once (1 Gyr, 50 Myr, as the pipeline) and kept, then scored on every grid exactly as the pipeline: the data
fitted at 27 deg on that grid, the library (sfratio, K = 10) rotated to 27 deg, scheme B (17 terms per block, rings
inside 6 kpc, the run's data covariance).
  sigma_1  std of chi2 - chi2(the point itself);   sigma_2  std of the chi2 difference between two perturbed runs.

Usage: python chi2_jitter_per_grid.py --case osarg_s22 [--n 40]
Output: results/chi2_jitter_per_grid_<case>.csv (one row per run and grid) and a printed summary
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test2b_qmc_distance_error_run import load_design
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from scan_omega_at_fixed_bar import CASES
from chi2_jitter_small_perturbations import RANDOM_WIDTH

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
GRIDS = [(0.3, 3.0, 1.0), (0.5, 2.0, 1.0), (0.5, 2.0, 2.0), (1.0, 2.0, 2.0), (0.5, 3.0, 3.0), (0.3, 2.0, 2.0)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="osarg_s22", choices=list(CASES))
    ap.add_argument("--n", type=int, default=40)
    args = ap.parse_args()
    c = CASES[args.case]
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    ev = pd.read_csv(RESULTS_DIR / c["run"] / "evals.csv").set_index("idx")
    best = design.set_index("idx").loc[int(ev[f"logLs_{B}_sfratio"].idxmax())]
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    tdir = Path.home() / "mira_dr4_runs" / f"jitter_traj_{args.case}"
    tdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2123)
    runs = []
    ssl.ACCURACY = old.ACCURACY
    for r in range(args.n + 1):
        q = best.copy()
        if r:
            for k, wdt in RANDOM_WIDTH.items():
                q[k] = q[k] + wdt * rng.standard_normal()
        path = tdir / f"traj_{r:03d}.npy"
        if not path.exists():
            t0 = time.time()
            pot = ssl.bar_aligned_potential(float(q.bar_frac), float(q.S), c["galaxy"])
            _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(float(q.phi_b_deg))), -abs(float(q.Omega_abs)))
            np.save(path, np.ascontiguousarray(traj[:, ::2]))
            print(f"integrated run {r} ({time.time() - t0:.0f} s)", flush=True)
        runs.append((r, float(q.phi_b_deg), path))
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    rows = []
    for g in GRIDS:
        tag = fs.set_grid(*g) or "default"
        cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
        cd0 = cd0.reshape(fs.NB, 6, fs.NK)
        t0 = time.time()
        for r, phi, path in runs:
            traj = np.load(path)
            cl, _ = library_fit_empirical(traj, phi, sigma_dm, dat, ok)
            cl = (cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - phi).T).reshape(fs.NB, 6, fs.NK)
            valid = ok & np.all(np.isfinite(cd0), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            x2, nd = fs.scheme_chi2(cd0, cl, cov0, valid, fs.SCHEMES[B])
            rows.append(dict(grid=tag, ring_width=g[0], z_max=g[1], layer_width=g[2], run=r, chi2=x2, ndata=nd))
        t = pd.DataFrame([x for x in rows if x["grid"] == tag])
        x0, x = t.chi2.iloc[0], t.chi2.values[1:]
        pairs = (x[:, None] - x[None, :])[np.triu_indices(len(x), 1)]
        nd = int(t.ndata.mode().iloc[0])
        print(f"grid {tag:14s}: N_data {nd:4d}, chi2 of the point {x0:7.1f} (chi2/N {x0 / nd:.3f}); sigma_1 "
              f"{np.std(x - x0, ddof=1):5.1f}, sigma_2 {np.std(pairs, ddof=1):5.1f}, sqrt(2N) {np.sqrt(2 * nd):5.1f}  "
              f"({time.time() - t0:.0f} s)", flush=True)
    fs.set_grid()
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"chi2_jitter_per_grid_{args.case}.csv", index=False)


if __name__ == "__main__":
    main()

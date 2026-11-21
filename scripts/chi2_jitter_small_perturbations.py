"""How much does χ² jitter between models that are essentially the same? The model-to-model scatter that a Δχ²
comparison actually faces.

Every design point is built from the same stars with the same distance-draw seed, so noise common to all models
(which stars there are, their re-drawn distances, the K = 10 draws) cancels in a χ² DIFFERENCE. What does not cancel
is how each potential scrambles the orbital phases at the 15 snapshots. This perturbs the best-fit point by amounts
far below the statistical resolution (|Ω_b| +-0.1, 0.25, 0.5 km/s/kpc; φ_b +-0.1, 0.25 deg; S and bar_frac +-0.005,
0.01), integrates each, scores it exactly as the pipeline does, and takes the scatter of χ² against the data.
A repeat of the best point itself (must give the same χ²) is the control.
χ² is computed with C_data (the pipeline) and with C_data + C_model (model_uncertainty_best_fit.py), Hartlap 0.94.

Usage: python chi2_jitter_small_perturbations.py --case osarg_h24
Output: results/chi2_jitter_<case>.csv
"""
import argparse
import json
import time

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

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
RANDOM_WIDTH = dict(Omega_abs=0.25, phi_b_deg=0.15, S=0.007, bar_frac=0.007)   # far below the statistical resolution
PERTURB = [("Omega_abs", d) for d in (-0.5, -0.25, -0.1, 0.1, 0.25, 0.5)] + \
          [("phi_b_deg", d) for d in (-0.25, -0.1, 0.1, 0.25)] + \
          [("S", d) for d in (-0.01, -0.005, 0.005, 0.01)] + [("bar_frac", d) for d in (-0.01, -0.005, 0.005, 0.01)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="osarg_h24", choices=list(CASES))
    ap.add_argument("--random", type=int, default=0, metavar="N",
                    help="instead of the fixed list, N runs with all four parameters perturbed together, drawn from "
                         "Gaussians of width RANDOM_WIDTH")
    args = ap.parse_args()
    c = CASES[args.case]
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    ev = pd.read_csv(RESULTS_DIR / c["run"] / "evals.csv").set_index("idx")
    ibest = int(ev[f"logLs_{B}_sfratio"].idxmax())
    best = design.set_index("idx").loc[ibest]
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    mu = np.load(RESULTS_DIR / f"model_uncertainty_{args.case}.npz")
    cd0, Cd, ok = mu["coef_data"].reshape(fs.NB, -1), mu["cov_data"], mu["ok"]
    nr, spec = fs.SCHEMES[B]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    h = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)
    Ws = {}
    for tag, C in (("C_data", Cd), ("C_data + C_model", Cd + mu["cov_model"])):
        Ws[tag] = {b: np.linalg.inv(np.linalg.cholesky(C[b][np.ix_(idx, idx)]))      # scheme B: rings inside 6 kpc
                   for b in (r * fs.NL + z for r in range(nr) for z in range(fs.NL)) if ok[b]}
    ssl.ACCURACY = old.ACCURACY
    rows = []
    if args.random:
        rng = np.random.default_rng(2123)
        todo = [("repeat", 0.0)] + [("random", r) for r in range(args.random)]
    else:
        todo = [("repeat", 0.0)] + PERTURB
    for par, dlt in todo:
        q = best.copy()
        if par == "random":
            for k, wdt in RANDOM_WIDTH.items():
                q[k] = q[k] + wdt * rng.standard_normal()
        elif par != "repeat":
            q[par] = q[par] + dlt
        t0 = time.time()
        pot = ssl.bar_aligned_potential(float(q.bar_frac), float(q.S), c["galaxy"])
        _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(float(q.phi_b_deg))), -abs(float(q.Omega_abs)))
        cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, ::2]), float(q.phi_b_deg), sigma_dm, dat, ok)
        cl = cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - float(q.phi_b_deg)).T
        row = dict(param=par, delta=dlt, Omega_abs=float(q.Omega_abs), phi_b_deg=float(q.phi_b_deg), S=float(q.S),
                   bar_frac=float(q.bar_frac))
        for tag, W in Ws.items():
            blocks = [b for b in W if np.all(np.isfinite(cl[b][idx]))]
            row[tag] = h * sum(float(np.sum((W[b] @ (cd0[b][idx] - cl[b][idx])) ** 2)) for b in blocks)
            row["n_blocks"] = len(blocks)
        rows.append(row)
        print(f"  {par:10s} {dlt:+7.3f}: chi2 {row['C_data']:8.1f} | {row['C_data + C_model']:8.1f}  "
              f"({time.time() - t0:.0f} s)", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / f"chi2_jitter_{args.case}{'_random' if args.random else ''}.csv", index=False)
    N = int(t.n_blocks.iloc[0]) * len(idx)
    pert = t[t.param != "repeat"]
    if args.random:
        for tag in Ws:
            x = pert[tag].values
            d = x - t[tag].iloc[0]
            pairs = (x[:, None] - x[None, :])[np.triu_indices(len(x), 1)]
            print(f"[{tag}] {len(x)} random perturbations: dchi2 = chi2 - chi2(best run): mean {d.mean():+.1f}, "
                  f"std {d.std(ddof=1):.1f}, 16/50/84% {np.percentile(d, 16):+.1f}/{np.percentile(d, 50):+.1f}/"
                  f"{np.percentile(d, 84):+.1f}; between two perturbed runs: std {pairs.std(ddof=1):.1f} "
                  f"(sqrt(2) x {pairs.std(ddof=1) / np.sqrt(2):.1f}); best-run chi2/N {t[tag].iloc[0] / N:.3f}")
        return
    for tag in Ws:
        print(f"[{tag}] best-point chi2 {t[tag].iloc[0]:.1f} (chi2/N {t[tag].iloc[0] / N:.3f}); across the "
              f"{len(pert)} tiny perturbations: std {pert[tag].std():.1f}, range {pert[tag].min() - t[tag].iloc[0]:+.1f} "
              f"to {pert[tag].max() - t[tag].iloc[0]:+.1f}; by parameter: " +
              ", ".join(f"{p} {g[tag].std():.1f}" for p, g in pert.groupby("param")))


if __name__ == "__main__":
    main()

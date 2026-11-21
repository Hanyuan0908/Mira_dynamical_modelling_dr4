"""Re-score every design point with one fixed covariance C_data + C_model (model_uncertainty_best_fit.py).

The saved library coefficients of each design point are rotated from the trial bar frame to the data frame (where the
data fit, C_data and C_model live); χ² is unchanged by that rotation, which is orthogonal within each (cos, sin)
pair. Scheme B terms, block-diagonal. The Hartlap factor uses the smaller of the two sample counts (300 bootstraps,
N_REAL realisations), the conservative choice for a sum of two estimated covariances. No integration.

Reported for each covariance (data only; data + model; data + model with star bootstrap): χ²/N at the best point,
and the profile of Δχ² in |Ω_b| (smallest per bin over the other parameters), on the RAW χ² scale (no tempering).

Usage: python rescore_with_model_covariance.py --case osarg_h24
Output: results/<run>/evals_modelcov.csv (chi2 per design point for each covariance)
"""
import argparse
import json

import numpy as np
import pandas as pd

from common import ROOT, RESULTS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test2b_qmc_distance_error_run import load_design
from scan_omega_at_fixed_bar import CASES

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
SDIR = {"osarg_h24": "/Users/hanyuan/mira_dr4_runs/osarg_zhang2024_mw_hunter24",
        "t7_own": "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs/test2_t7_sf_12p5pct",
        "t7_h24": "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs/"
                  "test2_t7_sf_12p5pct_mwH24_osargdesign"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="osarg_h24", choices=list(SDIR))
    args = ap.parse_args()
    c = CASES[args.case]
    mu = np.load(RESULTS_DIR / f"model_uncertainty_{args.case}.npz")
    cd0, Cd, ok, n_real = mu["coef_data"].reshape(fs.NB, -1), mu["cov_data"], mu["ok"], int(mu["n_real"])
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    ev = pd.read_csv(RESULTS_DIR / c["run"] / "evals.csv").set_index("idx")
    ev = ev[ev.kind == "qmc"]
    nr, spec = fs.SCHEMES[B]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    p = len(idx)
    libs = {}
    for i in ev.index:
        cl = np.load(f"{SDIR[args.case]}/sfratio_{i:03d}.npz")["coef_lib"].reshape(fs.NB, -1)
        libs[i] = cl @ fs.rotation(BAR_ANGLE_DEG - float(design.loc[i, "phi_b_deg"])).T
    blocks = [r * fs.NL + z for r in range(nr) for z in range(fs.NL)
              if ok[r * fs.NL + z] and all(np.all(np.isfinite(libs[i][r * fs.NL + z][idx])) for i in libs)]
    N = len(blocks) * p
    covs = {"data only": (Cd, fs.N_BOOT), "data + model": (Cd + mu["cov_model"], min(fs.N_BOOT, n_real)),
            "data + model (star bootstrap)": (Cd + mu["cov_model_boot"], min(fs.N_BOOT, n_real))}
    om = design.loc[ev.index, "Omega_abs"].values
    out = pd.DataFrame(index=ev.index)
    edges = np.arange(15.0, 45.01, 2.5)
    print(f"{args.case}: {len(ev)} design points, {len(blocks)} blocks x {p} terms = N {N}; C_model from {n_real} realisations")
    for tag, (C, nsamp) in covs.items():
        h = (nsamp - p - 2) / (nsamp - 1)
        W = {b: np.linalg.inv(np.linalg.cholesky(C[b][np.ix_(idx, idx)])) for b in blocks}
        chi2 = np.array([h * sum(float(np.sum((W[b] @ (cd0[b][idx] - libs[i][b][idx])) ** 2)) for b in blocks)
                         for i in ev.index])
        out[tag] = chi2
        ib = int(np.argmin(chi2))
        prof = [np.min(chi2[(om >= a) & (om < b)]) - chi2.min() if ((om >= a) & (om < b)).any() else np.nan
                for a, b in zip(edges[:-1], edges[1:])]
        print(f"  [{tag}] Hartlap {h:.3f}; best idx {ev.index[ib]} (Omega {om[ib]:.2f}); best chi2/N {chi2[ib] / N:.3f}; "
              f"median chi2/N {np.median(chi2) / N:.3f}")
        print("     Δχ² profile per 2.5 km/s/kpc, Omega 15 -> 45 (raw):", " ".join(f"{x:6.0f}" for x in prof))
    out.to_csv(RESULTS_DIR / c["run"] / "evals_modelcov.csv")


if __name__ == "__main__":
    main()

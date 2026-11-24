"""Scheme B extended from R < 6 to R < 8 kpc, re-scored from the saved library coefficients (no integration).

Why: in the Milky Way family the R < 6 kpc kinematics hardly react to the pattern speed (PROGRESS 2026-11-24), while
corotation sweeps 5.8-8.8 kpc for |Omega_b| 32.5 +-20%. The pipeline already computes the Fourier coefficients of
every block out to 8 kpc (24 rings: 20 of 0.3 kpc, then 4 of 0.5 kpc); scheme B only sums the first 20. This adds the
6-8 kpc rings with the same terms (<v_R> m2,4; <v_phi> m0,2,4; sigma_R^2 m2,4; sigma_phi^2 m2,4), the same joint
block covariance (300 bootstraps, Hartlap) and the same per-point block validity, in the data frame (library rotated
from the trial angle to 27 deg; the data fitted once at 27 deg). The R < 6 scheme B chi2 is recomputed the same way
as a control against evals.csv.

Writes to results/<run>/evals.csv the columns chi2_/ndata_/logLs_fourier_stars_R0to8_vR_sigR_sigphi_m24_vphi_m024_sfratio,
logLs tempered as the pipeline does for data without a truth row: -chi2 / (2 sqrt(2 N)), N the most common N_data.

Usage: python rescore_scheme_b_to_8kpc.py --run mw_hunter24_mock_from_osarg_ics --mock mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5
       then test1_gp_emulator.py --run <run> --column logLs_fourier_stars_R0to8_vR_sigR_sigphi_m24_vphi_m024_sfratio
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test2b_qmc_distance_error_run import load_design

B6 = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
B8 = "fourier_stars_R0to8_vR_sigR_sigphi_m24_vphi_m024"
SCHEME8 = (len(fs.RING_EDGES) - 1, fs.SCHEMES[B6][1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--mock", required=True)
    ap.add_argument("--score-dir", default=None, help="default: $MIRA_DR4_RUNS (or ~/mira_dr4_runs)/<run>")
    args = ap.parse_args()
    sdir = args.score_dir or os.path.join(os.environ.get("MIRA_DR4_RUNS", os.path.expanduser("~/mira_dr4_runs")), args.run)
    meta = json.load(open(RESULTS_DIR / args.run / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    with h5py.File(MOCK_DIR / args.mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    cd0 = cd0.reshape(fs.NB, 6, fs.NK)
    ev = pd.read_csv(RESULTS_DIR / args.run / "evals.csv")
    chi6, chi8, nd8 = [], [], []
    for i in ev["idx"]:
        cl = np.load(f"{sdir}/sfratio_{int(i):03d}.npz")["coef_lib"].reshape(fs.NB, -1)
        cl = (cl @ fs.rotation(BAR_ANGLE_DEG - float(design.loc[i, "phi_b_deg"])).T).reshape(fs.NB, 6, fs.NK)
        valid = ok & np.all(np.isfinite(cd0), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        chi6.append(fs.scheme_chi2(cd0, cl, cov0, valid, fs.SCHEMES[B6])[0])
        c, n = fs.scheme_chi2(cd0, cl, cov0, valid, SCHEME8)
        chi8.append(c)
        nd8.append(n)
    chi6, chi8, nd8 = np.array(chi6), np.array(chi8), np.array(nd8)
    saved = ev[f"chi2_{B6}_sfratio"].values
    n_mode = int(pd.Series(nd8).mode().iloc[0])
    print(f"{args.run}: {len(ev)} points; control, R < 6 scheme B chi2 recomputed vs evals.csv: max |diff| "
          f"{np.max(np.abs(chi6 - saved)):.3g} (median chi2 {np.median(saved):.0f})")
    print(f"   R < 8: N_data mode {n_mode} (range {nd8.min()}-{nd8.max()}), best chi2/N {chi8.min() / n_mode:.3f}; "
          f"R < 6: best chi2/N {saved.min() / int(ev[f'ndata_{B6}_sfratio'].mode().iloc[0]):.3f}")
    ev[f"chi2_{B8}_sfratio"] = chi8
    ev[f"ndata_{B8}_sfratio"] = nd8
    ev[f"logLs_{B8}_sfratio"] = -chi8 / (2 * np.sqrt(2 * n_mode))
    ev.to_csv(RESULTS_DIR / args.run / "evals.csv", index=False)
    print(f"   wrote chi2_/ndata_/logLs_{B8}_sfratio to {RESULTS_DIR / args.run / 'evals.csv'}")


if __name__ == "__main__":
    main()

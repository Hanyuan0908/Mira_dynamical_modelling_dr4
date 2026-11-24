"""Scheme B without the <v_phi> m = 0 term, re-scored from the saved library coefficients (no integration).

Why (user, 2026-11-24): the m = 0 term of <v_phi> is the azimuthally averaged rotation, essentially the rotation
curve and the asymmetric drift. In the mocks' own potentials it changes strongly with the pattern speed; in the Milky
Way rotation curve it hardly does (cross_response_stars_vs_potential.py). If it carries much of the mocks' Omega_b
information, the mock constraints are partly a rotation-curve measurement. This removes it:
  scheme B          <v_R> m2,4; <v_phi> m0,2,4; sigma_R^2 m2,4; sigma_phi^2 m2,4   (17 terms per block)
  without m = 0     <v_R> m2,4; <v_phi> m2,4;   sigma_R^2 m2,4; sigma_phi^2 m2,4   (16 terms per block)
Same rings (R < 6 kpc), same joint block covariance (300 bootstraps, Hartlap for 16 terms), same per-point block
validity, in the data frame (library rotated from the trial angle to 27 deg). The data fit comes from the run's own
fourier_000.npz (coef_data, cov_data_joint, n_stars >= 50), or, for a --keep sub-sample, is recomputed on the masked
mock. Control: scheme B recomputed the same way against evals.csv.

Writes chi2_/ndata_/logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m24_<variant> into results/<run>/evals.csv,
logLs = -chi2 / (2 sqrt(2 N)) with N at the truth row for a mock, else the most common N (as the pipeline).

Usage: python rescore_without_vphi_m0.py --run test2_t7_sf_12p5pct --score-dir <USB>/runs/test2_t7_sf_12p5pct
       [--variant sfratio_osarg_cover --results-run test2_t7_sf_12p5pct_osarg_cover --mock ... --keep ...]
"""
import argparse
import json

import numpy as np
import pandas as pd

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test2b_qmc_distance_error_run import load_design

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
B_NOM0 = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m24"
NR, SPEC = fs.SCHEMES[B]
SCHEME_NOM0 = (NR, {f: [k for k in terms if not (f == 1 and k == 0)] for f, terms in SPEC.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="run whose design and trajectories/scores are used")
    ap.add_argument("--score-dir", required=True)
    ap.add_argument("--variant", default="sfratio")
    ap.add_argument("--results-run", default=None, help="results/<this>/evals.csv to read and extend (default --run)")
    ap.add_argument("--mock", default=None, help="recompute the data fit on this mock (needed with --keep)")
    ap.add_argument("--keep", default=None, help="boolean mask .npy of the sub-sample")
    args = ap.parse_args()
    res_run = args.results_run or args.run
    meta = json.load(open(RESULTS_DIR / args.run / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    if args.mock:
        import h5py
        with h5py.File(MOCK_DIR / args.mock, "r") as f:
            w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
        if args.keep:
            w = w[np.load(ROOT / args.keep)]
        boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
        cd, cov, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    else:
        z = np.load(f"{args.score_dir}/fourier_000.npz")
        cd, cov, ok = z["coef_data"], z["cov_data_joint"], z["n_stars"] >= fs.N_MIN_BLOCK
    cd = cd.reshape(fs.NB, 6, fs.NK)
    ev = pd.read_csv(RESULTS_DIR / res_run / "evals.csv")
    c6, c0, n0 = [], [], []
    for i, kind in zip(ev["idx"], ev["kind"]):
        if kind not in ("qmc", "truth"):          # e.g. a corrupt row in test2_au22's evals.csv
            c6.append(np.nan); c0.append(np.nan); n0.append(0)
            continue
        cl = np.load(f"{args.score_dir}/{args.variant}_{int(i):03d}.npz")["coef_lib"].reshape(fs.NB, -1)
        cl = (cl @ fs.rotation(BAR_ANGLE_DEG - float(design.loc[i, "phi_b_deg"])).T).reshape(fs.NB, 6, fs.NK)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        c6.append(fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[B])[0])
        c, n = fs.scheme_chi2(cd, cl, cov, valid, SCHEME_NOM0)
        c0.append(c)
        n0.append(n)
    c6, c0, n0 = np.array(c6), np.array(c0), np.array(n0)
    saved = ev[f"chi2_{B}_{args.variant}"].values
    truth = ev["kind"].values == "truth"
    use = np.isfinite(c6)
    n_ref = int(n0[truth][0]) if truth.any() else int(pd.Series(n0[use]).mode().iloc[0])
    print(f"{res_run}: {len(ev)} rows ({(~use).sum()} skipped); control, scheme B recomputed vs evals.csv: max |diff| "
          f"{np.max(np.abs(c6[use] - saved[use])):.3g} "
          f"(median chi2 {np.median(saved):.0f}); without vphi m=0: N {n_ref} (scheme B {n_ref * 17 // 16}), "
          f"best chi2/N {np.nanmin(c0) / n_ref:.3f}")
    col = f"{B_NOM0}_{args.variant}"
    ev[f"chi2_{col}"], ev[f"ndata_{col}"], ev[f"logLs_{col}"] = c0, n0, -c0 / (2 * np.sqrt(2 * n_ref))
    ev.to_csv(RESULTS_DIR / res_run / "evals.csv", index=False)
    print(f"   wrote chi2_/ndata_/logLs_{col}")


if __name__ == "__main__":
    main()

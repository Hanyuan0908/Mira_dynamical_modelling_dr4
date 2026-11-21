"""One model covariance for all models, from the best-fit model (user's design, 2026-11-22).

A covariance computed separately for every model would differ between models, and χ² would then partly reward the
noisier models instead of measuring goodness of fit. Instead: take the best-fit design point, re-draw every star's
distance from its error N_REAL times, recompute its Galactocentric position and velocity, integrate in the best-fit
potential (1 Gyr, 50 Myr) and score the library exactly as the pipeline does (sfratio, K = 10, against the OBSERVED
data, which stay fixed). The scatter of the library coefficients over the realisations, per (ring, layer) block, is
C_model. It is used unchanged for every design point, rotated to each trial bar angle like the data covariance:
    χ² = Σ_blocks h · d_S^T (C_data + C_model)_SS^-1 d_S.
A second covariance, C_model_boot, is computed from the same integrations with multinomial star weights added to the
orbits (finite-sample noise, which re-drawing distances alone does not vary). It is reported, not used by default.

Distances: d_r = d_obs · 10^(0.2 σ_i N(0,1)), with σ_i = per-star DM_err where the file has it (OSARG: the paper's
M_JKs_sigma), else the file's sigma_DM. Block-diagonal: model noise correlated ACROSS blocks (orbits cross blocks)
is not captured.

Usage: python model_uncertainty_best_fit.py --case osarg_h24 [--n-real 100]
Output: results/model_uncertainty_<case>.npz; prints C_model against C_data and the re-scored χ²
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import observe
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test2b_qmc_distance_error_run import load_design
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from scan_omega_at_fixed_bar import CASES

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
SEED = 2122


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="osarg_h24", choices=list(CASES))
    ap.add_argument("--n-real", type=int, default=100)
    args = ap.parse_args()
    c = CASES[args.case]
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    ev = pd.read_csv(RESULTS_DIR / c["run"] / "evals.csv").set_index("idx")
    ev = ev[ev.kind == "qmc"]
    ibest = int(ev[f"logLs_{B}_sfratio"].idxmax())
    p = design.loc[ibest]
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        obs = {k: f[k][:].astype(np.float64) for k in ("l", "b", "d", "pml", "pmb", "vlos", "x", "y", "z", "vx", "vy", "vz")}
        sig_i = f["DM_err"][:].astype(np.float64) if "DM_err" in f else np.full(len(obs["d"]), sigma_dm)
    if np.allclose(sig_i, 0):
        sig_i = np.full(len(obs["d"]), sigma_dm)
    w_obs = np.column_stack([obs[k] for k in ("x", "y", "z", "vx", "vy", "vz")])
    n = len(w_obs)
    boot = np.random.default_rng(fs.SEED).multinomial(n, np.full(n, 1 / n), size=fs.N_BOOT).astype(float)
    cd0, cov0, _, ok = fs.data_fit(w_obs, BAR_ANGLE_DEG, boot)
    phi = float(p.phi_b_deg)
    Tinv = fs.rotation(BAR_ANGLE_DEG - phi)            # trial frame -> data frame
    print(f"{args.case}: best point idx {ibest}: Omega {p.Omega_abs:.2f}, phi {phi:.2f}, S {p.S:.3f}, f {p.bar_frac:.3f}; "
          f"{n:,} stars; per-star sigma_DM median {np.median(sig_i):.4f}; {args.n_real} realisations", flush=True)

    ssl.ACCURACY = old.ACCURACY
    pot = ssl.bar_aligned_potential(float(p.bar_frac), float(p.S), c["galaxy"])
    rng = np.random.default_rng(SEED)
    libs, libs_boot = [], []
    t0 = time.time()
    for r in range(args.n_real):
        d_r = obs["d"] * 10.0 ** (0.2 * sig_i * rng.standard_normal(n))
        w_r = np.column_stack(observe.to_galactocentric(obs["l"], obs["b"], d_r, obs["pml"], obs["pmb"], obs["vlos"]))
        _, traj = ssl.integrate(pot, ssl.rotate_z(w_r, np.radians(phi)), -abs(float(p.Omega_abs)))
        traj = np.ascontiguousarray(traj[:, ::2])
        cl, _ = library_fit_empirical(traj, phi, sigma_dm, dat, ok, seed=11 + r)
        cb, _ = library_fit_empirical(traj, phi, sigma_dm, dat, ok, seed=11 + r,
                                      star_w=rng.multinomial(n, np.full(n, 1 / n)).astype(float))
        libs.append((cl.reshape(fs.NB, -1) @ Tinv.T))
        libs_boot.append((cb.reshape(fs.NB, -1) @ Tinv.T))
        if r % 10 == 0:
            print(f"  realisation {r + 1}/{args.n_real}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
    libs, libs_boot = np.array(libs), np.array(libs_boot)
    covm = {}
    for tag, L in (("model", libs), ("model_boot", libs_boot)):
        dev = L - np.nanmean(L, 0)
        covm[tag] = np.einsum("rbi,rbj->bij", dev, dev) / (len(L) - 1)
    np.savez(RESULTS_DIR / f"model_uncertainty_{args.case}.npz", coef_lib_real=libs, coef_lib_real_boot=libs_boot,
             cov_model=covm["model"], cov_model_boot=covm["model_boot"], cov_data=cov0, coef_data=cd0, ok=ok,
             best=np.array([p.Omega_abs, phi, p.S, p.bar_frac]), ibest=ibest, n_real=args.n_real)

    # size of the model noise against the data noise, scheme B terms, per field
    nr, spec = fs.SCHEMES[B]
    names = {0: "<v_R> m2,4", 1: "<v_phi> m0,2,4", 3: "sig_R^2 m2,4", 4: "sig_phi^2 m2,4"}
    blocks = [r * fs.NL + z for r in range(nr) for z in range(fs.NL) if ok[r * fs.NL + z]]
    print("median over usable blocks of sqrt(C_model / C_data) on the diagonal:")
    for f, terms in spec.items():
        ii = [f * fs.NK + k for k in terms]
        rat = [np.sqrt(np.diag(covm["model"][b])[ii] / np.diag(cov0[b])[ii]) for b in blocks]
        ratb = [np.sqrt(np.diag(covm["model_boot"][b])[ii] / np.diag(cov0[b])[ii]) for b in blocks]
        print(f"   {names[f]:16s} distance re-draws {np.nanmedian(rat):.2f}   + star bootstrap {np.nanmedian(ratb):.2f}")


if __name__ == "__main__":
    main()

"""Emulate the library Fourier coefficients instead of logL, then evaluate chi2 exactly.

`test1_gp_emulator.py` fits a GP to logL itself. That is the hard way round: logL is a quadratic form in
~750 coefficients, spans thousands in range, and the GP has to represent all of that from 256 points - which
is why its leave-one-out error (about 1 in tempered logL units) is twice the 0.5 that defines the 1 sigma
contour.

The alternative, and what the cosmology emulators do (Coyote/CosmicEmu, the Lyman-alpha emulators): emulate the
*data vector* - here the library's Fourier coefficients c_lib(theta) - and then compute chi2 exactly, with the
real covariance. The coefficients are smooth, low-dynamic-range functions of theta, so they are far easier to
interpolate, and the quadratic structure of chi2 is then exact rather than learned.

Many correlated outputs are handled the standard way: PCA the coefficient vectors over the design, keep the
leading components, and fit one GP per component. ~750 outputs collapse to ~10-20.

The data side needs no emulation: at a trial angle phi_b the data coefficients and their covariance follow
exactly from the fit at 27 deg by the rotation c' = T c, C' = T C T^T (`fs.rotation`).

Validation is like-for-like with the logL emulator: hold out points, predict their chi2 through the emulated
coefficients, and compare with the chi2 actually computed from the saved library.

Usage: python emulate_fourier_coefficients.py --run test2_t7_sf_10pct --variant sf [--n-pc 16]
"""
import argparse
import json

import numpy as np
import pandas as pd
import h5py
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from common import ROOT, MOCK_DIR, RUNS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test1_nbody_qmc_focused_fourier_stars import rotation
from test2b_qmc_distance_error_run import scaled_design

SCHEMES = {"A": "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024",
           "B": "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"}


def scheme_blocks(scheme, valid):
    """(block index, coefficient indices) actually entering chi2, in scheme_chi2's own order."""
    n_ring, spec = fs.SCHEMES[scheme]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    return [(r * fs.NL + z, idx) for r in range(n_ring) for z in range(fs.NL) if valid[r * fs.NL + z]]


def chi2_from_coef(cl_flat, blocks, cd, cov):
    """Exactly fs.scheme_chi2, but taking the library coefficients as a flat vector."""
    chi, nd, pos = 0.0, 0, 0
    for i, idx in blocks:
        c = cl_flat[pos:pos + len(idx)]
        pos += len(idx)
        d = cd[i].reshape(-1)[idx] - c
        h = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)
        chi += h * float(d @ np.linalg.solve(cov[i][np.ix_(idx, idx)], d))
        nd += len(idx)
    return chi, nd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--variant", default="sf")
    ap.add_argument("--scheme", default="B", choices=list(SCHEMES))
    ap.add_argument("--n-pc", type=int, default=16)
    ap.add_argument("--n-test", type=int, default=40, help="held-out points for validation")
    args = ap.parse_args()
    sch = SCHEMES[args.scheme]

    meta = json.load(open(ROOT / "results" / args.run / "run_meta.json"))
    design, _ = scaled_design(meta["galaxy"])
    out_dir = RUNS_DIR / args.run
    pref = {"sf": "fourier", "asis": "fourier", "observed": "fourier",
            "sfratio": "sfratio", "sfemp": "sfemp"}[args.variant]

    s0 = np.load(out_dir / "fourier_000.npz")
    cd0 = s0["coef_data"].reshape(fs.NB, 6, fs.NK)
    cov0 = s0["cov_data_joint"].reshape(fs.NB, 6 * fs.NK, 6 * fs.NK)
    ok = s0["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK
    key = f"coef_lib_{args.variant}" if pref == "fourier" else "coef_lib"

    X, Y, chi_true, kinds = [], [], [], []
    for i in design.idx:
        i = int(i)
        f = out_dir / (f"fourier_{i:03d}.npz" if pref == "fourier" else f"{pref}_{i:03d}.npz")
        if not f.exists():
            continue
        s = np.load(f)
        cl = s[key].reshape(fs.NB, 6, fs.NK)
        valid = ok & np.all(np.isfinite(cl), axis=(1, 2))
        blocks = scheme_blocks(sch, valid)
        p = design.set_index("idx").loc[i]
        T = rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        cd = (cd0.reshape(fs.NB, -1) @ T.T).reshape(fs.NB, 6, fs.NK)
        cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
        # whiten by the data covariance: chi2 = h * |e|^2 with e = L^-1 (cd - cl), L the block Cholesky.
        # PCA truncation error in e then maps straight onto chi2 error, which raw coefficients do not.
        ev = []
        for b, idx in blocks:
            d = cd[b].reshape(-1)[idx] - cl[b].reshape(-1)[idx]
            L = np.linalg.cholesky(cov[b][np.ix_(idx, idx)])
            ev.append(np.linalg.solve(L, d))
        vec = np.concatenate(ev)
        h = (fs.N_BOOT - len(blocks[0][1]) - 2) / (fs.N_BOOT - 1)
        c2, nd = h * float(vec @ vec), sum(len(i) for _, i in blocks)
        X.append([p.Omega_abs, p.phi_b_deg, np.log(p.S), np.log(p.bar_frac)])
        Y.append(vec); chi_true.append(c2); kinds.append(p.kind)
        if len(X) == 1:
            blocks0, nd0 = blocks, nd
    X = np.array(X); Y = np.array(Y); chi_true = np.array(chi_true); kinds = np.array(kinds)
    print(f"{args.run} / {args.variant} / scheme {args.scheme}: {len(X)} points, "
          f"{Y.shape[1]} coefficients (N_data {nd0})")

    lo, hi = X.min(0), X.max(0)
    U = (X - lo) / (hi - lo)
    mu, sd = Y.mean(0), Y.std(0) + 1e-12
    Z = (Y - mu) / sd
    Uc, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    npc = min(args.n_pc, len(S))
    W = Uc[:, :npc] * S[:npc]
    print(f"  PCA: {npc} components explain {100 * (S[:npc]**2).sum() / (S**2).sum():.2f}% of the variance")

    rng = np.random.default_rng(0)
    qmc = np.flatnonzero(kinds == "qmc")
    test = rng.choice(qmc, min(args.n_test, len(qmc)), replace=False)
    train = np.setdiff1d(np.arange(len(X)), test)

    def fit_predict(tr, te):
        pred = np.zeros((len(te), npc))
        for k in range(npc):
            y = W[tr, k]
            kern = (ConstantKernel(np.var(y), (1e-3 * np.var(y) + 1e-12, 1e3 * np.var(y) + 1)) *
                    Matern(length_scale=[0.3] * 4, length_scale_bounds=(0.02, 10), nu=2.5) +
                    WhiteKernel(1e-3 * np.var(y) + 1e-12, (1e-12, 1e2)))
            g = GaussianProcessRegressor(kern, normalize_y=True, n_restarts_optimizer=2,
                                         random_state=0).fit(U[tr], y)
            pred[:, k] = g.predict(U[te])
        return (pred @ Vt[:npc] + Z.mean(0)) * sd + mu

    Yp = fit_predict(train, test)
    h0 = (fs.N_BOOT - len(blocks0[0][1]) - 2) / (fs.N_BOOT - 1)
    chi_pred = h0 * np.einsum("ij,ij->i", Yp, Yp)
    temper = 2 * np.sqrt(2 * nd0)
    err = (chi_true[test] - chi_pred) / temper
    print(f"  held-out {len(test)} points: logLs error rms {np.sqrt(np.mean(err**2)):.3f}, "
          f"median |err| {np.median(np.abs(err)):.3f}")
    j = json.load(open(ROOT / "results" / args.run / f"emulator_results_logLs_{sch}_{args.variant}.json"))
    print(f"  direct logL emulator, same run: LOO rms {j['loo_rms']:.3f}")
    print(f"  -> coefficient emulator is {j['loo_rms'] / np.sqrt(np.mean(err**2)):.1f}x more accurate "
          f"(1 sigma of the tempered logL = 0.5)")


if __name__ == "__main__":
    main()

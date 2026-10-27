"""Test 1, N-body, focused global search: Fourier moments fitted directly to the individual stars and library samples
(no azimuthal sectors).

Grid: R rings 0-6 kpc in 20 rings of 0.3 kpc (requested), extended by 4 rings of 0.5 kpc to 8 kpc; |z| layers 0-1, 1-2,
2-3 kpc (folded: |z|, v_z sign z). In each (ring, layer) block, for v = v_R, v_phi, v_z sign z:
  stage 1: least squares over the stars of v_i = sum_k a_k f_k(phi_i), f = 1, cos 2phi, sin 2phi, cos 4phi, sin 4phi,
           giving the Fourier coefficients of the mean field <v>(phi). This does not depend on how the stars are spread
           in phi as long as these modes capture the field. (A plain sum sum_i v_i cos m phi_i / N would instead give
           the density-weighted field n <v>, mixing in the bar's density and, later, the selection function.)
  stage 2: least squares of r_i^2 / (1 - h_ii) on the same basis (r_i = residual of stage 1, h_ii = leverage; for a
           constant-only basis this is the N / (N - 1) correction), giving the coefficients of sigma^2(phi).
The normal equations are sums over stars (weights: 1 for the data, bootstrap counts for resamplings, 1/tau for library
samples). Library: all samples with t >= 0.25 Gyr from the saved trajectories (no leverage correction needed: ~10^3
effective orbits per block). Data errors: full 5x5 covariance of each block's coefficients from N_BOOT bootstrap
resamplings of the stars. A block needs >= N_MIN_BLOCK data stars.
Fields (index): 0 <v_R>, 1 <v_phi>, 2 <v_z>, 3 sigma_R^2, 4 sigma_phi^2, 5 sigma_z^2.
Saved per point: ~/mira_dr4_runs/test1_nbody_qmc_focused/fourier_stars_NNN.npz (coef_data, coef_lib, cov_data, valid,
n_stars). Columns chi2_/ndata_/logLs_fourier_stars_<scheme> are merged into results/test1_nbody_qmc_focused/evals.csv;
chi2 over kept terms = dc_S^T (C_SS)^-1 dc_S (dropping a term = projecting it out of the fit); tempered
logLs = -chi2 / (2 sqrt(2 N_data)), N_data = number of terms used at the truth.
"""
import json

import numpy as np
import pandas as pd
import h5py

import steady_state_likelihood as ssl
from test1_nbody_qmc_focused_run import OUT_DIR, DESIGN_CSV
import test1_nbody_qmc as old

RING_EDGES = np.round(np.r_[np.arange(0.0, 6.0 + 1e-9, 0.3), [6.5, 7.0, 7.5, 8.0]], 6)
N_RING_DEFAULT = 20
LAYER_EDGES = np.array([0.0, 1.0, 2.0, 3.0])
NR, NL = len(RING_EDGES) - 1, len(LAYER_EDGES) - 1
NB = NR * NL
TERMS = ["m0", "cos2", "sin2", "cos4", "sin4"]
FIELD_NAMES = ["<v_R>", "<v_phi>", "<v_z>", "sigma_R^2", "sigma_phi^2", "sigma_z^2"]
NK = len(TERMS)
PAIRS = [(k, l) for k in range(NK) for l in range(k, NK)]
N_MIN_BLOCK = 50
N_BOOT = 200
SEED = 0
ALL, NO_M0 = [0, 1, 2, 3, 4], [1, 2, 3, 4]
SCHEMES = {   # name: (rings used, {field: kept terms})
    "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024": (N_RING_DEFAULT, {0: NO_M0, 1: ALL, 3: ALL, 4: ALL}),
    "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024": (N_RING_DEFAULT, {0: NO_M0, 1: ALL, 3: NO_M0, 4: NO_M0}),
    "fourier_stars_R0to6_all6_m024": (N_RING_DEFAULT, {f: ALL for f in range(6)}),
    "fourier_stars_R0to8_vR_m24_vphi_sigR_sigphi_m024": (NR, {0: NO_M0, 1: ALL, 3: ALL, 4: ALL}),
}


def prepare(w):
    """Block index, basis (n, 5), velocities (v_R, v_phi, v_z sign z) and the in-grid mask for bar-frame points."""
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.arctan2(w[:, 1], w[:, 0])                  # even m only, so no folding in phi is needed
    ring = np.digitize(R, RING_EDGES) - 1
    layer = np.digitize(np.abs(w[:, 2]), LAYER_EDGES) - 1
    m = (ring >= 0) & (ring < NR) & (layer >= 0) & (layer < NL)
    R, phi, w = R[m], phi[m], w[m]
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5] * np.sign(w[:, 2])])
    X = np.column_stack([np.ones_like(phi), np.cos(2 * phi), np.sin(2 * phi), np.cos(4 * phi), np.sin(4 * phi)])
    return ring[m] * NL + layer[m], X, v, m


def normal_eq(b, X, wt, Y):
    """Per block: G = X^T W X (NB, 5, 5) and H = X^T W Y (NB, ny, 5)."""
    G = np.zeros((NB, NK, NK))
    for k, l in PAIRS:
        G[:, k, l] = G[:, l, k] = np.bincount(b, weights=wt * X[:, k] * X[:, l], minlength=NB)
    H = np.stack([np.stack([np.bincount(b, weights=wt * X[:, k] * Y[:, a], minlength=NB) for k in range(NK)], -1)
                  for a in range(Y.shape[1])], 1)
    return G, H


def fit(b, X, v, wt, ok, leverage=True):
    """Two-stage least squares; coefficients (NB, 6, 5): means of v (3 fields), then variances (3); NaN where not ok."""
    G, H = normal_eq(b, X, wt, v)
    Ginv = np.full_like(G, np.nan)
    Ginv[ok] = np.linalg.inv(G[ok])
    beta = np.einsum("bkl,bal->bak", Ginv, H)
    r2 = (v - np.einsum("ik,iak->ia", X, beta[b]))**2
    if leverage:
        h = wt * np.einsum("ik,ikl,il->i", X, Ginv[b], X)
        r2 = r2 / np.clip(1 - h, 1e-3, None)[:, None]
    _, H2 = normal_eq(b, X, wt, np.nan_to_num(r2))
    gamma = np.einsum("bkl,bal->bak", Ginv, H2)
    return np.concatenate([beta, gamma], axis=1)


def data_fit(stars_mw, phi_deg, boot_w):
    """Data coefficients, bootstrap covariance (NB, 6, 5, 5), star counts and usable blocks at a trial bar angle."""
    b, X, v, m = prepare(ssl.rotate_z(stars_mw, np.radians(phi_deg)))
    n = np.bincount(b, minlength=NB)
    ok = n >= N_MIN_BLOCK
    c = fit(b, X, v, np.ones(len(b)), ok)
    boots = np.array([fit(b, X, v, bw[m], ok) for bw in boot_w])
    dev = boots - boots.mean(0)
    cov = np.einsum("sbfk,sbfl->bfkl", dev, dev) / (len(boot_w) - 1)
    return c, cov, n, ok


def library_fit(traj, times, ok):
    R = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])
    tau = np.clip(np.mean(R < ssl.R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
    keep = times >= ssl.T_MIN_GYR - 1e-9
    w = traj[:, keep].reshape(-1, 6).astype(np.float64)
    wt = np.repeat(1.0 / tau, keep.sum())
    b, X, v, m = prepare(w)
    return fit(b, X, v, wt[m], ok, leverage=False)


def scheme_chi2(cd, cl, cov, valid, scheme):
    n_ring, spec = scheme
    chi, nd = 0.0, 0
    for r in range(n_ring):
        for z in range(NL):
            i = r * NL + z
            if not valid[i]:
                continue
            for f, terms in spec.items():
                d = (cd[i, f] - cl[i, f])[terms]
                chi += float(d @ np.linalg.solve(cov[i, f][np.ix_(terms, terms)], d))
                nd += len(terms)
    return chi, nd


def main():
    design = pd.read_csv(DESIGN_CSV).sort_values("idx").reset_index(drop=True)
    with h5py.File(old.MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    rng = np.random.default_rng(SEED)
    boot_w = rng.multinomial(len(stars), np.full(len(stars), 1 / len(stars)), size=N_BOOT).astype(float)
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    rows, cache = [], {}
    for _, p in design.iterrows():
        i = int(p.idx)
        phi = round(float(p.phi_b_deg), 6)
        if phi not in cache:
            cache = {phi: data_fit(stars, phi, boot_w)}
        cd, cov, n, ok = cache[phi]
        cl = library_fit(np.load(p.traj_file), times, ok)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        np.savez(OUT_DIR / f"fourier_stars_{i:03d}.npz", ring_edges=RING_EDGES, layer_edges=LAYER_EDGES,
                 terms=TERMS, fields=FIELD_NAMES, coef_data=cd.reshape(NR, NL, 6, NK), coef_lib=cl.reshape(NR, NL, 6, NK),
                 cov_data=cov.reshape(NR, NL, 6, NK, NK), valid=valid.reshape(NR, NL), n_stars=n.reshape(NR, NL))
        row = dict(idx=i)
        for name, sc in SCHEMES.items():
            row[f"chi2_{name}"], row[f"ndata_{name}"] = scheme_chi2(cd, cl, cov, valid, sc)
        rows.append(row)
        print(f"{i:3d} {p.batch} {p.kind:5s} " + "  ".join(f"{k[19:]} {v:.0f}" for k, v in row.items() if k.startswith("chi2_")),
              flush=True)
    new = pd.DataFrame(rows)
    out_csv = DESIGN_CSV.parent / "evals.csv"
    base = pd.read_csv(out_csv) if out_csv.exists() else design[["idx", "kind", "batch", "Omega_abs", "phi_b_deg", "S", "bar_frac"]]
    base = base.drop(columns=[c for c in base.columns if "fourier_stars_" in c]).merge(new, on="idx")
    t = base.kind == "truth"
    summary = {}
    for name in SCHEMES:
        nd = int(base.loc[t, f"ndata_{name}"].iloc[0])
        base[f"logLs_{name}"] = -0.5 * base[f"chi2_{name}"] / np.sqrt(2 * nd)
        dl = (base[f"logLs_{name}"] - base.loc[t, f"logLs_{name}"].iloc[0])[~t]
        summary[name] = dict(N_data=nd, truth_chi2_per_N=float(base.loc[t, f"chi2_{name}"].iloc[0] / nd),
                             truth_rank=int((dl > 0).sum()) + 1, best_other=float(dl.max()),
                             best_other_idx=int(base.idx[dl.idxmax()]), n_within_2=int((dl > -2).sum()),
                             n_within_4p5=int((dl > -4.5).sum()))
    base.to_csv(out_csv, index=False)
    with open(DESIGN_CSV.parent / "summary_fourier_stars.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    for k, v in summary.items():
        print(k, v)


if __name__ == "__main__":
    main()

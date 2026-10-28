"""Test 1, N-body, focused global search: Fourier moment likelihood on an R-|z| grid (all folded to |z| >= 0, phi mod pi).

Grid: R rings 0-6 kpc in 20 rings of 0.3 kpc (the requested grid), extended by 4 rings of 0.5 kpc to 8 kpc; |z| layers
0-1, 1-2, 2-3 kpc. Azimuthal sectors per (ring, layer), fixed once from the data at the true bar angle: 12 sectors of
15 deg if >= 6 of them hold >= N_MIN stars, else 6 of 30 deg (same rule), else 4 of 45 deg (then only m = 0, 2 are fitted).
The |z| 2-3 layer (820 stars) is too sparse for 0.3 kpc rings and drops out. All edges lie on the fine grid of
test1_nbody_qmc_focused_run.py, so the library sector sums are exact.
Per (ring, layer) and field (<v_R>, <v_phi>, <v_z sign z>, sigma_R, sigma_phi, sigma_z): generalised least-squares fit
of the sector values with the data's bootstrap errors (sectors with >= N_MIN stars; at least one more than the number
of terms) on TERMS = 1, cos 2phi, sin 2phi, cos 4phi, sin 4phi (45-deg sectors: first three); the same operator gives the library coefficients. Dispersions are
unbiased (data N/(N-1); library W^2/(W^2 - sum w^2)).
Saved per point, ~/mira_dr4_runs/test1_nbody_qmc_focused/fourier_NNN.npz: coef_data, coef_lib (ring, layer, field, term),
cov_data (ring, layer, field, term, term), valid (ring, layer, field), n_sectors. Any subset of fields / terms can be
scored later exactly: chi2 = dc_S^T (C_SS)^-1 dc_S, which equals projecting the dropped terms out of the GLS fit.
Table: results/test1_nbody_qmc_focused/evals.csv with chi2_<scheme>, tempered logLs_<scheme> = -chi2 / (2 sqrt(2 N_data))
(N_data = number of coefficients used at the truth), and the current sector likelihood (sectors_means,
sectors_meansdisp; test1_nbody_moments bootstrap errors) for comparison.
"""
import json

import numpy as np
import pandas as pd
import h5py

import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, chi2 as sector_chi2, N_MIN, SEED
from test1_nbody_qmc_focused_run import (OUT_DIR, DESIGN_CSV, FINE_R, FINE_PHI, FINE_Z, FINE_SHAPE, TRUTH, fine_sums)
import test1_nbody_qmc as old

RING_EDGES = np.round(np.r_[np.arange(0.0, 6.0 + 1e-9, 0.3), [6.5, 7.0, 7.5, 8.0]], 6)
N_RING_DEFAULT = 20                                   # rings with R < 6 kpc
LAYER_EDGES = np.array([0.0, 1.0, 2.0, 3.0])
N_SEC = 12
TERMS = ["m0", "cos2", "sin2", "cos4", "sin4"]
FIELD_NAMES = ["<v_R>", "<v_phi>", "<v_z>", "sigma_R", "sigma_phi", "sigma_z"]
N_BOOT = 200
NR, NL = len(RING_EDGES) - 1, len(LAYER_EDGES) - 1
ALL = [0, 1, 2, 3, 4]
NO_M0 = [1, 2, 3, 4]
# scheme: (number of rings used, {field: term indices kept})
SCHEMES = {
    "fourier_R0to6_vR_m24_vphi_sigR_sigphi_m024": (N_RING_DEFAULT, {0: NO_M0, 1: ALL, 3: ALL, 4: ALL}),
    "fourier_R0to6_vR_sigR_sigphi_m24_vphi_m024": (N_RING_DEFAULT, {0: NO_M0, 1: ALL, 3: NO_M0, 4: NO_M0}),
    "fourier_R0to6_all6_m024": (N_RING_DEFAULT, {f: ALL for f in range(6)}),
    "fourier_R0to8_vR_m24_vphi_sigR_sigphi_m024": (NR, {0: NO_M0, 1: ALL, 3: ALL, 4: ALL}),
}


def target_map():
    """Fine (R, phi, |z|) cell -> target (ring, layer, sector) index, -1 outside the grid."""
    rc = 0.5 * (FINE_R[1:] + FINE_R[:-1])
    zc = 0.5 * (FINE_Z[:-1] + np.r_[FINE_Z[1:-1], FINE_Z[-2] + 1.0])
    ring = np.digitize(rc, RING_EDGES) - 1
    layer = np.digitize(zc, LAYER_EDGES) - 1
    sec = np.arange(FINE_SHAPE[1]) // 2
    I, J, K = np.meshgrid(ring, sec, layer, indexing="ij")
    ok = (I >= 0) & (I < NR) & (K >= 0) & (K < NL)
    t = np.full(FINE_SHAPE, -1)
    t[ok] = np.ravel_multi_index((I[ok], K[ok], J[ok]), (NR, NL, N_SEC))
    return t.ravel()


TMAP = target_map()
NT = NR * NL * N_SEC


def sector_stats(W, W2, S1, S2d):
    """Means and unbiased dispersions (last axis: <v_R>, <v_phi>, <v_z>, sigma_R, sigma_phi, sigma_z)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        m = S1 / W[..., None]
        corr = W**2 / (W**2 - W2)
        s = np.sqrt(np.maximum((S2d / W[..., None] - m**2) * corr[..., None], 0))
    return np.concatenate([m, s], axis=-1)


NSEC_OPTIONS = (12, 6, 4)
CENTRES = {n: np.radians((180.0 / n) * (np.arange(n) + 0.5)) for n in NSEC_OPTIONS}
NTERMS = {12: 5, 6: 5, 4: 3}


def merge(a, n):
    """Sum 15-deg target sums (..., NT, k) into n sectors per (ring, layer): shape (..., NR, NL, n, k)."""
    a = a.reshape(a.shape[:-2] + (NR, NL, N_SEC) + a.shape[-1:])
    return a.reshape(a.shape[:-2] + (n, N_SEC // n) + a.shape[-1:]).sum(-2)


def library_sectors(npz):
    """Library sector values (means, unbiased dispersions) for every sector option: {n: (NR, NL, n, 6)}."""
    ok = TMAP >= 0
    t = TMAP[ok]
    cols = [npz["fine_W"].ravel()[ok], npz["fine_W2"].ravel()[ok]]
    cols += [npz["fine_S1"].reshape(len(TMAP), 3)[ok, a] for a in range(3)]
    cols += [npz["fine_S2"].reshape(len(TMAP), 6)[ok, a] for a in range(3)]
    sums = np.column_stack([np.bincount(t, weights=c, minlength=NT) for c in cols])     # (NT, 8)
    out = {}
    for n in NSEC_OPTIONS:
        m = merge(sums, n)
        out[n] = sector_stats(m[..., 0], m[..., 1], m[..., 2:5], m[..., 5:8])
    return out


def star_targets(stars_mw, phi_deg):
    w = ssl.rotate_z(stars_mw, np.radians(phi_deg))
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.mod(np.arctan2(w[:, 1], w[:, 0]), np.pi)
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5] * np.sign(w[:, 2])])
    i = np.digitize(R, FINE_R) - 1
    j = np.minimum(np.digitize(phi, FINE_PHI) - 1, FINE_SHAPE[1] - 1)
    k = np.digitize(np.abs(w[:, 2]), FINE_Z) - 1
    ok = (i >= 0) & (i < FINE_SHAPE[0]) & (k >= 0) & (k < FINE_SHAPE[2])
    t = np.full(len(w), -1)
    t[ok] = TMAP[np.ravel_multi_index((i[ok], j[ok], k[ok]), FINE_SHAPE)]
    return t, v


def data_sums(t, v, weights):
    ok = t >= 0
    tt, w, vv = t[ok], weights[ok], v[ok]
    return np.column_stack([np.bincount(tt, weights=w, minlength=NT)]
                           + [np.bincount(tt, weights=w * vv[:, a], minlength=NT) for a in range(3)]
                           + [np.bincount(tt, weights=w * vv[:, a]**2, minlength=NT) for a in range(3)])   # (NT, 7)


def data_with_errors(stars_mw, phi_deg, boot_w):
    """Per sector option n: star counts, values and bootstrap errors, each (NR, NL, n[, 6]).
    Frequency weights (bootstrap counts): unbiased variance N / (N - 1)."""
    t, v = star_targets(stars_mw, phi_deg)
    main = data_sums(t, v, np.ones(len(v)))
    boots = np.array([data_sums(t, v, w) for w in boot_w])
    out = {}
    for n in NSEC_OPTIONS:
        m, b = merge(main, n), merge(boots, n)
        val = sector_stats(m[..., 0], m[..., 0], m[..., 1:4], m[..., 4:7])
        bv = sector_stats(b[..., 0], b[..., 0], b[..., 1:4], b[..., 4:7])
        out[n] = (m[..., 0], val, np.nanstd(bv, axis=0))
    return out


def choose_sectors(data_n):
    """Sector option per (ring, layer) from the data counts at the true angle (finest with enough usable sectors)."""
    nsec = np.zeros((NR, NL), int)
    for r in range(NR):
        for z in range(NL):
            for n in NSEC_OPTIONS:
                if (data_n[n][0][r, z] >= N_MIN).sum() >= NTERMS[n] + 1:
                    nsec[r, z] = n
                    break
    return nsec


def fourier(data, lib, nsec):
    """GLS coefficients of data and library per (ring, layer, field) with the data's errors; terms beyond the block's
    sector option (m = 4 for 45-deg sectors) are NaN."""
    cd = np.full((NR, NL, 6, 5), np.nan)
    cl = np.full_like(cd, np.nan)
    cov = np.full((NR, NL, 6, 5, 5), np.nan)
    valid = np.zeros((NR, NL, 6), bool)
    nused = np.zeros((NR, NL), int)
    for r in range(NR):
        for z in range(NL):
            n = nsec[r, z]
            if n == 0:
                continue
            N, d_val, d_err = data[n][0][r, z], data[n][1][r, z], data[n][2][r, z]
            l_val = lib[n][r, z]
            use = (N >= N_MIN) & np.all(np.isfinite(l_val), 1) & np.all(d_err > 0, 1)
            nused[r, z] = int(use.sum())
            k = NTERMS[n]
            if use.sum() < k + 1:
                continue
            ph = CENTRES[n][use]
            F = np.column_stack([np.ones_like(ph), np.cos(2 * ph), np.sin(2 * ph), np.cos(4 * ph), np.sin(4 * ph)])[:, :k]
            for f in range(6):
                w = 1.0 / d_err[use, f]**2
                C = np.linalg.inv(F.T @ (w[:, None] * F))
                cd[r, z, f, :k] = C @ (F.T @ (w * d_val[use, f]))
                cl[r, z, f, :k] = C @ (F.T @ (w * l_val[use, f]))
                cov[r, z, f, :k, :k] = C
                valid[r, z, f] = True
    return dict(coef_data=cd, coef_lib=cl, cov_data=cov, valid=valid, n_sectors_used=nused, n_sectors=nsec)


def scheme_chi2(fr, scheme):
    n_ring, spec = scheme
    chi, nd = 0.0, 0
    for r in range(n_ring):
        for z in range(NL):
            for f, terms in spec.items():
                if not fr["valid"][r, z, f]:
                    continue
                tt = [k for k in terms if np.isfinite(fr["coef_data"][r, z, f, k])]
                d = (fr["coef_data"][r, z, f] - fr["coef_lib"][r, z, f])[tt]
                C = fr["cov_data"][r, z, f][np.ix_(tt, tt)]
                chi += float(d @ np.linalg.solve(C, d))
                nd += len(tt)
    return chi, nd


def main():
    design = pd.read_csv(DESIGN_CSV).sort_values("idx").reset_index(drop=True)
    with h5py.File(old.MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    rng = np.random.default_rng(SEED)
    boot_w = rng.multinomial(len(stars), np.full(len(stars), 1 / len(stars)), size=N_BOOT).astype(float)
    rows = []
    nsec = choose_sectors(data_with_errors(stars, TRUTH["phi_b_deg"], boot_w[:2]))
    print("sectors per (ring, layer) [rows: |z| 0-1, 1-2, 2-3]:\n", nsec.T)
    cache = {}
    out_csv = DESIGN_CSV.parent / "evals.csv"
    prev = pd.read_csv(out_csv).set_index("idx") if out_csv.exists() else None
    for _, p in design.iterrows():
        i = int(p.idx)
        cols = [f"{a}_{n}" for n in list(SCHEMES) + ["sectors_means", "sectors_meansdisp"] for a in ("chi2", "ndata")]
        if prev is not None and i in prev.index and (OUT_DIR / f"fourier_{i:03d}.npz").exists() \
                and all(c in prev.columns for c in cols) and prev.loc[i, cols].notna().all():
            # re-use: the design is fixed by its seed and the evaluation is deterministic (seeded bootstraps)
            rows.append({**{k: p[k] for k in ("idx", "kind", "batch", "Omega_abs", "phi_b_deg", "S", "bar_frac")},
                         **prev.loc[i, cols].to_dict()})
            continue
        lib = np.load(OUT_DIR / f"lib_{i:03d}.npz")
        phi = round(float(p.phi_b_deg), 6)
        if phi not in cache:
            cache = {phi: (data_with_errors(stars, phi, boot_w),
                           boot_errors(ssl.rotate_z(stars, np.radians(phi)), np.random.default_rng(SEED)))}
        data, obs_sec = cache[phi]
        fr = fourier(data, library_sectors(lib), nsec)
        np.savez(OUT_DIR / f"fourier_{i:03d}.npz", ring_edges=RING_EDGES, layer_edges=LAYER_EDGES, terms=TERMS,
                 fields=FIELD_NAMES, **fr)
        row = {k: p[k] for k in ("idx", "kind", "batch", "Omega_abs", "phi_b_deg", "S", "bar_frac")}
        for name, sc in SCHEMES.items():
            row[f"chi2_{name}"], row[f"ndata_{name}"] = scheme_chi2(fr, sc)
        row["chi2_sectors_means"], row["ndata_sectors_means"] = sector_chi2(obs_sec, lib["lib_m1"], lib["lib_s"], 0.0, False)
        row["chi2_sectors_meansdisp"], row["ndata_sectors_meansdisp"] = sector_chi2(obs_sec, lib["lib_m1"], lib["lib_s"],
                                                                                   0.0, True)
        row["ndata_sectors_means"] *= 3
        row["ndata_sectors_meansdisp"] *= 6
        rows.append(row)
        print(f"{i:3d} {p.batch} {p.kind:5s} " + "  ".join(f"{k[5:]} {v:.0f}" for k, v in row.items() if k.startswith("chi2_")),
              flush=True)
    df = pd.DataFrame(rows)
    t = df.kind == "truth"
    summary = {}
    for name in list(SCHEMES) + ["sectors_means", "sectors_meansdisp"]:
        nd = int(df.loc[t, f"ndata_{name}"].iloc[0])
        df[f"logLs_{name}"] = -0.5 * df[f"chi2_{name}"] / np.sqrt(2 * nd)
        dl = (df[f"logLs_{name}"] - df.loc[t, f"logLs_{name}"].iloc[0])[~t]
        summary[name] = dict(N_data=nd, sqrt_2N=float(np.sqrt(2 * nd)), truth_chi2_per_N=float(df.loc[t, f"chi2_{name}"].iloc[0] / nd),
                             truth_rank=int((dl > 0).sum()) + 1, best_other=float(dl.max()), best_other_idx=int(df.idx[dl.idxmax()]),
                             n_within_2=int((dl > -2).sum()), n_within_4p5=int((dl > -4.5).sum()))
    df.to_csv(out_csv, index=False)
    with open(DESIGN_CSV.parent / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    for k, v in summary.items():
        print(k, v)


if __name__ == "__main__":
    main()

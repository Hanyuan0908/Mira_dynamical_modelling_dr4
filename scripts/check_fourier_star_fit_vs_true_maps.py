"""Check of the per-star Fourier fit (test1_nbody_qmc_focused_fourier_stars.py) against the ground truth.

Ground truth: all star particles of the N-body snapshot with R < 8 kpc (1,493,610; the pool the 50k mock was drawn from
at random), in the bar frame at the true angle (27 deg), number-weighted like the mock. Three checks:
  1. Reconstruction: the m <= 4 maps of <v_R>, <v_phi>, sigma_R, sigma_phi rebuilt from the 50k-mock coefficients (the
     ones in the likelihood; truth point 000) against the true maps, in units of km/s and of the combined error.
  2. Truncation: the same fit applied to the whole pool against the pool's own binned maps (can m <= 4 per 0.3 kpc ring
     and |z| layer represent the true field at all?).
  3. Coefficients: chi2 of (mock - pool) coefficients with the mock's joint bootstrap covariance and the Hartlap factor,
     exactly as in the likelihood. If the estimator and its errors are right this is (1 - 50000 / 1493610) N = 0.967 N
     (the mock is a subset of the pool). The truth library (global search point 000) is scored against the pool the same
     way, which splits the truth's chi2 / N = 1.50 against the mock into estimator/noise and model parts.
  4. Error model: N_INDEP independent 50k mocks drawn from the pool (without replacement, like the mock) give the true
     sampling covariance of the coefficients. Compared with the mock's bootstrap covariance: variance ratios per layer and
     field, and chi2 / N of the real mock (and of the independent mocks) against the pool with either covariance.
  5. What m <= 4 leaves out of the true maps: amplitudes of m = 1, 3, 6, 8 per ring (|z| < 1, R < 6) fitted to the pool
     together with the even m <= 4 terms.
Maps: face-on square pixels, bar frame, one |z| layer. A pixel's reconstruction is evaluated at the pool particles in that
pixel (mean of <v>(phi_i) of each particle's ring; sigma^2 = mean of sigma^2(phi_i) + the scatter of <v>(phi_i)), so it
is compared like for like with the binned truth. Pixel errors: reconstruction from the joint bootstrap covariance; truth
s / sqrt(N) and s / sqrt(2N).
Outputs: figures/test1_nbody_fourier_star_fit_vs_true_maps_absz_0to1.png, ..._absz_1to2.png,
results/test1_nbody_qmc_focused/fourier_star_fit_vs_true_maps.json.
"""
import json

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, set_style, savefig, INK_2, CMAP_SEQ_W
from plot_au18_velocity_maps import CMAP_DIV
from plot_bar_length_strength_square_maps import grey_bad, VALUE_LIMS
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
from test1_nbody_qmc_focused_run import OUT_DIR, DESIGN_CSV
from test1_nbody_qmc import TRUTH, MOCK

POOL = MOCK_DIR / "nbody_stars_mwframe.hdf5"
N_PIX_MIN = 20
FIELDS = [("<v_R>", 0, "mean"), ("<v_phi>", 1, "mean"), ("sigma_R", 0, "disp"), ("sigma_phi", 1, "disp")]
LAYERS = [(0, 0.25, "absz_0to1", "|z| < 1 kpc"), (1, 0.5, "absz_1to2", "1 < |z| < 2 kpc")]
OUT_JSON = ROOT / "results" / "test1_nbody_qmc_focused" / "fourier_star_fit_vs_true_maps.json"
NK = fs.NK
N_INDEP = 200
PROPOSED = [0 * NK + k for k in fs.NO_M0] + [f * NK + k for f in (1, 3, 4) for k in fs.ALL]   # 19 terms per block


def chi2_per_N(d, cov, n_samples, blocks, idx=PROPOSED):
    """chi2 / N over blocks for coefficient differences d (NB, 30), covariance (NB, 30, 30) estimated from n_samples."""
    h = (n_samples - len(idx) - 2) / (n_samples - 1)
    tot = sum(h * d[i][idx] @ np.linalg.solve(cov[i][np.ix_(idx, idx)], d[i][idx]) for i in blocks)
    return float(tot / (len(idx) * len(blocks)))


def error_model_check(pool, bp, Xp, vp, mp, ok, valid, nstars, c_mock, c_lib, c_pool, cov):
    """Check 4: bootstrap covariance vs the true sampling covariance from independent mocks."""
    in_grid = np.full(len(pool), -1)
    in_grid[np.flatnonzero(mp)] = np.arange(mp.sum())
    rng = np.random.default_rng(123)
    ens = np.empty((N_INDEP, fs.NB, 6 * NK))
    for k in range(N_INDEP):
        g = in_grid[rng.choice(len(pool), 50000, replace=False)]
        g = g[g >= 0]
        ens[k] = fs.fit(bp[g], Xp[g], vp[g], np.ones(len(g)), ok).reshape(fs.NB, 6 * NK)
    p = 50000 / len(pool)
    dev = ens - ens.mean(0)
    C_ens = np.einsum("mbi,mbj->bij", dev, dev) / (N_INDEP - 1) / (1 - p)    # covariance of one 50k mock
    cm, cp = c_mock.reshape(fs.NB, -1), c_pool.reshape(fs.NB, -1)
    ring, layer = np.divmod(np.arange(fs.NB), fs.NL)
    use = valid & (ring < fs.N_RING_DEFAULT)
    out = dict(n_independent_mocks=N_INDEP, layers={})
    for z, lab in enumerate(["absz_0to1", "absz_1to2", "absz_2to3"]):
        blocks = np.flatnonzero(use & (layer == z))
        ratio = np.diagonal(cov[blocks], axis1=1, axis2=2) / np.diagonal(C_ens[blocks], axis1=1, axis2=2)
        sims = [chi2_per_N(ens[k] - cp, cov, 300, blocks) for k in range(N_INDEP)]
        out["layers"][lab] = dict(
            n_blocks=int(len(blocks)), stars_per_block=[int(x) for x in (np.min(nstars[blocks]), np.max(nstars[blocks]))],
            variance_ratio_bootstrap_over_true={nm: round(float(np.median(ratio[:, f * NK:(f + 1) * NK])), 3)
                                                for f, nm in enumerate(fs.FIELD_NAMES)},
            real_mock_chi2_per_N_bootstrap_cov=round(chi2_per_N(cm - cp, cov, 300, blocks), 3),
            real_mock_chi2_per_N_true_cov=round(chi2_per_N((cm - cp) / np.sqrt(1 - p), C_ens, N_INDEP, blocks), 3),
            independent_mocks_chi2_per_N_bootstrap_cov=[round(float(np.mean(sims)), 3),
                                                        round(float(np.percentile(sims, 16)), 3),
                                                        round(float(np.percentile(sims, 84)), 3)])
    blocks = np.flatnonzero(use)
    out["all_R0to6"] = dict(real_mock_chi2_per_N_bootstrap_cov=round(chi2_per_N(cm - cp, cov, 300, blocks), 3),
                            real_mock_chi2_per_N_true_cov=round(chi2_per_N((cm - cp) / np.sqrt(1 - p), C_ens,
                                                                           N_INDEP, blocks), 3),
                            truth_library_chi2_per_N_true_cov_vs_pool=round(chi2_per_N(
                                c_lib.reshape(fs.NB, -1) - cp, C_ens, N_INDEP, blocks), 3))
    return out


def missing_modes(pool, b, X, v, mp):
    """Check 5: amplitudes (km/s) of modes outside the basis in the true maps, |z| < 1, R < 6 kpc, per ring."""
    w = pool[mp]
    phi = np.arctan2(w[:, 1], w[:, 0])
    out = {}
    for a, nm in ((0, "<v_R>"), (1, "<v_phi>")):
        amps = []
        for r in range(fs.N_RING_DEFAULT):
            sel = ((b % fs.NL) == 0) & (b // fs.NL == r)
            F = np.column_stack([X[sel]] + [g(m * phi[sel]) for m in (1, 3, 6, 8) for g in (np.cos, np.sin)])
            coef = np.linalg.lstsq(F, v[sel, a], rcond=None)[0][NK:]
            amps.append(np.hypot(coef[0::2], coef[1::2]))
        amps = np.array(amps)
        out[nm] = {f"m{m}": [round(float(np.median(amps[:, j])), 2), round(float(amps[:, j].max()), 2)]
                   for j, m in enumerate((1, 3, 6, 8))}
    n = np.median([np.sum(((b % fs.NL) == 0) & (b // fs.NL == r)) for r in range(fs.N_RING_DEFAULT)])
    out["note"] = (f"[median, max] over the 20 rings; pool noise of an amplitude ~ sigma sqrt(2/n) = "
                   f"{100 * np.sqrt(2 / n):.2f} km/s for sigma = 100 km/s (median {int(n)} particles per ring)")
    return out


def load(path, sel_R=None):
    with h5py.File(path, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    if sel_R is not None:
        w = w[np.hypot(w[:, 0], w[:, 1]) < sel_R]
    return ssl.rotate_z(w, np.radians(TRUTH["phi_b_deg"]))


def chi2_terms(c1, c2, cov, valid, scheme):
    chi, nd = fs.scheme_chi2(c1, c2, cov, valid, scheme)
    return chi, nd


def pixel_maps(pos, b, X, v, edges, coefs, cov):
    """Truth and reconstruction maps for particles (already restricted to one layer and to valid blocks)."""
    n = len(edges) - 1
    i, j = np.digitize(pos[:, 0], edges) - 1, np.digitize(pos[:, 1], edges) - 1
    inside = (i >= 0) & (i < n) & (j >= 0) & (j < n)
    p = (j * n + i)[inside]
    b, X, v = b[inside], X[inside], v[inside]
    npix = n * n
    N = np.bincount(p, minlength=npix).astype(float)
    out = {"N": N.reshape(n, n)}
    with np.errstate(invalid="ignore", divide="ignore"):
        # sum of the basis per (pixel, block), for the error propagation
        pb = p * fs.NB + b
        Sx = np.stack([np.bincount(pb, weights=X[:, k], minlength=npix * fs.NB) for k in range(NK)], -1)
        Sx = Sx.reshape(npix, fs.NB, NK) / N[:, None, None]
        for name, a, kind in FIELDS:
            s1 = np.bincount(p, weights=v[:, a], minlength=npix)
            s2 = np.bincount(p, weights=v[:, a]**2, minlength=npix)
            mean_t = s1 / N
            var_t = (s2 / N - mean_t**2) * N / (N - 1)
            res = {}
            for tag, c in coefs.items():
                mu = np.einsum("ik,ik->i", X, c[b, a])
                s2f = np.einsum("ik,ik->i", X, c[b, 3 + a])
                m_mu = np.bincount(p, weights=mu, minlength=npix) / N
                v_mu = np.bincount(p, weights=mu**2, minlength=npix) / N - m_mu**2
                m_s2 = np.bincount(p, weights=s2f, minlength=npix) / N
                res[tag] = m_mu if kind == "mean" else np.sqrt(np.maximum(m_s2 + v_mu, 0))
                if tag == "mock":
                    f = a if kind == "mean" else 3 + a
                    C = cov[:, f * NK:(f + 1) * NK, f * NK:(f + 1) * NK]
                    var_c = np.nansum(np.einsum("pbk,bkl,pbl->pb", Sx, np.nan_to_num(C), Sx), axis=1)
                    res["mock_err"] = np.sqrt(var_c) if kind == "mean" else np.sqrt(var_c) / (2 * res["mock"])
            truth = mean_t if kind == "mean" else np.sqrt(var_t)
            terr = np.sqrt(var_t / N) if kind == "mean" else np.sqrt(var_t / (2 * N))
            bad = N < N_PIX_MIN
            out[name] = {k: np.where(bad, np.nan, x).reshape(n, n)
                         for k, x in dict(truth=truth, truth_err=terr, **res).items()}
    return out


def panel_stats(t, te, rec, re=None):
    ok = np.isfinite(t) & np.isfinite(rec)
    d = (rec - t)[ok]
    s = dict(n_pix=int(ok.sum()), rms_kms=float(np.sqrt(np.mean(d**2))), mean_kms=float(np.mean(d)),
             rms_truth_noise_kms=float(np.sqrt(np.mean(te[ok]**2))))
    if re is not None:
        pull = d / np.sqrt(re[ok]**2 + te[ok]**2)
        s.update(pull_mean=float(np.mean(pull)), pull_std=float(np.std(pull)),
                 frac_abs_pull_gt3=float(np.mean(np.abs(pull) > 3)))
    return s


def figure(maps, edges, label, tag, stats):
    cols = ["true map (all 1.49 M particles)", "Fourier m ≤ 4 from the 50k mock",
            "mock Fourier − true [km/s]", "(mock Fourier − true) / error",
            "Fourier m ≤ 4 from all particles − true [km/s]"]
    fig, axes = plt.subplots(len(FIELDS), len(cols), figsize=(3.3 * len(cols), 3.0 * len(FIELDS)),
                             gridspec_kw=dict(hspace=0.32, wspace=0.42))
    th = np.linspace(0, 2 * np.pi, 200)
    for r, (name, _, kind) in enumerate(FIELDS):
        m = maps[name]
        cmap, lo, hi = VALUE_LIMS[name]
        dlim = 6.0 if kind == "mean" else 4.0
        panels = [(m["truth"], cmap, lo, hi), (m["mock"], cmap, lo, hi),
                  (m["mock"] - m["truth"], CMAP_DIV, -3 * dlim, 3 * dlim),
                  ((m["mock"] - m["truth"]) / np.sqrt(m["mock_err"]**2 + m["truth_err"]**2), CMAP_DIV, -3, 3),
                  (m["pool"] - m["truth"], CMAP_DIV, -dlim, dlim)]
        st = stats[name]
        notes = [None, None,
                 f"rms {st['mock']['rms_kms']:.1f} km/s",
                 f"pulls: mean {st['mock']['pull_mean']:+.2f}, std {st['mock']['pull_std']:.2f}",
                 f"rms {st['pool']['rms_kms']:.1f} (truth-pixel noise {st['pool']['rms_truth_noise_kms']:.1f})"]
        for c, ((img, cm, vlo, vhi), note) in enumerate(zip(panels, notes)):
            ax = axes[r, c]
            im = ax.pcolormesh(edges, edges, img, cmap=grey_bad(cm), vmin=vlo, vmax=vhi, rasterized=True)
            ax.plot(6 * np.cos(th), 6 * np.sin(th), color=INK_2, lw=0.6, ls="--")
            ax.set_aspect("equal")
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(edges[0], edges[-1])
            ax.tick_params(labelsize=7)
            if r == 0:
                ax.set_title(cols[c], fontsize=8.5)
            if c == 0:
                ax.set_ylabel(f"{name}\ny [kpc]", fontsize=9)
            if r == len(FIELDS) - 1:
                ax.set_xlabel("x [kpc] (bar along x)", fontsize=8)
            if note:
                ax.text(0.03, 0.03, note, transform=ax.transAxes, fontsize=6.5, color=INK_2,
                        bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
            cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
            cb.ax.tick_params(labelsize=6.5)
    fig.suptitle(f"Per-star Fourier fit (m ≤ 4 per 0.3 kpc ring) vs the true N-body maps, bar frame, {label}, "
                 f"{edges[1] - edges[0]:g} kpc pixels; dashed: R = 6 kpc (edge of the proposed scheme); grey: < "
                 f"{N_PIX_MIN} particles or no fit", fontsize=10, y=0.995)
    savefig(fig, f"test1_nbody_fourier_star_fit_vs_true_maps_{tag}")


def main():
    set_style()
    design = pd.read_csv(DESIGN_CSV)
    p0 = design[design.idx == 0].iloc[0]
    assert p0.kind == "truth" and abs(p0.phi_b_deg - TRUTH["phi_b_deg"]) < 1e-9
    s = np.load(OUT_DIR / "fourier_stars_000.npz")
    c_mock = s["coef_data"].reshape(fs.NB, 6, NK)
    c_lib = s["coef_lib"].reshape(fs.NB, 6, NK)
    cov = s["cov_data_joint"].reshape(fs.NB, 6 * NK, 6 * NK)
    valid = s["valid"].reshape(fs.NB)
    ok = s["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK

    # sanity: refit the mock here and compare with the saved coefficients
    bm, Xm, vm, _ = fs.prepare(load(MOCK))
    c_check = fs.fit(bm, Xm, vm, np.ones(len(bm)), ok)
    refit_diff = float(np.nanmax(np.abs(c_check - c_mock)))

    pool = load(POOL, sel_R=8.0)
    with h5py.File(MOCK, "r") as f:
        attrs = dict(pool_size=int(f.attrs["pool_size"]))
    bp, Xp, vp, mp = fs.prepare(pool)
    c_pool = fs.fit(bp, Xp, vp, np.ones(len(bp)), ok)
    res = dict(n_pool=len(pool), pool_size_in_mock_file=attrs["pool_size"], refit_max_abs_diff=refit_diff,
               expected_chi2_per_N_mock_vs_pool=1 - 50000 / len(pool))

    # 3. coefficient chi2 with the likelihood's covariance
    v_all = valid & np.all(np.isfinite(c_pool), axis=(1, 2))
    schemes = dict(fs.SCHEMES)
    for f, fname in enumerate(fs.FIELD_NAMES):
        schemes[f"field_{fname}_m024"] = (fs.N_RING_DEFAULT, {f: fs.ALL})
        schemes[f"field_{fname}_m0"] = (fs.N_RING_DEFAULT, {f: [0]})
        schemes[f"field_{fname}_m24"] = (fs.N_RING_DEFAULT, {f: fs.NO_M0})
    res["chi2_per_N"] = {}
    for name, sc in schemes.items():
        out = {}
        for tag, (a, bcoef) in dict(mock_vs_pool=(c_mock, c_pool), truth_library_vs_pool=(c_lib, c_pool),
                                    truth_library_vs_mock=(c_lib, c_mock)).items():
            chi, nd = fs.scheme_chi2(a, bcoef, cov, v_all, sc)
            out[tag] = round(chi / nd, 3)
        out["N_data"] = nd
        res["chi2_per_N"][name] = out
    # bias check of the m = 0 terms: mean pull over blocks (R < 6 kpc)
    res["m0_mean_pull_mock_minus_pool"] = {}
    rings = np.repeat(np.arange(fs.NR), fs.NL)
    use = v_all & (rings < fs.N_RING_DEFAULT)
    for f, fname in enumerate(fs.FIELD_NAMES):
        pull = (c_mock[use, f, 0] - c_pool[use, f, 0]) / np.sqrt(cov[use, f * NK, f * NK])
        res["m0_mean_pull_mock_minus_pool"][fname] = [round(float(pull.mean()), 3),
                                                     round(float(pull.std() / np.sqrt(len(pull))), 3), int(len(pull))]

    # 4, 5. error model and missing modes
    res["error_model"] = error_model_check(pool, bp, Xp, vp, mp, ok, v_all, s["n_stars"].reshape(fs.NB), c_mock,
                                            c_lib, c_pool, cov)
    res["modes_outside_basis_kms"] = missing_modes(pool, bp, Xp, vp, mp)

    # 1, 2. maps
    pos = pool[mp]
    res["maps"] = {}
    for layer, pix, tag, label in LAYERS:
        sel = ((bp % fs.NL) == layer) & valid[bp]
        edges = np.round(np.arange(-8.0, 8.0 + 1e-9, pix), 6)
        maps = pixel_maps(pos[sel], bp[sel], Xp[sel], vp[sel], edges, dict(mock=c_mock, pool=c_pool), cov)
        stats = {}
        for name, _, _ in FIELDS:
            m = maps[name]
            stats[name] = dict(mock=panel_stats(m["truth"], m["truth_err"], m["mock"], m["mock_err"]),
                               pool=panel_stats(m["truth"], m["truth_err"], m["pool"]))
        res["maps"][tag] = stats
        figure(maps, edges, label, tag, stats)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

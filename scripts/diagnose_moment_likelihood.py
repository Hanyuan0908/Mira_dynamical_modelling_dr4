"""Diagnostic: KDE likelihood vs moment-based chi^2 for the true and two wrong pattern speeds.

Moment chi^2: bar-frame cells in (R, phi mod pi, |z|) (the 4-fold folding of the library applied to the
data too); per cell the observed means <v_R>, <v_phi>, <v_z sign z> and dispersions sigma_R, sigma_phi,
sigma_z are compared with the library's (1/tau-weighted) values, with observational errors sigma/sqrt(N)
for means and sigma/sqrt(2N) for dispersions; cells with N_obs >= N_MIN. Library noise is neglected.
Also reports how much of the per-star KDE Delta log p comes from the most extreme stars.
"""
import json

import numpy as np
import h5py

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl

OMEGAS = (36.27, 45.0, 30.0)
R_EDGES = np.arange(0.0, 8.0 + 1e-9, 0.5)
PHI_EDGES = np.radians(np.arange(0, 181, 15))
Z_EDGES = np.array([0.0, 0.5, 1.0, 2.0, np.inf])
N_MIN = 20


def cell_moments(w, weight=None):
    """Per-cell weighted N, means and dispersions of (v_R, v_phi, v_z sign z) in folded polar cells."""
    weight = np.ones(len(w)) if weight is None else weight
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.mod(np.arctan2(w[:, 1], w[:, 0]), np.pi)
    az = np.abs(w[:, 2])
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5] * np.sign(w[:, 2])])
    i = np.digitize(R, R_EDGES) - 1
    j = np.digitize(phi, PHI_EDGES) - 1
    k = np.digitize(az, Z_EDGES) - 1
    ok = (i >= 0) & (i < len(R_EDGES) - 1) & (j >= 0) & (j < len(PHI_EDGES) - 1)
    shape = (len(R_EDGES) - 1, len(PHI_EDGES) - 1, len(Z_EDGES) - 1)
    flat = np.ravel_multi_index((i[ok], j[ok], k[ok]), shape)
    n = np.prod(shape)
    W = np.bincount(flat, weights=weight[ok], minlength=n)
    N = np.bincount(flat, minlength=n)
    m1 = np.column_stack([np.bincount(flat, weights=weight[ok] * v[ok, a], minlength=n) for a in range(3)]) / W[:, None]
    m2 = np.column_stack([np.bincount(flat, weights=weight[ok] * v[ok, a]**2, minlength=n) for a in range(3)]) / W[:, None]
    return N, m1, np.sqrt(np.maximum((m2 - m1**2) * N[:, None] / (N[:, None] - 1), 0)), shape   # unbiased variance


def chi2(obs, lib, use_disp, rmin=0.0):
    N, mo, so, shape = obs
    _, ml, sl, _ = lib
    Rc = np.repeat(0.5 * (R_EDGES[1:] + R_EDGES[:-1]), shape[1] * shape[2])
    good = (N >= N_MIN) & (Rc >= rmin) & np.all(np.isfinite(ml), 1)
    c = ((mo - ml)**2 / (so**2 / N[:, None]))[good]
    out = dict(means=float(c.sum()))
    if use_disp:
        cd = ((so - sl)**2 / (so**2 / (2 * N[:, None])))[good]
        out["disp"] = float(cd.sum())
    out["n_cells"] = int(good.sum())
    out["by_R"] = {f"{a:g}-{b:g}": float((c.sum(1) + (cd.sum(1) if use_disp else 0))[(Rc[good] >= a) & (Rc[good] < b)].sum())
                   for a, b in ((0, 1.5), (1.5, 3), (3, 5), (5, 8))}
    return out


def main():
    with h5py.File(MOCK_DIR / "mock_perfect_R8_N50k.hdf5", "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    stars = ssl.rotate_z(stars_mw, np.radians(27.0))
    obs = cell_moments(stars)
    pot = ssl.bar_aligned_potential(1.0, 1.0)
    res = {}
    for Om in OMEGAS:
        times, traj = ssl.integrate(pot, stars, -Om)
        lib, _, lib_w, _ = ssl.build_library(times, traj)
        lm = cell_moments(lib, lib_w)
        res[Om] = {f"means, R>={rmin:g}": chi2(obs, lm, False, rmin) for rmin in (0.0, 1.5)}
        res[Om].update({f"means+disp, R>={rmin:g}": chi2(obs, lm, True, rmin) for rmin in (0.0, 1.5)})
        del lib, lib_w, traj

    print(f"cells with N_obs >= {N_MIN}: {res[OMEGAS[0]]['means, R>=0']['n_cells']} (all R)")
    for key in res[OMEGAS[0]]:
        ref = res[OMEGAS[0]][key]
        tot = lambda r: r["means"] + r.get("disp", 0.0)
        line = f"{key:20s} chi2(truth) = {tot(ref):8.1f} (dof~{ref['n_cells'] * (6 if 'disp' in key else 3)})"
        for Om in OMEGAS[1:]:
            r = res[Om][key]
            line += f" | Om={Om:g}: d(logL) = {-0.5 * (tot(r) - tot(ref)):+8.1f}"
        print(line)
        for Om in OMEGAS[1:]:
            d = {k: -0.5 * (res[Om][key]["by_R"][k] - ref["by_R"][k]) for k in ref["by_R"]}
            print(f"      Om={Om:g} by R: " + "  ".join(f"{k}: {v:+7.1f}" for k, v in d.items()))

    # tail dominance of the per-star KDE differences
    lp = np.load(ROOT / "results" / "test1_qmc" / "per_star_logp_omega.npz")
    tails = {}
    for Om in OMEGAS[1:]:
        d = lp[f"Om{Om:g}"] - lp["Om36.27"]
        o = np.argsort(-np.abs(d))
        tails[Om] = {f"top{p}%": float(d[o[: int(len(d) * p / 100)]].sum()) for p in (0.1, 1, 5)}
        tails[Om]["total"] = float(d.sum())
        print(f"KDE dlogL (Om={Om:g} - truth) = {d.sum():+.1f}; from the 0.1/1/5% most extreme stars: "
              f"{tails[Om]['top0.1%']:+.1f} / {tails[Om]['top1%']:+.1f} / {tails[Om]['top5%']:+.1f}; "
              f"rms per-star d = {d.std():.3f}")
    with open(ROOT / "results" / "test1_qmc" / "moment_vs_kde_diagnostic.json", "w") as f:
        json.dump(dict(moment_chi2={str(k): v for k, v in res.items()}, kde_tails={str(k): v for k, v in tails.items()}),
                  f, indent=2)


if __name__ == "__main__":
    main()

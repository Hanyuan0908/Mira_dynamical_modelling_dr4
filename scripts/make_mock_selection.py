"""Test 2a/2b mock: N-body stars observed with Mira-like distance errors and selected like the Gaia DR3 Miras.

Recipe (user, 2026-10-29):
  1. every N-body star particle (2 M) is observed from the Sun with a Gaussian distance-modulus error
     (DM_obs = DM_true + sigma_DM N(0, 1); l, b, proper motions and v_los exact; observe.py), and the observed
     Galactocentric phase space is recomputed;
  2. the observed particles are selected so that they follow the sky distribution (l, b) of the DR3 Miras, their
     distribution in heliocentric distance and their distribution in Galactocentric radius.
Target: DR3 Miras (data/gaia_dr3_lpv_mira.csv) with DM_err < 0.4 and Galactocentric R < 8 kpc (catalogue x, y, z),
the modelling volume. Particles enter only if their observed R < 8 kpc, as the Miras do.
Selection function, defined on observed quantities:  S(l, b, d) = S_sky(HEALPix pixel of (l, b)) * S_dist(d bin)
* S_R(R bin), R = R(l, b, d), fitted by iterative proportional fitting (cycling through the three corrections) so that
the expected selected counts reproduce the target's sky map (HEALPix nside NSIDE), distance histogram (D_BIN kpc) and
Galactocentric-R histogram (R_BIN kpc), scaled to N expected stars. Each particle is then kept independently with probability min(S, 1) (Bernoulli; seed SEED).
Writes mock/nbody_mock_selection_sigmaDM<sigma>_N<n>k.hdf5:
  observed x..vz, d, l, b, pml, pmb, vlos under the usual names; true values as <name>_true; DM_true, DM_obs, DM_err;
  star_index (into nbody_stars_mwframe.hdf5), comp, mass; and the selection function: S_sky (HEALPix, nest=False,
  Galactic l, b), dist_edges and S_dist, R_edges and S_R, so that S can be evaluated anywhere (e.g. to weight the
  orbit library).
Diagnostic figure: figures/test2_nbody_mock_selection_vs_dr3_miras_sigmaDM<sigma>.png
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import healpy as hp
import matplotlib.pyplot as plt

from common import MOCK_DIR, MIRA_DR3_FILE, R0_KPC, Z_SUN_KPC, V_SUN_KMS, set_style, savefig, INK, INK_2, ACCENT, \
    BLUE_RAMP, CMAP_SEQ_W
from galaxies import galaxy
import observe

NSIDE = 32
D_BIN = 0.25
D_MAX = 20.0
R_BIN = 0.25
R_MAX = 8.0
SIGMA_DM_CUT = 0.4
N_ITER = 10
SEED_ERR = 2031
SEED_SEL = 2032


def mira_target():
    df = pd.read_csv(MIRA_DR3_FILE, usecols=["l", "b", "DM_err", "luminosity_distance", "x", "y", "z"]).dropna()
    df = df[(df.DM_err < SIGMA_DM_CUT) & (np.hypot(df.x, df.y) < R_MAX)]
    return df


def fit_selection(index, targets, n_target):
    """Iterative proportional fitting of S = prod_k s_k[index_k] to the target margins (each scaled to n_target)."""
    targets = [t * n_target / t.sum() for t in targets]
    factors = [np.ones(len(t)) for t in targets]
    for _ in range(N_ITER):
        for k, (ix, t) in enumerate(zip(index, targets)):
            others = np.prod([f[i] for j, (f, i) in enumerate(zip(factors, index)) if j != k], axis=0)
            have = np.bincount(ix, weights=others, minlength=len(t))
            factors[k] = np.divide(t, have, out=np.zeros(len(t)), where=have > 0)
    return factors


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sigma-dm", type=float, default=0.075 / (np.log(10) / 5), help="default 7.5% (0.1629 mag)")
    ap.add_argument("--n", type=int, default=50000)
    args = ap.parse_args()
    set_style()
    out = MOCK_DIR / f"nbody_mock_selection_sigmaDM{args.sigma_dm:.3f}_N{args.n // 1000}k.hdf5"

    with h5py.File(galaxy("nbody")["stars"], "r") as f:
        cols = {k: f[k][:] for k in ("x", "y", "z", "vx", "vy", "vz", "mass", "comp")}
    obs_true = observe.to_galactic(*(cols[k] for k in ("x", "y", "z", "vx", "vy", "vz")))
    rng = np.random.default_rng(SEED_ERR)
    d_obs = observe.scatter_distance_modulus(obs_true["d"], args.sigma_dm, rng)
    x, y, z, vx, vy, vz = observe.to_galactocentric(obs_true["l"], obs_true["b"], d_obs, obs_true["pml"],
                                                    obs_true["pmb"], obs_true["vlos"])
    R_obs = np.hypot(x, y)
    inside = R_obs < R_MAX

    edges = np.arange(0.0, D_MAX + 1e-9, D_BIN)
    pix_all = hp.ang2pix(NSIDE, obs_true["l"], obs_true["b"], lonlat=True)
    dbin_all = np.clip(np.digitize(d_obs, edges) - 1, 0, len(edges) - 2)
    tgt = mira_target()
    t_sky = np.bincount(hp.ang2pix(NSIDE, tgt.l.values, tgt.b.values, lonlat=True), minlength=hp.nside2npix(NSIDE))
    t_dist = np.bincount(np.clip(np.digitize(tgt.luminosity_distance.values, edges) - 1, 0, len(edges) - 2),
                         minlength=len(edges) - 1)
    r_edges = np.arange(0.0, R_MAX + 1e-9, R_BIN)
    rbin_all = np.clip(np.digitize(R_obs, r_edges) - 1, 0, len(r_edges) - 2)
    t_R = np.bincount(np.clip(np.digitize(np.hypot(tgt.x.values, tgt.y.values), r_edges) - 1, 0, len(r_edges) - 2),
                      minlength=len(r_edges) - 1)
    s_sky, s_dist, s_R = fit_selection([pix_all[inside], dbin_all[inside], rbin_all[inside]],
                                       [t_sky.astype(float), t_dist.astype(float), t_R.astype(float)], args.n)
    S = np.where(inside, s_sky[pix_all] * s_dist[dbin_all] * s_R[rbin_all], 0.0)
    capped = S > 1
    keep = np.random.default_rng(SEED_SEL).random(len(S)) < np.minimum(S, 1.0)   # independent of the error draw
    idx = np.flatnonzero(keep)

    with h5py.File(out, "w") as f:
        for k, v in zip(("x", "y", "z", "vx", "vy", "vz"), (x, y, z, vx, vy, vz)):
            f[k] = v[idx]
            f[f"{k}_true"] = cols[k][idx]
        f["d"], f["d_true"] = d_obs[idx], obs_true["d"][idx]
        for k in ("l", "b", "pml", "pmb", "vlos"):
            f[k] = obs_true[k][idx]
        f["DM_true"] = 5 * np.log10(obs_true["d"][idx] * 100.0)
        f["DM_obs"] = 5 * np.log10(d_obs[idx] * 100.0)
        f["DM_err"] = np.full(len(idx), args.sigma_dm)
        f["star_index"], f["comp"], f["mass"] = idx, cols["comp"][idx], cols["mass"][idx]
        f["S_sky"], f["dist_edges"], f["S_dist"], f["R_edges"], f["S_R"] = s_sky, edges, s_dist, r_edges, s_R
        f.attrs.update(dict(
            description="N-body star particles observed with Mira-like distance errors (Gaussian in DM; l, b, pm, "
                        "v_los exact) and selected like the DR3 Miras (DM_err < 0.4, R < 8 kpc): S(l, b, d_obs) = "
                        "S_sky[HEALPix(l, b)] * S_dist[d_obs bin] * S_R[R_obs bin] for observed R < 8 kpc, kept with probability "
                        "min(S, 1). Observed values under the usual names, true values as <name>_true.",
            source=galaxy("nbody")["stars"].name, sigma_DM=args.sigma_dm,
            frac_distance_error=np.log(10) / 5 * args.sigma_dm, n_expected=args.n, n=len(idx),
            healpix_nside=NSIDE, healpix_nest=False, d_bin_kpc=D_BIN, R_bin_kpc=R_BIN, R_max_obs=R_MAX, sigma_DM_cut_target=SIGMA_DM_CUT,
            n_target_miras=len(tgt), seed_distance_error=SEED_ERR, seed_selection=SEED_SEL,
            frac_particles_with_S_above_1=float(capped[inside].mean()), R0=R0_KPC, z_sun=Z_SUN_KPC, v_sun=V_SUN_KMS))

    # diagnostics: margins of mock vs target, and what the selection does to the face-on distribution
    l_wrap = lambda l: ((l + 180) % 360) - 180
    lm, bm = l_wrap(obs_true["l"][idx]), obs_true["b"][idx]
    lt, bt = l_wrap(tgt.l.values), tgt.b.values
    sky_mock = np.bincount(pix_all[idx], minlength=len(t_sky))
    occupied = t_sky > 0
    res = dict(n=len(idx), frac_capped=float(capped[inside].mean()),
               sky_corr=float(np.corrcoef(sky_mock[occupied], t_sky[occupied] * len(idx) / t_sky.sum())[0, 1]))
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), gridspec_kw=dict(hspace=0.32, wspace=0.28))
    le, be = np.arange(-90, 90.01, 2.0), np.arange(-20, 20.01, 1.0)
    for ax, (lv, bv, title) in zip(axes[0, :2], ((lt, bt, f"DR3 Miras (σ_DM < 0.4, R < 8 kpc), {len(tgt):,}"),
                                                (lm, bm, f"mock, {len(idx):,} stars"))):
        h = np.histogram2d(lv, bv, bins=[le, be])[0].T / len(lv)
        im = ax.pcolormesh(le, be, h, cmap=CMAP_SEQ_W, vmin=0, vmax=0.004, rasterized=True)
        ax.set_xlim(90, -90)
        ax.set_xlabel("l [deg]")
        ax.set_ylabel("b [deg]")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, label="fraction per 2° × 1° pixel")
    ax = axes[0, 2]
    db = np.arange(0, 16.01, 0.5)
    for v, lab, c, ls in ((tgt.luminosity_distance.values, "DR3 Miras", INK, "-"), (d_obs[idx], "mock (observed d)", ACCENT, "-"),
                          (obs_true["d"][idx], "mock (true d)", ACCENT, ":"),
                          (d_obs[inside], "all particles, observed R < 8", BLUE_RAMP[4], "--")):
        ax.hist(v, bins=db, density=True, histtype="step", color=c, ls=ls, lw=1.6, label=lab)
    ax.set_xlabel("heliocentric distance [kpc]")
    ax.set_ylabel("density")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("distance distribution", fontsize=10)
    xe = np.arange(-8, 8.01, 0.25)
    for ax, (xx, yy, title) in zip(axes[1, :2], ((tgt.x.values, tgt.y.values, "DR3 Miras, face-on"),
                                                (x[idx], y[idx], "mock, face-on (observed positions)"))):
        h = np.histogram2d(xx, yy, bins=[xe, xe])[0].T
        ax.pcolormesh(xe, xe, np.log10(np.maximum(h, 0.5)), cmap=CMAP_SEQ_W, rasterized=True)
        ax.plot(-R0_KPC, 0, marker="*", color=ACCENT, ms=12)
        ax.set_aspect("equal")
        ax.set_xlim(-8.5, 8)
        ax.set_ylim(-8, 8)
        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("y [kpc]")
        ax.set_title(title + " (log counts)", fontsize=10)
    ax = axes[1, 2]
    rb = np.arange(0, 8.01, 0.25)
    for v, lab, c, ls in ((np.hypot(tgt.x, tgt.y), "DR3 Miras", INK, "-"), (np.hypot(x[idx], y[idx]), "mock (observed)", ACCENT, "-"),
                          (np.hypot(cols["x"][idx], cols["y"][idx]), "mock (true)", ACCENT, ":")):
        ax.hist(v, bins=rb, density=True, histtype="step", color=c, ls=ls, lw=1.6, label=lab)
    ax.set_xlabel("Galactocentric R [kpc]")
    ax.set_ylabel("density")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("radial distribution", fontsize=10)
    fig.suptitle(f"Mock with distance errors (σ_DM = {args.sigma_dm:.3f} mag, {100 * np.log(10) / 5 * args.sigma_dm:.1f}%) "
                 f"and a DR3-Mira-like selection in (l, b), distance and R; HEALPix nside {NSIDE}, distance and R "
                 f"bins {D_BIN} kpc", fontsize=11, y=0.995)
    savefig(fig, f"test2_nbody_mock_selection_vs_dr3_miras_sigmaDM{args.sigma_dm:.3f}")
    print(f"wrote {out.name}: {len(idx):,} stars (expected {args.n:,}); particles with S > 1 (capped): "
          f"{100 * res['frac_capped']:.2f}%; sky-map correlation mock vs target {res['sky_corr']:.3f}")
    print(f"  distance median: Miras {np.median(tgt.luminosity_distance):.2f}, mock {np.median(d_obs[idx]):.2f} kpc; "
          f"R median: Miras {np.median(np.hypot(tgt.x, tgt.y)):.2f}, mock {np.median(np.hypot(x[idx], y[idx])):.2f} kpc; "
          f"|b| median: Miras {np.median(np.abs(bt)):.2f}, mock {np.median(np.abs(bm)):.2f} deg")
    print(f"  true |z| > 1 kpc: {np.mean(np.abs(cols['z'][idx]) > 1):.1%}; comp: disc {np.mean(cols['comp'][idx] == 0):.1%}")


if __name__ == "__main__":
    main()

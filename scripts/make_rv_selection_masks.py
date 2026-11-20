"""Thin an existing mock to Gaia-like RV sub-samples, without new integrations (user, 2026-11-21).

The mock already mimics the spatial distribution of the DR3 Miras (DM_err < 0.4). Each scenario keeps a mock
star with a probability measured on the REAL catalogues at that star's observed position: among the K_NN real stars
nearest to it in observed (x, y, z),
  g14, g15, g16   the fraction of DR3 Miras with G (phot_g_mean_mag) < 14, 15, 16;
  grvs14          the fraction with G_RVS < 14 (measured, else G_RP - 0.642; plot_presentation_mira_grvs_distribution);
                  G_RVS < 15 and < 16 keep 98% and ~100% of the Miras, so they are not run;
  osarg_cover     the OSARG-to-Mira density ratio (the fraction of OSARGs among the K_NN nearest of the pooled
                  Mira + OSARG catalogue, turned into a ratio), scaled so the expected count equals the OSARGs at
                  observed R < 8 kpc (26,630);
  random_osargN   a constant, giving the same expected count: sample size alone, coverage unchanged.
The last two separate the two explanations of the weak OSARG constraint: spatial coverage, or sample size.
Magnitude cuts depend on position only through distance and extinction; their period dependence is not modelled.

Usage: python make_rv_selection_masks.py [--mock nbody_mock_selection_sigmaDM0.271_N50k.hdf5]
Writes results/rv_selection_masks/<mock stem>__<scenario>.npy (bool, one per mock star) and a summary table, and
figures/test2_t7_mock_thinned_to_gaia_rv_scenarios_face_on.png.
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree

from common import ROOT, MOCK_DIR, RESULTS_DIR, R0_KPC, set_style, savefig, rel, INK, INK_2
from make_osarg_data import final_sample

K_NN = 32
SEED = 2121
OUT = RESULTS_DIR / "rv_selection_masks"
GRVS_OFFSET = -0.642            # G_RVS - G_RP on the 2,958 Miras with a measured grvs_mag


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", default="nbody_mock_selection_sigmaDM0.271_N50k.hdf5")
    args = ap.parse_args()
    with h5py.File(MOCK_DIR / args.mock, "r") as f:
        pos = np.column_stack([f[k][:] for k in ("x", "y", "z")])
    n = len(pos)
    m = pd.read_csv(ROOT / "data" / "gaia_dr3_lpv_mira.csv",
                    usecols=["x", "y", "z", "DM_err", "phot_g_mean_mag", "phot_rp_mean_mag", "grvs_mag"])
    m = m[m.DM_err < 0.4]
    grvs = np.where(np.isfinite(m.grvs_mag), m.grvs_mag, m.phot_rp_mean_mag + GRVS_OFFSET)
    mira = m[["x", "y", "z"]].values
    o = final_sample()
    osarg = np.column_stack([o["x"], o["y"], o["z"]])
    n_osarg_R8 = int((np.hypot(o["x"], o["y"]) < 8).sum())

    _, nn = cKDTree(mira).query(pos, k=K_NN)
    prob = {f"g{c}": (m.phot_g_mean_mag.values[nn] < c).mean(1) for c in (14, 15, 16)}
    prob["grvs14"] = (grvs[nn] < 14).mean(1)
    pooled = np.vstack([mira, osarg])
    is_osarg = np.r_[np.zeros(len(mira), bool), np.ones(len(osarg), bool)]
    _, nn2 = cKDTree(pooled).query(pos, k=2 * K_NN)
    q = is_osarg[nn2].mean(1)
    ratio = np.where(q < 1, q / np.clip(1 - q, 1e-9, None), np.inf) * len(mira) / len(osarg)
    c = n_osarg_R8 / n
    for _ in range(50):                  # scale so that sum(min(c * ratio, 1)) = n_osarg_R8
        c *= n_osarg_R8 / np.minimum(c * ratio, 1).sum()
    prob["osarg_cover"] = np.minimum(c * ratio, 1)
    prob["random_osargN"] = np.full(n, n_osarg_R8 / n)

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    stem = args.mock.removesuffix(".hdf5")
    R, az = np.hypot(pos[:, 0], pos[:, 1]), np.abs(pos[:, 2])
    rows = [dict(scenario="full mock", n=n, R_lt_6=int((R < 6).sum()),
                 R_lt_6_z0to1=int(((R < 6) & (az < 1)).sum()), R_lt_6_z1to2=int(((R < 6) & (az >= 1) & (az < 2)).sum()),
                 beyond_centre=int((pos[:, 0] > 0).sum()))]
    masks = {}
    for name, p in prob.items():
        keep = rng.random(n) < p
        masks[name] = keep
        np.save(OUT / f"{stem}__{name}.npy", keep)
        rows.append(dict(scenario=name, n=int(keep.sum()), R_lt_6=int((keep & (R < 6)).sum()),
                         R_lt_6_z0to1=int((keep & (R < 6) & (az < 1)).sum()),
                         R_lt_6_z1to2=int((keep & (R < 6) & (az >= 1) & (az < 2)).sum()),
                         beyond_centre=int((keep & (pos[:, 0] > 0)).sum())))
    tab = pd.DataFrame(rows)
    tab["frac"] = (tab.n / n).round(3)
    real = {f"G < {c}": int((m.phot_g_mean_mag < c).sum()) for c in (14, 15, 16)}
    real.update({f"G_RVS < {c}": int((grvs < c).sum()) for c in (14, 15, 16)})
    print(f"real DR3 Miras (DM_err < 0.4): {len(m):,}; " + ", ".join(f"{k}: {v:,} ({100 * v / len(m):.0f}%)"
                                                                   for k, v in real.items()))
    print(f"OSARGs: {len(osarg):,} final sample, {n_osarg_R8:,} at observed R < 8 kpc; OSARG scale factor c = {c:.3f}, "
          f"{(prob['osarg_cover'] >= 1).sum():,} mock stars capped at p = 1")
    print(tab.to_string(index=False))
    tab.to_csv(OUT / f"{stem}__summary.csv", index=False)

    set_style()
    order = ["full mock", "g16", "g15", "g14", "grvs14", "osarg_cover", "random_osargN"]
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6), sharex=True, sharey=True)
    edges = np.arange(-12, 12.01, 0.4)
    norm = LogNorm(vmin=1, vmax=400)
    for ax, name in zip(axes.ravel(), order):
        sel = np.ones(n, bool) if name == "full mock" else masks[name]
        h = np.histogram2d(pos[sel, 0], pos[sel, 1], bins=[edges, edges])[0].T
        im = ax.imshow(np.where(h > 0, h, np.nan), origin="lower", extent=[edges[0], edges[-1]] * 2, norm=norm,
                       cmap="viridis", interpolation="nearest")
        ax.plot(-R0_KPC, 0, "*", color="w", ms=9, mec=INK)
        ax.set_title(f"{name}: {int(sel.sum()):,} stars, {int((sel & (R < 6)).sum()):,} at R < 6",
                     color=INK, fontsize=10)
        ax.set_aspect("equal")
    axes.ravel()[-1].axis("off")
    for ax in axes[1]:
        ax.set_xlabel("x [kpc]")
    for ax in axes[:, 0]:
        ax.set_ylabel("y [kpc]")
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="stars per 0.4 kpc pixel")
    fig.suptitle(f"t = 7 mock (12.5% distance error, Mira-like selection) thinned to Gaia RV scenarios "
                 f"(probabilities from the real catalogues, {K_NN} nearest neighbours)", color=INK_2, y=0.98)
    savefig(fig, "test2_t7_mock_thinned_to_gaia_rv_scenarios_face_on")


if __name__ == "__main__":
    main()

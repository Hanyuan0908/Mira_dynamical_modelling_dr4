"""The Milky Way bar's strength profile A2/A0(R) from the OSARG posterior, with a bar-angle prior (user, 2026-11-24).

Posterior: OSARG data in the adopted Sormani+22 family (run osarg_zhang2024_sormani22, scheme B, sfratio, K = 10, the
default grid), at the adopted temperature T = sigma_2 (the chi2 jitter between two near-identical models at the best
fit; results/chi2_jitter_osarg_s22_random.csv), i.e. logL = -dchi2 / (2 T); emulator and sampler as the pipeline.
A phi_b prior N(25, 2 deg) is applied by importance weighting. N_DRAW draws are taken with those weights; for each, the
trial potential at (S, bar_frac) is rebuilt and the m = 2 amplitude of its midplane density measured
(plot_bar_strength_profile_recovery.a2_over_a0; phi_b and Omega_b do not change the potential's shape). Shown: the
median and 16-84% / 2.5-97.5% bands, against the published Sormani+22 bar and Hunter+24's (both S = bar_frac = 1).
Also: the peak A2/A0, its radius, and the half-peak radius beyond it (a bar-length proxy), per draw.
The family scales the non-axisymmetric density by bar_frac, so at large bar_frac the midplane density goes NEGATIVE
where the bar is thin (A2/A0 > 1 there). Such draws are not physical models: they are reported and left out of the
bands (a physical prior on bar_frac in effect).

Usage: python plot_bar_strength_profile_osarg.py [--n-draw 300] [--workers 6]
Figure: figures/osarg_zhang2024_sormani22_bar_strength_profile_phib25pm2.png
Table:  results/osarg_zhang2024_sormani22/bar_strength_profile_phib25pm2.npz
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge
from compare_grid_strategies import posterior
from plot_bar_strength_profile_recovery import a2_over_a0, R_GRID, _PTS, _RR, _PP
import steady_state_likelihood as ssl

RUN = "osarg_zhang2024_sormani22"
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
PHI0, SIG_PHI = 25.0, 2.0


def curve(args):
    """A2/A0(R) of the midplane density, and the fraction of the midplane grid where the density is negative."""
    f, S = args
    rho = ssl.bar_aligned_potential(f, S, "sormani22").density(_PTS).reshape(_RR.shape)
    a0 = rho.mean(axis=1)
    a2 = 2 * np.abs((rho * np.exp(-2j * _PP)).mean(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(a0 > 0, a2 / a0, np.nan), float(np.mean(rho < 0))


def r_half(prof):
    k = int(np.nanargmax(prof))
    half = 0.5 * prof[k]
    for i in range(k, len(R_GRID) - 1):
        if prof[i] >= half > prof[i + 1]:
            return R_GRID[i] + (prof[i] - half) / (prof[i] - prof[i + 1]) * (R_GRID[i + 1] - R_GRID[i])
    return np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draw", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    set_style()
    jit = pd.read_csv(RESULTS_DIR / "chi2_jitter_osarg_s22_random.csv")
    x = jit.loc[jit.param == "random", "C_data"].values
    T = float(np.std((x[:, None] - x[None, :])[np.triu_indices(len(x), 1)], ddof=1))
    ge.use_run(RUN)
    ev = pd.read_csv(RESULTS_DIR / RUN / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    chi2 = tr[f"chi2_{B}_sfratio"].values
    P = posterior(ge.to_unit(tr), -(chi2 - chi2.min()) / (2 * T))
    w = np.exp(-0.5 * ((P[:, 1] - PHI0) / SIG_PHI) ** 2)
    ess = w.sum() ** 2 / np.sum(w ** 2) / len(w)
    rng = np.random.default_rng(7)
    idx = rng.choice(len(P), args.n_draw, replace=True, p=w / w.sum())
    D = P[idx]
    q = lambda a: np.percentile(a, [16, 50, 84])
    print(f"T = sigma_2 = {T:.1f}; phi_b prior N({PHI0}, {SIG_PHI}): ESS fraction {ess:.2f}")
    for k, i in (("|Omega_b|", 0), ("phi_b", 1), ("S", 2), ("bar_frac", 3)):
        a = q(D[:, i])
        print(f"   {k:9s} {a[1]:.2f} ({a[0]:.2f}-{a[2]:.2f})")
    with ProcessPoolExecutor(args.workers) as ex:
        res = list(ex.map(curve, [(float(f), float(S)) for f, S in zip(D[:, 3], D[:, 2])]))
    curves_all = np.array([r[0] for r in res])
    negfrac = np.array([r[1] for r in res])
    phys = negfrac == 0
    print(f"   draws with negative midplane density somewhere: {np.sum(~phys)} of {len(D)} ({np.mean(~phys):.0%}); "
          f"their bar_frac median {np.median(D[~phys, 3]) if (~phys).any() else np.nan:.2f}, min "
          f"{D[~phys, 3].min() if (~phys).any() else np.nan:.2f}")
    with ProcessPoolExecutor(args.workers) as ex:
        onset = list(ex.map(curve, [(f, 1.0) for f in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5)]))
    print("   at S = 1, fraction of the midplane grid with negative density for bar_frac 1.0 ... 1.5: "
          + " ".join(f"{o[1]:.3f}" for o in onset))
    curves = curves_all[phys]
    D_all, D = D, D[phys]
    for k, i in (("|Omega_b|", 0), ("phi_b", 1), ("S", 2), ("bar_frac", 3)):
        a = q(D[:, i])
        print(f"   physical draws: {k:9s} {a[1]:.2f} ({a[0]:.2f}-{a[2]:.2f})")
    ref_s22 = a2_over_a0("sormani22", 1.0, 1.0)
    ref_h24 = a2_over_a0("mw", 1.0, 1.0)
    band = np.nanpercentile(curves, [2.5, 16, 50, 84, 97.5], axis=0)
    peak = np.nanmax(curves, axis=1)
    rpk = R_GRID[np.nanargmax(curves, axis=1)]
    rh = np.array([r_half(c) for c in curves])
    for lab, a in (("peak A2/A0", peak), ("R of the peak [kpc]", rpk), ("half-peak radius beyond it [kpc]", rh)):
        b = np.nanpercentile(a, [16, 50, 84])
        print(f"   {lab:33s} {b[1]:.3f} ({b[0]:.3f}-{b[2]:.3f})")
    print(f"   published Sormani+22 bar: peak {np.nanmax(ref_s22):.3f} at R {R_GRID[np.nanargmax(ref_s22)]:.1f}, half-peak "
          f"radius {r_half(ref_s22):.2f}; Hunter+24: peak {np.nanmax(ref_h24):.3f} at R {R_GRID[np.nanargmax(ref_h24)]:.1f}, "
          f"half-peak radius {r_half(ref_h24):.2f}")
    np.savez(RESULTS_DIR / RUN / "bar_strength_profile_phib25pm2.npz", R=R_GRID, curves_all=curves_all, draws_all=D_all,
             negative_density_fraction=negfrac, band_physical=band, ref_sormani22=ref_s22, ref_hunter24=ref_h24, T=T, ess=ess)

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.fill_between(R_GRID, band[0], band[4], color=BLUE_RAMP[2], alpha=0.35, lw=0, label="OSARG posterior, 95%")
    ax.fill_between(R_GRID, band[1], band[3], color=BLUE_RAMP[4], alpha=0.55, lw=0, label="OSARG posterior, 68%")
    ax.plot(R_GRID, band[2], color=BLUE_RAMP[6], lw=2.2, label="OSARG posterior, median")
    ax.plot(R_GRID, ref_s22, color=INK, lw=1.6, ls="--", label="Sormani+22 bar as published (S = bar_frac = 1)")
    ax.plot(R_GRID, ref_h24, color=ACCENT, lw=1.4, ls=":", label="Hunter+24 (same bar, own disc and centre)")
    ax.axvspan(6, R_GRID[-1], color=INK_MUTED, alpha=0.12, lw=0)
    ax.text(6.12, 0.55, "not fitted\n(R > 6 kpc)", fontsize=8, color=INK_2)
    ax.set(xlabel="R [kpc]", ylabel=r"$A_2/A_0$ of the midplane density", xlim=(R_GRID[0], R_GRID[-1]), ylim=(0, None))
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.set_title("Milky Way bar strength profile from the OSARGs (Sormani+22 family, T = σ₂ = "
                 f"{T:.0f}, φ_b prior {PHI0:g} ± {SIG_PHI:g}°)\n{len(D)} of {len(D_all)} draws; the {len(D_all) - len(D)} "
                 "with negative midplane density somewhere are left out", fontsize=10, color=INK)
    savefig(fig, "osarg_zhang2024_sormani22_bar_strength_profile_phib25pm2")


if __name__ == "__main__":
    main()

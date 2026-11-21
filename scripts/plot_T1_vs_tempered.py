"""The 25 Test 2 galaxies at T = 1 (no tempering) against the tempered posterior (rerun_emulators_T1.py).

1. Pattern speed: recovered |Omega_b| (median, 16-84%) against the truth for both, with the fractional bias, the
   median width and the coverage of the 16-84% interval. Truth convention as everywhere: the simulation table's
   finite-difference pattern speed for the N-body snapshots, our own measurement for Auriga.
2. Bar strength profile A2/A0(R) of the midplane density (plot_bar_strength_profile_recovery.a2_over_a0): the trial
   potential rebuilt at 300 draws of each chain, against the truth (bar_frac = S = 1), for the galaxies whose
   potential files are on disk (all 12 N-body snapshots and Au18; the other 12 Auriga potentials are not).

Usage: python plot_T1_vs_tempered.py [--skip-profiles]
Figures: figures/test2_pattern_speed_recovery_12p5pct_T1_vs_tempered.png,
         figures/test2_bar_strength_profile_recovery_T1_vs_tempered.png
Table:   results/T1_vs_tempered_summary.csv, results/bar_strength_profile_T1_vs_tempered.npz
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from plot_presentation_figures import NBODY, AURIGA, TRUTH_TABLE
import plot_test2_pattern_speed_recovery as psr

N_DRAW = 300
PROFILE_GALAXIES = NBODY + ["au18"]


def truth_omega(g, meta, tab):
    if g.startswith("au"):
        return float(meta["omega_truth"])
    t = 7.0 if meta["galaxy"] == "nbody" else float(meta["galaxy"].split("_t")[1])
    return float(tab.loc[np.isclose(tab.snap_time, t), "pattern_speed_fd"].iloc[0])


def profile_job(g):
    from plot_bar_strength_profile_recovery import a2_over_a0
    run = f"test2_{g}_sf_12p5pct"
    gal = json.load(open(RESULTS_DIR / run / "run_meta.json"))["galaxy"]
    ch = np.load(RESULTS_DIR / run / "posterior_chains_T1.npz")
    rng = np.random.default_rng(7)
    out = {"truth": a2_over_a0(gal, 1.0, 1.0)}
    for tag in ("T1", "tempered"):
        post = ch[tag]
        idx = rng.choice(len(post), N_DRAW, replace=True)
        curves = np.array([a2_over_a0(gal, float(post[i, 3]), float(post[i, 2])) for i in idx])
        out[tag] = np.nanpercentile(curves, [16, 50, 84], axis=0)
    return g, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-profiles", action="store_true")
    args = ap.parse_args()
    set_style()
    tab = pd.read_csv(TRUTH_TABLE)
    rows = []
    for g in NBODY + AURIGA:
        run = f"test2_{g}_sf_12p5pct"
        meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
        r = json.load(open(RESULTS_DIR / run / "emulator_results_T1.json"))
        tr = truth_omega(g, meta, tab)
        for tag in ("T1", "tempered"):
            p = r[tag]["posterior"]
            o = p["Omega_abs"]
            side = (o["p84"] - o["p50"]) if tr > o["p50"] else (o["p50"] - o["p16"])
            rows.append(dict(galaxy=g, auriga=g.startswith("au"), T=tag, truth=tr, med=o["p50"], lo=o["p16"],
                             hi=o["p84"], bias_pct=100 * (o["p50"] / tr - 1), width_pct=100 * (o["p84"] - o["p16"]) / tr,
                             z=(tr - o["p50"]) / side, loo=r[tag]["loo_rms"], phi=p["phi_b_deg"]["p50"],
                             S=p["S"]["p50"], bar_frac=p["bar_frac"]["p50"]))
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "T1_vs_tempered_summary.csv", index=False)
    for tag in ("tempered", "T1"):
        s = t[t["T"] == tag]
        for lab, ss in (("N-body", s[~s.auriga]), ("Auriga", s[s.auriga]), ("all 25", s)):
            print(f"[{tag:8s}] {lab:7s}: mean |bias| {ss.bias_pct.abs().mean():.2f}%, mean bias {ss.bias_pct.mean():+.2f}%, "
                  f"median 16-84% width {ss.width_pct.median():.2f}%, coverage {np.mean(ss.z.abs() < 1):.2f}, "
                  f"median |z| {ss.z.abs().median():.2f}, median LOO {ss.loo.median():.2f}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw=dict(wspace=0.28))
    lim = (12, 72)
    a1.plot(lim, lim, ls="--", color=INK_MUTED, lw=1)
    a2.axhline(0, ls="--", color=INK_MUTED, lw=1)
    for tag, col, off, lab in (("tempered", INK_MUTED, -0.35, "tempered (√(2N))"), ("T1", BLUE_RAMP[5], 0.35, "T = 1")):
        s = t[t["T"] == tag]
        for au, mk in ((False, "o"), (True, "D")):
            q = s[s.auriga == au]
            a1.errorbar(q.truth + off, q.med, yerr=[q.med - q.lo, q.hi - q.med], fmt=mk, color=col, ms=5.5,
                        capsize=0, lw=1.2, label=f"{lab}, {'Auriga' if au else 'N-body'}")
            a2.errorbar(q.truth + off, q.bias_pct, yerr=[100 * (q.med - q.lo) / q.truth, 100 * (q.hi - q.med) / q.truth],
                        fmt=mk, color=col, ms=5.5, capsize=0, lw=1.2)
    a1.set(xlim=lim, ylim=lim, xlabel=r"true $|\Omega_b|$", ylabel=r"recovered $|\Omega_b|$ (median, 16–84%)")
    a1.legend(frameon=False, fontsize=8.5, loc="upper left")
    a2.set(xlim=lim, ylim=(-16, 17), xlabel=r"true $|\Omega_b|$", ylabel="bias [%]")
    cov = {tag: np.mean(t[t["T"] == tag].z.abs() < 1) for tag in ("T1", "tempered")}
    a2.set_title(f"coverage of the 16–84% interval: T = 1 {cov['T1']:.2f}, tempered {cov['tempered']:.2f}",
                 color=INK, fontsize=10)
    fig.suptitle("Test 2, 25 galaxies at 12.5% distance error: pattern speed at T = 1 against the tempered posterior "
                 "(points offset ±0.35 for legibility)", color=INK_2, y=0.99, fontsize=11)
    savefig(fig, "test2_pattern_speed_recovery_12p5pct_T1_vs_tempered")

    if args.skip_profiles:
        return
    from plot_bar_strength_profile_recovery import R_GRID
    with ProcessPoolExecutor(6) as ex:
        res = dict(ex.map(profile_job, PROFILE_GALAXIES))
    np.savez_compressed(RESULTS_DIR / "bar_strength_profile_T1_vs_tempered.npz", R=R_GRID,
                        **{f"{g}_{k}": v for g, o in res.items() for k, v in o.items()})
    nc = 5
    nr = int(np.ceil(len(PROFILE_GALAXIES) / nc))
    fig, AX = plt.subplots(nr, nc, figsize=(3.3 * nc, 2.8 * nr), sharex=True, sharey=True)
    for i, g in enumerate(PROFILE_GALAXIES):
        ax = AX.ravel()[i]
        o = res[g]
        ax.fill_between(R_GRID, o["tempered"][0], o["tempered"][2], color=INK_MUTED, alpha=0.3, lw=0,
                        label="tempered, 16–84%")
        ax.fill_between(R_GRID, o["T1"][0], o["T1"][2], color=BLUE_RAMP[3], alpha=0.6, lw=0, label="T = 1, 16–84%")
        ax.plot(R_GRID, o["T1"][1], color=BLUE_RAMP[6], lw=1.4, label="T = 1 median")
        ax.plot(R_GRID, o["truth"], color=INK, lw=1.8, label="truth")
        ax.set_title(g, color=INK, fontsize=9)
        ax.set_ylim(0, 1.0)
        ax.set_xlim(R_GRID[0], R_GRID[-1])
        pk = int(np.nanargmax(o["truth"]))
        print(f"{g:7s} A2/A0 at the truth's peak (R {R_GRID[pk]:.1f}): truth {o['truth'][pk]:.3f}; "
              f"tempered {o['tempered'][1][pk]:.3f} ({o['tempered'][0][pk]:.3f}-{o['tempered'][2][pk]:.3f}); "
              f"T=1 {o['T1'][1][pk]:.3f} ({o['T1'][0][pk]:.3f}-{o['T1'][2][pk]:.3f})")
    for j in range(len(PROFILE_GALAXIES), nr * nc):
        AX.ravel()[j].axis("off")
    AX.ravel()[-1].legend(*AX.ravel()[0].get_legend_handles_labels(), frameon=False, fontsize=10, loc="center")
    for ax in AX[-1]:
        ax.set_xlabel("R [kpc]")
    for ax in AX[:, 0]:
        ax.set_ylabel(r"$A_2/A_0$")
    fig.suptitle("Bar strength profile, T = 1 against the tempered posterior (300 draws each); 12 N-body + Au18 "
                 "(other Auriga potentials not on disk)", color=INK, y=0.995, fontsize=11)
    fig.tight_layout()
    savefig(fig, "test2_bar_strength_profile_recovery_T1_vs_tempered")


if __name__ == "__main__":
    main()

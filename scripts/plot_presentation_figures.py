"""Presentation versions (black background) of two Test 2 results at 12.5% distance error.

1. Pattern speed recovered against the truth for the 12 N-body snapshots (blue circles) and the 13 Auriga
   galaxies (orange diamonds), with the fractional bias beside it. Same runs, column and truth convention as
   figures/test2_pattern_speed_recovery_12p5pct_12pt_plus_auriga.png (plot_test2_pattern_speed_recovery.py):
   scheme B with sfratio; N-body truth = the simulation table's finite-difference column, Auriga truth = our own
   Dehnen measurement (run_meta omega_truth). One difference: a run that carries an explicit K = 10 column
   uses it. Only t = 7 does, where the plain `sfratio` column is neither K = 1 nor K = 10 (PROGRESS
   2026-11-05) and Omega_b moves by 0.02 km/s/kpc.
2. The bar-strength profile A2/A0(R) of the midplane density for t = 7: the truth (the galaxy's own potential,
   bar_frac = S = 1) and the trial potential rebuilt at 300 draws of the flat-prior posterior (median and
   16-84%), as in plot_bar_strength_profile_recovery.py (same chain seeds and draw seed), but on the K = 10
   column.

Colours: the reference palette's dark steps, blue #3987e5 and orange #d95926, validated on #000000
(dataviz validate_palette.js --mode dark --pairs all: all checks pass, CVD dE 26.8); shapes differ too.

Usage: python plot_presentation_figures.py [--truth-table PATH] [--n-draw 300]
Figures: figures/presentation/{pattern_speed_recovered_vs_true_nbody_and_auriga_12p5pct_black,
         bar_strength_profile_A2_over_A0_truth_vs_recovered_nbody_t7_black}.{png,pdf}
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from common import RESULTS_DIR, FIG_DIR, rel
import plot_test2_pattern_speed_recovery as psr

SCHEME_B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
NBODY = ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25", "t7"]
AURIGA = ["au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24", "au26", "au27", "au28"]
TRUTH_TABLE = Path.home() / "Desktop/PhD_projects/SchwarMAX_data/Bar_model_TG21/model/bar_strength_and_pattern_speed.csv"
OUT = FIG_DIR / "presentation"

BG = "#000000"
TEXT = "#f0efec"           # primary text
TEXT_2 = "#c3c2b7"         # secondary text
AXIS = "#6b6a66"           # spines, ticks
GRID = "#262624"
BLUE, ORANGE = "#3987e5", "#d95926"


def dark_style():
    mpl.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "savefig.edgecolor": BG,
        "axes.edgecolor": AXIS, "axes.labelcolor": TEXT, "axes.titlecolor": TEXT, "text.color": TEXT,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2, "xtick.labelcolor": TEXT_2, "ytick.labelcolor": TEXT_2,
        "font.size": 14, "axes.labelsize": 16, "axes.titlesize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
        "legend.fontsize": 13, "axes.linewidth": 1.0, "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.8, "axes.axisbelow": True, "font.family": "sans-serif",
        "mathtext.fontset": "dejavusans", "legend.frameon": False,
    })


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext, dpi in (("png", 300), ("pdf", None)):
        fig.savefig(OUT / f"{name}.{ext}", dpi=dpi, facecolor=BG, bbox_inches="tight", pad_inches=0.15)
    print(f"saved {rel(OUT / name)}.png/.pdf")


def recovery_points(table):
    """One point per galaxy: truth, posterior median and 16-84% of |Omega_b|."""
    import pandas as pd
    tab = pd.read_csv(table)
    pts = []
    for g in NBODY + AURIGA:
        run = RESULTS_DIR / f"test2_{g}_sf_12p5pct"
        meta = json.load(open(run / "run_meta.json"))
        f = run / f"emulator_results_logLs_{SCHEME_B}_sfratioK10.json"
        f = f if f.exists() else run / f"emulator_results_logLs_{SCHEME_B}_sfratio.json"
        p = json.load(open(f))["posterior"]["Omega_abs"]
        auriga = g.startswith("au")
        if auriga:
            truth = float(meta["omega_truth"])
        else:
            t = 7.0 if meta["galaxy"] == "nbody" else float(meta["galaxy"].split("_t")[1])
            truth = float(tab.loc[np.isclose(tab.snap_time, t), "pattern_speed_fd"].iloc[0])
        pts.append(dict(galaxy=g, auriga=auriga, truth=truth, med=p["p50"], lo=p["p16"], hi=p["p84"],
                        column=f.stem[len("emulator_results_"):]))
    return pts


def plot_recovery(pts):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.0, 5.6), gridspec_kw=dict(wspace=0.28))
    lim = (12.0, 72.0)
    a1.plot(lim, lim, ls=(0, (5, 4)), color=TEXT_2, lw=1.2, zorder=1)
    a2.axhspan(-5, 5, color="#141413", lw=0, zorder=0)
    a2.axhline(0, ls=(0, (5, 4)), color=TEXT_2, lw=1.2, zorder=1)
    a2.text(lim[0] + 0.8, 5.4, "±5%", color=TEXT_2, fontsize=12, ha="left", va="bottom")
    for auriga, color, marker, ms, name in ((False, BLUE, "o", 8.5, "N-body snapshots"),
                                            (True, ORANGE, "D", 8.0, "Auriga galaxies")):
        s = [p for p in pts if p["auriga"] == auriga]
        x = np.array([p["truth"] for p in s])
        m = np.array([p["med"] for p in s])
        lo, hi = np.array([p["lo"] for p in s]), np.array([p["hi"] for p in s])
        bias = 100 * (m / x - 1)
        lab = f"{name} ({len(s)}): mean |bias| {np.mean(np.abs(bias)):.1f}%"
        kw = dict(fmt=marker, ms=ms, color=color, mfc=color, mec=BG, mew=1.4, elinewidth=1.8, capsize=0,
                  zorder=3)
        a1.errorbar(x, m, yerr=[m - lo, hi - m], label=lab, **kw)
        a2.errorbar(x, bias, yerr=[100 * (m - lo) / x, 100 * (hi - m) / x], **kw)
    a1.set(xlim=lim, ylim=lim, xlabel=r"true $|\Omega_b|$  [km s$^{-1}$ kpc$^{-1}$]",
           ylabel=r"recovered $|\Omega_b|$  [km s$^{-1}$ kpc$^{-1}$]")
    a1.set_aspect("equal")
    h, lab = a1.get_legend_handles_labels()      # above the panels, where no data can sit under it
    fig.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, 0.925), ncol=2, handletextpad=0.4,
               columnspacing=2.5)
    fig.subplots_adjust(top=0.80)
    a2.set(xlim=lim, ylim=(-16, 17), xlabel=r"true $|\Omega_b|$  [km s$^{-1}$ kpc$^{-1}$]",
           ylabel=r"bias in $|\Omega_b|$  [%]")
    for a in (a1, a2):
        a.tick_params(length=4, width=1.0)
    fig.text(0.5, 0.995, "Bar pattern speed recovered from Mira-like mocks with 12.5% distance errors "
             "and a Gaia-like selection", ha="center", va="top", color=TEXT, fontsize=16)
    save(fig, "pattern_speed_recovered_vs_true_nbody_and_auriga_12p5pct_black")


def t7_profile(n_draw):
    """Truth and the flat-posterior 16/50/84% of A2/A0(R) for t = 7 (K = 10 column)."""
    import emcee
    from scipy.optimize import minimize
    import test1_gp_emulator as ge
    from plot_bar_strength_profile_recovery import a2_over_a0, R_GRID
    run, column = "test2_t7_sf_12p5pct", f"logLs_{SCHEME_B}_sfratioK10"
    ge.use_run(run)
    df, _ = ge.load_table(column, None, None)
    tr = df[df.kind == "qmc"]
    u, y = ge.to_unit(tr), tr["logL"].values
    emu = ge.Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    sam = emcee.EnsembleSampler(32, 4, lambda q: -np.inf if np.any(q < 0) | np.any(q > 1)
                                else emu.predict(q)[0])
    sam.random_state = np.random.RandomState(2).get_state()
    sam.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1),
                 3000, progress=False)
    post = ge.from_unit(sam.get_chain(discard=1000, thin=5, flat=True))
    idx = np.random.default_rng(7).choice(len(post), n_draw, replace=True)
    curves = np.array([a2_over_a0("nbody", post[i, 3], post[i, 2]) for i in idx])
    q = np.nanpercentile(curves, [16, 50, 84], axis=0)
    truth = a2_over_a0("nbody", 1.0, 1.0)
    np.savez(RESULTS_DIR / "test2_t7_sf_12p5pct" / f"bar_strength_profile_t7_{column}.npz", R=R_GRID,
             truth=truth, p16=q[0], p50=q[1], p84=q[2], n_draw=n_draw)
    return R_GRID, truth, q, np.percentile(post, [16, 50, 84], axis=0)


def plot_profile(R, truth, q):
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.fill_between(R, q[0], q[2], color=BLUE, alpha=0.32, lw=0, label="recovered, 16–84%", zorder=2)
    ax.plot(R, q[1], color=BLUE, lw=2.4, label="recovered, median", zorder=3)
    ax.plot(R, truth, color=TEXT, lw=2.4, label="truth", zorder=4)
    ax.axvline(6.0, color=AXIS, lw=1.0, ls=(0, (2, 3)), zorder=1)
    ax.text(6.1, 0.30, "← fitted range\n    R < 6 kpc", color=TEXT_2, fontsize=12, ha="left", va="center")
    ax.set(xlim=(R[0], R[-1]), ylim=(0, 0.75), xlabel=r"$R$  [kpc]",
           ylabel=r"bar strength  $A_2/A_0$  (midplane density)")
    ax.legend(loc="upper right", handlelength=1.8)
    ax.tick_params(length=4, width=1.0)
    ax.set_title("N-body snapshot t = 7: bar strength profile, truth vs recovered", color=TEXT, fontsize=16,
                 pad=12)
    save(fig, "bar_strength_profile_A2_over_A0_truth_vs_recovered_nbody_t7_black")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth-table", type=Path, default=psr.TABLE if psr.TABLE.exists() else TRUTH_TABLE)
    ap.add_argument("--n-draw", type=int, default=300)
    args = ap.parse_args()
    dark_style()
    pts = recovery_points(args.truth_table)
    print(f"{'galaxy':8s} {'truth':>7s} {'median':>7s} {'16-84%':>13s} {'bias':>7s}  column")
    for p in pts:
        print(f"{p['galaxy']:8s} {p['truth']:7.2f} {p['med']:7.2f} {p['lo']:6.2f}-{p['hi']:6.2f} "
              f"{100 * (p['med'] / p['truth'] - 1):+6.1f}%  {p['column'][len('logLs_' + SCHEME_B) + 1:]}")
    plot_recovery(pts)
    saved = RESULTS_DIR / "test2_t7_sf_12p5pct" / f"bar_strength_profile_t7_logLs_{SCHEME_B}_sfratioK10.npz"
    if saved.exists() and int(np.load(saved)["n_draw"]) == args.n_draw:     # the draws are seeded: re-use
        s = np.load(saved)
        R, truth, q = s["R"], s["truth"], np.array([s["p16"], s["p50"], s["p84"]])
    else:
        R, truth, q, pq = t7_profile(args.n_draw)
        print(f"t = 7 posterior (K = 10): S {pq[1, 2]:.3f} ({pq[0, 2]:.3f}-{pq[2, 2]:.3f}), "
              f"bar_frac {pq[1, 3]:.3f} ({pq[0, 3]:.3f}-{pq[2, 3]:.3f})")
    i = int(np.nanargmax(truth))
    j = int(np.nanargmax(q[1]))
    print(f"A2/A0 peak: truth {truth[i]:.3f} at R = {R[i]:.1f} kpc; recovered median {q[1, j]:.3f} "
          f"({q[0, j]:.3f}-{q[2, j]:.3f}) at R = {R[j]:.1f} kpc")
    plot_profile(R, truth, q)


if __name__ == "__main__":
    main()

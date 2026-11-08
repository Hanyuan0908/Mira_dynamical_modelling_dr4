"""Dynamical bar length from x1 orbits, recovered against the truth.

The photometric length is not what a dynamicist means by the bar: Petersen, Weinberg & Katz (2021)
argue for the radial extent of the orbits actually TRAPPED by it, which is a property of the
potential and pattern speed rather than of the light. Zhang et al. (2024, arXiv:2406.06678) make
that operational with a frequency-ratio cut in the bar frame,

    Omega_R / Omega_x = 2 +- 0.2    and    r_apo < 7 kpc,

and read the length off the extent of the orbits that pass (R_x1 ~ 3.8 kpc for the Milky Way).

This reads the tables written by `integrate_orbit_frequencies.py` - one per truth potential and one
per posterior draw - applies that cut, and compares. The truth is the galaxy's own potential at
(bar_frac, S) = (1, 1) with the tabulated Omega_b, integrating ORIGINAL SNAPSHOT PARTICLES (no
distance error, no selection function); the model is N draws from the posterior re-weighted by the
kinematic phi_b prior, integrating the mock stars AS OBSERVED.

WHICH STATISTIC IS THE LENGTH is a choice, and the maximum is not it: the largest R_apo in the x1
set lands on the r_apo < 7 kpc cut itself (6.95 kpc for Au26), so it measures the cut, not the bar.
Every candidate is computed and tabulated; `--stat` picks the one that is plotted. p95 of R_apo is
the default as the steadiest analogue of "the extent of the x1 orbits".

Usage: python bar_dynamical_length.py [--stat p95_R_apo] [--min-x1 50]
Figure: figures/test2_bar_dynamical_length.png
Table:  results/bar_dynamical_length.csv
"""
import argparse
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

IN_DIR = RESULTS_DIR / "bar_dynamics"
RATIO, TOL, R_APO_MAX = 2.0, 0.2, 7.0          # Zhang et al. (2024)
STATS = ["max_R_apo", "p99.5_R_apo", "p99_R_apo", "p95_R_apo", "p90_R_apo",
         "p99.5_R_apo_mean", "p95_R_apo_mean", "p90_R_apo_mean", "med_R_apo"]


def lengths(path, min_x1):
    """Every candidate length statistic for one frequency table, plus the x1 count."""
    d = np.load(path, allow_pickle=True)
    c, T = list(d["cols"]), d["table"]
    g = lambda k: T[:, c.index(k)]
    with np.errstate(invalid="ignore", divide="ignore"):
        r = g("Omega_R_withv") / g("Omega_x_withv")
    sel = (np.abs(r - RATIO) <= TOL) & (g("r_apo") < R_APO_MAX)
    if sel.sum() < min_x1:
        return None
    Ra, Rm = g("R_apo")[sel], g("R_apo_mean")[sel]
    out = dict(n_x1=int(sel.sum()), f_x1=float(sel.mean()),
               max_R_apo=float(Ra.max()),
               # p99.5 is the usual literature choice; at n_x1 ~ 1400 it is the top ~7 orbits, so
               # it is far noisier than p95 and sits much closer to the r_apo < 7 kpc cut.
               **{"p99.5_R_apo": float(np.percentile(Ra, 99.5)),
                  "p99.5_R_apo_mean": float(np.percentile(Rm, 99.5))},
               p99_R_apo=float(np.percentile(Ra, 99)),
               p95_R_apo=float(np.percentile(Ra, 95)),
               p90_R_apo=float(np.percentile(Ra, 90)),
               p95_R_apo_mean=float(np.percentile(Rm, 95)),
               p90_R_apo_mean=float(np.percentile(Rm, 90)),
               med_R_apo=float(np.median(Ra)),
               axis_ratio=float(np.median(g("y_max")[sel] / g("x_max")[sel])))
    # integrate_orbit_frequencies.py stores the four parameters by name; the pipeline's own
    # freq_NNN.npz (test2b --freq-snap) stores them as a `params` array in the same order.
    if "params" in d.files:
        for k, v in zip(("Omega_abs", "phi_b_deg", "S", "bar_frac"), np.asarray(d["params"])):
            out[k] = float(v)
    else:
        for k in ("Omega_abs", "phi_b_deg", "S", "bar_frac"):
            out[k] = float(d[k])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stat", default="p95_R_apo", choices=STATS)
    ap.add_argument("--min-x1", type=int, default=50,
                    help="skip a table with fewer trapped orbits than this")
    args = ap.parse_args()
    set_style()

    rows = []
    for td in sorted(IN_DIR.glob("truth_*.npz")):
        gal = td.name[len("truth_"):-len(".npz")]
        t = lengths(td, args.min_x1)
        if t is None:
            print(f"  skip {gal}: truth has too few x1 orbits")
            continue
        run = next((d for d in IN_DIR.glob("test2_*")
                    if json.load(open(RESULTS_DIR / d.name / "run_meta.json"))["galaxy"] == gal),
                   None)
        draws = [lengths(f, args.min_x1) for f in sorted(run.glob("draw_*.npz"))] if run else []
        draws = [x for x in draws if x is not None]
        if not draws:
            print(f"  {gal}: truth only, no draws yet")
            continue
        v = np.array([x[args.stat] for x in draws])
        q = np.percentile(v, [16, 50, 84])
        rows.append(dict(galaxy=gal, run=run.name, n_draw=len(draws),
                         L_true=t[args.stat], L16=q[0], L50=q[1], L84=q[2],
                         bias=100 * (q[1] / t[args.stat] - 1),
                         n_x1_true=t["n_x1"], f_x1_true=t["f_x1"],
                         f_x1_model=float(np.mean([x["f_x1"] for x in draws])),
                         axis_ratio_true=t["axis_ratio"],
                         **{f"true_{s}": t[s] for s in STATS}))
    if not rows:
        raise SystemExit("no galaxy has both a truth and draws yet")
    D = pd.DataFrame(rows).sort_values("L_true")
    D.to_csv(RESULTS_DIR / "bar_dynamical_length.csv", index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 5.5), gridspec_kw=dict(wspace=0.26))
    au = D.galaxy.str.startswith("au").values
    lim = (0.85 * min(D.L_true.min(), D.L16.min()), 1.15 * max(D.L_true.max(), D.L84.max()))
    a1.plot(lim, lim, "--", color=INK, lw=1.1, zorder=1, label="perfect recovery")
    for m, c, mk, lab in ((au, ACCENT, "s", "Auriga"), (~au, BLUE_RAMP[5], "o", "N-body")):
        if not m.any():
            continue
        a1.errorbar(D.L_true[m], D.L50[m], yerr=[(D.L50 - D.L16)[m], (D.L84 - D.L50)[m]],
                    fmt=mk, ms=7, color=c, mec=INK, mew=0.8, capsize=3, lw=1.2, zorder=4, label=lab)
    a1.set_xlim(*lim); a1.set_ylim(*lim); a1.set_aspect("equal")
    a1.set_xlabel(f"true dynamical bar length [kpc]  ({args.stat})")
    a1.set_ylabel("recovered dynamical bar length (16-84%) [kpc]")
    a1.legend(loc="upper left", fontsize=9)
    a1.set_title(f"{len(D)} galaxies, {int(D.n_draw.min())}-{int(D.n_draw.max())} draws each",
                 color=INK_2, fontsize=10)

    a2.axhline(0, color=INK, lw=1.1, ls="--", zorder=1)
    rms = float(np.sqrt((D.bias ** 2).mean()))
    a2.axhspan(-rms, rms, color=BLUE_RAMP[1], alpha=0.5, lw=0, zorder=0, label=f"rms {rms:.1f}%")
    for _, r in D.iterrows():
        c = ACCENT if r.galaxy.startswith("au") else BLUE_RAMP[5]
        a2.plot(r.L_true, r.bias, "s" if r.galaxy.startswith("au") else "o", ms=8, color=c,
                mec=INK, mew=0.8, zorder=4)
        a2.annotate(r.galaxy.replace("nbody_", "").replace("nbody", "t7"), (r.L_true, r.bias),
                    xytext=(0, 9), textcoords="offset points", ha="center", fontsize=7.5,
                    color=INK_2)
    a2.set_xlabel(f"true dynamical bar length [kpc]  ({args.stat})")
    a2.set_ylabel("bias in dynamical length [%]")
    a2.legend(loc="lower right", fontsize=9)
    a2.set_title(f"mean {D.bias.mean():+.1f}%, mean |bias| {D.bias.abs().mean():.1f}%",
                 color=INK_2, fontsize=10)
    fig.suptitle("Dynamical bar length from trapped x1 orbits "
                 r"($\Omega_R/\Omega_x = 2 \pm 0.2$, $r_{\rm apo} < 7$ kpc; Zhang et al. 2024), "
                 "recovered against truth", color=INK, y=0.98, fontsize=11.5)
    fig.tight_layout()
    savefig(fig, "test2_bar_dynamical_length")

    print(f"  {'galaxy':8s} {'n_x1':>6s} {'truth':>6s} {'model':>6s} {'bias':>7s}  "
          f"{'f_x1 true/model':>16s}")
    for _, r in D.iterrows():
        print(f"  {r.galaxy:8s} {int(r.n_x1_true):6d} {r.L_true:6.2f} {r.L50:6.2f} "
              f"{r.bias:+6.1f}%  {100*r.f_x1_true:7.1f}% /{100*r.f_x1_model:6.1f}%")
    print(f"\n  mean bias {D.bias.mean():+.1f}%, mean |bias| {D.bias.abs().mean():.1f}%, rms {rms:.1f}%")
    print("\n  how the choice of statistic moves the truth (kpc):")
    print("  " + " ".join(f"{s:>15s}" for s in ["galaxy"] + STATS))
    for _, r in D.iterrows():
        print("  " + f"{r.galaxy:>15s}" + " ".join(f"{r['true_' + s]:15.2f}" for s in STATS))
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_dynamical_length.csv')}")


if __name__ == "__main__":
    main()

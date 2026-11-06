"""Angular bias of the bar-angle estimators against distance error, at a fixed true bar angle.

How much the measured m = 2 phase is dragged away from the truth depends on how badly the distances
are known, and the two estimators respond very differently: <v_R> is a signed streaming signal whose
inflow and outflow lobes cancel asymmetrically when smearing mixes azimuths, while <|v_R|/v_tot> is
unsigned and only loses amplitude. This measures both, on one galaxy, with everything else held
fixed.

Held fixed deliberately:
  - the true bar angle (default 27 deg, the Wegg & Gerhard value the mocks are built at);
  - the selection function, taken from the 12.5% mock and applied at every error level, so the
    curve isolates the distance error rather than mixing in a changing selection;
  - the window, default self:0.3 - rings whose own m = 2 amplitude exceeds 0.3 of its peak, adopted
    2026-11-05 because it gives the best blind transfer (1.80 deg rms) and is the only rule whose
    formal band and blind error agree.

Usage: python plot_bar_angle_bias_vs_distance_error.py [--galaxy au18] [--window kin0.3]
Figure: figures/test2_bar_angle_bias_vs_distance_error.png
Table:  results/bar_angle_bias_vs_distance_error.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import (MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel,
                    INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP)
from galaxies import galaxy
from bar_angle_vecsum_verification import phase
import observe, selection, steady_state_likelihood as ssl

EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}
FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
SEL_SIGMA = 0.271                      # the selection tables to use, held fixed across the scan


def one(job):
    """Measured phase of both estimators at one (distance error, realisation)."""
    gname, phi_true, frac, seed, window = job
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{SEL_SIGMA:.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - phi_true)).T)
    sig = frac * 5 / np.log(10)
    d = observe.scatter_distance_modulus(o["d"], sig, np.random.default_rng(seed)) if frac > 0 \
        else o["d"]
    w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
    S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
    sel = w[np.random.default_rng(seed + 500_000).random(len(S)) < np.minimum(S, 1.0)]
    win = window if window.startswith("kin") else tuple(float(v) for v in window.split("-"))
    # phase() returns the bar AXIS in the input frame, which is MINUS the bar angle measured from
    # the Sun-GC line (test_bar_angle_transfer_12p5.phases negates it for exactly this reason).
    # So phi_measured = -phase(), and the bias is -phase() - phi_true. The zero-error row is the
    # check: it must come out near 0, not near -27.
    return frac, {e: -phase(sel, e, win) - phi_true for e in EST}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", default="au18")
    ap.add_argument("--phi-true", type=float, default=27.0)
    ap.add_argument("--window", default="kin0.3", help="kin<F> for self:<F>, or 'LO-HI' kpc")
    ap.add_argument("--n-real", type=int, default=8)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    jobs = [(args.galaxy, args.phi_true, f, 101 + s, args.window)
            for f in FRACS for s in range(args.n_real)]
    print(f"{len(jobs)} measurements: {args.galaxy}, phi_true {args.phi_true}, "
          f"window {args.window}", flush=True)
    acc = {}
    with ProcessPoolExecutor(args.workers) as ex:
        for k, (frac, ph) in enumerate(ex.map(one, jobs), 1):
            for e, v in ph.items():
                acc.setdefault((frac, e), []).append(v)
            if k % 10 == 0:
                print(f"  {k}/{len(jobs)}", flush=True)

    rows = []
    for f in FRACS:
        r = dict(frac=f, pct=100 * f)
        for e in EST:
            v = np.array(acc[(f, e)])
            r[f"{e}_bias"] = float(v.mean())
            r[f"{e}_sem"] = float(v.std(ddof=1) / np.sqrt(len(v)))
            r[f"{e}_scatter"] = float(v.std(ddof=1))
        rows.append(r)
    D = pd.DataFrame(rows)
    out = RESULTS_DIR / f"bar_angle_bias_vs_distance_error{args.tag}.csv"
    D.to_csv(out, index=False)

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    col = {"vR": ACCENT, "absvR": BLUE_RAMP[5]}
    ax.axhline(0, color=INK_MUTED, lw=1.1, ls="--", zorder=1)
    for e, lab in EST.items():
        ax.errorbar(D.pct, D[f"{e}_bias"], yerr=D[f"{e}_sem"], marker="o", ms=6, lw=1.9,
                    capsize=3, color=col[e], label=lab, zorder=4)
        ax.fill_between(D.pct, D[f"{e}_bias"] - D[f"{e}_scatter"],
                        D[f"{e}_bias"] + D[f"{e}_scatter"], color=col[e], alpha=0.16, lw=0,
                        zorder=2)
    ax.axvline(12.5, color=INK_MUTED, lw=0.9, ls=":", zorder=1)
    ax.annotate("DR3 Miras\n12.5%", (12.5, ax.get_ylim()[0]), xytext=(4, 6),
                textcoords="offset points", fontsize=8.5, color=INK_2, va="bottom")
    ax.set_xlabel("fractional distance error [%]")
    ax.set_ylabel(r"angular bias  $\phi_{\rm measured} - \phi_{\rm true}$  [deg]")
    g = galaxy(args.galaxy)
    wl = (f"rings above {args.window[3:]} of the peak $m=2$ amplitude"
          if args.window.startswith("kin") else f"fixed {args.window} kpc")
    ax.set_title(f"{g['label']}, true $\\phi_b$ = {args.phi_true:g}$^\\circ$   |   "
                 f"window {args.window} ({wl})\n"
                 f"shading = scatter over {args.n_real} distance-error realisations; "
                 f"selection fixed at the 12.5% tables", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=10, loc="lower left")
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_bias_vs_distance_error{args.tag}")

    print(f"\n  {'error':>7s} | " + " | ".join(f"{EST[e]:>26s}" for e in EST))
    print("  " + "-" * 68)
    for _, r in D.iterrows():
        print(f"  {r.pct:6.1f}% | " +
              " | ".join(f"{r[e + '_bias']:+9.2f} +- {r[e + '_sem']:.2f} "
                         f"(sd {r[e + '_scatter']:.2f})" for e in EST))
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()

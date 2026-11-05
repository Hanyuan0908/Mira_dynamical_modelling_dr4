"""Bar-angle bias against distance error, for four snapshots, both kinematic estimators.

No orbit integration and no potential: the bar angle is read straight off the kinematics. A bar drives
an m = 2 pattern in the velocity field whose phase tracks the bar's major axis, so the angle is
measurable from the stars alone. Two estimators:

  <v_R>            the quadrupole in mean radial velocity; its node lies along the bar axis
  <|v_R|/v_tot>    Zhang et al. (2024); unsigned, so inflow and outflow lobes do not cancel

Per-ring m = 2 coefficients are combined by **vector sum** over 1-4 kpc (each ring weights itself by
its own amplitude), which PROGRESS 2026-11-01 showed halves the snapshot-to-snapshot spread against a
plain median. The galaxy is rotated so the bar sits at TRUE_ANGLE to the Sun-GC line, observed from the
fixed Sun with a distance-modulus error, and the DR3-Mira-like selection applied.

The quantity plotted is (measured phase - true angle): the bias a real measurement would carry. What
matters for calibration is not its size but its **spread across snapshots**, since a bias common to all
galaxies can be subtracted and one that varies cannot.

Usage: python plot_bar_angle_drift_multi.py [--snaps nbody_t2 nbody_t3.625 nbody_t5 nbody]
                                            [--n-real 10] [--sigma-ref 0.125]
Figure: figures/test2_bar_angle_drift_multi.png   Table: results/bar_angle_drift_multi.csv
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, \
    INK_MUTED, ACCENT, BLUE_RAMP
from galaxies import galaxy
from nbody import SNAPSHOT_DIR
import observe, selection, steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
from bar_angle_vecsum_verification import phase, EST

FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
SNAPS = ["nbody_t2", "nbody_t3.625", "nbody_t5", "nbody"]      # Omega_b 43.9, 34.8, 30.6, 25.2
TRUE_ANGLE = BAR_ANGLE_DEG                                      # 27 deg

# Radius at which the disc's density m = 2 amplitude peaks, i.e. the bar's own scale. Measured in
# scripts/diagnose_bar_evolution.py: the bar is SHORTEST when it is strongest, so a window fixed in
# kpc samples a different fraction of the bar at each epoch and lets spiral arms in early on.
R_PEAK = {"nbody_t2": 2.38, "nbody_t3.625": 2.62, "nbody_t5": 3.12, "nbody": 3.38}


def t_of(gname):
    return 7.0 if gname == "nbody" else float(gname.split("_t")[1])


def measure(gname, phi_true, frac, n_real, seed0, sigma_ref, window=None):
    """Mean and SEM of both estimators at one (galaxy, bar angle, distance error)."""
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{sigma_ref * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - phi_true)).T)
    sig = frac * 5 / np.log(10)
    got = {k: [] for k in EST}
    for s in range(n_real):
        rng = np.random.default_rng(seed0 + s)
        d = observe.scatter_distance_modulus(o["d"], sig, rng) if frac > 0 else o["d"]
        w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
        S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
        sel = w[np.random.default_rng(seed0 + 500_000 + s).random(len(S)) < np.minimum(S, 1.0)]
        for k in got:
            got[k].append(phase(sel, k, window) if window else phase(sel, k))
    # phase() returns MINUS the measured bar angle: phi_measured = -phase(). (phase() + BAR_ANGLE_DEG
    # is NOT the bias - it is 27 deg minus the measured angle, which only coincides with the bias
    # when phi_true happens to be 27. Verified against the density m = 2 axis measured directly in
    # the galactocentric frame.) The bias a real measurement carries is phi_measured - phi_true, and
    # it is NEGATIVE: line-of-sight smearing stretches the galaxy along the Sun-GC direction, pulling
    # the apparent major axis towards that line, so the bar always looks MORE END-ON than it is.
    return {k: (-float(np.mean(v)) - phi_true, float(np.std(v, ddof=1) / np.sqrt(n_real)))
            for k, v in got.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snaps", nargs="+", default=SNAPS)
    ap.add_argument("--fracs", nargs="+", type=float, default=FRACS)
    ap.add_argument("--n-real", type=int, default=10)
    ap.add_argument("--angle", type=float, default=TRUE_ANGLE)
    ap.add_argument("--sigma-ref", type=float, default=0.125,
                    help="which mock's selection tables to use (the survey selection, not the error)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--replot", action="store_true", help="re-draw from the saved csv, no re-measuring")
    ap.add_argument("--scaled-window", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                    help="fit over LO-HI times each galaxy's own R_peak instead of a fixed 1-4 kpc")
    ap.add_argument("--kin-window", type=float, default=None, metavar="F",
                    help="define the bar region from the data's own kinematic m=2 amplitude: keep "
                         "rings above F times its peak (e.g. 0.7). Needs no prior bar length.")
    args = ap.parse_args()
    set_style()
    tb = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
    om = {g: float(tb[np.isclose(tb.snap_time, t_of(g))]["pattern_speed_fd"].iloc[0]) for g in args.snaps}

    out = RESULTS_DIR / f"bar_angle_drift_multi{args.tag}.csv"
    if args.replot and out.exists():
        D = pd.read_csv(out)
        rows = None
    else:
        rows = []
    for gi, gn in enumerate(args.snaps if rows is not None else []):
        for i, fr in enumerate(args.fracs):
            win = (f"kin{args.kin_window:g}" if args.kin_window
                   else (tuple(np.array(args.scaled_window) * R_PEAK[gn])
                         if args.scaled_window else None))
            r = measure(gn, args.angle, fr, args.n_real,
                        seed0=10_000 + 7919 * i + 104_729 * gi, sigma_ref=args.sigma_ref,
                        window=win)
            rows.append(dict(galaxy=gn, t=t_of(gn), Omega_b=om[gn], frac=fr,
                             **{f"{k}_{s}": r[k][j] for k in EST for j, s in enumerate(("mean", "err"))}))
            print(f"  {gn:14s} frac {fr:5.3f}: " + "  ".join(
                f"{k} {r[k][0]:+7.2f}+/-{r[k][1]:.2f}" for k in EST), flush=True)
    if rows is not None:
        D = pd.DataFrame(rows)
        D.to_csv(out, index=False)
        print(f"\nwrote {rel(out)}")

    cols = [BLUE_RAMP[5], "#1baf7a", "#d4a017", ACCENT]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)
    for ax, k in zip(axes, EST):
        ax.axhline(0, color=INK_MUTED, lw=1, ls="--")
        ax.axvline(12.5, color=INK_MUTED, lw=1, ls=":")
        for c, gn in zip(cols, args.snaps):
            d = D[D.galaxy == gn]
            ax.errorbar(d.frac * 100, d[f"{k}_mean"], yerr=d[f"{k}_err"], marker="o", ms=5, lw=1.5,
                        capsize=3, color=c, label=f"t = {t_of(gn):g}  ($\\Omega_b$ {om[gn]:.1f})")
        at = D[np.isclose(D.frac, 0.125)][f"{k}_mean"].values
        ax.set_title(f"{EST[k]}\nat 12.5%: {at.mean():+.2f}$^\\circ$ $\\pm$ {at.std(ddof=1):.2f}$^\\circ$ "
                     f"across snapshots", color=INK, fontsize=10.5)
        ax.set_xlabel("distance error [%]")
    axes[0].set_ylabel(f"measured $-$ true bar angle [deg]   (truth {args.angle:g}$^\\circ$)"
                       "\n(negative = bar looks more end-on)")
    axes[0].legend(frameon=False, fontsize=9, loc="upper left")
    axes[1].text(12.7, axes[1].get_ylim()[0] + 0.5, "our Mira error", fontsize=8, color=INK_2)
    fig.suptitle("Bar angle from the kinematics alone — bias against distance error, "
                 f"{len(args.snaps)} snapshots, amplitude-weighted m = 2 phase over "
                 + (f"rings above {args.kin_window:.0%} of the kinematic $m=2$ peak"
                    if args.kin_window else
                    (f"{args.scaled_window[0]:g}–{args.scaled_window[1]:g} $R_{{\\rm peak}}$"
                     if args.scaled_window else "1–4 kpc")),
                 color=INK, y=0.99, fontsize=11)
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_drift_multi{args.tag}")

    print(f"\n{'estimator':28s} {'bias at 12.5%':>16s} {'spread':>8s}")
    for k in EST:
        at = D[np.isclose(D.frac, 0.125)][f"{k}_mean"].values
        print(f"{EST[k]:28s} {at.mean():+9.2f} deg    {at.std(ddof=1):6.2f} deg")


if __name__ == "__main__":
    main()

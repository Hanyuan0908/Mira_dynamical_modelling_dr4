"""Calibrate the bar angle on one snapshot, apply it to others: which estimator transfers better?

The bar angle is not directly observable once distances are uncertain - the measured m = 2 phase is
a biased, compressed version of the truth. The practical question is therefore not "how big is the
bias" but "does a relation calibrated on one galaxy work on a different one", because that is what
applying a simulation calibration to the Milky Way amounts to.

Procedure, at a 12.5% distance error throughout:
  1. rotate each galaxy so its bar sits at a chosen angle to the Sun-GC line, observe it with the
     distance error and the DR3-Mira-like selection, and measure the m = 2 phase;
  2. fit phi_obs = a phi_true + b on the CALIBRATOR snapshot;
  3. invert that relation on every other snapshot and compare the inferred angle with the truth.

Four combinations are run: {<v_R>, <|v_R|/v_tot>} x {fixed 1-4 kpc, rings above a fraction of the
kinematic m = 2 peak}. The kinematic window needs no prior knowledge of the bar length, which is
what makes it usable on real data.

Usage: python test_bar_angle_transfer_12p5.py [--calib nbody] [--n-real 8] [--kin 0.7]
Table: results/bar_angle_transfer_12p5.csv;  figure: figures/test2_bar_angle_transfer_12p5.png
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
from bar_angle_vecsum_verification import ring_coeffs, Rc

from scan_bar_angle_windows import NBODY, AURIGA

# 12 snapshots of the one isolated N-body galaxy, plus the 13 non-interacting barred Auriga
# haloes. The N-body sequence alone spans a narrow range of bar strength, which is the axis the
# residual bias runs along (plot_bar_angle_amplitude_correction); Auriga widens it.
SNAPS = NBODY + AURIGA
ANGLES = [15.0, 20.0, 25.0, 27.0, 30.0, 35.0]
FRAC = 0.125
EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}
WINS = {"fixed": "fixed 1–4 kpc", "kin": "kinematic window"}


def is_auriga(g):
    return str(g).startswith("au")


def t_of(g):
    """Snapshot model time; each Auriga halo is a separate galaxy, flagged with nan."""
    if is_auriga(g):
        return float("nan")
    return 7.0 if g == "nbody" else float(g.split("_t")[1])


def label(g):
    return f"Au{str(g)[2:]}" if is_auriga(g) else f"t = {t_of(g):g}"


def phases(w, kin):
    """m = 2 phase for both estimators and both windows, from one pass of ring_coeffs each."""
    out = {}
    for which in EST:
        a2, b2 = ring_coeffs(w, which)
        ok = np.isfinite(a2) & np.isfinite(b2)
        amp = np.where(ok, np.hypot(a2, b2), np.nan)
        off = 45.0 if which == "vR" else 0.0
        for wn, m in (("fixed", ok & (Rc > 1) & (Rc < 4)),
                      ("kin", ok & (amp > kin * np.nanmax(amp)))):
            if not m.any():
                out[(which, wn)] = np.nan
                continue
            z = np.nansum(a2[m] + 1j * b2[m])
            # phase() returns the bar axis as an angle in the input frame, which is MINUS the bar
            # angle measured from the Sun-GC line. So phi_measured = -phase(), with no further
            # offset. Verified against the density m = 2 axis measured directly (bar at -26.9 deg
            # in the mw frame -> phi_measured = +26.9 for a perfect measurement).
            out[(which, wn)] = -float((0.5 * np.degrees(np.angle(z)) - off - BAR_ANGLE_DEG
                                       + 90) % 180 - 90)
    return out


def measure(gname, phi_true, n_real, seed0, kin):
    """Mean observed phase over n_real distance-error realisations."""
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{FRAC * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - phi_true)).T)
    sig = FRAC * 5 / np.log(10)
    acc = {k: [] for k in [(e, w) for e in EST for w in WINS]}
    for s in range(n_real):
        d = observe.scatter_distance_modulus(o["d"], sig, np.random.default_rng(seed0 + s))
        w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
        S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
        sel = w[np.random.default_rng(seed0 + 500_000 + s).random(len(S)) < np.minimum(S, 1.0)]
        p = phases(sel, kin)
        for k in acc:
            acc[k].append(p[k])
    return {k: (float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(n_real))) for k, v in acc.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snaps", nargs="+", default=SNAPS)
    ap.add_argument("--calib", default="nbody", help="snapshot the relation is fitted on")
    ap.add_argument("--angles", nargs="+", type=float, default=ANGLES)
    ap.add_argument("--n-real", type=int, default=8)
    ap.add_argument("--kin", type=float, default=0.7)
    ap.add_argument("--replot", action="store_true", help="re-draw from the saved csv")
    args = ap.parse_args()
    set_style()
    tb = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
    au_tb = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")

    out = RESULTS_DIR / "bar_angle_transfer_12p5.csv"
    if args.replot and out.exists():
        D = pd.read_csv(out)
        rows = None
    else:
        rows = []
    for gi, gn in enumerate(args.snaps if rows is not None else []):
        for ai, a in enumerate(args.angles):
            r = measure(gn, a, args.n_real, 20_000 + 7919 * ai + 104_729 * gi, args.kin)
            if is_auriga(gn):
                h = int(str(gn)[2:])
                om = float(au_tb.loc[h, "Omega_p"]) if h in au_tb.index else float("nan")
            else:
                om = float(tb[np.isclose(tb.snap_time, t_of(gn))]["pattern_speed_fd"].iloc[0])
            rec = dict(galaxy=gn, t=t_of(gn), phi_true=a, Omega_b=om)
            for (e, w), (m, s) in r.items():
                rec[f"{e}_{w}"] = m
                rec[f"{e}_{w}_err"] = s
            rows.append(rec)
            print(f"  {gn:14s} phi_true {a:5.1f}: " +
                  "  ".join(f"{e}/{w} {r[(e, w)][0]:6.2f}" for e in EST for w in WINS), flush=True)
    if rows is not None:
        D = pd.DataFrame(rows)
        if out.exists():                       # keep galaxies measured on an earlier run
            old = pd.read_csv(out)
            D = pd.concat([old[~old.galaxy.isin(D.galaxy.unique())], D], ignore_index=True)
            D = D.sort_values(["galaxy", "phi_true"]).reset_index(drop=True)
        D.to_csv(out, index=False)
        print(f"\nwrote {rel(out)}\n")

    # Leave-one-out over the calibrator: fit on each snapshot in turn, transfer to the other 11.
    # This measures the systematic a user would incur without knowing which galaxy to calibrate on.
    print("leave-one-out over the calibrator: every snapshot used as the calibrator in turn\n")
    print(f"{'estimator / window':34s} {'best calib':>14s} {'worst calib':>14s} "
          f"{'mean':>7s} {'worst':>7s}")
    print("-" * 82)
    loo = {}
    for e in EST:
        for w in WINS:
            c = f"{e}_{w}"
            per = []
            for cg in args.snaps:
                cc = D[D.galaxy == cg]
                A = np.polyfit(cc.phi_true, cc[c], 1)
                err = []
                for gn in args.snaps:
                    if gn == cg:
                        continue
                    d = D[D.galaxy == gn]
                    err.extend((d[c].values - A[1]) / A[0] - d.phi_true.values)
                per.append(np.sqrt(np.mean(np.array(err) ** 2)))
            per = np.array(per)
            loo[(e, w)] = per
            print(f"{EST[e] + ' / ' + WINS[w]:34s} "
                  f"{label(args.snaps[int(per.argmin())]):>14s} "
                  f"{label(args.snaps[int(per.argmax())]):>14s} "
                  f"{per.mean():7.2f} {per.max():7.2f}")
    print("-" * 82)
    bl = min(loo, key=lambda k: loo[k].mean())
    print(f"  best on the leave-one-out mean: {EST[bl[0]]} / {WINS[bl[1]]} -> "
          f"{loo[bl].mean():.2f} deg systematic\n")

    cal = D[D.galaxy == args.calib]
    print(f"calibration fitted on {label(args.calib)}, "
          f"applied to the other {len(args.snaps) - 1} snapshots\n")
    print(f"{'estimator / window':34s} {'slope':>7s} {'intcpt':>7s} | "
          f"{'mean err':>9s} {'rms':>6s} {'max|err|':>9s}")
    print("-" * 82)
    summary = {}
    for e in EST:
        for w in WINS:
            c = f"{e}_{w}"
            A = np.polyfit(cal.phi_true, cal[c], 1)
            err = []
            for gn in args.snaps:
                if gn == args.calib:
                    continue
                d = D[D.galaxy == gn]
                inferred = (d[c].values - A[1]) / A[0]
                err.extend(inferred - d.phi_true.values)
            err = np.array(err)
            summary[(e, w)] = (A, err)
            print(f"{EST[e] + ' / ' + WINS[w]:34s} {A[0]:7.3f} {A[1]:7.2f} | "
                  f"{err.mean():+9.2f} {np.sqrt((err**2).mean()):6.2f} {np.abs(err).max():9.2f}")
    print("-" * 82)
    best = min(summary, key=lambda k: np.sqrt((summary[k][1]**2).mean()))
    print(f"  best: {EST[best[0]]} with the {WINS[best[1]]}, "
          f"rms {np.sqrt((summary[best][1]**2).mean()):.2f} deg")

    fig, AX = plt.subplots(2, 2, figsize=(11.4, 8.4))
    # two suites, two colour families: viridis along the N-body time sequence (it IS ordered),
    # a warm ramp across the Auriga galaxies (they are not ordered by anything). Au18 stays black
    # because it is the galaxy the earlier blind test held out.
    nb = [g for g in args.snaps if not is_auriga(g)]
    au = [g for g in args.snaps if is_auriga(g)]
    cmap, amap = plt.get_cmap("viridis"), plt.get_cmap("autumn")
    cols = {g: cmap(i / max(len(nb) - 1, 1)) for i, g in enumerate(nb)}
    cols.update({g: amap(0.82 * i / max(len(au) - 1, 1)) for i, g in enumerate(au)})
    if "au18" in cols:
        cols["au18"] = INK
    for i, e in enumerate(EST):
        for j, w in enumerate(WINS):
            ax = AX[i, j]
            A, err = summary[(e, w)]
            x = np.linspace(min(args.angles) - 2, max(args.angles) + 2, 50)
            ax.plot(x, x, color=INK_MUTED, lw=1.0, ls=":", label="1:1 (no bias)")
            ax.plot(x, np.polyval(A, x), color=INK, lw=1.5, ls="--",
                    label=f"fit on {label(args.calib)}:  "
                          f"{A[0]:.2f}$\\phi$ + {A[1]:.1f}")
            for gn in args.snaps:
                d = D[D.galaxy == gn]
                isau, ref = is_auriga(gn), gn == "au18"
                # N-body: thin circles. Auriga: squares, so the two suites read apart even
                # where the colours overlap. Au18 heavy black.
                # Each suite is labelled in one panel only: N-body in (0,0), Auriga in (0,1).
                legend_panel = (i == 0) and (j == (1 if isau else 0))
                ax.errorbar(d.phi_true, d[f"{e}_{w}"], yerr=d[f"{e}_{w}_err"],
                            marker="s" if isau else "o", capsize=2, color=cols[gn],
                            ms=7 if ref else (5 if isau else 4),
                            lw=2.4 if ref else (1.4 if isau else 1.0),
                            zorder=6 if ref else (4 if isau else 2),
                            label=label(gn) if legend_panel else None)
            ax.set_ylim(0, max(args.angles) + 3)
            ax.set_title(f"{EST[e]} / {WINS[w]}\ntransfer rms "
                         f"{np.sqrt((err**2).mean()):.2f}$^\\circ$, max "
                         f"{np.abs(err).max():.2f}$^\\circ$", color=INK, fontsize=10)
            ax.set_xlabel(r"true bar angle $\phi_b$ [deg]")
            ax.set_ylabel("measured bar angle [deg]")
            if i == 0 and j in (0, 1):
                ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left",
                          title="N-body sequence" if j == 0 else "Auriga galaxies",
                          title_fontsize=8)
    fig.suptitle(f"Does a bar-angle calibration transfer between galaxies? "
                 f"{FRAC*100:g}% distance error, fitted on {label(args.calib)} — "
                 f"{len(nb)} N-body snapshots + {len(au)} Auriga galaxies",
                 color=INK, y=0.995)
    fig.tight_layout()
    savefig(fig, "test2_bar_angle_transfer_12p5")


if __name__ == "__main__":
    main()

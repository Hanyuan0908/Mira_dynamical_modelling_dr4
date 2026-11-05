"""Face-on <v_R> and <|v_R|/v_tot> before and after a distance error, for two snapshots.

The bar angle is read from the m = 2 phase of these maps. A distance error displaces every star
along its own line of sight, which mixes azimuths asymmetrically about the Sun-GC line and so
rotates the pattern. This shows what that does, side by side with the unsmeared field, at the two
snapshots that behave most differently: t = 2 (the youngest, strongest, shortest bar, which drifts
+19 deg at 12.5%) and t = 5 (a settled bar, +13 deg).

Everything is in the bar frame, with the bar along x, so the quadrupole node should lie along x.
The overlaid line is the measured amplitude-weighted m = 2 phase over the window used for the
calibration, which is what a real measurement would report.

Usage: python plot_smeared_velocity_maps.py [--snaps nbody_t2 nbody_t5] [--frac 0.125]
Figure: figures/test2_smeared_velocity_maps.png
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, INK, INK_2, INK_MUTED
from galaxies import galaxy
from nbody import SNAPSHOT_DIR
import observe, selection, steady_state_likelihood as ssl
from bar_angle_vecsum_verification import phase

R_PEAK = {"nbody_t2": 2.38, "nbody_t3.625": 2.62, "nbody_t5": 3.12, "nbody": 3.38}


def is_auriga(g):
    return str(g).startswith("au")


def t_of(g):
    return float("nan") if is_auriga(g) else (7.0 if g == "nbody" else float(g.split("_t")[1]))


def galaxy_tag(gname):
    """The bold caption on each galaxy's first panel: what distinguishes this galaxy.

    For the N-body sequence that is the epoch; for an Auriga halo it is the quantities the
    outlier analysis turns on - pattern speed, corotation radius and bar strength.
    """
    if not is_auriga(gname):
        tb = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
        om = float(tb[np.isclose(tb.snap_time, t_of(gname))]["pattern_speed_fd"].iloc[0])
        return (f"t = {t_of(gname):g}\n$\\Omega_b$ = {om:.1f}   "
                f"$R_{{\\rm peak}}$ = {R_PEAK[gname]:.2f} kpc")
    h = int(str(gname)[2:])
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    S = pd.read_csv(RESULTS_DIR / "auriga_bar_survey.csv").set_index("halo")
    return (f"Au{h}\n$\\Omega_b$ = {Q.loc[h, 'Omega_p']:.1f}   "
            f"$R_{{\\rm CR}}$ = {Q.loc[h, 'R_CR']:.2f}\n"
            f"$R_{{\\rm bar}}$ = {S.loc[h, 'R_bar']:.2f}   "
            f"$A_2^{{\\rm max}}$ = {S.loc[h, 'A2_max']:.2f}")


def sample(gname, frac, seed, sigma_ref=0.125):
    """Observed, selection-applied stars in the bar frame, with and without the distance error."""
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{sigma_ref * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*W.T)
    rng = np.random.default_rng(seed)
    d = observe.scatter_distance_modulus(o["d"], frac * 5 / np.log(10), rng) if frac > 0 else o["d"]
    w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
    S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
    keep = np.random.default_rng(seed + 500_000).random(len(S)) < np.minimum(S, 1.0)
    return w[keep]                     # MW frame: phase() rotates to the bar frame itself


def fields(w):
    R = np.hypot(w[:, 0], w[:, 1])
    vR = (w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / np.maximum(R, 1e-9)
    q = np.abs(vR) / np.maximum(np.linalg.norm(w[:, 3:], axis=1), 1e-9)
    return vR, q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snaps", nargs="+", default=["nbody_t2", "nbody_t5"])
    ap.add_argument("--frac", type=float, default=0.125)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--tag", default="")
    ap.add_argument("--window", default="1-4",
                    help="m=2 window for the measured phase: 'LO-HI' kpc, or 'kin<F>' for the "
                         "self-scaled window the calibration actually uses (e.g. kin0.4)")
    args = ap.parse_args()
    set_style()
    win = (args.window if args.window.startswith("kin")
           else tuple(float(v) for v in args.window.split("-")))
    wlab = (f"self-scaled kinematic window, f = {args.window[3:]}"
            if args.window.startswith("kin") else f"fixed {args.window} kpc window")
    b = np.linspace(-6, 6, 73)

    nrow = 2 * len(args.snaps)
    fig, AX = plt.subplots(nrow, 2, figsize=(9.6, 4.9 * len(args.snaps)))
    for si, gn in enumerate(args.snaps):
        for fi, fr in enumerate((0.0, args.frac)):
            w_mw = sample(gn, fr, args.seed)
            w = ssl.rotate_z(w_mw, np.radians(BAR_ANGLE_DEG))    # for plotting: bar along x
            vR, q = fields(w)
            for ki, (val, lab, cm, lim) in enumerate(
                    ((vR, r"$\langle v_R\rangle$", "RdBu_r", 40),
                     (q, r"$\langle |v_R|/v_{\rm tot}\rangle$", "viridis", None))):
                ax = AX[2 * si + ki, fi]
                n, *_ = np.histogram2d(w[:, 0], w[:, 1], bins=[b, b])
                s, *_ = np.histogram2d(w[:, 0], w[:, 1], bins=[b, b], weights=val)
                M = np.where(n > 8, s / np.maximum(n, 1), np.nan)
                kw = dict(norm=TwoSlopeNorm(0, -lim, lim)) if lim else dict(vmin=0.05, vmax=0.55)
                im = ax.imshow(M.T, origin="lower", extent=[-6, 6, -6, 6], cmap=cm, **kw)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
                # the measured m = 2 phase, which is what the calibration uses
                ph = phase(w_mw, "vR" if ki == 0 else "absvR", win) + BAR_ANGLE_DEG  # measured - true
                t = np.radians(ph)
                ax.plot([-5.5 * np.cos(t), 5.5 * np.cos(t)], [-5.5 * np.sin(t), 5.5 * np.sin(t)],
                        color="k", lw=1.6, ls="--")
                ax.axhline(0, color=INK_MUTED, lw=0.8)
                ax.set_xticks([-6, -3, 0, 3, 6]); ax.set_yticks([-6, -3, 0, 3, 6])
                ax.set_title(f"{lab}   {'no error' if fr == 0 else f'{fr*100:g}% distance error'}"
                             f"\nmeasured phase {ph:+.1f}$^\\circ$ (truth 0$^\\circ$)",
                             fontsize=9.5, color=INK)
                if fi == 0:
                    ax.set_ylabel("y [kpc]  (bar along x)", fontsize=9)
                if 2 * si + ki == nrow - 1:
                    ax.set_xlabel("x [kpc]", fontsize=9)
    # One caption per galaxy, inside its first panel on a white patch. Above the panel it
    # collides with the two-line titles, and at the old (-5.6, 5.0) it runs into the colourbar
    # as soon as the caption carries more than two quantities.
    for si, gn in enumerate(args.snaps):
        AX[2 * si, 0].text(0.025, 0.975, galaxy_tag(gn), transform=AX[2 * si, 0].transAxes,
                           fontsize=8.5, color=INK, fontweight="bold", ha="left", va="top",
                           bbox=dict(facecolor="white", alpha=0.82, edgecolor="none", pad=2.5))
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.suptitle("What the distance error does to the bar's kinematic signature "
                 f"(bar frame, selection applied; {wlab})", color=INK, y=0.999, fontsize=11)
    savefig(fig, f"test2_smeared_velocity_maps{args.tag}")


if __name__ == "__main__":
    main()

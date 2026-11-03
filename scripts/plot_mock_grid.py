"""Face-on view of every Test 2 mock, in pattern-speed order, against the DR3 Mira target.

One panel per test galaxy: the observed (distance-smeared, selection-applied) stars seen face-on in
the bar frame, with the snapshot's Omega_b and bar strength. The bottom row summarises how well each
mock reproduces the DR3 Miras - sky-map correlation and the medians of distance, Galactocentric R and
|b| - which is what `make_mock_selection.py` prints one mock at a time.

Usage: python plot_mock_grid.py [--frac 0.125] [--tag NAME]
Figure: figures/test2_mock_grid_<frac>pct<tag>.png
"""
import argparse
import math

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, MIRA_DR3_FILE, R0_KPC, set_style, savefig, rel, INK, INK_2, CMAP_SEQ_W
from nbody import SNAPSHOT_DIR
from galaxies import nbody_name

TIMES = [2.0, 2.375, 2.625, 3.0, 3.375, 3.625, 4.0, 4.5, 5.0, 5.5, 6.25, 7.0]
NEW = {2.375, 2.625, 3.375, 3.625, 4.5, 5.0, 6.25}
BAR_ANGLE_DEG = 27.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frac", type=float, default=0.125)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()
    sigma = args.frac * 5 / math.log(10)
    tab = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")

    def prop(t, col):
        return float(tab.loc[np.isclose(tab.snap_time, t), col].iloc[0])

    rows, found = [], []
    for t in TIMES:
        g = nbody_name(t)
        pre = "nbody" if g == "nbody" else g
        f = MOCK_DIR / f"{pre}_mock_selection_sigmaDM{sigma:.3f}_N50k.hdf5"
        if not f.exists():
            print(f"  missing: {rel(f)}")
            continue
        with h5py.File(f, "r") as h:
            x, y = h["x"][:], h["y"][:]
            d, b = h["d"][:], h["b"][:]
            R = np.hypot(x, y)
        found.append(t)
        rows.append(dict(t=t, om=prop(t, "pattern_speed_fd"), a2=prop(t, "bar_strength"),
                         n=len(x), d_med=np.median(d), R_med=np.median(R), b_med=np.median(np.abs(b)),
                         x=x, y=y))
    if not rows:
        raise SystemExit("no mocks found; build them with make_mock_selection.py")

    # the same target cut make_mock_selection.py uses: DM_err < 0.4 and R < 8 kpc
    mira = pd.read_csv(MIRA_DR3_FILE,
                       usecols=["l", "b", "DM_err", "luminosity_distance", "x", "y", "z"]).dropna()
    mira = mira[(mira.DM_err < 0.4) & (np.hypot(mira.x, mira.y) < 8.0)]
    tgt_d = float(np.median(mira.luminosity_distance))
    tgt_R = float(np.median(np.hypot(mira.x, mira.y)))
    tgt_b = float(np.median(np.abs(mira.b)))

    ncol = 6
    nrow = int(np.ceil(len(rows) / ncol))
    fig = plt.figure(figsize=(2.35 * ncol, 2.65 * nrow + 2.9))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[1] * nrow + [1.25], hspace=0.42, wspace=0.12)
    bins = np.linspace(-8, 8, 90)
    for k, r in enumerate(rows):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        c = math.cos(math.radians(-BAR_ANGLE_DEG)); s = math.sin(math.radians(-BAR_ANGLE_DEG))
        xb, yb = c * r["x"] - s * r["y"], s * r["x"] + c * r["y"]
        H, *_ = np.histogram2d(xb, yb, bins=[bins, bins])
        ax.imshow(np.log10(H.T + 1), origin="lower", extent=[-8, 8, -8, 8], cmap=CMAP_SEQ_W,
                  aspect="equal", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        new = r["t"] in NEW
        ax.set_title(f"t = {r['t']:g}   $\\Omega_b$ = {r['om']:.1f}", fontsize=9.5,
                     color="#b4480f" if new else INK, fontweight="bold" if new else "normal")
        for sp in ax.spines.values():
            sp.set_color("#eb6834" if new else "#cccccc")
            sp.set_linewidth(1.8 if new else 0.8)

    bx = fig.add_subplot(gs[nrow, :])
    t = [r["t"] for r in rows]
    xx = np.arange(len(rows))
    for key, tgt, lbl, col in (("d_med", tgt_d, "median distance [kpc]", "#256abf"),
                               ("R_med", tgt_R, "median R [kpc]", "#1baf7a"),
                               ("b_med", tgt_b, "median |b| [deg]", "#d4a017")):
        v = np.array([r[key] for r in rows])
        bx.plot(xx, v / tgt, "o-", ms=5, lw=1.1, color=col, label=f"{lbl}  (DR3 {tgt:.2f})")
    bx.axhline(1, color=INK_2, lw=1, ls="--")
    bx.axhspan(0.98, 1.02, color=INK_2, alpha=0.1, lw=0)
    bx.set_xticks(xx)
    bx.set_xticklabels([f"{v:g}" for v in t])
    bx.set_xlabel("snapshot model time (orange = new)")
    bx.set_ylabel("mock / DR3 Miras")
    bx.set_ylim(0.9, 1.1)
    bx.legend(frameon=False, fontsize=9, ncol=3, loc="upper center")
    for lab, tt in zip(bx.get_xticklabels(), t):
        if tt in NEW:
            lab.set_color("#b4480f"); lab.set_fontweight("bold")

    fig.suptitle(f"Test 2 mocks at a {args.frac*100:g}% distance error "
                 f"($\\sigma_{{DM}}$ = {sigma:.3f} mag), face-on in the bar frame, "
                 f"against the DR3 Miras", color=INK, y=0.995)
    savefig(fig, f"mock_validation/test2_mock_grid_{args.frac*100:g}pct{args.tag}")

    print(f"\n{'t':>6s} {'Omega_b':>8s} {'A2':>6s} {'N':>8s} {'d_med':>7s} {'R_med':>7s} {'|b|_med':>8s}  new")
    for r in rows:
        print(f"{r['t']:6g} {r['om']:8.2f} {r['a2']:6.3f} {r['n']:8,d} {r['d_med']:7.2f} "
              f"{r['R_med']:7.2f} {r['b_med']:8.2f}  {'yes' if r['t'] in NEW else ''}")
    print(f"{'DR3':>6s} {'':8s} {'':6s} {len(mira):8,d} {tgt_d:7.2f} {tgt_R:7.2f} {tgt_b:8.2f}")


if __name__ == "__main__":
    main()

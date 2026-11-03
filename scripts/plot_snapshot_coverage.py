"""Where the test galaxies sit along the N-body run's pattern-speed track.

Top: Omega_b against snapshot time for every snapshot in 2 <= t <= 8, with a smooth fit. The spread
of the points about that fit is the finite-difference column's own noise (~0.66 km/s/kpc), which sets
a floor on how well any recovery bias can be measured snapshot by snapshot.
Bottom: the same galaxies projected onto the Omega_b axis alone - what the 1:1 recovery plot samples -
before and after the six fill-in runs.

Usage: python plot_snapshot_coverage.py [--tag NAME]
Figure: figures/test2_snapshot_coverage<tag>.png
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from nbody import SNAPSHOT_DIR, available_times

DONE = [2.0, 3.0, 4.0, 5.5, 7.0]
FILL = [2.375, 2.625, 3.375, 3.625, 4.5, 5.0, 6.25]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    tab = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv").sort_values("snap_time")
    have = set(np.round(available_times(), 4))
    d = tab[tab.snap_time.round(4).isin(have) & (tab.snap_time >= 2) & (tab.snap_time <= 8)
            & (tab.pattern_speed_fd > 0)]
    x, y = d.snap_time.values, d.pattern_speed_fd.values
    fit = np.poly1d(np.polyfit(x, y, 3))
    rms = (y - fit(x)).std(ddof=4)
    om = lambda t: float(d.loc[np.isclose(d.snap_time, t), "pattern_speed_fd"].iloc[0])

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.5, 7.4), height_ratios=[2.4, 1],
                                 gridspec_kw=dict(hspace=0.42))

    xs = np.linspace(2, 8, 300)
    ax.fill_between(xs, fit(xs) - rms, fit(xs) + rms, color=INK_MUTED, alpha=0.16, lw=0,
                    label=f"smooth fit $\\pm$ {rms:.2f} (the fd column's own scatter)")
    ax.plot(xs, fit(xs), color=INK_MUTED, lw=1.2, ls="--")
    ax.plot(x, y, "o", ms=4, color=INK_MUTED, mfc="none", lw=0, label="all 49 snapshots, 2 ≤ t ≤ 8")
    ax.plot(DONE, [om(t) for t in DONE], "o", ms=11, color=BLUE_RAMP[5], lw=0,
            label="run (5)", zorder=5)
    ax.plot(FILL, [om(t) for t in FILL], "D", ms=10, color=ACCENT, lw=0,
            label="fill-in (7)", zorder=5)
    for t in DONE + FILL:
        ax.annotate(f"{t:g}", (t, om(t)), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=8.5, color=INK_2)
    ax.set_xlabel("snapshot model time  (1 unit = 0.978 Gyr)")
    ax.set_ylabel(r"$\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    ax.set_title("the N-body run's pattern-speed track, and which snapshots are tested", color=INK)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_xlim(1.85, 8.15)

    # the Omega_b axis alone - what the recovery plot actually samples
    for j, (lbl, pts, colour, marker) in enumerate(
            [("before: 5 runs", DONE, BLUE_RAMP[5], "o"),
             ("after: 12 runs", DONE + FILL, None, None)]):
        yy = 1 - j
        v = sorted(om(t) for t in pts)
        gaps = np.diff(v)
        bx.hlines(yy, 24, 45, color=INK_MUTED, lw=0.8, zorder=1)
        if j == 0:
            bx.plot(v, [yy] * len(v), "o", ms=10, color=BLUE_RAMP[5], zorder=3)
        else:
            bx.plot([om(t) for t in DONE], [yy] * len(DONE), "o", ms=10, color=BLUE_RAMP[5], zorder=3)
            bx.plot([om(t) for t in FILL], [yy] * len(FILL), "D", ms=9, color=ACCENT, zorder=3)
        bx.annotate(f"{lbl}\nlargest gap {gaps.max():.1f}", (23.4, yy), ha="right", va="center",
                    fontsize=9, color=INK_2)
    bx.set_ylim(-0.6, 1.6)
    bx.set_xlim(18, 45.6)
    bx.set_yticks([])
    for s in ("left", "right", "top"):
        bx.spines[s].set_visible(False)
    bx.set_xlabel(r"true $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]  —  coverage of the 1:1 recovery plot")
    savefig(fig, f"test2_snapshot_coverage{args.tag}")

    print(f"{'t':>5s} {'Omega_b':>8s} {'A2':>6s}   set")
    for t in sorted(DONE + FILL):
        a = float(d.loc[np.isclose(d.snap_time, t), "bar_strength"].iloc[0])
        print(f"{t:5g} {om(t):8.2f} {a:6.3f}   {'run' if t in DONE else 'fill-in'}")
    v0 = sorted(om(t) for t in DONE)
    v1 = sorted(om(t) for t in DONE + FILL)
    print(f"\nlargest gap: {np.diff(v0).max():.2f} -> {np.diff(v1).max():.2f} km/s/kpc")
    print(f"range: {v1[0]:.2f}-{v1[-1]:.2f} (factor {v1[-1]/v1[0]:.2f}), unchanged by the fill-in")
    print(f"fd column scatter about a smooth fit: {rms:.2f} km/s/kpc ({rms/y.mean()*100:.1f}%)")


if __name__ == "__main__":
    main()

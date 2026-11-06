"""A2/A0(R) of the Portail+17 / Sormani+22 Milky Way bar, against the 25 modelled galaxies.

The trial-potential family in this project is built by rescaling the m != 0 part of the
Portail+17 / Sormani+22 AGAMA model (README), so it is worth knowing where that model's own bar
strength profile sits relative to the simulated galaxies the method has been tested on.

`Portail17.ini` ships with AGAMA. It is a FULL Milky Way model, not a bar-only component -
v_c(8 kpc) = 241 km/s and M(<20 kpc) = 1.7e11 Msun - so its midplane density includes the halo
and disc exactly as the mock potentials do, and A2/A0 is diluted the same way. That makes the
comparison like-for-like without any rescaling.

A2/A0 is measured identically to plot_bar_strength_profile_recovery: the m = 2 Fourier amplitude
of the midplane density on the same radial grid and the same 64 azimuths. Only the amplitude is
used, so the model's bar orientation does not matter.

Usage: python plot_portail17_bar_strength.py
Figure: figures/test2_portail17_bar_strength.png
Table:  results/portail17_bar_strength.csv
"""
import os

import numpy as np
import pandas as pd
import agama
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl        # noqa: F401  (sets AGAMA units)
from plot_bar_strength_profile_recovery import R_GRID, _RR, _PP, _PTS


def a2_over_a0(pot):
    rho = pot.density(_PTS).reshape(_RR.shape)
    a0 = rho.mean(axis=1)
    a2 = 2 * np.abs((rho * np.exp(-2j * _PP)).mean(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(a0 > 0, a2 / a0, np.nan)


def main():
    set_style()
    data = os.path.join(os.path.dirname(agama.__file__), "data")
    mw = agama.Potential(os.path.join(data, "Portail17.ini"))
    mw_a2 = a2_over_a0(mw)
    vc8 = float(np.sqrt(-mw.force([[8.0, 0, 0]])[0][0] * 8))

    d = np.load(RESULTS_DIR / "bar_strength_profile_recovery.npz")
    runs = sorted({k.rsplit("_", 1)[0] for k in d.files if k != "R"})
    truth = {r: d[f"{r}_truth"] for r in runs}
    au = np.array([v for r, v in truth.items() if "_au" in r])
    nb = np.array([v for r, v in truth.items() if "_au" not in r])
    allg = np.vstack([au, nb])

    pd.DataFrame(dict(R=R_GRID, A2_over_A0=mw_a2)).to_csv(
        RESULTS_DIR / "portail17_bar_strength.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    for i, c in enumerate(allg):
        ax.plot(R_GRID, c, "-", color=INK_MUTED, lw=0.7, alpha=0.5, zorder=2,
                label="25 modelled galaxies (truth)" if i == 0 else None)
    ax.fill_between(R_GRID, np.nanpercentile(allg, 16, 0), np.nanpercentile(allg, 84, 0),
                    color=BLUE_RAMP[2], alpha=0.5, lw=0, zorder=1,
                    label="16–84% of the 25")
    ax.plot(R_GRID, np.nanmedian(allg, 0), "-", color=BLUE_RAMP[6], lw=2.0, zorder=4,
            label="median of the 25")
    ax.plot(R_GRID, np.nanmedian(au, 0), "--", color=BLUE_RAMP[6], lw=1.5, zorder=4,
            label="median, 13 Auriga")
    ax.plot(R_GRID, mw_a2, "-", color=ACCENT, lw=3.0, zorder=6,
            label="Portail+17 / Sormani+22 Milky Way")

    pk = np.nanargmax(mw_a2)
    ax.plot(R_GRID[pk], mw_a2[pk], "o", color=ACCENT, ms=9, mec=INK, mew=1.2, zorder=7)
    ax.annotate(f"peak {mw_a2[pk]:.2f} at {R_GRID[pk]:.1f} kpc", (R_GRID[pk], mw_a2[pk]),
                xytext=(10, 8), textcoords="offset points", fontsize=9.5, color=ACCENT)
    ax.set_xlabel("R [kpc]")
    ax.set_ylabel(r"$A_2/A_0$ (midplane density, halo included)")
    ax.set_xlim(R_GRID[0], R_GRID[-1])
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title(f"Bar strength profile of the Portail+17 / Sormani+22 Milky Way "
                 f"($v_c$(8 kpc) = {vc8:.0f} km/s)\nagainst the 25 galaxies the method is "
                 f"tested on, measured identically", color=INK, fontsize=11)
    fig.tight_layout()
    savefig(fig, "test2_portail17_bar_strength")

    med = np.nanmedian(allg, 0)
    print(f"  Portail17 v_c(8 kpc) = {vc8:.1f} km/s")
    print(f"  {'R':>5s} {'MW A2/A0':>9s} {'median of 25':>13s} {'16-84% of 25':>18s} {'MW pctile':>10s}")
    for i in range(0, len(R_GRID), 4):
        lo, hi = np.nanpercentile(allg[:, i], [16, 84])
        pc = 100 * np.mean(allg[:, i] < mw_a2[i])
        print(f"  {R_GRID[i]:5.1f} {mw_a2[i]:9.3f} {med[i]:13.3f} "
              f"{f'{lo:.3f} - {hi:.3f}':>18s} {pc:9.0f}%")
    print(f"\n  MW peak A2/A0 = {mw_a2[pk]:.3f} at R = {R_GRID[pk]:.1f} kpc")
    pks = np.nanmax(allg, axis=1)
    print(f"  peak of the 25: median {np.median(pks):.3f}, range {pks.min():.3f}-{pks.max():.3f}")
    print(f"  the MW's peak is at the {100*np.mean(pks < mw_a2[pk]):.0f}th percentile of the 25")
    print(f"\nwrote {rel(RESULTS_DIR / 'portail17_bar_strength.csv')}")


if __name__ == "__main__":
    main()

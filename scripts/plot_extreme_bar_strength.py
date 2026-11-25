"""Figure for extreme_bar_strength_test.py: chi2 and the library <v_R> quadrupole against bar_frac, all else fixed.
Left: Delta chi2 / T1 (default grid and geo20 |z|<2), OSARG in Sormani+22 at its best point and the steady-state mock in
Hunter+24 at its truth (bar_frac 1). Right: the library's <v_R> m = 2 amplitude (1 < R < 4 kpc, |z| < 1) against
bar_frac, with the data's as a horizontal line.
Usage: python plot_extreme_bar_strength.py   ->  figures/extreme_bar_strength_test_osarg_and_mock.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import RESULTS_DIR, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP

T1 = {"osarg_s22": {"default": 32.1, "geo20 |z|<2": 11.6}, "h24mock": {"default": 25.2, "geo20 |z|<2": 10.6}}
LAB = {"osarg_s22": "OSARG, Sormani+22 (Ω_b 32.06, φ_b 20.5, S 0.95 fixed)", "h24mock": "steady-state mock, Hunter+24 (truth 32.5, 25, 1; true bar_frac 1)"}
COL = {"osarg_s22": BLUE_RAMP[5], "h24mock": ACCENT}


def main():
    set_style()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), layout="constrained")
    for case in ("osarg_s22", "h24mock"):
        d = pd.read_csv(RESULTS_DIR / f"extreme_bar_strength_{case}.csv")
        for g, ls in (("default", "-"), ("geo20 |z|<2", "--")):
            x = d[f"chi2_{g}"].values
            ax[0].plot(d.bar_frac, (x - x.min()) / T1[case][g], ls, marker="o", ms=4, color=COL[case], lw=1.8,
                       label=f"{LAB[case].split(' (')[0]}, {g} grid")
        ax[1].plot(d.bar_frac, d.vR_m2_lib, marker="o", ms=4, color=COL[case], lw=1.8, label=f"library: {LAB[case]}")
        ax[1].axhline(d.vR_m2_data.iloc[0], color=COL[case], lw=1, ls=":", label=f"data ({d.vR_m2_data.iloc[0]:.1f} km/s)")
    ax[0].axvline(1.0, color=INK_2, lw=0.8, ls=":")
    ax[0].set(xlabel="bar_frac (everything else fixed)", ylabel="Δχ² / T1 (from each curve's minimum)", ylim=(0, 50),
              title="χ² against bar strength: no bar is strongly rejected; the mock's minimum is at the truth (1)")
    ax[1].set(xlabel="bar_frac", ylabel="⟨v_R⟩ m=2 amplitude, 1 < R < 4 kpc, |z| < 1 [km/s]",
              title="The library quadrupole saturates near bar_frac 1.5 and then falls")
    for a in ax:
        a.legend(frameon=False, fontsize=7.5)
        a.tick_params(top=False, right=False)
    savefig(fig, "extreme_bar_strength_test_osarg_and_mock")


if __name__ == "__main__":
    main()

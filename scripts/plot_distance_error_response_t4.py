"""Figure: how the pattern-speed signal responds to the distance error (t = 4 mock) compared with OSARG.

Cells of cross_response_stars_vs_potential.py, all with N = 26,630 stars and the bar held at the truth (or the best
point), Omega_b stepped by -30 ... +40%:
  t = 4 mock, 12.5% data & 12.5% library;  20% data & 20% library;  20% data, library smeared as if 10%;
  OSARG in Sormani+22;  the steady-state OSARG-like mock in Hunter+24.
Panel 1: how much worse the DATA fit (Delta chi2, tempered by 2 sqrt(2N)), relative to the mean of the +-0.8% steps
         (removes each cell's reference-library draw). Higher = the data reject that pattern speed more strongly.
Panel 2: how much the PREDICTION changes (D, same units); dotted: each cell's noise floor (D at +-0.8%).
Panel 3: chi2/N at the true (reference) pattern speed: the goodness of fit.
Usage: python plot_distance_error_response_t4.py  ->  figures/distance_error_response_t4_vs_osarg.png
"""
import numpy as np
import matplotlib.pyplot as plt
from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

CELLS = [("t = 4 mock: 12.5% errors, modelled as 12.5%", "t4m12p5_nbody_t4_clusteraxifix", BLUE_RAMP[3], "-", "o"),
         ("t = 4 mock: 20% errors, modelled as 20%", "t4m20_nbody_t4_clusteraxifix", BLUE_RAMP[6], "-", "s"),
         ("t = 4 mock: 20% errors, modelled as 10%", "t4m20as10_nbody_t4_clusteraxifix", "#8e44ad", "--", "D"),
         ("OSARG data, Sormani+22 (9.6%)", "osarg_sormani22", ACCENT, "-", "*"),
         ("steady-state OSARG-like mock, Hunter+24 (9.6%)", "h24mock_mw", INK_MUTED, ":", "^")]


def main():
    set_style()
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.4), gridspec_kw=dict(width_ratios=[1.2, 1.2, 0.8]), layout="constrained")
    chi = []
    for lab, k, c, ls, mk in CELLS:
        d = np.load(RESULTS_DIR / f"cross_response_{k}.npz")
        fr, t = d["fracs"], 2 * np.sqrt(2 * d["ndata"])
        x = d["chi2"] / t
        nb = 0.5 * (x[np.isclose(fr, -0.008)][0] + x[np.isclose(fr, 0.008)][0])
        sel = ~np.isclose(np.abs(fr), 0.008)
        ax[0].plot(100 * fr[sel], (x - nb)[sel], ls, marker=mk, color=c, lw=2, ms=6, label=lab)
        D = d["D"] / t
        fl = 0.5 * (D[np.isclose(fr, -0.008)][0] + D[np.isclose(fr, 0.008)][0])
        ax[1].plot(100 * fr[sel], D[sel], ls, marker=mk, color=c, lw=2, ms=6)
        ax[1].axhline(fl, color=c, lw=0.8, ls=":")
        chi.append((lab, d["chi2"][fr == 0][0] / d["ndata"], c))
    ax[0].axhline(0, color=INK_2, lw=0.8)
    ax[0].set(xlabel="change of Ω_b from the truth / reference [%]", ylabel="how much worse the DATA fit (Δχ², tempered)",
              title="1. Do the data reject a wrong pattern speed?\n(higher = stronger rejection; slower bar on the left)")
    ax[0].legend(frameon=False, fontsize=8, loc="upper center")
    ax[1].set(xlabel="change of Ω_b from the truth / reference [%]", ylabel="how much the PREDICTION changes (D, tempered)",
              title="2. Does the model react to the pattern speed?\n(dotted: each cell's noise floor)")
    y = np.arange(len(chi))
    ax[2].barh(y, [v for _, v, _ in chi], color=[c for _, _, c in chi], height=0.6)
    for i, (_, v, _) in enumerate(chi):
        ax[2].text(v + 0.03, i, f"{v:.2f}", va="center", fontsize=9, color=INK)
    ax[2].set_yticks(y, [l.replace(": ", ":\n", 1) for l, _, _ in chi], fontsize=8)
    ax[2].invert_yaxis()
    ax[2].set(xlabel="χ²/N at the true pattern speed", xlim=(0, 3.8), title="3. Goodness of fit")
    for a in ax:
        a.tick_params(top=False, right=False)
    fig.suptitle("Distance error vs the pattern-speed signal: t = 4 mock (N 26,630, its own corrected family, bar at the "
                 "truth) against OSARG", fontsize=11, x=0.01, ha="left")
    savefig(fig, "distance_error_response_t4_vs_osarg")


if __name__ == "__main__":
    main()

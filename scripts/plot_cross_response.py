"""Figure for cross_response_stars_vs_potential.py: how much the predicted kinematics change with the pattern speed,
for every (star set, family) cell, against the fractional change of |Omega_b| from each family's reference.

Left: the change of the PREDICTION, D (tempered), one colour per family, one line style per star set.
Right: the same for the <v_phi> terms alone (m = 0, 2, 4), which carry most of the pattern-speed information.
The library noise floor (D at +-0.8%) is the grey band.

Usage: python plot_cross_response.py
Figure: figures/cross_response_stars_vs_potential.png
"""
import numpy as np
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

FAMILIES = {"nbody": ("t = 7 family (own rotation curve + own bar)", BLUE_RAMP[6]),
            "hyb_t7axi_mwbar": ("t = 7 rotation curve + Milky Way bar", BLUE_RAMP[3]),
            "hyb_mwaxi_t7bar": ("Milky Way rotation curve + t = 7 bar", "#E8A33D"),
            "mw": ("Milky Way family (Hunter+24 rotation curve + Portail+17 bar)", ACCENT)}
STARS = {"t7cov": ("t = 7 stars (OSARG coverage + N)", "-", "o"), "osarg": ("OSARG stars", "--", "s"),
         "h24mock": ("steady-state OSARG-like mock", ":", "^")}


def main():
    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw=dict(wspace=0.22))
    floors = []
    for fam, (flab, col) in FAMILIES.items():
        for st, (slab, ls, mk) in STARS.items():
            p = RESULTS_DIR / f"cross_response_{st}_{fam}.npz"
            if not p.exists():
                continue
            s = np.load(p)
            temp = 2 * np.sqrt(2 * int(s["ndata"]))
            fr, D, Dv = s["fracs"], s["D"] / temp, s["D_field"][:, 1] / temp
            floors.append(0.5 * (D[np.argmin(np.abs(fr + 0.008))] + D[np.argmin(np.abs(fr - 0.008))]))
            for ax, y in ((axes[0], D), (axes[1], Dv)):
                ax.plot(100 * fr, y, ls=ls, marker=mk, color=col, lw=1.8, ms=4.5)
    lo, hi = min(floors), max(floors)
    for ax in axes:
        ax.axvline(0, color=INK_MUTED, lw=0.8)
        ax.set_xlabel(r"$|\Omega_b| / |\Omega_b|_{\rm ref} - 1$ [%]")
    axes[0].axhspan(lo, hi, color=INK_MUTED, alpha=0.3, lw=0)
    axes[0].text(8, lo - 1.1, f"library noise floor (D at ±0.8%): {lo:.1f}–{hi:.1f}", fontsize=8.5, color=INK_2)
    axes[0].set(ylabel="change of the prediction D (tempered)", ylim=(0, 28),
                title="all scheme B terms")
    axes[1].set(ylabel=r"change of the $\langle v_\phi\rangle$ terms alone (tempered)", ylim=(0, 8),
                title=r"$\langle v_\phi\rangle$ (m = 0, 2, 4)")
    hf = [plt.Line2D([], [], color=c, lw=2.5, label=l) for l, c in FAMILIES.values()]
    hs = [plt.Line2D([], [], color=INK, ls=ls, marker=mk, lw=1.5, label=l) for l, ls, mk in STARS.values()]
    axes[0].legend(handles=hf, frameon=False, fontsize=8.5, loc="upper center", title="family (colour)",
                   title_fontsize=8.5)
    axes[1].legend(handles=hs, frameon=False, fontsize=8.5, loc="upper center", title="stars (line style)",
                   title_fontsize=8.5)
    fig.suptitle("Pattern speed changed at a fixed bar (S = bar_frac = 1): does the predicted kinematics react? "
                 "Same stars, different potentials", color=INK, fontsize=11)
    savefig(fig, "cross_response_stars_vs_potential")


if __name__ == "__main__":
    main()

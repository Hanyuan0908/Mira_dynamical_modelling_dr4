"""Are any of the 25 simulated galaxies Milky-Way-like? Resonance geometry against the pattern-speed constraint.
(user, 2026-11-25)

Numbers, all at the sqrt(2N) tempering so that the widths are comparable:
  - 13 Auriga: ILR and CR radii at each galaxy's true pattern speed, from the cluster ground-truth files
    (results/au*_ground_truth.json, tabulated on 2026-11-23/24), and the 16-84% |Omega_b| width (PROGRESS 2026-11-24).
    Their rotation curves are NOT on this machine (potentials built on the cluster), so they cannot be drawn.
  - 12 N-body snapshots: rotation curves from their trial-family potentials (portable drive), ILR radius at their
    finite-difference pattern speed (rotation_curves_nbody_mocks_vs_milky_way.csv), ILR/CR and width as tabulated with
    the Auriga on 2026-11-24.
  - Milky Way: OSARG in Hunter+24 (width 43.3%) and in Sormani+22 (38.1%), and the steady-state mock in Hunter+24
    (46.0%); ILR, CR at Omega_b = 30.53 (milky_way_rotation_curves_and_resonances.csv). McMillan17, Cautun+20 and
    gala MWP2022 have no posterior; their ILR and ILR/CR are marked on the axis.
Panels: (a) rotation curves that can be drawn; (b) ILR/CR against the width; (c) ILR radius against the width.
Usage: python plot_resonances_vs_constraint_simulations_and_milky_way.py
Figure: figures/resonances_vs_pattern_speed_constraint_simulations_and_milky_way.png
"""
import os
import sys
import csv

import numpy as np
import agama
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_milky_way_rotation_curves as mw

REPO = "/Users/hanyuan/Library/CloudStorage/Dropbox/python_script/mira_dynamics_dr4"
USB = "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/mock"
# Auriga: ILR, CR [kpc] at the true Omega_b (cluster ground-truth files), ILR/CR, width % (sqrt(2N)); 2026-11-23/24
AURIGA = {"Au5": (1.92, 5.39, 0.356, 21.6), "Au9": (2.46, 6.67, 0.369, 12.8), "Au10": (1.79, 4.77, 0.375, 8.1),
          "Au12": (1.72, 5.60, 0.307, 18.4), "Au14": (1.85, 5.47, 0.338, 15.4), "Au17": (2.35, 5.52, 0.426, 8.5),
          "Au18": (2.15, 7.19, 0.299, 14.3), "Au20": (3.46, 12.31, 0.281, 17.0), "Au22": (1.97, 5.33, 0.370, 12.3),
          "Au24": (2.81, 6.54, 0.430, 16.4), "Au26": (1.93, 5.15, 0.375, 11.8), "Au27": (2.10, 6.68, 0.314, 18.4),
          "Au28": (2.16, 7.10, 0.304, 10.5)}
# N-body: ILR/CR, width % (sqrt(2N)), as tabulated on 2026-11-24
NBODY = {2: (0.382, 10.0), 2.375: (0.391, 12.4), 2.625: (0.399, 11.4), 3: (0.399, 10.7), 3.375: (0.404, 11.2),
         3.625: (0.406, 12.0), 4: (0.409, 14.1), 4.5: (0.410, 13.8), 5: (0.409, 14.4), 5.5: (0.401, 7.3),
         6.25: (0.401, 9.6), 7: (0.396, 6.8)}
C_AUR, C_NB = "#e07b28", "#7a7a7a"


def read_csv(name):
    p = os.path.join(HERE, name)
    return {r["model"]: r for r in csv.DictReader(open(p))}


def main():
    nbc = read_csv("rotation_curves_nbody_mocks_vs_milky_way.csv")
    mwc = read_csv("milky_way_rotation_curves_and_resonances.csv")
    fig, ax = plt.subplots(1, 3, figsize=(18.5, 5.8), layout="constrained")
    R = mw.R
    for t in sorted(NBODY):
        d = f"{USB}/nbody_potential" if t == 7 else f"{USB}/nbody_t{t:g}_potential"
        pot = agama.Potential(agama.Potential(f"{d}/dm_spherical.ini"), agama.Potential(f"{d}/Nbody_baryon_axi.ini"))
        ax[0].plot(R, mw.vc_phi_avg(pot, R), color=C_NB, lw=1, alpha=0.8)
    for lab, faxi, _, col, ls in mw.MODELS:
        ax[0].plot(R, mw.vc_axi(agama.Potential(faxi), R), color=col, ls=ls, lw=2.2)
    ax[0].set(xlim=(0, 10), ylim=(0, 320), xlabel="R [kpc]", ylabel=r"$v_c$ [km s$^{-1}$]",
              title="(a) Rotation curves available on this machine")
    ax[0].text(0.3, 18, "13 Auriga: potentials on the cluster, not here (only Au18 measured before:\n"
               "flat, v_c(2)/v_c(5) = 0.99). Their resonance radii are in (b) and (c).", fontsize=8, color="#444444")
    h = [Line2D([], [], color=C_NB, lw=1.2, label="12 N-body snapshots")]
    h += [Line2D([], [], color=c, ls=ls, lw=2.2, label=l) for l, _, _, c, ls in mw.MODELS]
    ax[0].legend(handles=h, frameon=False, fontsize=8, loc="center right")

    # (b) ILR/CR vs width, (c) ILR radius vs width
    for t, (ic, w) in NBODY.items():
        ilr = float(nbc[f"N-body t = {t:g}"]["ILR"])
        ax[1].plot(ic, w, "o", color=C_NB, ms=7, mec="white", mew=0.6)
        ax[2].plot(ilr, w, "o", color=C_NB, ms=7, mec="white", mew=0.6)
    for g, (ilr, cr, ic, w) in AURIGA.items():
        ax[1].plot(ic, w, "^", color=C_AUR, ms=8, mec="white", mew=0.6)
        ax[2].plot(ilr, w, "^", color=C_AUR, ms=8, mec="white", mew=0.6)
        ax[1].annotate(g, (ic, w), xytext=(4, 3), textcoords="offset points", fontsize=7.5, color="#8a4a12")
        ax[2].annotate(g, (ilr, w), xytext=(4, 3), textcoords="offset points", fontsize=7.5, color="#8a4a12")
    mwpts = [("OSARG, Hunter+24", "Hunter et al. (2024)", 43.3, "*", "#e0002a"),
             ("OSARG, Sormani+22", "Sormani et al. (2022) [adopted family]", 38.1, "*", "#1f5fbf"),
             ("steady-state mock, Hunter+24 (truth known)", "Hunter et al. (2024)", 46.0, "*", "white")]
    for lab, key, w, mk, fc in mwpts:
        r = mwc[key]
        ic, ilr = float(r["ILR"]) / float(r["CR"]), float(r["ILR"])
        for a, x in ((ax[1], ic), (ax[2], ilr)):
            a.plot(x, w, mk, ms=17, color=fc, mec="#222222", mew=1.0)
            a.annotate(lab, (x, w), xytext=(9, -3), textcoords="offset points", fontsize=8, color="#222222")
    for lab, faxi, _, col, ls in mw.MODELS:
        r = mwc[lab]
        ic, ilr = float(r["ILR"]) / float(r["CR"]), float(r["ILR"])
        ax[1].plot([ic, ic], [0, 2.5], color=col, lw=3)
        ax[2].plot([ilr, ilr], [0, 2.5], color=col, lw=3)
    for a, xl, xr, ttl in ((ax[1], "ILR / CR (at each galaxy's true Ω_b; Milky Way at 30.53)", (0.12, 0.46),
                            "(b) Resonance geometry vs Ω_b constraint"),
                           (ax[2], "ILR radius [kpc]", (0.9, 3.8), "(c) ILR radius vs Ω_b constraint")):
        a.set(xlabel=xl, ylabel="Ω_b posterior width, 16–84% / median [%] (√(2N) tempering)", xlim=xr, ylim=(0, 50),
              title=ttl)
        a.tick_params(top=False, right=False)
        a.text(xr[0] + 0.01 * (xr[1] - xr[0]), 3.2, "bars on the axis: Milky Way models (colours as in a)", fontsize=7.5,
               color="#555555")
    ax[2].axvspan(2, 3.8, color="#f0f0f0", zorder=0)
    ax[2].text(2.05, 47.5, "inside the fitted disc (R > 2 kpc)", fontsize=8, color="#666666")
    hb = [Line2D([], [], marker="o", ls="", color=C_NB, ms=7, label="12 N-body (own family)"),
          Line2D([], [], marker="^", ls="", color=C_AUR, ms=8, label="13 Auriga (own family)"),
          Line2D([], [], marker="*", ls="", color="#888888", mec="#222222", ms=13, label="Milky Way family runs")]
    ax[1].legend(handles=hb, frameon=False, fontsize=8, loc="center right")
    fig.suptitle("Is any simulated galaxy Milky-Way-like? Resonances of the 25 test galaxies and the Milky Way potentials, "
                 "against how well Ω_b is constrained", fontsize=10.5, x=0.01, ha="left")
    out = os.path.join(REPO, "figures", "resonances_vs_pattern_speed_constraint_simulations_and_milky_way.png")
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()

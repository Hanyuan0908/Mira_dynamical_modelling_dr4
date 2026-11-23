"""Where are the Lindblad resonances of each trial family, and how far do they move when the pattern speed changes?

A change of |Omega_b| changes the predicted kinematics mostly through the resonances it moves: orbits switch family
(x1/x2 near the ILR, trapped/untrapped near corotation) where Omega - kappa/2, Omega or Omega + kappa/2 crosses
|Omega_b|. From the azimuthally averaged midplane force of each family at S = bar_frac = 1:
  Omega(R) = v_c / R,  kappa^2 = R dOmega^2/dR + 4 Omega^2
and the ILR (Omega - kappa/2 = Omega_b), CR (Omega = Omega_b) and OLR (Omega + kappa/2 = Omega_b) radii at the
reference pattern speed and at -20% / +20%. Also max(Omega - kappa/2) over 0.3-6 kpc: if it is below Omega_b there is
no ILR at all.

Usage: python resonances_vs_pattern_speed.py
Output: results/resonances_vs_pattern_speed.csv, figures/resonances_vs_pattern_speed_mw_vs_nbody.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
from plot_presentation_figures import NBODY

R = np.arange(0.2, 25.0 + 1e-9, 0.02)


def curves(gal):
    pot = ssl.bar_aligned_potential(1.0, 1.0, gal)
    phi = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    RR, PP = np.meshgrid(R, phi, indexing="ij")
    F = pot.force(np.column_stack([(RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel(), np.zeros(RR.size)]))
    F = F.reshape(len(R), len(phi), 3)
    FR = (F[..., 0] * np.cos(PP) + F[..., 1] * np.sin(PP)).mean(axis=1)
    om2 = -FR / R
    om = np.sqrt(om2)
    kap = np.sqrt(np.clip(R * np.gradient(om2, R) + 4 * om2, 0, None))
    return om, kap


def roots(f, level):
    s = np.where((f[:-1] - level) * (f[1:] - level) <= 0)[0]
    return [float(R[i] + (f[i] - level) / (f[i] - f[i + 1]) * (R[i + 1] - R[i])) for i in s]


def main():
    set_style()
    rows, prof = [], {}
    cases = [(g, json.load(open(RESULTS_DIR / f"test2_{g}_sf_12p5pct" / "run_meta.json"))) for g in NBODY]
    cases = [(g, m["galaxy"], float(m["omega_truth"])) for g, m in cases] + [("mw", "mw", 32.5)]
    for g, gal, om_b in cases:
        om, kap = curves(gal)
        prof[g] = (om, kap, om_b)
        inner = (R > 0.3) & (R < 6)
        for f in (-0.2, 0.0, 0.2):
            ob = om_b * (1 + f)
            ilr = [r for r in roots(om - kap / 2, ob) if r < 8]
            rows.append(dict(galaxy=g, Omega_ref=om_b, frac=f, Omega_b=ob, max_om_minus_kap2_inner=float(np.max((om - kap / 2)[inner])),
                             ILR=";".join(f"{r:.2f}" for r in ilr) or "none", CR=(roots(om, ob) or [np.nan])[0],
                             OLR=(roots(om + kap / 2, ob) or [np.nan])[0]))
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "resonances_vs_pattern_speed.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 100, "display.precision", 2):
        print(t)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, g, col in ((axes[0], "t7", BLUE_RAMP[5]), (axes[1], "mw", ACCENT)):
        om, kap, om_b = prof[g]
        ax.plot(R, om, color=col, lw=2, label=r"$\Omega$")
        ax.plot(R, om - kap / 2, color=col, lw=2, ls="--", label=r"$\Omega - \kappa/2$")
        ax.plot(R, om + kap / 2, color=col, lw=2, ls=":", label=r"$\Omega + \kappa/2$")
        ax.axhspan(0.8 * om_b, 1.2 * om_b, color=INK_MUTED, alpha=0.25, lw=0, label=r"$\Omega_b$ ± 20%")
        ax.axhline(om_b, color=INK, lw=1)
        ax.axvspan(6, 15, color=INK_MUTED, alpha=0.10, lw=0)
        ax.set(xlim=(0, 12), ylim=(0, 160), xlabel="R [kpc]",
               title=("t = 7 N-body family" if g == "t7" else "Milky Way family (Hunter+24)") + f", $\\Omega_b$ = {om_b:g}")
        ax.legend(frameon=False, fontsize=9, loc="upper right")
    axes[0].set_ylabel(r"frequency [km s$^{-1}$ kpc$^{-1}$]")
    fig.suptitle("Resonance curves of the two families; grey band = the pattern speeds a ±20% change spans; shaded R > 6 kpc "
                 "is outside the fitted data", color=INK, fontsize=10.5)
    fig.tight_layout()
    savefig(fig, "resonances_vs_pattern_speed_mw_vs_nbody")


if __name__ == "__main__":
    main()

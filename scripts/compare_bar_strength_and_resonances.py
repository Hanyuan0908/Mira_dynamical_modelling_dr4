"""How do the trial families' bars compare: the Milky Way (Hunter+24, i.e. the Portail+17 bar) against the mock galaxies?

A bar whose non-axisymmetric force is weak where the data are, or whose corotation lies far outside the data, changes
the predicted kinematics little when only the pattern speed changes. This measures, for every family whose potential
files are on this machine (the 12 N-body snapshots, Au18, and the Milky Way family "mw"), at S = bar_frac = 1 (the
galaxy's own bar), in the midplane:
  v_c(R)        from the azimuthally averaged radial force
  Q_T(R)        max over phi of |F_phi| / |<F_R>_phi|, the bar's tangential force relative to the radial one
  A2/A0(R)      of the midplane density (plot_bar_strength_profile_recovery.a2_over_a0)
  R_bar         the radius beyond the A2/A0 peak where it falls to half the peak (same definition for every family)
  R_CR          corotation at the pattern speed: Omega_circ(R_CR) = |Omega_b|, with |Omega_b| the run's truth for a
                mock and 32.5 (the steady-state mock's truth) or 30 (the OSARG best) for the Milky Way
For the other 12 Auriga galaxies only their ground-truth files are here: R_CR and the particle-based A2 half-length.

Usage: python compare_bar_strength_and_resonances.py
Output: results/bar_strength_and_resonances.csv, figures/bar_strength_and_resonances_mw_vs_mocks.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
from plot_bar_strength_profile_recovery import a2_over_a0, R_GRID
from plot_presentation_figures import NBODY, AURIGA

R = np.arange(0.3, 10.0 + 1e-9, 0.1)
NPHI = 72
PHI = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)


def forces(pot):
    RR, PP = np.meshgrid(R, PHI, indexing="ij")
    pts = np.column_stack([(RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel(), np.zeros(RR.size)])
    F = pot.force(pts).reshape(len(R), NPHI, 3)
    c, s = np.cos(PP), np.sin(PP)
    FR = F[..., 0] * c + F[..., 1] * s
    Fphi = -F[..., 0] * s + F[..., 1] * c
    FR0 = FR.mean(axis=1)
    return np.sqrt(-FR0 * R), np.max(np.abs(Fphi), axis=1) / np.abs(FR0)


def r_half(prof, r):
    k = int(np.nanargmax(prof))
    half = 0.5 * prof[k]
    for i in range(k, len(r) - 1):
        if prof[i] >= half > prof[i + 1]:
            return float(r[i] + (prof[i] - half) / (prof[i] - prof[i + 1]) * (r[i + 1] - r[i]))
    return np.nan


def r_cr(vc, omega):
    om = vc / R
    i = np.where((om[:-1] - omega) * (om[1:] - omega) <= 0)[0]
    return float(R[i[0]] + (om[i[0]] - omega) / (om[i[0]] - om[i[0] + 1]) * (R[1] - R[0])) if len(i) else np.nan


def main():
    set_style()
    rows, prof = [], {}
    cases = [(g, f"test2_{g}_sf_12p5pct") for g in NBODY] + [("au18", "test2_au18_sf_12p5pct"), ("mw", None)]
    for g, run in cases:
        gal = json.load(open(RESULTS_DIR / run / "run_meta.json"))["galaxy"] if run else "mw"
        pot = ssl.bar_aligned_potential(1.0, 1.0, gal)
        vc, qt = forces(pot)
        a2 = a2_over_a0(gal, 1.0, 1.0)
        prof[g] = (vc, qt, a2)
        omegas = [float(json.load(open(RESULTS_DIR / run / "run_meta.json"))["omega_truth"])] if run else [32.5, 30.0]
        rb = r_half(a2, R_GRID)
        for om in omegas:
            rc = r_cr(vc, om)
            rows.append(dict(galaxy=g if run else f"mw (Omega {om:g})", Omega=om, Q_b=float(qt[R < 8].max()),
                             R_Qb=float(R[R < 8][np.argmax(qt[R < 8])]),
                             Q_T_3to5=float(qt[(R >= 3) & (R <= 5)].mean()), Q_T_4to6=float(qt[(R >= 4) & (R <= 6)].mean()),
                             A2_peak=float(np.nanmax(a2)), R_bar=rb, R_CR=rc, R_CR_over_R_bar=rc / rb,
                             vc_2=float(np.interp(2, R, vc)), vc_5=float(np.interp(5, R, vc)), source="potential"))
        print(rows[-1], flush=True)
    for g in AURIGA:
        if g == "au18":
            continue
        t = json.load(open(RESULTS_DIR / f"{g}_ground_truth.json"))
        rows.append(dict(galaxy=g, Omega=t["Omega_p_kms_kpc"], R_bar=t["bar_half_length_A2_kpc"], R_CR=t["R_CR_kpc"],
                         R_CR_over_R_bar=t["R_CR_over_Rbar_A2"], source="ground-truth file (particle A2 half-length)"))
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "bar_strength_and_resonances.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 3):
        print(t)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for g, (vc, qt, a2) in prof.items():
        mw = g == "mw"
        kw = dict(color=ACCENT if mw else (BLUE_RAMP[5] if g == "t7" else INK_MUTED), lw=2.4 if mw or g == "t7" else 1.0,
                  alpha=1 if mw or g == "t7" else 0.6, label="Milky Way (Hunter+24 / Portail+17 bar)" if mw else
                  ("t = 7 N-body" if g == "t7" else ("other N-body snapshots, Au18" if g == "t2" else None)))
        axes[0].plot(R, vc, **kw)
        axes[1].plot(R, qt, **kw)
        axes[2].plot(R_GRID, a2, **kw)
    axes[0].set(xlabel="R [kpc]", ylabel=r"$v_c$ [km/s]", xlim=(0, 10))
    axes[1].set(xlabel="R [kpc]", ylabel=r"$Q_T = \max|F_\phi| / |\langle F_R\rangle|$", xlim=(0, 10))
    axes[2].set(xlabel="R [kpc]", ylabel=r"midplane $A_2/A_0$", xlim=(0, 8))
    for ax in axes:
        ax.axvspan(6, 10, color=INK_MUTED, alpha=0.12, lw=0)
    axes[1].text(6.2, axes[1].get_ylim()[1] * 0.92, "outside the fitted\nR < 6 kpc", fontsize=8, color=INK_2, va="top")
    axes[0].legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle("Trial-family bars at S = bar_frac = 1: the Milky Way family against the N-body snapshots and Au18",
                 color=INK, fontsize=11)
    fig.tight_layout()
    savefig(fig, "bar_strength_and_resonances_mw_vs_mocks")


if __name__ == "__main__":
    main()

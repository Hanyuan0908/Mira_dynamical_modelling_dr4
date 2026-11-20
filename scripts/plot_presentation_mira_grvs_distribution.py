"""G_RVS distribution of the DR3 O-rich Miras (data/gaia_dr3_lpv_mira.csv), black background.

O-rich: the file as it is. Its XP-spectrum classification flags no carbon star (is_cstar_xpspectrum: 52,055 not
carbon, 3,019 unclassified); the high-amplitude XP classifier flags 12. Gaia DR3's own is_cstar would flag 8,764 -
almost all of them classified O-rich by their XP spectra - and is not used.

Only 2,958 stars (5.4%) have a measured grvs_mag, all brighter than G_RVS = 14 (DR3 publishes it no fainter). For
the rest G_RVS is estimated from G_RP with the offset measured on those 2,958: G_RVS - G_RP = OFFSET, scatter
(MAD-sigma) 0.26 mag. A quadratic in G - RP or BP - RP leaves the same scatter (0.255), so no colour term is used.
The calibration stars are bluer than the full sample (G - RP 5-95%: 1.18-1.84 against 1.43-2.00), so the
estimate is an extrapolation for the reddest stars, untested there.

Usage: python plot_presentation_mira_grvs_distribution.py
Figure: figures/presentation/dr3_mira_grvs_distribution_measured_and_estimated_black.{png,pdf}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import ROOT
from plot_presentation_figures import dark_style, save, TEXT, TEXT_2, BLUE, ORANGE

EDGES = np.arange(6.0, 18.0 + 1e-9, 0.25)


def main():
    dark_style()
    m = pd.read_csv(ROOT / "data" / "gaia_dr3_lpv_mira.csv",
                    usecols=["grvs_mag", "phot_rp_mean_mag", "phot_g_mean_mag", "DM_err", "radial_velocity"])
    has = np.isfinite(m.grvs_mag).values
    d = (m.grvs_mag - m.phot_rp_mean_mag).values[has]
    offset = float(np.median(d))
    scatter = float(1.4826 * np.median(np.abs(d - offset)))
    est = np.where(has, m.grvs_mag, m.phot_rp_mean_mag + offset)        # measured where there is one
    good = (m.DM_err < 0.4).values
    print(f"{len(m):,} Miras; measured G_RVS {has.sum():,} ({100 * has.mean():.1f}%), RV {np.isfinite(m.radial_velocity).sum():,}")
    print(f"G_RVS - G_RP on the measured: median {offset:+.3f}, MAD-sigma {scatter:.3f}")
    print(f"measured G_RVS 5/16/50/84/95%: {np.round(np.percentile(m.grvs_mag[has], [5, 16, 50, 84, 95]), 2)}")
    for lab, sel in (("all", np.ones(len(m), bool)), ("DM_err < 0.4", good)):
        e = est[sel]
        print(f"{lab:13s} N {sel.sum():,}: G_RVS 5/16/50/84/95% {np.round(np.percentile(e, [5, 16, 50, 84, 95]), 2)}; "
              + ", ".join(f"< {t}: {(e < t).sum():,} ({100 * (e < t).mean():.0f}%)" for t in (12, 13, 14, 15, 16)))

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    n, _, _ = ax.hist(est, bins=EDGES, color=BLUE, alpha=0.85, lw=0,
                      label=f"all O-rich Miras ({len(m):,}): measured, else $G_{{RP}}$ {offset:+.2f} ± {scatter:.2f}")
    ax.hist(m.grvs_mag[has], bins=EDGES, histtype="step", color=ORANGE, lw=2.2,
            label=f"measured $G_{{RVS}}$ in DR3 ({has.sum():,})")
    ax.set(xlim=(EDGES[0], EDGES[-1]), ylim=(0, 1.32 * n.max()), xlabel=r"$G_{\rm RVS}$  [mag]",
           ylabel="stars per 0.25 mag")
    ax.axvline(14.0, ymax=0.76, color=TEXT_2, lw=1.2, ls=(0, (4, 3)))     # stops below the legend
    ax.text(14.12, 1.02 * n.max(), "$G_{RVS}$ = 14\n(DR3 limit)", color=TEXT_2, fontsize=12, va="bottom")
    ax.legend(loc="upper left")
    ax.set_title(r"Gaia DR3 O-rich Miras: $G_{\rm RVS}$ distribution", color=TEXT, fontsize=16, pad=10)
    save(fig, "dr3_mira_grvs_distribution_measured_and_estimated_black")


if __name__ == "__main__":
    main()

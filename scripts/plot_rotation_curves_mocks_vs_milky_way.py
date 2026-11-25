"""Rotation curves of the 12 N-body test galaxies against the Milky Way potentials (user, 2026-11-25).

N-body snapshot t: the axisymmetric part its trial family uses, i.e. the spherical DM halo (dm_spherical.ini) plus the
axisymmetrised baryons (Nbody_baryon_axi.ini), from the portable drive's copy of mock/nbody[_t<T>]_potential/.
Pattern speeds: the finite-difference truths of the 25-galaxy table (PROGRESS 2026-11-21). Milky Way models as in
plot_milky_way_rotation_curves.py, at the OSARG pattern speed 30.53 (Sormani+22, geo20 |z|<2, T1).
Panels: v_c(R); v_c(R)/v_c(5 kpc) (shape only); (Omega - kappa/2)/Omega_b, which crosses 1 at the ILR for every galaxy at
its own pattern speed. Printed: v_c at 1/2/5 kpc, v_c(2)/v_c(5), radius of the inner peak, ILR, CR, ILR/CR, and the
pattern-speed width of each N-body mock's posterior (sqrt(2N) tempered, PROGRESS 2026-11-21 table) for reference.
Usage: python plot_rotation_curves_mocks_vs_milky_way.py OUTDIR
"""
import os
import sys
import csv

import numpy as np
import agama
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_milky_way_rotation_curves as mw

U = "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/mock"
# t: (truth Omega_b, tempered posterior p16, p50, p84)  -- PROGRESS 2026-11-21 (25 galaxies at T = 1 vs tempered)
NB = {2: (43.90, 40.12, 42.23, 44.52), 2.375: (40.87, 39.21, 41.53, 44.28), 2.625: (40.21, 38.10, 40.26, 42.70),
      3: (38.67, 37.14, 39.12, 41.27), 3.375: (36.83, 35.91, 37.77, 40.04), 3.625: (34.79, 33.90, 35.81, 38.08),
      4: (33.28, 32.64, 35.14, 37.34), 4.5: (32.07, 30.79, 33.25, 35.23), 5: (30.60, 27.70, 29.32, 32.11),
      5.5: (28.21, 25.75, 26.70, 27.82), 6.25: (26.58, 25.11, 26.28, 27.67), 7: (25.20, 23.98, 24.81, 25.71)}
R = mw.R


def resonances(vc, ob):
    om, kap = mw.freqs(vc)
    cr = mw.crossings(om, ob)
    il = [r for r in mw.crossings(om - kap / 2, ob) if r < (cr or [99])[0]]
    return om, kap, (il[-1] if il else np.nan), (cr[0] if cr else np.nan)


def main():
    out = sys.argv[1]
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), layout="constrained")
    cmap = plt.get_cmap("Greys")
    ts = sorted(NB)
    for k, t in enumerate(ts):
        d = f"{U}/nbody_potential" if t == 7 else f"{U}/nbody_t{t:g}_potential"
        pot = agama.Potential(agama.Potential(f"{d}/dm_spherical.ini"), agama.Potential(f"{d}/Nbody_baryon_axi.ini"))
        vc = mw.vc_phi_avg(pot, R)   # azimuthal average: the cluster-built 'axi' baryon files are not axisymmetric
        ob, p16, p50, p84 = NB[t]
        om, kap, ilr, cr = resonances(vc, ob)
        col = cmap(0.35 + 0.6 * k / (len(ts) - 1))
        ax[0].plot(R, vc, color=col, lw=1.2)
        ax[1].plot(R, vc / np.interp(5, R, vc), color=col, lw=1.2)
        ax[2].plot(R, (om - kap / 2) / ob, color=col, lw=1.2)
        if np.isfinite(ilr):
            ax[2].plot(ilr, 1, "o", ms=5, color=col, mec="white", mew=0.6)
        rows.append(dict(model=f"N-body t = {t:g}", Omega_b=ob, vc_1=np.interp(1, R, vc), vc_2=np.interp(2, R, vc),
                         vc_5=np.interp(5, R, vc), ratio_2_5=np.interp(2, R, vc) / np.interp(5, R, vc),
                         R_peak_inner=R[np.argmax(np.where(R < 6, vc, 0))], ILR=ilr, CR=cr, ILR_over_CR=ilr / cr,
                         width_pct_tempered=100 * (p84 - p16) / p50))
    for lab, faxi, ffull, col, ls in mw.MODELS:
        vc = mw.vc_axi(agama.Potential(faxi), R)
        om, kap, ilr, cr = resonances(vc, mw.OMB)
        ax[0].plot(R, vc, color=col, ls=ls, lw=2.2)
        ax[1].plot(R, vc / np.interp(5, R, vc), color=col, ls=ls, lw=2.2)
        ax[2].plot(R, (om - kap / 2) / mw.OMB, color=col, ls=ls, lw=2.2)
        if np.isfinite(ilr):
            ax[2].plot(ilr, 1, "s", ms=6, color=col, mec="white", mew=0.6)
        rows.append(dict(model=lab, Omega_b=mw.OMB, vc_1=np.interp(1, R, vc), vc_2=np.interp(2, R, vc),
                         vc_5=np.interp(5, R, vc), ratio_2_5=np.interp(2, R, vc) / np.interp(5, R, vc),
                         R_peak_inner=R[np.argmax(np.where(R < 6, vc, 0))], ILR=ilr, CR=cr, ILR_over_CR=ilr / cr,
                         width_pct_tempered=np.nan))
    for r in rows:
        print(f"{r['model']:42s} Om_b {r['Omega_b']:5.2f}  v_c 1/2/5 kpc {r['vc_1']:3.0f} {r['vc_2']:3.0f} {r['vc_5']:3.0f}  "
              f"v_c(2)/v_c(5) {r['ratio_2_5']:.3f}  peak (R<6) at {r['R_peak_inner']:.2f}  ILR {r['ILR']:.2f}  CR {r['CR']:.2f}  "
              f"ILR/CR {r['ILR_over_CR']:.2f}" + (f"  Omega_b width {r['width_pct_tempered']:.1f}%" if np.isfinite(r['width_pct_tempered']) else ""))
    ax[0].set(xlim=(0, 10), ylim=(0, 320), xlabel="R [kpc]", ylabel=r"$v_c$ [km s$^{-1}$]",
              title="Circular velocity (axisymmetric part of each trial family)")
    ax[1].set(xlim=(0, 8), ylim=(0.5, 1.35), xlabel="R [kpc]", ylabel=r"$v_c(R)\,/\,v_c(5\,\mathrm{kpc})$",
              title="Shape: normalised at 5 kpc")
    ax[1].axhline(1, color="#999999", lw=0.8, ls=":")
    ax[1].axvspan(0, 6, color="#eeeeee", zorder=0)
    ax[1].text(0.15, 0.53, "fitted region R < 6 kpc", fontsize=8, color="#666666")
    ax[2].axhline(1, color="#555555", lw=1, ls="--")
    ax[2].axvspan(2, 6, color="#eeeeee", zorder=0)
    ax[2].set(xlim=(0, 8), ylim=(0, 3), xlabel="R [kpc]", ylabel=r"$(\Omega - \kappa/2)\,/\,\Omega_b$",
              title="ILR where each curve crosses 1 (own Ω_b; dots / squares)")
    for a in ax:
        a.tick_params(top=False, right=False)
    handles = [Line2D([], [], color=cmap(0.65), lw=1.5, label="12 N-body snapshots t = 2–7 (light → dark)")]
    handles += [Line2D([], [], color=c, ls=ls, lw=2.2, label=l) for l, _, _, c, ls in mw.MODELS]
    ax[0].legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("The 12 N-body test galaxies against the Milky Way potentials (Auriga: potentials not on this "
                 "machine outside the locked Dropbox folder)", fontsize=10, x=0.01, ha="left")
    fig.savefig(os.path.join(out, "rotation_curves_nbody_mocks_vs_milky_way.png"), dpi=150)
    with open(os.path.join(out, "rotation_curves_nbody_mocks_vs_milky_way.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()

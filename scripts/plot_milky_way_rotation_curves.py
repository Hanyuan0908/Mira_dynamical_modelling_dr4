"""Rotation curves of published Milky Way potentials, and where their inner Lindblad resonance falls (user, 2026-11-25).

Models (AGAMA files, Msun / kpc / km/s):
  McMillan (2017)                agama/data/McMillan17.ini               pre-Gaia-DR2 constraints
  Sormani et al. (2022)          Portail17_axi.ini (AGAMA example_mw_bar_potential.py; bar from Portail+17 M2M)
  Hunter et al. (2024)           MWPotentialHunter24_axi.ini (AGAMA example_mw_potential_hunter24.py)
  Cautun et al. (2020)           agama/data/Cautun20.ini, contracted NFW halo fitted to the Gaia DR2 rotation curve
  Price-Whelan / gala (2022)     agama/data/PriceWhelan22.ini, gala MilkyWayPotential2022 (disc fitted to the
                                 Eilers+2019 Gaia-DR2-based circular velocities inside ~15 kpc)
The barred models are taken in their axisymmetrised form (the *_axi files); the azimuthal average of the full barred
potential is checked against it. v_c = sqrt(R dPhi/dR) in the midplane; Omega = v_c/R; kappa^2 = R dOmega^2/dR + 4 Omega^2.
Panels: v_c to 20 kpc; v_c in the fitted region (R < 6 kpc); Omega and Omega - kappa/2 with the OSARG pattern-speed
posterior (Sormani+22, geo20 |z|<2, T1: 30.53, 16-84% 28.10-36.07) as a band. ILR: Omega - kappa/2 = Omega_b
(outermost crossing inside corotation); CR: Omega = Omega_b.
Usage: python plot_milky_way_rotation_curves.py
Output: figures/milky_way_rotation_curves_and_resonances.png, results/milky_way_rotation_curves_and_resonances.csv
"""
import os
import sys

import numpy as np
import agama
import matplotlib.pyplot as plt

agama.setUnits(mass=1, length=1, velocity=1)
AD = os.path.join(os.path.dirname(agama.__file__), "data")
OC = os.path.expanduser("~/Desktop/PhD_projects/oCen_bar/agama_potentials")
MODELS = [  # label, axisymmetric file, full file (or None), colour, line style
    ("McMillan (2017)", f"{AD}/McMillan17.ini", None, "#2b2b2b", "-"),
    ("Sormani et al. (2022) [adopted family]", f"{OC}/Portail17_axi.ini", f"{OC}/Portail17.ini", "#1f5fbf", "-"),
    ("Hunter et al. (2024)", f"{OC}/MWPotentialHunter24_axi.ini", f"{OC}/MWPotentialHunter24_full.ini", "#e07b28", "-"),
    ("Cautun et al. (2020), Gaia DR2", f"{AD}/Cautun20.ini", None, "#2e8b57", "--"),
    ("Price-Whelan / gala MWP2022, Gaia DR2", f"{AD}/PriceWhelan22.ini", None, "#8e44ad", "--"),
]
OMB, OMB_LO, OMB_HI = 30.53, 28.10, 36.07
R = np.linspace(0.05, 20.0, 2000)


def vc_axi(pot, R):
    F = pot.force(np.column_stack([R, np.zeros_like(R), np.zeros_like(R)]))
    return np.sqrt(np.clip(-R * F[:, 0], 0, None))


def vc_phi_avg(pot, R, nphi=72):
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    RR, PP = np.meshgrid(R, ph, indexing="ij")
    x, y = RR * np.cos(PP), RR * np.sin(PP)
    F = pot.force(np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)]))
    FR = (F[:, 0] * np.cos(PP.ravel()) + F[:, 1] * np.sin(PP.ravel())).reshape(RR.shape)
    return np.sqrt(np.clip(-R * FR.mean(axis=1), 0, None))


def freqs(vc):
    om = vc / R
    kap2 = R * np.gradient(om ** 2, R) + 4 * om ** 2
    return om, np.sqrt(np.clip(kap2, 0, None))


def crossings(y, level):
    s = np.sign(y - level)
    i = np.where(np.diff(s) != 0)[0]
    return [R[k] + (level - y[k]) * (R[k + 1] - R[k]) / (y[k + 1] - y[k]) for k in i]


def main():
    from common import FIG_DIR, RESULTS_DIR
    out, res = (sys.argv[1], sys.argv[1]) if len(sys.argv) > 1 else (str(FIG_DIR), str(RESULTS_DIR))
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4), layout="constrained")
    for lab, faxi, ffull, col, ls in MODELS:
        pot = agama.Potential(faxi)
        vc = vc_axi(pot, R)
        chk = ""
        if ffull:
            vf = vc_phi_avg(agama.Potential(ffull), R[::50])
            chk = f"; full barred model, azimuthal average, max |dv_c| {np.max(np.abs(vf - vc[::50])):.2f} km/s"
        om, kap = freqs(vc)
        ilr = lambda ob: [r for r in crossings(om - kap / 2, ob) if r < (crossings(om, ob) or [99])[0]]
        cr = crossings(om, OMB)
        il = ilr(OMB)
        r = dict(model=lab, vc_1=np.interp(1, R, vc), vc_2=np.interp(2, R, vc), vc_5=np.interp(5, R, vc),
                 vc_8=np.interp(8.122, R, vc), vc_15=np.interp(15, R, vc), ratio_2_5=np.interp(2, R, vc) / np.interp(5, R, vc),
                 ILR=il[-1] if il else np.nan, ILR_lo=(ilr(OMB_HI) or [np.nan])[-1], ILR_hi=(ilr(OMB_LO) or [np.nan])[-1],
                 CR=cr[0] if cr else np.nan, peak_inner=R[np.argmax(np.where(R < 4, vc, 0))])
        rows.append(r)
        print(f"{lab:42s} v_c at 1/2/5/8.12/15 kpc: {r['vc_1']:.0f} {r['vc_2']:.0f} {r['vc_5']:.0f} {r['vc_8']:.0f} "
              f"{r['vc_15']:.0f}; v_c(2)/v_c(5) {r['ratio_2_5']:.3f}; at Omega_b {OMB}: ILR {r['ILR']:.2f} kpc "
              f"(band {r['ILR_lo']:.2f}-{r['ILR_hi']:.2f}), CR {r['CR']:.2f} kpc{chk}")
        for a in ax[:2]:
            a.plot(R, vc, color=col, ls=ls, lw=2, label=lab)
        ax[2].plot(R, om - kap / 2, color=col, ls=ls, lw=2)
        ax[2].plot(R, om, color=col, ls=ls, lw=1, alpha=0.5)
    ax[0].set(xlim=(0, 20), ylim=(0, 300), xlabel="R [kpc]", ylabel=r"$v_c$ [km s$^{-1}$]",
              title="Circular velocity (axisymmetrised), midplane")
    ax[0].axvline(8.122, color="#999999", lw=0.8, ls=":")
    ax[0].text(8.3, 15, "Sun", fontsize=8, color="#666666")
    ax[0].legend(frameon=False, fontsize=8.5, loc="lower right")
    ax[1].set(xlim=(0, 6), ylim=(120, 290), xlabel="R [kpc]", ylabel=r"$v_c$ [km s$^{-1}$]",
              title="Zoom: the fitted region (R < 6 kpc)")
    ax[2].axhspan(OMB_LO, OMB_HI, color="#cccccc", alpha=0.5, lw=0)
    ax[2].axhline(OMB, color="#555555", lw=1, ls="--")
    ax[2].text(3.0, OMB_HI + 1.5, f"OSARG Ω_b posterior (Sormani+22, geo20, T1):\n{OMB} (16–84%: {OMB_LO}–{OMB_HI})", fontsize=8, color="#444444")
    ax[2].set(xlim=(0, 9), ylim=(0, 80), xlabel="R [kpc]", ylabel=r"[km s$^{-1}$ kpc$^{-1}$]",
              title="Ω − κ/2 (thick) and Ω (thin): ILR and corotation")
    for a in ax:
        a.tick_params(top=False, right=False)
    fig.suptitle("Milky Way potentials: rotation curves and resonances (solid: our trial families' sources and McMillan17; "
                 "dashed: post-2017 models fitted to Gaia-DR2-based rotation curves)", fontsize=10, x=0.01, ha="left")
    fig.savefig(os.path.join(out, "milky_way_rotation_curves_and_resonances.png"), dpi=150)
    import csv
    with open(os.path.join(res, "milky_way_rotation_curves_and_resonances.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()

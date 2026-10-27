"""Bar pattern speed of Au18 with the Fourier method of Dehnen, Semczuk & Schönrich (2023, MNRAS
518, 2712; code in python_script/patternSpeed-main), and the bar resonances in the ground-truth potential.

Usage follows python_script/get_bar_property_ForSingleSnapshot.py and Nexus/get_pot_nexus.py:
the disc is mirrored to rotate anticlockwise, and FourierMethod is given centred x, y, vx, vy
in kpc and km/s, so Omega comes out in km/s/kpc. The automatic bar finder
(FourierMethod.patternSpeed) is the fiducial result, and measureOmega on fixed radial ranges is
the cross-check. No extra unit factor is applied.
"""
import argparse
import json
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from build_au18_potential import vcirc
from galaxies import galaxy, add_argument as add_galaxy_argument, load_potential

PATTERN_SPEED_DIR = "/Users/hanyuan/Library/CloudStorage/Dropbox/python_script/patternSpeed-main"

if not hasattr(np, "asfarray"):   # removed in NumPy 2.0; patternSpeed.py still calls it
    np.asfarray = lambda a: np.asarray(a, dtype=np.float64)
sys.path.insert(0, PATTERN_SPEED_DIR)
import patternSpeed  # noqa: E402


def resonances(pot, Omega_p, R=np.linspace(0.2, 20, 2000)):
    """Omega(R), kappa(R) of the axisymmetrised potential and the CR / OLR / ILR radii."""
    Om = vcirc(pot, R) / R
    kappa = np.sqrt(R * np.gradient(Om**2, R) + 4 * Om**2)

    def roots(f):
        s = np.flatnonzero(np.diff(np.sign(f)) != 0)
        return [float(R[i] - f[i] * (R[i + 1] - R[i]) / (f[i + 1] - f[i])) for i in s]
    return R, Om, kappa, dict(CR=roots(Om - Omega_p), OLR=roots(Om + kappa / 2 - Omega_p),
                              ILR=roots(Om - kappa / 2 - Omega_p))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    add_galaxy_argument(ap)
    name = ap.parse_args().galaxy
    g = galaxy(name)
    RESULTS_FILE = g["truth"]
    set_style()
    with h5py.File(g["stars"], "r") as f:
        x, y, vx, vy, m = (f[k][:].astype(np.float64) for k in ("x", "y", "vx", "vy", "mass"))
        R_bar_A2 = float(f.attrs["bar_half_length_kpc"])

    # MW-like frame rotates clockwise: mirror y -> -y to give the method an anticlockwise disc
    pat, disc = patternSpeed.FourierMethod.patternSpeed(x, -y, vx, -vy, mu=m)
    tool = patternSpeed.FourierMethod(x, -y, vx, -vy, mu=m)
    checks = {f"{r0:g}-{r1:g}": tool.measureOmega(region=(r0, r1)).to_dict()
              for r0, r1 in ((0.1, R_bar_A2), (0.1, 5.0), (0.5, 4.0))}
    Om_p, Om_e = float(pat["Om"]), float(pat["Om_e"])
    print(f"fiducial: Omega_p = {Om_p:.2f} +/- {Om_e:.2f} km/s/kpc, bar region "
          f"{pat['R0']:.2f}-{pat['R1']:.2f} kpc (Rm = {pat['Rm']:.2f}), "
          f"bar phase {np.degrees(pat['psi']):.2f} deg (mirrored frame)")
    for k, v in checks.items():
        print(f"  measureOmega R = {k} kpc: {v['Om']:.2f} +/- {v['Om_e']:.2f}")

    comps, total = load_potential(name)
    R, Om, kappa, res = resonances(total, Om_p)
    R_CR = res["CR"][0]
    print(f"corotation R_CR = {R_CR:.2f} kpc; OLR = {res['OLR']}; ILR = {res['ILR']}")
    print(f"R_CR / R_bar = {R_CR / R_bar_A2:.2f} (R_bar = {R_bar_A2:.2f} kpc, A2 half-max), "
          f"{R_CR / pat['R1']:.2f} (R1 = {pat['R1']:.2f} kpc, Dehnen bar region)")

    truth = dict(
        Omega_p_kms_kpc=Om_p, Omega_p_err=Om_e,
        sense="clockwise in the MW-like frame: bar angle decreases with time, "
              "phi_bar(t) = -27 deg - Omega_p t",
        bar_angle_deg_at_t0=-27.0, bar_region_dehnen_kpc=[float(pat["R0"]), float(pat["R1"])],
        bar_Rm_kpc=float(pat["Rm"]), bar_half_length_A2_kpc=R_bar_A2,
        R_CR_kpc=R_CR, R_OLR_kpc=res["OLR"], R_ILR_kpc=res["ILR"],
        R_CR_over_Rbar_A2=R_CR / R_bar_A2, cross_checks={k: {"Om": v["Om"], "Om_e": v["Om_e"]}
                                                       for k, v in checks.items()},
        potential=f"{g['pot_dir'].relative_to(ROOT)}/{{{','.join(g['components'])}}}.ini",
        method="Dehnen, Semczuk & Schoenrich (2023) Fourier method, patternSpeed.py v" +
               patternSpeed.version)
    RESULTS_FILE.parent.mkdir(exist_ok=True)
    with open(RESULTS_FILE, "w") as f:
        json.dump(truth, f, indent=2)
    print(f"wrote {RESULTS_FILE.relative_to(ROOT)}")

    # figure: Fourier analysis of the disc and the resonance diagram
    fig, (a1, a3) = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw=dict(wspace=0.32))
    Rd = disc["R"].values
    a1.plot(Rd, disc["A2"], color=BLUE_RAMP[4], lw=2, label=r"$A_2 = \Sigma_2/\Sigma_0$")
    a1.axvspan(pat["R0"], pat["R1"], color=BLUE_RAMP[0], alpha=0.6, lw=0, label="bar region (method)")
    a1.axvline(R_bar_A2, color=INK_MUTED, ls="-.", lw=0.8)
    a1.text(R_bar_A2 * 1.05, 0.55, r"$R_{\rm bar}$ ($A_2$ half max)", color=INK_2, fontsize=8, rotation=90,
            va="bottom", transform=a1.get_xaxis_transform())
    a1.set_xscale("log")
    a1.set_xlim(0.1, 20)
    a1.set_ylim(0, None)
    a1.set_xlabel(r"$R$ [kpc]")
    a1.set_ylabel(r"$A_2$")
    a2 = a1.twinx()
    a2.plot(Rd, np.degrees(disc["ps2"]), color=ACCENT, lw=1.2, label=r"$\psi_2$")
    a2.set_ylabel(r"$\psi_2$ [deg] (mirrored frame)", color=ACCENT)
    a2.set_ylim(-90, 180)
    a1.legend(loc="upper right", fontsize=8)
    a1.set_title(rf"Fourier method: $\Omega_p$ = {Om_p:.2f} $\pm$ {Om_e:.2f} km s$^{{-1}}$ kpc$^{{-1}}$",
                 color=INK_2)

    a3.plot(R, Om, color=INK, lw=2, label=r"$\Omega$")
    a3.plot(R, Om - kappa / 2, color=BLUE_RAMP[4], lw=1.5, label=r"$\Omega \pm \kappa/2$")
    a3.plot(R, Om + kappa / 2, color=BLUE_RAMP[4], lw=1.5)
    a3.axhspan(Om_p - Om_e, Om_p + Om_e, color=ACCENT, alpha=0.35, lw=0)
    a3.axhline(Om_p, color=ACCENT, lw=1.2, label=r"$\Omega_p$")
    for name, rr in (("CR", res["CR"]), ("OLR", res["OLR"])):
        for r in rr:
            a3.plot(r, Om_p, "o", color=ACCENT, ms=6, mec="white", zorder=5)
            a3.annotate(f"{name}\n{r:.2f} kpc", (r, Om_p), xytext=(4, 8), textcoords="offset points",
                        fontsize=8, color=INK_2)
    a3.axvline(R_bar_A2, color=INK_MUTED, ls="-.", lw=0.8)
    a3.set_xlim(0, 15)
    a3.set_ylim(0, 150)
    a3.set_xlabel(r"$R$ [kpc]")
    a3.set_ylabel(r"frequency [km s$^{-1}$ kpc$^{-1}$]")
    a3.legend(loc="upper right", fontsize=8)
    a3.set_title(rf"Resonances: $R_{{\rm CR}}/R_{{\rm bar}}$ = {R_CR / R_bar_A2:.2f}", color=INK_2)
    savefig(fig, f"{g['prefix']}_pattern_speed")


if __name__ == "__main__":
    main()

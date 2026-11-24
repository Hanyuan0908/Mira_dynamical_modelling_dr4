"""Trial-potential family for the Milky Way from Sormani et al. (2022): AGAMA's example_mw_bar_potential.py, the
analytic fit to the Portail et al. (2017) M2M model, as exported to Portail17.ini / Portail17_axi.ini.

Why a second Milky Way family (user, 2026-11-24): the OSARG pattern-speed constraint is weak because the Hunter+24
rotation curve makes the R < 6 kpc kinematics insensitive to Omega_b (PROGRESS 2026-11-24). Sormani+22 shares the
bar but has its own disc, central mass concentration and halo, hence its own rotation curve.

The exported model is two components:
  CylSpline  X-shaped bar + two long bars + disc + central mass concentration, mmax 12, triaxial, bar along x
  Multipole  the dark halo, FLATTENED (axisymmetric, lmax 8; q = 0.8 in the script)
The project's rule is a spherical halo in every trial potential (memory: spherical DM halo), so fixed.ini is the
halo re-expanded with lmax = 0 (its spherical average). The effect on v_c is printed. Everything else mirrors
build_mw_potential.py:
  fixed.ini           the spherical halo, never rescaled
  MW_baryon_full.ini  the triaxial CylSpline
  MW_baryon_axi.ini   its axisymmetrised version (from Portail17_axi.ini)
  MW_full.ini, MW_axi.ini   fixed + each of the above
so steady_state_likelihood.bar_aligned_potential(bar_frac, S, "sormani22") = fixed + baryon_axi
+ bar_frac [baryon_full(x/S) - baryon_axi(x/S)], no rotation.

Usage: python build_sormani22_potential.py [--source DIR]
Output: mock/sormani22_potential/{fixed,MW_baryon_full,MW_baryon_axi,MW_full,MW_axi}.ini
"""
import argparse
from pathlib import Path

import numpy as np
import agama

from common import MOCK_DIR
from build_mw_potential import sections, write

SOURCE = Path.home() / "Desktop/PhD_projects/oCen_bar/agama_potentials"
OUT = MOCK_DIR / "sormani22_potential"


def vc_curve(pot, R):
    phi = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    RR, PP = np.meshgrid(R, phi, indexing="ij")
    F = pot.force(np.column_stack([(RR * np.cos(PP)).ravel(), (RR * np.sin(PP)).ravel(), np.zeros(RR.size)]))
    F = F.reshape(len(R), len(phi), 3)
    return np.sqrt(-R * (F[..., 0] * np.cos(PP) + F[..., 1] * np.sin(PP)).mean(axis=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=SOURCE)
    args = ap.parse_args()
    agama.setUnits(mass=1, length=1, velocity=1)
    full, axi = sections(args.source / "Portail17.ini"), sections(args.source / "Portail17_axi.ini")
    assert len(full) == 2 and len(axi) == 2, "expected CylSpline + Multipole in each file"
    assert "type=CylSpline" in full[0] and "type=Multipole" in full[1] and "type=CylSpline" in axi[0]
    assert full[1] == axi[1], "the halo sections of the full and axi files differ"
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = OUT / "_halo_flattened.ini"
    write([full[1]], tmp)
    halo = agama.Potential(str(tmp))
    radii = [float(l.split()[0]) for l in full[1].split("#radius", 1)[1].splitlines()[1:] if l.strip()[:1].isdigit()]
    sph = agama.Potential(type="Multipole", potential=halo, lmax=0, symmetry="spherical", gridSizeR=len(radii),
                          rmin=min(radii), rmax=max(radii))
    sph.export(str(OUT / "fixed.ini"))
    tmp.unlink()
    write([full[0]], OUT / "MW_baryon_full.ini")
    write([axi[0]], OUT / "MW_baryon_axi.ini")
    fixed_block = sections(OUT / "fixed.ini")
    write(fixed_block + [full[0]], OUT / "MW_full.ini")
    write(fixed_block + [axi[0]], OUT / "MW_axi.ini")
    print(f"wrote {OUT} (halo re-expanded spherically: {len(radii)} radii, {min(radii):g}-{max(radii):g} kpc)")

    import steady_state_likelihood as ssl
    R = np.array([1, 2, 3, 4, 5, 6, 8.0])
    src = agama.Potential(str(args.source / "Portail17.ini"))
    fam = ssl.bar_aligned_potential(1.0, 1.0, "sormani22")
    h24 = ssl.bar_aligned_potential(1.0, 1.0, "mw")
    for lab, p in (("Sormani+22 as published (flattened halo)", src), ("sormani22 family, S = bar_frac = 1 (spherical halo)", fam),
                   ("Hunter+24 family, S = bar_frac = 1", h24)):
        print(f"v_c at R = 1, 2, 3, 4, 5, 6, 8 kpc, {lab}: " + " ".join(f"{v:.0f}" for v in vc_curve(p, R)))
    rng = np.random.default_rng(0)
    Rr, ph = rng.uniform(0.2, 8.0, 4000), rng.uniform(0, 2 * np.pi, 4000)
    pts = np.column_stack([Rr * np.cos(ph), Rr * np.sin(ph), rng.uniform(-3, 3, 4000)])
    rel = np.linalg.norm(fam.force(pts) - src.force(pts), axis=1) / np.linalg.norm(src.force(pts), axis=1)
    print(f"family (spherical halo, mmax 8 re-expansion) against the published model, 4000 points in R 0.2-8, |z| < 3: "
          f"relative force difference median {np.median(rel):.1e}, 99th pct {np.percentile(rel, 99):.1e}")
    phi = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    for r in (1.0, 2.0, 3.0, 4.0, 5.0):
        rho = fam.density(np.column_stack([r * np.cos(phi), r * np.sin(phi), np.zeros_like(phi)]))
        a2 = np.sum(rho * np.exp(-2j * phi))
        print(f"  R = {r:.0f} kpc: m = 2 density phase {np.degrees(np.angle(a2)) / 2:+.2f} deg (0 = bar along x), "
              f"A2/A0 {2 * abs(a2) / rho.sum():.3f}")


if __name__ == "__main__":
    main()

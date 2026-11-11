"""Trial-potential family for the Milky Way: the Hunter et al. (2024) model, which carries the Sormani et al. (2022)
analytic bar (an approximation of the Portail et al. 2017 M2M bar), as exported by AGAMA's
example_mw_potential_hunter24.py.

That model is two components (AGAMA's own split):
  Multipole  every spheroidal piece: SgrA*, NSC, NSD and a SPHERICAL Einasto dark halo (no axisRatioZ), so the
             project's spherical-halo rule holds without changes;
  CylSpline  bar + thin/thick stellar discs + HI and H2 gas discs, mmax 8, bar along the x axis, non-rotating.
The family files mirror the N-body ones (potential_family.files):
  fixed.ini           the Multipole, never rescaled
  MW_baryon_full.ini  the triaxial CylSpline
  MW_baryon_axi.ini   its axisymmetrised (mmax 0) version
  MW_full.ini, MW_axi.ini   fixed + each of the above
so steady_state_likelihood.bar_aligned_potential(bar_frac, S, "mw") = fixed + baryon_axi
+ bar_frac [baryon_full(x/S) - baryon_axi(x/S)], with no rotation (the bar is already along x).

The source files are the ones AGAMA writes (MWPotentialHunter24_full.ini / _axi.ini); --source points at them.
Checks: the family at bar_frac = S = 1 against the source model (forces), and the m = 2 density phase (bar on x).

Usage: python build_mw_potential.py [--source DIR]
Output: mock/mw_potential/{fixed,MW_baryon_full,MW_baryon_axi,MW_full,MW_axi}.ini
"""
import argparse
from pathlib import Path

import numpy as np
import agama

from common import MOCK_DIR

SOURCE = Path.home() / "Library/CloudStorage/Dropbox/python_script/oCen_bar/agama_potentials"
OUT = MOCK_DIR / "mw_potential"


def sections(path):
    """The [Potential...] blocks of an AGAMA .ini file, header lines stripped."""
    return [b.split("]", 1)[1] for b in open(path).read().split("[Potential")[1:]]


def write(blocks, out):
    with open(out, "w") as f:
        for i, body in enumerate(blocks):
            f.write(f"[Potential{i if i else ''}]" + body.rstrip("\n") + "\n\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=SOURCE)
    args = ap.parse_args()
    agama.setUnits(mass=1, length=1, velocity=1)
    full, axi = sections(args.source / "MWPotentialHunter24_full.ini"), sections(args.source / "MWPotentialHunter24_axi.ini")
    assert len(full) == 2 and len(axi) == 2, "expected Multipole + CylSpline in each file"
    assert full[0] == axi[0], "the Multipole sections of the full and axi files differ"
    assert "type=Multipole" in full[0] and "type=CylSpline" in full[1] and "type=CylSpline" in axi[1]
    OUT.mkdir(parents=True, exist_ok=True)
    write([full[0]], OUT / "fixed.ini")
    write([full[1]], OUT / "MW_baryon_full.ini")
    write([axi[1]], OUT / "MW_baryon_axi.ini")
    write(full, OUT / "MW_full.ini")
    write(axi, OUT / "MW_axi.ini")
    print(f"wrote {OUT}")

    import steady_state_likelihood as ssl
    src = agama.Potential(str(args.source / "MWPotentialHunter24_full.ini"))
    rng = np.random.default_rng(0)
    R = rng.uniform(0.2, 8.0, 4000)
    ph = rng.uniform(0, 2 * np.pi, 4000)
    pts = np.column_stack([R * np.cos(ph), R * np.sin(ph), rng.uniform(-3, 3, 4000)])
    trial = ssl.bar_aligned_potential(1.0, 1.0, "mw")
    f0, f1 = src.force(pts), trial.force(pts)
    rel = np.linalg.norm(f1 - f0, axis=1) / np.linalg.norm(f0, axis=1)
    print(f"family at bar_frac = S = 1 against the source model, 4000 points in R 0.2-8, |z| < 3: "
          f"relative force error median {np.median(rel):.1e}, 99th pct {np.percentile(rel, 99):.1e}, "
          f"max {rel.max():.1e}")
    vc = lambda p, r: float(np.sqrt(-r * agama.Potential(str(OUT / "MW_axi.ini")).force([[r, 0, 0]])[0][0]))
    print(f"v_c (axisymmetrised): R = 4 kpc {vc(None, 4.0):.1f}, R = 8.122 kpc {vc(None, 8.122):.1f} km/s")
    phi = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    for r in (1.0, 2.0, 3.0, 4.0):
        rho = trial.density(np.column_stack([r * np.cos(phi), r * np.sin(phi), np.zeros_like(phi)]))
        a2 = np.sum(rho * np.exp(-2j * phi))
        print(f"  R = {r:.0f} kpc: m = 2 density phase {np.degrees(np.angle(a2)) / 2:+.2f} deg "
              f"(0 = bar along x), A2/A0 {2 * abs(a2) / rho.sum():.3f}")
    for f, s in ((0.5, 1.0), (1.0, 1.5), (1.5, 0.5)):
        p = ssl.bar_aligned_potential(f, s, "mw")
        rho = p.density(np.column_stack([2.0 * np.cos(phi), 2.0 * np.sin(phi), np.zeros_like(phi)]))
        a2 = np.sum(rho * np.exp(-2j * phi))
        print(f"  bar_frac {f}, S {s}: R = 2 kpc A2/A0 {2 * abs(a2) / rho.sum():.3f}, "
              f"phase {np.degrees(np.angle(a2)) / 2:+.2f} deg")


if __name__ == "__main__":
    main()

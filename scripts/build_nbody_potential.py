"""Ground-truth potential of the isolated N-body galaxy in the MW-like frame (AGAMA; Msun, kpc, km/s).

Components (split as in AGAMA's py/example_actions_nbody.py: spheroids as Multipole, disc as CylSpline):
  dm_spherical  Multipole, spherical: the ground-truth halo (the intended MW set-up has a spherical halo)
  alternatives: dm.ini (no symmetry, lmax 8, keeps the halo's bar response) and dm_axisymmetric.ini
  bulge  Multipole, no symmetry, lmax 8
  disc   CylSpline, no symmetry, mmax 8, with the grids used for Au18 (python_script/Nexus/get_potential*.py),
         zmax = Rmax so that the grid encloses the sphere r < Rmax
Radial ranges follow the particle extents: rmin about the 0.1% radius, rmax just beyond the outermost particle.
lmax = 8 is where the AGAMA reference documentation says N-body noise saturates the multipole expansion.
"""
import json
import time

import numpy as np
import h5py
import agama

from galaxies import galaxy, load_potential
from build_au18_potential import vcirc

agama.setUnits(mass=1, length=1, velocity=1)
G = galaxy("nbody")
POT_DIR = G["pot_dir"]
PARAMS = {
    "dm": dict(type="Multipole", symmetry="none", lmax=8, gridSizeR=30, rmin=0.5, rmax=600.0),
    "dm_spherical": dict(type="Multipole", symmetry="spherical", gridSizeR=30, rmin=0.5, rmax=600.0),
    "dm_axisymmetric": dict(type="Multipole", symmetry="axisymmetric", lmax=8, gridSizeR=30, rmin=0.5, rmax=600.0),
    "bulge": dict(type="Multipole", symmetry="none", lmax=8, gridSizeR=30, rmin=0.03, rmax=20.0),
    "disc": dict(type="CylSpline", symmetry="none", mmax=8, gridSizeR=30, gridSizez=30,
                 Rmin=0.1, Rmax=50.0, zmin=0.05, zmax=50.0),
}
SOURCE = {"dm": ("dm", 2), "dm_spherical": ("dm", 2), "dm_axisymmetric": ("dm", 2), "bulge": ("stars", 1),
          "disc": ("stars", 0)}


def main():
    POT_DIR.mkdir(exist_ok=True)
    data = {}
    for key in ("stars", "dm"):
        with h5py.File(G[key], "r") as f:
            data[key] = (np.column_stack([f[k][:].astype(np.float64) for k in "xyz"]),
                         f["mass"][:].astype(np.float64), f["comp"][:])
    summary = dict(units="Msun, kpc, km/s", frame="MW-like frame of nbody_stars_mwframe.hdf5", params=PARAMS)
    for name, kw in PARAMS.items():
        pos, m, comp = data[SOURCE[name][0]]
        s = comp == SOURCE[name][1]
        t0 = time.time()
        pot = agama.Potential(particles=(pos[s], m[s]), **kw)
        pot.export(str(POT_DIR / f"{name}.ini"))
        summary[name] = dict(N=int(s.sum()), mass=float(m[s].sum()), build_s=round(time.time() - t0, 1))
        print(f"{name:13s}: {s.sum():9,} particles, {m[s].sum():.3e} Msun, built in {time.time() - t0:.1f} s")
    comps, total = load_potential("nbody")
    R = np.array([2.0, 4.0, 8.0, 12.0])
    summary["vc"] = {k: dict(zip(map(str, R), map(float, vcirc(p, R)))) for k, p in {**comps, "total": total}.items()}
    with open(POT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    for k, v in summary["vc"].items():
        print(f"v_c {k:6s}: " + "  ".join(f"R={r}: {x:6.1f}" for r, x in v.items()))


if __name__ == "__main__":
    main()

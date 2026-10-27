"""Registry of the mock galaxies, so that the analysis scripts can run on any of them (--galaxy NAME).

Every galaxy provides a star file in the MW-like frame (x, y, z, vx, vy, vz, mass; kpc, km/s, Msun;
bar at -27 deg, disc spinning clockwise seen from +z) and a directory of AGAMA potential components.
"""
import numpy as np
import agama

from common import ROOT, MOCK_DIR

GALAXIES = {
    "au18": dict(label="Auriga Au18", prefix="au18",
                 stars=MOCK_DIR / "au18_stars_mwframe.hdf5",
                 pot_dir=MOCK_DIR / "au18_potential", components=("dm", "stars", "gas", "bh"),
                 truth=ROOT / "results" / "au18_ground_truth.json",
                 checks=ROOT / "results" / "au18_potential_checks.json"),
    "nbody": dict(label="N-body (t = 7)", prefix="nbody",
                  stars=MOCK_DIR / "nbody_stars_mwframe.hdf5", dm=MOCK_DIR / "nbody_dm_mwframe.hdf5",
                  pot_dir=MOCK_DIR / "nbody_potential", components=("dm", "bulge", "disc"),
                  truth=ROOT / "results" / "nbody_ground_truth.json",
                  checks=ROOT / "results" / "nbody_potential_checks.json"),
}


def galaxy(name):
    if name not in GALAXIES:
        raise ValueError(f"unknown galaxy {name!r}; choose from {list(GALAXIES)}")
    return GALAXIES[name]


def add_argument(parser):
    parser.add_argument("--galaxy", default="au18", choices=list(GALAXIES), help="mock galaxy (default au18)")


def load_potential(name="au18"):
    """Components (dict) and their sum, from the galaxy's exported coefficient files."""
    agama.setUnits(mass=1, length=1, velocity=1)
    g = galaxy(name)
    comps = {k: agama.Potential(str(g["pot_dir"] / f"{k}.ini")) for k in g["components"]}
    return comps, agama.Potential(*comps.values())


def stellar_components(name):
    """Names of the potential components that represent the star particles (for density checks)."""
    return ("stars",) if name == "au18" else ("bulge", "disc")

"""Registry of the mock galaxies, so that the analysis scripts can run on any of them (--galaxy NAME).

Every galaxy provides a star file in the MW-like frame (x, y, z, vx, vy, vz, mass; kpc, km/s, Msun;
bar at -27 deg, disc spinning clockwise seen from +z) and a directory of AGAMA potential components.
"""
import re

import numpy as np
import agama

from common import ROOT, MOCK_DIR, RESULTS_DIR
import nbody as _nbody

GALAXIES = {
    "au18": dict(label="Auriga Au18", prefix="au18",
                 stars=MOCK_DIR / "au18_stars_mwframe.hdf5",
                 pot_dir=MOCK_DIR / "au18_potential", components=("dm", "stars", "gas", "bh"),
                 truth=ROOT / "results" / "au18_ground_truth.json",
                 checks=ROOT / "results" / "au18_potential_checks.json"),
}

# The N-body run is one secularly evolving galaxy: each snapshot time is its own test galaxy, with its
# own pattern speed (see measure_nbody_pattern_speed_evolution.py). "nbody" is the fiducial snapshot
# t = DEFAULT_TIME and keeps the original un-tagged file names, so earlier results stay valid;
# "nbody_t<T>" (e.g. nbody_t5.5) tags every file with its time.
NBODY_DEFAULT_TIME = _nbody.DEFAULT_TIME


def nbody_name(t=NBODY_DEFAULT_TIME):
    """Galaxy name of the N-body snapshot at model time t."""
    return "nbody" if abs(float(t) - NBODY_DEFAULT_TIME) < 1e-9 else f"nbody_t{float(t):g}"


def _nbody_entry(t):
    tag = "" if abs(float(t) - NBODY_DEFAULT_TIME) < 1e-9 else f"_t{float(t):g}"
    return dict(label=f"N-body (t = {float(t):g})", prefix=f"nbody{tag}", time=float(t),
                stars=MOCK_DIR / f"nbody{tag}_stars_mwframe.hdf5",
                dm=MOCK_DIR / f"nbody{tag}_dm_mwframe.hdf5",
                # spherical DM halo, as in the intended MW set-up (dm.ini / dm_axisymmetric.ini are alternatives)
                pot_dir=MOCK_DIR / f"nbody{tag}_potential", components=("dm_spherical", "bulge", "disc"),
                truth=RESULTS_DIR / f"nbody{tag}_ground_truth.json",
                checks=RESULTS_DIR / f"nbody{tag}_potential_checks.json")


GALAXIES["nbody"] = _nbody_entry(NBODY_DEFAULT_TIME)

# The real Milky Way: no star file and no truth, only the trial-potential family (build_mw_potential.py).
GALAXIES["mw"] = dict(label="Milky Way (Hunter+24)", prefix="mw", pot_dir=MOCK_DIR / "mw_potential",
                      components=("MW_full",), truth=RESULTS_DIR / "mw_ground_truth_DOES_NOT_EXIST.json")


def galaxy(name):
    """Galaxy entry by name: "au18", "nbody", or "nbody_t<T>" for any snapshot time T."""
    if name in GALAXIES:
        return GALAXIES[name]
    if re.fullmatch(r"au\d+", str(name)):          # any Auriga halo, same layout as au18
        n = str(name)[2:]
        GALAXIES[name] = dict(label=f"Auriga Au{n}", prefix=name,
                              stars=MOCK_DIR / f"{name}_stars_mwframe.hdf5",
                              pot_dir=MOCK_DIR / f"{name}_potential",
                              components=("dm", "stars", "gas", "bh"),
                              truth=RESULTS_DIR / f"{name}_ground_truth.json",
                              checks=RESULTS_DIR / f"{name}_potential_checks.json")
        return GALAXIES[name]
    if str(name).startswith("nbody_t"):
        try:
            t = float(str(name)[len("nbody_t"):])
        except ValueError:
            raise ValueError(f"cannot read a snapshot time from galaxy name {name!r} (expected e.g. nbody_t5.5)")
        GALAXIES[name] = _nbody_entry(t)       # cache, so --galaxy choices stay informative
        return GALAXIES[name]
    raise ValueError(f"unknown galaxy {name!r}; choose from {list(GALAXIES)} or nbody_t<T>")


def _galaxy_name(value):
    galaxy(value)        # argparse type: validates, and registers an nbody_t<T> entry
    return value


def add_argument(parser):
    parser.add_argument("--galaxy", default="au18", type=_galaxy_name, metavar="NAME",
                        help="mock galaxy: au18, nbody, or nbody_t<T> for an N-body snapshot (default au18)")


def add_time_argument(parser, default=NBODY_DEFAULT_TIME):
    """--time T for the scripts that work on one N-body snapshot (model time, not Gyr)."""
    parser.add_argument("--time", type=float, default=default, metavar="T",
                        help=f"N-body snapshot model time (default {default:g}; 1 unit = 0.9778 Gyr)")


def load_potential(name="au18"):
    """Components (dict) and their sum, from the galaxy's exported coefficient files."""
    agama.setUnits(mass=1, length=1, velocity=1)
    g = galaxy(name)
    comps = {k: agama.Potential(str(g["pot_dir"] / f"{k}.ini")) for k in g["components"]}
    return comps, agama.Potential(*comps.values())


def stellar_components(name):
    """Names of the potential components that represent the star particles (for density checks)."""
    return ("stars",) if name == "au18" else ("bulge", "disc")

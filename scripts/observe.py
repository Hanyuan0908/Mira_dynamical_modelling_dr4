"""Mock observations: Galactocentric <-> Galactic (l, b, d, pm_l*, pm_b, v_los) with AGAMA,
and observational noise. Frame and solar parameters as in common.py (astropy-default
Galactocentric axes); proper motions in mas/yr."""
import numpy as np
import agama

from common import R0_KPC, Z_SUN_KPC, V_SUN_KMS

K_PM = 4.740470446          # km/s/kpc per mas/yr
FRAME = dict(galcen_distance=R0_KPC, z_sun=Z_SUN_KPC, galcen_v_sun=V_SUN_KMS)


def to_galactic(x, y, z, vx, vy, vz):
    """Galactocentric (kpc, km/s) -> dict of l, b [deg], d [kpc], pml, pmb [mas/yr], vlos [km/s]."""
    l, b, d, pml, pmb, vlos = agama.getGalacticFromGalactocentric(
        *(np.asarray(a, dtype=np.float64) for a in (x, y, z, vx, vy, vz)), **FRAME)
    return dict(l=np.degrees(l) % 360.0, b=np.degrees(b), d=d, pml=pml / K_PM, pmb=pmb / K_PM, vlos=vlos)


def to_galactocentric(l, b, d, pml, pmb, vlos):
    """Inverse of to_galactic; returns x, y, z [kpc], vx, vy, vz [km/s]."""
    return agama.getGalactocentricFromGalactic(
        np.radians(l), np.radians(b), np.asarray(d, dtype=np.float64),
        np.asarray(pml, dtype=np.float64) * K_PM, np.asarray(pmb, dtype=np.float64) * K_PM,
        np.asarray(vlos, dtype=np.float64), **FRAME)


def scatter_distance(d, frac_err, rng):
    """Gaussian fractional distance errors, sigma_d = frac_err * d; non-positive draws are redrawn."""
    d_obs = d + frac_err * d * rng.standard_normal(len(d))
    bad = d_obs <= 0
    while bad.any():
        d_obs[bad] = d[bad] + frac_err * d[bad] * rng.standard_normal(bad.sum())
        bad = d_obs <= 0
    return d_obs

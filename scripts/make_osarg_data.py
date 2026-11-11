"""The real data: the Zhang et al. (2024) low-amplitude LPV (OSARG) sample, written in the mock-file format so that
every pipeline script reads it unchanged.

Input: data/osarg_zhang2024/gaia_dr3_lpv_withdist.fits (Zenodo 13285574; 107,845 LA-LPV candidates with the paper's
PL distances). The paper's final cuts are applied here: M_JKs_sigma < 0.4 (45,819 stars) and a Gaia DR3 RV
(33,704), both reproduced exactly. M_JKs_sigma is the distance-modulus error (the paper's sigma_M < 0.4 is its
~20% distance cut).

Representative distance error (user, 2026-11-11): the MEDIAN M_JKs_sigma of that 33,704-star sample, stored as
the file's sigma_DM attribute, which is the single number the pipeline smears the library with. Per-star values
go in DM_err, as in the Test 3 mocks.

Positions and velocities: (l, b, luminosity_dist, pm_l*, pm_b, RV) -> Galactocentric with observe.to_galactocentric
(R0 8.122 kpc, z_sun 20.8 pc, v_sun (12.9, 245.6, 7.78) km/s, the frame of every mock and of the paper). Proper
motions are rotated to Galactic with astropy. As in the mocks, only stars at observed R < R_MAX = 8 kpc are kept.

Usage: python make_osarg_data.py
Output: mock/osarg_zhang2024_observed_R8.hdf5
"""
import numpy as np
import h5py
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from common import ROOT, MOCK_DIR, rel
import observe

SRC = ROOT / "data" / "osarg_zhang2024" / "gaia_dr3_lpv_withdist.fits"
OUT = MOCK_DIR / "osarg_zhang2024_observed_R8.hdf5"
SIGMA_M_MAX = 0.4
R_MAX = 8.0


def main():
    d = fits.open(SRC)[1].data
    sig = np.asarray(d["M_JKs_sigma"], float)
    rv = np.asarray(d["radial_velocity"], float)
    cut1 = sig < SIGMA_M_MAX
    final = cut1 & np.isfinite(rv)
    print(f"{rel(SRC)}: {len(d):,} stars; M_JKs_sigma < {SIGMA_M_MAX}: {cut1.sum():,}; with RV: {final.sum():,} "
          f"(Zhang et al. 2024: 45,819 and 33,704)")
    sigma_dm = float(np.median(sig[final]))
    frac = np.log(10) / 5 * sig[final]
    print(f"distance-modulus error of the final sample: median {sigma_dm:.4f} mag "
          f"(fractional {100 * np.log(10) / 5 * sigma_dm:.2f}%), mean fractional {100 * frac.mean():.2f}%, "
          f"16-84% {100 * np.percentile(frac, 16):.2f}-{100 * np.percentile(frac, 84):.2f}%")

    s = d[final]
    c = SkyCoord(ra=s["ra"] * u.deg, dec=s["dec"] * u.deg, pm_ra_cosdec=s["pmra"] * u.mas / u.yr,
                 pm_dec=s["pmdec"] * u.mas / u.yr, frame="icrs").galactic
    l, b = c.l.deg, c.b.deg
    pml, pmb = c.pm_l_cosb.to_value(u.mas / u.yr), c.pm_b.to_value(u.mas / u.yr)
    dist = np.asarray(s["luminosity_dist"], float)
    vlos = np.asarray(s["radial_velocity"], float)
    x, y, z, vx, vy, vz = observe.to_galactocentric(l, b, dist, pml, pmb, vlos)
    R = np.hypot(x, y)
    keep = R < R_MAX
    print(f"observed R < {R_MAX:g} kpc: {keep.sum():,} stars (R < 6: {(R < 6).sum():,})")

    # check against astropy's own transformation (same solar parameters)
    from astropy.coordinates import Galactocentric, CartesianDifferential
    from common import R0_KPC, Z_SUN_KPC, V_SUN_KMS
    ga = SkyCoord(l=l[:2000] * u.deg, b=b[:2000] * u.deg, distance=dist[:2000] * u.kpc,
                  pm_l_cosb=pml[:2000] * u.mas / u.yr, pm_b=pmb[:2000] * u.mas / u.yr,
                  radial_velocity=vlos[:2000] * u.km / u.s, frame="galactic").transform_to(
        Galactocentric(galcen_distance=R0_KPC * u.kpc, z_sun=Z_SUN_KPC * u.kpc,
                       galcen_v_sun=CartesianDifferential(V_SUN_KMS * u.km / u.s)))
    dv = np.max(np.abs(np.column_stack([ga.v_x.value, ga.v_y.value, ga.v_z.value])
                       - np.column_stack([vx[:2000], vy[:2000], vz[:2000]])))
    dx = np.max(np.abs(np.column_stack([ga.x.value, ga.y.value, ga.z.value])
                       - np.column_stack([x[:2000], y[:2000], z[:2000]])))
    print(f"AGAMA vs astropy Galactocentric, 2000 stars: max |dx| {dx:.1e} kpc, max |dv| {dv:.1e} km/s")

    k = keep
    with h5py.File(OUT, "w") as f:
        for name, a in (("x", x), ("y", y), ("z", z), ("vx", vx), ("vy", vy), ("vz", vz), ("l", l), ("b", b),
                        ("d", dist), ("pml", pml), ("pmb", pmb), ("vlos", vlos), ("DM_err", sig[final]),
                        ("DM_obs", np.asarray(s["DM"], float))):
            f[name] = np.asarray(a, float)[k]
        f["source_id"] = np.asarray(s["source_id"])[k]
        f.attrs["sigma_DM"] = sigma_dm
        f.attrs["R_max_obs"] = R_MAX
        f.attrs["source"] = "Zhang et al. 2024 (MNRAS 534, 1175), Zenodo 13285574"
        f.attrs["cuts"] = f"M_JKs_sigma < {SIGMA_M_MAX}, Gaia DR3 RV, observed R < {R_MAX:g} kpc"
        f.attrs["sigma_DM_definition"] = "median M_JKs_sigma of the 33,704-star final sample (before the R cut)"
    print(f"wrote {rel(OUT)}: {k.sum():,} stars, sigma_DM = {sigma_dm:.4f}")
    Rk, vphi = R[k], -(x * vy - y * vx)[k] / R[k]
    for lo in (2, 4, 6):
        m = (Rk > lo) & (Rk < lo + 2) & (np.abs(z[k]) < 1)
        print(f"  {lo}-{lo + 2} kpc, |z| < 1: {m.sum():6d} stars, median v_phi {np.median(vphi[m]):6.1f} km/s "
              f"(clockwise positive)")


if __name__ == "__main__":
    main()

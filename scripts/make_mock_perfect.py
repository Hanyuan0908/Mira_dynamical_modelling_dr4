"""Test 1 mock ('perfect' scenario): N random Au18 star particles with Galactocentric R < R_MAX,
no observational errors and no selection function.

Writes true Galactocentric phase space plus the corresponding error-free observables
(l, b, d, pm_l*, pm_b, v_los) seen from the Sun in the MW-like frame."""
import argparse

import numpy as np
import h5py

from common import MOCK_DIR, R0_KPC, Z_SUN_KPC, V_SUN_KMS, BAR_ANGLE_DEG
from prepare_au18_mwframe import OUT_FILE as AU18_FILE
import observe

N_STARS = 50_000
R_MAX = 8.0      # kpc, cylindrical Galactocentric radius
SEED = 2026


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, default=N_STARS)
    p.add_argument("--rmax", type=float, default=R_MAX)
    p.add_argument("--seed", type=int, default=SEED)
    args = p.parse_args()
    out = MOCK_DIR / f"mock_perfect_R{args.rmax:g}_N{args.n // 1000}k.hdf5"

    with h5py.File(AU18_FILE, "r") as f:
        R = np.hypot(f["x"][:], f["y"][:])
        pool = np.flatnonzero(R < args.rmax)
        rng = np.random.default_rng(args.seed)
        idx = np.sort(rng.choice(pool, size=args.n, replace=False))
        cols = {k: f[k][:][idx] for k in ("x", "y", "z", "vx", "vy", "vz", "mass", "age", "Z", "ParticleIDs")}
        bar_half_length = f.attrs["bar_half_length_kpc"]

    obs = observe.to_galactic(*(cols[k] for k in ("x", "y", "z", "vx", "vy", "vz")))
    with h5py.File(out, "w") as f:
        for k, v in {**cols, **obs}.items():
            f[k] = v
        f.attrs.update(dict(
            description="Perfect mock: random Au18 star particles, no errors, no selection. "
                        "Galactocentric x..vz [kpc, km/s]; l, b [deg]; d [kpc]; pml, pmb [mas/yr]; "
                        "vlos [km/s]; mass [Msun]; age [Gyr].",
            source=AU18_FILE.name, n=args.n, R_max=args.rmax, seed=args.seed, pool_size=len(pool),
            R0=R0_KPC, z_sun=Z_SUN_KPC, v_sun=V_SUN_KMS, bar_angle_deg=BAR_ANGLE_DEG,
            bar_half_length_kpc=bar_half_length))
    z = np.abs(cols["z"])
    print(f"wrote {out.name}: {args.n:,} of {len(pool):,} stars with R < {args.rmax:g} kpc; "
          f"|z| > 1 kpc: {np.mean(z > 1):.1%}, |z| > 3 kpc: {np.mean(z > 3):.1%}; "
          f"d range {obs['d'].min():.2f}-{obs['d'].max():.2f} kpc")


if __name__ == "__main__":
    main()

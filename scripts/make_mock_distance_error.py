"""Test 2b mock: the perfect N-body mock with Mira-like distance errors.

The error is the median distance-modulus error of the Gaia DR3 Mira catalogue (data/gaia_dr3_lpv_mira.csv) after the
usual cut sigma_DM < SIGMA_DM_CUT, applied to every star as a Gaussian in distance modulus (log-normal in distance):
DM_obs = DM_true + sigma_DM * N(0, 1). Angles (l, b), proper motions and v_los are kept exact; the observed
Galactocentric phase space is recomputed from (l, b, d_obs, pm_l*, pm_b, v_los), so positions and tangential velocities
are both smeared along the line of sight. The stars are the same 50k as in the perfect mock (selected on their true
R < 8 kpc); some now have an observed R > 8 kpc.

Writes mock/nbody_mock_distance_error_sigmaDM<sigma>_R8_N50k.hdf5 with the observed values under the usual names
(x, y, z, vx, vy, vz, d) so existing readers work unchanged, the true values as <name>_true, and DM_true, DM_obs, DM_err.
"""
import argparse

import numpy as np
import pandas as pd
import h5py

from common import MOCK_DIR, MIRA_DR3_FILE, R0_KPC, Z_SUN_KPC, V_SUN_KMS
import observe

PERFECT = MOCK_DIR / "nbody_mock_perfect_R8_N50k.hdf5"
SIGMA_DM_CUT = 0.4
SEED = 2029
PHASE = ("x", "y", "z", "vx", "vy", "vz")


def catalogue_sigma_dm(cut=SIGMA_DM_CUT):
    """Median DM_err of the DR3 Miras with DM_err < cut, and how many stars that is."""
    s = pd.read_csv(MIRA_DR3_FILE, usecols=["DM", "DM_err"]).dropna()
    s = s[s.DM_err < cut]
    return float(np.median(s.DM_err)), len(s)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sigma-dm", type=float, default=None, help="default: catalogue median with DM_err < 0.4")
    p.add_argument("--seed", type=int, default=SEED)
    args = p.parse_args()
    sigma_dm, n_cat = catalogue_sigma_dm()
    if args.sigma_dm is not None:
        sigma_dm = args.sigma_dm
    out = MOCK_DIR / f"nbody_mock_distance_error_sigmaDM{sigma_dm:.3f}_R8_N50k.hdf5"

    with h5py.File(PERFECT, "r") as f:
        cols = {k: f[k][:] for k in f.keys()}
        attrs = dict(f.attrs)
    # the stored observables reproduce the stored phase space (float32) when converted back at the true distance
    back = observe.to_galactocentric(*(cols[k] for k in ("l", "b", "d", "pml", "pmb", "vlos")))
    roundtrip = max(float(np.max(np.abs(a - cols[k]))) for a, k in zip(back, PHASE))

    rng = np.random.default_rng(args.seed)
    d_obs = observe.scatter_distance_modulus(cols["d"], sigma_dm, rng)
    obs = observe.to_galactocentric(cols["l"], cols["b"], d_obs, cols["pml"], cols["pmb"], cols["vlos"])
    with h5py.File(out, "w") as f:
        for k, v in cols.items():
            if k in PHASE or k == "d":
                f[f"{k}_true"] = v
            else:
                f[k] = v
        for k, v in zip(PHASE, obs):
            f[k] = v
        f["d"] = d_obs
        f["DM_true"] = 5 * np.log10(cols["d"] * 100.0)
        f["DM_obs"] = 5 * np.log10(d_obs * 100.0)
        f["DM_err"] = np.full(len(d_obs), sigma_dm)
        attrs.update(dict(
            description="Perfect N-body mock with Mira-like distance errors: DM_obs = DM_true + sigma_DM N(0,1) for "
                        "every star; l, b, pm, v_los exact; x..vz, d observed (true values in <name>_true). Same stars "
                        "as the perfect mock (selected on true R < 8 kpc).",
            source=PERFECT.name, sigma_DM=sigma_dm, frac_distance_error=np.log(10) / 5 * sigma_dm,
            sigma_DM_from=f"median DM_err of {MIRA_DR3_FILE.name} with DM_err < {SIGMA_DM_CUT} ({n_cat} stars)",
            seed_distance_error=args.seed, R0=R0_KPC, z_sun=Z_SUN_KPC, v_sun=V_SUN_KMS))
        f.attrs.update(attrs)

    x, y, z, vx, vy, vz = obs
    R_obs, R_true = np.hypot(x, y), np.hypot(cols["x"], cols["y"])
    dpos = np.sqrt((x - cols["x"])**2 + (y - cols["y"])**2 + (z - cols["z"])**2)
    frac = d_obs / cols["d"] - 1
    print(f"wrote {out.name}: sigma_DM = {sigma_dm:.4f} mag (fractional distance error {np.log(10) / 5 * sigma_dm:.4f}); "
          f"round trip at true d: max |diff| {roundtrip:.2e}")
    print(f"  d_obs / d_true - 1: median {np.median(frac):+.4f}, std {np.std(frac):.4f}, 16-84% "
          f"{np.percentile(frac, 16):+.4f} {np.percentile(frac, 84):+.4f}")
    print(f"  position shift: median {np.median(dpos):.2f} kpc, 84% {np.percentile(dpos, 84):.2f} kpc; "
          f"observed R > 8 kpc: {np.mean(R_obs > 8):.1%}; |R_obs - R_true| median {np.median(np.abs(R_obs - R_true)):.2f} kpc")
    print(f"  stars with |z| > 3 kpc: true {np.mean(np.abs(cols['z']) > 3):.1%}, observed {np.mean(np.abs(z) > 3):.1%}")


if __name__ == "__main__":
    main()

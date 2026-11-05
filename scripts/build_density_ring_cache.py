"""Cache the OBSERVED density m = 2 coefficients per ring, to allow density-defined windows.

The kinematic windows tried so far all leave a bar-strength residual. A window defined by the
observed DENSITY profile is a different idea worth testing: the bar is a density structure, the
observed density is available in the same sample, and selecting rings by it may pick a region whose
compression under smearing is more nearly universal.

Same galaxies, angles, realisations and seeds as scan_bar_angle_windows, so the two caches align
element for element.

Usage: python build_density_ring_cache.py [--n-real 8]
Cache: results/bar_angle_density_cache.npz
"""
import argparse

import numpy as np
import h5py

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy
import observe, selection, steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
from scan_bar_angle_windows import SNAPS, ANGLES, FRAC

RING_C = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
OUT = RESULTS_DIR / "bar_angle_density_cache.npz"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-real", type=int, default=8)
    args = ap.parse_args()
    Z = np.full((len(SNAPS), len(ANGLES), args.n_real, len(RING_C)), np.nan, complex)
    for gi, gn in enumerate(SNAPS):
        g = galaxy(gn)
        tab = selection.load_tables(
            MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{FRAC * 5 / np.log(10):.3f}_N50k.hdf5")
        with h5py.File(g["stars"], "r") as f:
            W = np.column_stack([f[k][:].astype(np.float64)
                                 for k in ("x", "y", "z", "vx", "vy", "vz")])
        for ai, a in enumerate(ANGLES):
            o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - a)).T)
            seed0 = 20_000 + 7919 * ai + 104_729 * gi      # identical to the kinematic cache
            for s in range(args.n_real):
                d = observe.scatter_distance_modulus(o["d"], FRAC * 5 / np.log(10),
                                                     np.random.default_rng(seed0 + s))
                w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"],
                                                              o["pmb"], o["vlos"]))
                S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
                sel = w[np.random.default_rng(seed0 + 500_000 + s).random(len(S)) < np.minimum(S, 1)]
                wb = ssl.rotate_z(sel, np.radians(BAR_ANGLE_DEG))
                R = np.hypot(wb[:, 0], wb[:, 1])
                phi = np.arctan2(wb[:, 1], wb[:, 0])
                thin = np.abs(wb[:, 2]) < 1.5
                for i in range(len(RING_C)):
                    m = thin & (R >= fs.RING_EDGES[i]) & (R < fs.RING_EDGES[i + 1])
                    if m.sum() >= 50:
                        Z[gi, ai, s, i] = np.sum(np.exp(2j * phi[m])) / m.sum()
            print(f"  {gn:14s} phi {a:5.1f} done", flush=True)
    np.savez_compressed(OUT, Z=Z, snaps=SNAPS, angles=ANGLES, Rc=RING_C)
    print(f"wrote {rel(OUT)}")


if __name__ == "__main__":
    main()

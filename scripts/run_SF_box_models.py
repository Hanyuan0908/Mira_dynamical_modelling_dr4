"""Integrate the N-body mock in four bar variants at the true Omega_b and phi_b (S = 1.25, 0.8; bar_frac = 1.4, 0.7)
and save the library box sums (0.5 kpc Cartesian boxes, bar frame) for square-pixel maps. The truth model's box
sums come from the global search (point 000). Output: ~/mira_dr4_runs/test1_nbody/maps_S<S>_f<f>.npz
"""
import numpy as np
import h5py

import steady_state_likelihood as ssl
from test1_nbody_qmc import ARR_DIR, MOCK, TRUTH

VARIANTS = [(1.25, 1.0), (0.8, 1.0), (1.0, 1.4), (1.0, 0.7)]


def main():
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    ssl.ACCURACY = 1e-5
    for S, fr in VARIANTS:
        out = ARR_DIR / f"maps_S{S:g}_f{fr:g}.npz"
        if out.exists():
            continue
        _, arr = ssl.evaluate(stars, TRUTH["Omega_abs"], TRUTH["phi_b_deg"], S, fr, galaxy="nbody", extended=True)
        np.savez(out, **{k: v for k, v in arr.items() if k.startswith("box_")})
        print(f"saved {out.name}", flush=True)


if __name__ == "__main__":
    main()

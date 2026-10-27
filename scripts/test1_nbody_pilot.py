"""Test 1 on the N-body galaxy, pilot: the truth (twice, at two integrator accuracies: noise floor) and
single-parameter offsets, each scored with the per-star KDE, the per-star local Gaussian and the
moment chi^2 from the same integrations. Results: results/test1_nbody/pilot.csv and per-point arrays
results/test1_nbody/pilot_NN.npz (per-star log p, local mu and C, cell moments).
"""
import csv
import time

import numpy as np
import h5py

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl

OUT_DIR = ROOT / "results" / "test1_nbody"
MOCK = MOCK_DIR / "nbody_mock_perfect_R8_N50k.hdf5"
TRUTH = dict(Omega_abs=26.02, phi_b_deg=27.0, S=1.0, bar_frac=1.0)
ACC = 1e-5
POINTS = [("truth", {}, ACC), ("truth, acc 3e-6 (noise)", {}, 3e-6),
          ("Omega=20", dict(Omega_abs=20.0), ACC), ("Omega=32", dict(Omega_abs=32.0), ACC),
          ("phi=17", dict(phi_b_deg=17.0), ACC), ("phi=37", dict(phi_b_deg=37.0), ACC),
          ("f=0.7", dict(bar_frac=0.7), ACC), ("f=1.4", dict(bar_frac=1.4), ACC),
          ("S=0.8", dict(S=0.8), ACC), ("S=1.25", dict(S=1.25), ACC)]
FIELDS = ["i", "label", "Omega_abs", "phi_b_deg", "S", "bar_frac", "accuracy", "logL", "logL_gauss",
          "chi2_means", "chi2_means_R1", "chi2_meansdisp_R1", "chi2_means_ncells", "tau_median",
          "t_integrate_s", "t_likelihood_s"]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    t_start = time.time()
    with open(OUT_DIR / "pilot.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        wr.writeheader()
        for i, (label, change, acc) in enumerate(POINTS):
            p = {**TRUTH, **change}
            ssl.ACCURACY = acc
            res, arr = ssl.evaluate(stars, p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"],
                                    galaxy="nbody", extended=True)
            np.savez(OUT_DIR / f"pilot_{i:02d}.npz", **arr)
            wr.writerow(dict(i=i, label=label, accuracy=acc, **p, **res))
            fh.flush()
            print(f"[{i + 1}/{len(POINTS)}] {label:24s} logL_kde {res['logL']:.1f}  logL_gauss {res['logL_gauss']:.1f}  "
                  f"chi2_means {res['chi2_means']:.1f}  ({res['t_integrate_s'] + res['t_likelihood_s']:.0f} s, "
                  f"elapsed {(time.time() - t_start) / 60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()

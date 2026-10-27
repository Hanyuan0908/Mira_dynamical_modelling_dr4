"""Test 1 global search on the N-body galaxy: the truth plus N_QMC scrambled Sobol points in
(|Omega_b|, phi_b, ln S, ln bar_frac), every likelihood computed from the same integrations, and
everything needed to recompute other likelihoods later saved (no re-integration needed):

  results/test1_nbody/qmc_evals.csv           totals: KDE, local-Gaussian logL, moment chi^2 variants
  ~/mira_dr4_runs/test1_nbody/NNN.npz         per-star KDE / Gaussian log p, local mu and C, tau,
                                              polar cell moments (data and library), 0.5 kpc box sums
  ~/mira_dr4_runs/test1_nbody/traj_NNN.npy    full trajectories, bar frame, float32 (N, 41, 6),
                                              snapshots every 25 Myr over 1 Gyr (times in the npz)
Resumes if restarted (skips indices already in the CSV).
"""
import csv
import time
from pathlib import Path

import numpy as np
import h5py
from scipy.stats import qmc

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl

OUT_CSV = ROOT / "results" / "test1_nbody" / "qmc_evals.csv"
ARR_DIR = Path.home() / "mira_dr4_runs" / "test1_nbody"
MOCK = MOCK_DIR / "nbody_mock_perfect_R8_N50k.hdf5"
N_QMC = 64
SEED = 2027
ACCURACY = 1e-5
TRUTH = dict(Omega_abs=26.02, phi_b_deg=27.0, S=1.0, bar_frac=1.0)
BOX = dict(Omega_abs=(18.0, 38.0), phi_b_deg=(14.0, 44.0), lnS=(np.log(0.75), np.log(1.45)),
           lnf=(np.log(0.55), np.log(1.6)))
FIELDS = ["idx", "kind", "Omega_abs", "phi_b_deg", "S", "bar_frac", "logL", "logL_gauss", "chi2_means",
          "chi2_means_R1", "chi2_meansdisp_R1", "chi2_means_ncells", "chi2_means_R1_ncells", "logp_median",
          "tau_median", "tau_min", "t_integrate_s", "t_likelihood_s"]


def design():
    u = qmc.Sobol(d=4, scramble=True, seed=SEED).random(N_QMC)
    lo = np.array([v[0] for v in BOX.values()])
    hi = np.array([v[1] for v in BOX.values()])
    p = lo + u * (hi - lo)
    pts = [dict(kind="truth", **TRUTH)]
    pts += [dict(kind="qmc", Omega_abs=a, phi_b_deg=b, S=float(np.exp(c)), bar_frac=float(np.exp(d)))
            for a, b, c, d in p]
    return pts


def main():
    ARR_DIR.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            done = {int(r["idx"]) for r in csv.DictReader(f)}
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    ssl.ACCURACY = ACCURACY
    pts = design()
    t_start = time.time()
    with open(OUT_CSV, "a", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if not done:
            wr.writeheader()
        for i, p in enumerate(pts):
            if i in done:
                continue
            res, arr = ssl.evaluate(stars, p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"],
                                    galaxy="nbody", extended=True, keep_traj=True)
            np.save(ARR_DIR / f"traj_{i:03d}.npy", arr.pop("traj"))
            np.savez(ARR_DIR / f"{i:03d}.npz", params=np.array([p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"]]),
                     accuracy=ACCURACY, **arr)
            wr.writerow(dict(idx=i, **p, **res))
            fh.flush()
            el = time.time() - t_start
            print(f"[{i + 1}/{len(pts)}] {p['kind']:5s} Omega={p['Omega_abs']:.2f} phi={p['phi_b_deg']:.1f} "
                  f"S={p['S']:.3f} f={p['bar_frac']:.3f}  kde {res['logL']:.1f}  gauss {res['logL_gauss']:.1f}  "
                  f"chi2 {res['chi2_means']:.1f}  ({res['t_integrate_s'] + res['t_likelihood_s']:.0f} s; "
                  f"elapsed {el / 60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()

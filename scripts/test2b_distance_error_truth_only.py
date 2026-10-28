"""Test 2b, cheap tolerance check: extra misfit at the TRUE potential caused by distance errors, per error level.

For each fractional distance error (same deviates as test2b_distance_error_tolerance.py), the smeared stars are
integrated in the true potential (one integration per level) and the current method's chi2 (proposed scheme, joint
bootstrap covariance) is computed against the smeared data. Also the library from the true positions against the
smeared data (no integration). Compared with the perfect case (chi2 at the truth = 1.50 N): the extra chi2 against
sqrt(2 N), the 1-sigma step of the tempered likelihood.
Numbers: results/test2b/distance_error_truth_only.json
"""
import json

import numpy as np

from common import ROOT
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test2b_distance_error_tolerance import FRACS, TRUTH, PROPOSED, observed_stars, trajectory

RESULTS = ROOT / "results" / "test2b" / "distance_error_truth_only.json"


def main():
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(50000, np.full(50000, 1 / 50000), size=fs.N_BOOT).astype(float)
    traj_true = np.load(old.ARR_DIR / "traj_000.npy")
    res = {}
    for frac in FRACS:
        w = observed_stars(frac)
        cd, cov, n, ok = fs.data_fit(w, TRUTH["phi_b_deg"], boot_w)
        out = {}
        for key, traj in (("orbits_from_observed_stars", trajectory(frac, "truth", TRUTH, w, print)),
                          ("orbits_from_true_positions", traj_true)):
            cl = fs.library_fit(traj, times, ok)
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            chi, nd = fs.scheme_chi2(cd, cl, cov, valid, PROPOSED)
            out[key] = dict(chi2=round(chi, 1), N_data=nd, chi2_per_N=round(chi / nd, 3))
        res[f"{100 * frac:.1f}%"] = out
        print(f"{100 * frac:5.1f}%  observed-star orbits: chi2/N {out['orbits_from_observed_stars']['chi2_per_N']}"
              f"   true-position orbits: chi2/N {out['orbits_from_true_positions']['chi2_per_N']}", flush=True)
    ref = res["0.0%"]["orbits_from_observed_stars"]
    for v in res.values():
        o = v["orbits_from_observed_stars"]
        o["extra_chi2_vs_perfect"] = round(o["chi2"] - ref["chi2"], 1)
        o["extra_in_tempered_logL"] = round(-0.5 * (o["chi2"] - ref["chi2"]) / np.sqrt(2 * ref["N_data"]), 2)
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

"""Test 2b: chi2 / N of the true potential with the co-added Fourier error (user's recipe for a real observation).

Data: the 7.5% distance-error mock at the true bar angle (coefficients saved by test2b_qmc_distance_error_run.py).
Error = C_MC + C_boot:
  C_MC    covariance of the coefficients over N_MC re-draws of every star's distance (DM_obs + sigma_DM N(0, 1)),
          same stars, l, b, pm, v_los;
  C_boot  bootstrap covariance (300) of the unsmeared data, as in Test 1 (fourier_stars_000.npz of the focused run);
          also, for comparison, the bootstrap of the observed data.
Hartlap factor with 300 samples. Libraries at the true potential from the observed stars (the current method):
observed again with the same distance error, and as it is. Proposed scheme.
Numbers: results/test2b/coadded_errors_truth_chi2.json
"""
import json

import numpy as np
import h5py

from common import ROOT
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import observe
from test1_nbody_qmc_focused_run import OUT_DIR as FOCUSED_DIR
from test2b_qmc_distance_error_run import MOCK, OUT_DIR

N_MC = 300
PHI = 27.0
RESULTS = ROOT / "results" / "test2b" / "coadded_errors_truth_chi2.json"


def main():
    s = np.load(OUT_DIR / "fourier_000.npz")
    cd = s["coef_data"]
    ok = s["n_stars"] >= fs.N_MIN_BLOCK
    sigma_dm = float(h5py.File(MOCK, "r").attrs["sigma_DM"])
    with h5py.File(MOCK, "r") as f:
        m = {k: f[k][:].astype(np.float64) for k in ("l", "b", "pml", "pmb", "vlos", "d")}
    rng = np.random.default_rng(21)
    mc = np.empty((N_MC,) + cd.shape)
    for k in range(N_MC):
        d = observe.scatter_distance_modulus(m["d"], sigma_dm, rng)
        obs = np.column_stack(observe.to_galactocentric(m["l"], m["b"], d, m["pml"], m["pmb"], m["vlos"]))
        b, X, v, _ = fs.prepare(ssl.rotate_z(obs, np.radians(PHI)))
        mc[k] = fs.fit(b, X, v, np.ones(len(b)), ok)
    dev = (mc - mc.mean(0)).reshape(N_MC, fs.NB, -1)
    C_mc = np.einsum("mbi,mbj->bij", dev, dev) / (N_MC - 1)
    p = np.load(FOCUSED_DIR / "fourier_stars_000.npz")
    C_boot_unsmeared = p["cov_data_joint"].reshape(fs.NB, 30, 30)
    ok_both = ok & (p["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK)
    covs = {"C_MC + C_boot(unsmeared data)": C_mc + C_boot_unsmeared,
            "C_MC + C_boot(observed data)": C_mc + s["cov_data_joint"],
            "C_boot(observed data) only": s["cov_data_joint"]}
    res = dict(N_MC=N_MC, sigma_DM=sigma_dm, chi2_per_N_truth={})
    for cname, C in covs.items():
        row = {}
        for lib in ("observed", "asis"):
            cl = s[f"coef_lib_{lib}"]
            valid = ok_both & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            chi, nd = fs.scheme_chi2(cd, cl, C, valid, fs.SCHEMES["fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"])
            row[lib] = dict(chi2_per_N=round(chi / nd, 3), N_data=nd)
        res["chi2_per_N_truth"][cname] = row
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

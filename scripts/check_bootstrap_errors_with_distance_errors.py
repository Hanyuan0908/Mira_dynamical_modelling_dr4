"""Are the bootstrap errors of the per-star Fourier coefficients still right when the data have distance errors?

Truth for the comparison: the scatter of the coefficients over N_INDEP independent 50k mocks drawn from the 1.49 M N-body
star particles (R < 8 kpc), each with its own distance-modulus errors (sigma_DM of the 7.5% mock), fitted at the true bar
angle exactly like the data (test1_nbody_qmc_focused_fourier_stars.py). The population value of the smeared statistic is
the same fit to all 1.49 M particles, smeared once. Compared with the bootstrap covariance of the 7.5% mock (saved by
test2b_qmc_distance_error_run.py at the truth, 300 bootstraps): variance ratios per |z| layer and field, and chi2 / N of
the mock against the population with either covariance (proposed scheme; about 1 if the errors are right). Also how much
the distance errors inflate the bootstrap errors relative to the perfect mock.
Numbers: results/test2b/bootstrap_errors_with_distance_errors.json
"""
import json

import numpy as np
import h5py

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import observe
from test1_nbody_qmc_focused_run import OUT_DIR as FOCUSED_DIR
from check_fourier_star_fit_vs_true_maps import POOL, chi2_per_N
from test2b_qmc_distance_error_run import MOCK, OUT_DIR

N_INDEP = 200
PHI = 27.0
RESULTS = ROOT / "results" / "test2b" / "bootstrap_errors_with_distance_errors.json"


def smear_and_rotate(w_mw, sigma_dm, rng):
    g = observe.to_galactic(*w_mw.T)
    d = observe.scatter_distance_modulus(g["d"], sigma_dm, rng)
    obs = np.column_stack(observe.to_galactocentric(g["l"], g["b"], d, g["pml"], g["pmb"], g["vlos"]))
    return ssl.rotate_z(obs, np.radians(PHI))


def main():
    sigma_dm = float(h5py.File(MOCK, "r").attrs["sigma_DM"])
    with h5py.File(POOL, "r") as f:
        pool = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    pool = pool[np.hypot(pool[:, 0], pool[:, 1]) < 8.0]
    s = np.load(OUT_DIR / "fourier_000.npz")
    c_mock, cov = s["coef_data"].reshape(fs.NB, -1), s["cov_data_joint"]
    ok = s["n_stars"] >= fs.N_MIN_BLOCK
    rng = np.random.default_rng(5)
    b, X, v, _ = fs.prepare(smear_and_rotate(pool, sigma_dm, rng))
    c_pop = fs.fit(b, X, v, np.ones(len(b)), ok).reshape(fs.NB, -1)
    ens = np.empty((N_INDEP, fs.NB, c_mock.shape[1]))
    for k in range(N_INDEP):
        sub = pool[rng.choice(len(pool), 50000, replace=False)]
        b, X, v, _ = fs.prepare(smear_and_rotate(sub, sigma_dm, rng))
        ens[k] = fs.fit(b, X, v, np.ones(len(b)), ok).reshape(fs.NB, -1)
    dev = ens - ens.mean(0)
    C_true = np.einsum("mbi,mbj->bij", dev, dev) / (N_INDEP - 1)
    ring, layer = np.divmod(np.arange(fs.NB), fs.NL)
    valid = ok & np.all(np.isfinite(c_mock), 1) & np.all(np.isfinite(c_pop), 1)
    use = valid & (ring < fs.N_RING_DEFAULT)
    perfect = np.load(FOCUSED_DIR / "fourier_stars_000.npz")["cov_data_joint"].reshape(fs.NB, 30, 30)
    res = dict(sigma_DM=sigma_dm, n_independent_mocks=N_INDEP, layers={})
    for z, lab in enumerate(["absz_0to1", "absz_1to2", "absz_2to3"]):
        blocks = np.flatnonzero(use & (layer == z))
        if len(blocks) == 0:
            continue
        ratio = np.diagonal(cov[blocks], axis1=1, axis2=2) / np.diagonal(C_true[blocks], axis1=1, axis2=2)
        infl = np.sqrt(np.diagonal(cov[blocks], axis1=1, axis2=2) / np.diagonal(perfect[blocks], axis1=1, axis2=2))
        res["layers"][lab] = dict(
            n_blocks=int(len(blocks)),
            variance_ratio_bootstrap_over_true={nm: round(float(np.median(ratio[:, f * 5:(f + 1) * 5])), 3)
                                                for f, nm in enumerate(fs.FIELD_NAMES)},
            error_inflation_vs_perfect_data={nm: round(float(np.median(infl[:, f * 5:(f + 1) * 5])), 3)
                                             for f, nm in enumerate(fs.FIELD_NAMES)},
            chi2_per_N_mock_vs_population_bootstrap_cov=round(chi2_per_N(c_mock - c_pop, cov, fs.N_BOOT, blocks), 3),
            chi2_per_N_mock_vs_population_true_cov=round(chi2_per_N(c_mock - c_pop, C_true, N_INDEP, blocks), 3))
    blocks = np.flatnonzero(use)
    res["all_R0to6"] = dict(chi2_per_N_mock_vs_population_bootstrap_cov=round(chi2_per_N(c_mock - c_pop, cov, fs.N_BOOT, blocks), 3),
                            chi2_per_N_mock_vs_population_true_cov=round(chi2_per_N(c_mock - c_pop, C_true, N_INDEP, blocks), 3))
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

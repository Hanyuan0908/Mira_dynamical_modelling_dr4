"""Test 2b, best case for de-smeared initial conditions: several orbits per star, with distances drawn from each star's
distance posterior using the TRUE N-body distribution as the prior.

For star i (exact l, b, pm, v_los; DM_obs with sigma_DM = 0.265 mag), on a grid of N_GRID distance moduli
DM_obs + sigma_DM * [-4, 4], the posterior is
  A, density prior:  p(DM) ∝ N(DM_obs | DM, sigma_DM) rho(x(d)) d^3        (positions only; ignores pm and v_los)
  B, full DF prior:  p(DM) ∝ N(DM_obs | DM, sigma_DM) f(x(d), v(d)) d^5    (uses the measured pm and v_los)
d^2 and d^4 are the volume elements of (l, b, d) and of (l, b, d, pm, v_los); one more power of d for the uniform DM
grid. With the true f, draws from B are distributed exactly like the true stars; A ignores how f depends on velocity.
rho and f are k-nearest-neighbour estimates (k = K_NN) from the N-body star particles with R < 8 kpc without the 50k mock
stars; f uses (x, v / V_SCALE). K_DRAW distances per star are drawn (uniform jitter inside the grid cell), i.e. K_DRAW
orbits per star of weight 1 / K_DRAW, integrated in the true potential at the true pattern speed (as in
test2b_smeared_initial_conditions.py).
Checks, against the data the likelihood would use (per-star Fourier fit, proposed scheme, joint covariance):
  t = 0: the drawn initial conditions against the perfect data (do the draws reproduce the true distribution in the
         moments we fit?), with the smeared stars themselves for reference;
  library: each library observed with the distance errors (and as it is) against the smeared data, all of t >= 0.25 Gyr
           and in time windows; reference: the truth library observed, 1.41.
Files: ~/mira_dr4_runs/test2b/posterior_draws.npz, traj_posterior_draws_{A,B}.npy (re-used if present).
Figure: figures/test2b_posterior_distance_draws_fourier_coefficients_absz_0to1.png
Numbers: results/test2b/posterior_distance_draws.json
"""
import json
import time

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

from common import ROOT, MOCK_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import observe
from plot_distance_error_effect_on_fourier_coefficients import SMEARED, stars, forward_smeared_library_fit
from plot_fourier_coefficients_data_vs_models import as_kms, FIELDS, TERMS, RC
from test2b_smeared_initial_conditions import OUT_DIR, TRAJ as TRAJ_SMEARED, window_fit

POOL = MOCK_DIR / "nbody_stars_mwframe.hdf5"
N_GRID = 41
K_DRAW = 5
K_NN = 64
V_SCALE = 40.0           # km/s per kpc in the 6D distance of the f estimate
SEED = 7
DRAWS = OUT_DIR / "posterior_draws.npz"
RESULTS = ROOT / "results" / "test2b" / "posterior_distance_draws.json"
WINDOWS = ((0.0, 0.1), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0))
PROPOSED = fs.SCHEMES["fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"]


def pool_without_mock(mock_index):
    with h5py.File(POOL, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    keep = np.hypot(w[:, 0], w[:, 1]) < 8.0
    keep[mock_index] = False
    return w[keep]


def posterior_draws(sigma_dm):
    """K_DRAW Galactocentric (MW-frame) initial conditions per star for priors A and B, and diagnostics."""
    with h5py.File(SMEARED, "r") as f:
        m = {k: f[k][:].astype(np.float64) for k in ("l", "b", "pml", "pmb", "vlos", "DM_obs", "d_true", "star_index")}
    n = len(m["l"])
    pool = pool_without_mock(m["star_index"].astype(int))
    u = np.linspace(-4.0, 4.0, N_GRID)
    DM = m["DM_obs"][:, None] + sigma_dm * u[None, :]
    d = 10.0 ** (0.2 * DM - 2.0)
    rep = lambda a: np.repeat(a, N_GRID)
    pts = np.column_stack(observe.to_galactocentric(rep(m["l"]), rep(m["b"]), d.ravel(), rep(m["pml"]), rep(m["pmb"]),
                                                    rep(m["vlos"])))
    t0 = time.time()
    r3 = cKDTree(pool[:, :3]).query(pts[:, :3], k=K_NN, workers=-1)[0][:, -1].reshape(n, N_GRID)
    scale = np.r_[np.ones(3), np.full(3, 1.0 / V_SCALE)]
    r6 = cKDTree(pool * scale).query(pts * scale, k=K_NN, workers=-1)[0][:, -1].reshape(n, N_GRID)
    t_knn = time.time() - t0
    loglike = -0.5 * u**2
    logp = dict(A=loglike - 3 * np.log(r3) + 3 * np.log(d), B=loglike - 6 * np.log(r6) + 5 * np.log(d))
    rng = np.random.default_rng(SEED)
    out = dict(t_knn_s=t_knn)
    for key, lp in logp.items():
        p = np.exp(lp - lp.max(axis=1, keepdims=True))
        cdf = np.cumsum(p, axis=1)
        cdf /= cdf[:, -1:]
        idx = (cdf[:, None, :] < rng.random((n, K_DRAW))[:, :, None]).sum(-1)
        dm = DM[np.arange(n)[:, None], idx] + sigma_dm * (u[1] - u[0]) * (rng.random((n, K_DRAW)) - 0.5)
        dd = 10.0 ** (0.2 * dm - 2.0)
        out[f"d_{key}"] = dd
        out[f"ic_{key}"] = np.column_stack(observe.to_galactocentric(
            *(np.repeat(m[k], K_DRAW) for k in ("l", "b")), dd.ravel(),
            *(np.repeat(m[k], K_DRAW) for k in ("pml", "pmb", "vlos"))))
        # posterior mean distance error relative to the truth, and the spread of the posterior
        mean_d = np.sum(p / p.sum(1, keepdims=True) * d, axis=1)
        out[f"frac_err_posterior_mean_{key}"] = mean_d / m["d_true"] - 1
    out["frac_err_observed"] = 10.0 ** (0.2 * m["DM_obs"] - 2.0) / m["d_true"] - 1
    return out


def integrate(ic_mw, path):
    if path.exists():
        return np.load(path), 0.0
    ssl.ACCURACY = old.ACCURACY
    t0 = time.time()
    pot = ssl.bar_aligned_potential(old.TRUTH["bar_frac"], old.TRUTH["S"], "nbody")
    _, traj = ssl.integrate(pot, ssl.rotate_z(ic_mw, np.radians(old.TRUTH["phi_b_deg"])), -abs(old.TRUTH["Omega_abs"]))
    np.save(path, traj)
    return traj, time.time() - t0


def chi2N(cd, cl, cov, valid, scheme=PROPOSED):
    v = valid & np.all(np.isfinite(cl), axis=(1, 2))
    chi, nd = fs.scheme_chi2(cd, cl, cov, v, scheme)
    return round(chi / nd, 3)


def main():
    set_style()
    sigma_dm = float(h5py.File(SMEARED, "r").attrs["sigma_DM"])
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    if DRAWS.exists():
        dr = dict(np.load(DRAWS))
    else:
        dr = posterior_draws(sigma_dm)
        np.savez(DRAWS, **dr)
    res = dict(sigma_DM=sigma_dm, K_DRAW=K_DRAW, K_NN=K_NN, V_SCALE=V_SCALE, t_knn_s=float(dr["t_knn_s"]))
    res["distance_error_vs_truth"] = {
        "observed": dict(median_abs=round(float(np.median(np.abs(dr["frac_err_observed"]))), 4)),
        **{f"posterior_mean_{k}": dict(median_abs=round(float(np.median(np.abs(dr[f"frac_err_posterior_mean_{k}"]))), 4),
                                       median=round(float(np.median(dr[f"frac_err_posterior_mean_{k}"])), 4))
           for k in "AB"}}

    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(50000, np.full(50000, 1 / 50000), size=fs.N_BOOT).astype(float)
    phi = old.TRUTH["phi_b_deg"]
    cd_p, cov_p, _, ok_p = fs.data_fit(stars(old.MOCK), phi, boot_w)          # perfect data
    cd_s, cov_s, _, ok_s = fs.data_fit(stars(SMEARED), phi, boot_w)           # smeared data
    valid_p = ok_p & np.all(np.isfinite(cd_p), axis=(1, 2))
    valid_s = ok_s & np.all(np.isfinite(cd_s), axis=(1, 2))

    # t = 0: do the drawn initial conditions reproduce the true distribution in the fitted moments?
    res["t0_initial_conditions_vs_perfect_data"] = {}
    for key, ic in (("smeared stars", stars(SMEARED)), ("draws A (density prior)", dr["ic_A"]),
                    ("draws B (full DF prior)", dr["ic_B"])):
        b, X, v, mm = fs.prepare(ssl.rotate_z(ic, np.radians(phi)))
        res["t0_initial_conditions_vs_perfect_data"][key] = chi2N(cd_p, fs.fit(b, X, v, np.ones(len(b)), ok_p, leverage=False),
                                                                  cov_p, valid_p)

    # libraries
    trajs, res["integration_s"] = {}, {}
    for key in "AB":
        trajs[key], res["integration_s"][key] = integrate(dr[f"ic_{key}"], OUT_DIR / f"traj_posterior_draws_{key}.npy")
    trajs["smeared ICs"] = np.load(TRAJ_SMEARED)
    trajs["true ICs"] = np.load(old.ARR_DIR / "traj_000.npy")
    libs = {}
    res["library_vs_smeared_data"] = {}
    for key, tr in trajs.items():
        obs = forward_smeared_library_fit(tr, times, ok_s, sigma_dm, phi)
        libs[key] = obs
        res["library_vs_smeared_data"][key] = dict(
            observed=chi2N(cd_s, obs, cov_s, valid_s), as_is=chi2N(cd_s, fs.library_fit(tr, times, ok_s), cov_s, valid_s),
            observed_by_window={f"{a:g}-{b:g} Gyr": chi2N(cd_s, window_fit(tr, times, ok_s, a, b, sigma_dm), cov_s, valid_s)
                                for a, b in WINDOWS})
    res["library_vs_smeared_data_per_field_observed"] = {}
    for f, name in enumerate(fs.FIELD_NAMES):
        for terms, lab in (([0], "m0"), (fs.NO_M0, "m24")):
            res["library_vs_smeared_data_per_field_observed"][f"{name} {lab}"] = {
                key: chi2N(cd_s, libs[key], cov_s, valid_s, (fs.N_RING_DEFAULT, {f: terms})) for key in libs}

    # figure: |z| < 1, smeared data vs the observed libraries
    shape = (fs.NR, fs.NL, 6, fs.NK)
    lines = [("true ICs", "true initial conditions (target)", BLUE_RAMP[5], "--"),
             ("B", "posterior draws, full DF prior", ACCENT, "-"), ("A", "posterior draws, density prior", INK, ":"),
             ("smeared ICs", "smeared stars", INK_MUTED, "-.")]
    fig, axes = plt.subplots(len(FIELDS), len(TERMS), figsize=(17.5, 12.5), gridspec_kw=dict(hspace=0.42, wspace=0.32))
    layer = 0
    use = valid_s.reshape(fs.NR, fs.NL)[:, layer]
    r6 = use & (np.arange(fs.NR) < fs.N_RING_DEFAULT)
    for r, (name, f, kind) in enumerate(FIELDS):
        e = np.sqrt(np.diagonal(cov_s, axis1=1, axis2=2)).reshape(shape)[:, layer, f]
        dat, err = as_kms(cd_s.reshape(shape)[:, layer, f], e, kind)
        for k, term in enumerate(TERMS):
            ax = axes[r, k]
            ax.axvspan(6.0, 8.0, color=INK_MUTED, alpha=0.12, lw=0)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            notes = []
            for key, lab, colr, ls in lines:
                mod, _ = as_kms(libs[key].reshape(shape)[:, layer, f], kind=kind)
                ax.plot(RC[use], mod[use, k], color=colr, ls=ls, lw=1.6, label=f"{lab}, observed" if (r, k) == (0, 0) else None)
                notes.append(f"{np.mean(((dat[r6, k] - mod[r6, k]) / err[r6, k])**2):.2f}")
            ax.errorbar(RC[use], dat[use, k], yerr=err[use, k], fmt="o", ms=3.2, color=INK, ecolor=INK_2, elinewidth=0.9,
                        capsize=0, zorder=4, label="smeared data" if (r, k) == (0, 0) else None)
            ax.text(0.02, 0.97, "χ²/N  target " + notes[0] + ",  B " + notes[1] + ",  A " + notes[2] + ",  smeared "
                    + notes[3], transform=ax.transAxes, fontsize=6.4, va="top", color=INK_2,
                    bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
            ax.set_xlim(0, 8)
            ax.tick_params(labelsize=7.5)
            if r == 0:
                ax.set_title(term, fontsize=10)
            if k == 0:
                ax.set_ylabel(f"{name}\n" + ("km/s" if kind == "mean" else "σ [km/s]"), fontsize=9.5)
            if r == len(FIELDS) - 1:
                ax.set_xlabel("R [kpc] (observed)", fontsize=9)
    fig.legend(loc="upper center", ncol=5, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"Orbits from posterior distance draws ({K_DRAW} per star, true N-body as prior) in the true potential, "
                 f"observed with the distance errors; |z| < 1 kpc\nχ²/N per panel against the smeared data (diagonal "
                 f"errors, R < 6 kpc); σ terms: m = 0 is σ₀, others c_m / 2σ₀", fontsize=10.5, y=0.997)
    savefig(fig, "test2b_posterior_distance_draws_fourier_coefficients_absz_0to1")
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

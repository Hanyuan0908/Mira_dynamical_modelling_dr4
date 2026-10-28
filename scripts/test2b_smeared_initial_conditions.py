"""Test 2b: can the orbit library start from distance-smeared stars?

Real data give only smeared initial conditions. Here the smeared mock (sigma_DM = 0.265 mag;
make_mock_distance_error.py) is rotated by the true bar angle and integrated in the true potential at the true pattern
speed (accuracy 1e-5, 1 Gyr, 41 snapshots, as in the global search). Against the smeared data (per-star Fourier fit,
300 bootstraps, proposed scheme and per field), compare:
  (a) the library from smeared initial conditions, as it is;
  (b) the same library observed again with the same distance errors (the errors then enter twice);
  (c) the truth library (from the true positions) observed with the errors: the target, what a correct treatment of
      the initial conditions should reach;
  (d) the truth library as it is (no correction).
Also chi2 / N of (a) and (b) in time windows, to see how far the smeared stars drift from their starting distribution.
Trajectory: ~/mira_dr4_runs/test2b/traj_smeared_ics_true_potential.npy (re-used if present).
Figures: figures/test2b_smeared_initial_conditions_fourier_coefficients_absz_{0to1,1to2}.png
Numbers: results/test2b/smeared_initial_conditions.json
"""
import json
import time
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test1_nbody_qmc_focused_run import library_samples
from plot_distance_error_effect_on_fourier_coefficients import (SMEARED, stars, forward_smeared_library_fit,
                                                                observe_library_samples)
from plot_fourier_coefficients_data_vs_models import as_kms, FIELDS, TERMS, LAYERS, RC

OUT_DIR = Path.home() / "mira_dr4_runs" / "test2b"
TRAJ = OUT_DIR / "traj_smeared_ics_true_potential.npy"
RESULTS = ROOT / "results" / "test2b" / "smeared_initial_conditions.json"
WINDOWS = ((0.0, 0.1), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0))
NK = fs.NK
LINES = {"a": ("smeared ICs, as is", ACCENT, "-"), "b": ("smeared ICs, observed again", INK, ":"),
         "c": ("true ICs (truth library), observed", BLUE_RAMP[5], "--")}


def integrate_smeared(times):
    if TRAJ.exists():
        return np.load(TRAJ), 0.0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ssl.ACCURACY = old.ACCURACY
    t0 = time.time()
    pot = ssl.bar_aligned_potential(old.TRUTH["bar_frac"], old.TRUTH["S"], "nbody")
    tt, traj = ssl.integrate(pot, ssl.rotate_z(stars(SMEARED), np.radians(old.TRUTH["phi_b_deg"])),
                             -abs(old.TRUTH["Omega_abs"]))
    assert np.allclose(tt, times)
    np.save(TRAJ, traj)
    return traj, time.time() - t0


def window_fit(traj, times, ok, t_lo, t_hi, sigma_dm=None):
    """Library fit using only snapshots with t_lo <= t < t_hi (the last window includes t = 1 Gyr), 1/tau weights from
    the whole orbit as usual; with sigma_dm the samples are first observed with distance errors."""
    R = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])
    tau = np.clip(np.mean(R < ssl.R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
    keep = (times >= t_lo - 1e-9) & ((times < t_hi - 1e-9) | (t_hi >= times[-1] - 1e-9))
    w = traj[:, keep].reshape(-1, 6).astype(np.float64)
    wt = np.repeat(1.0 / tau, keep.sum())
    if sigma_dm is not None:
        w, wt = observe_library_samples(w, wt, sigma_dm, old.TRUTH["phi_b_deg"])
    b, X, v, m = fs.prepare(w)
    return fs.fit(b, X, v, wt[m], ok, leverage=False)


def main():
    set_style()
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    traj_s, t_int = integrate_smeared(times)
    traj_t = np.load(old.ARR_DIR / "traj_000.npy")
    sigma_dm = float(h5py.File(SMEARED, "r").attrs["sigma_DM"])
    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(50000, np.full(50000, 1 / 50000), size=fs.N_BOOT).astype(float)
    cd, cov, n, ok = fs.data_fit(stars(SMEARED), old.TRUTH["phi_b_deg"], boot_w)
    libs = dict(a=fs.library_fit(traj_s, times, ok),
                b=forward_smeared_library_fit(traj_s, times, ok, sigma_dm, old.TRUTH["phi_b_deg"]),
                c=forward_smeared_library_fit(traj_t, times, ok, sigma_dm, old.TRUTH["phi_b_deg"]),
                d=fs.library_fit(traj_t, times, ok))
    valid = ok & np.all(np.isfinite(cd), axis=(1, 2))
    for c in libs.values():
        valid &= np.all(np.isfinite(c), axis=(1, 2))
    _, _, tau_s = library_samples(times, traj_s)
    _, _, tau_t = library_samples(times, traj_t)
    res = dict(sigma_DM=sigma_dm, integration_s=round(t_int, 1),
               tau_median=dict(smeared_ics=float(np.median(tau_s)), true_ics=float(np.median(tau_t))),
               labels={k: v[0] for k, v in LINES.items()} | {"d": "true ICs (truth library), as is"},
               chi2_per_N={})
    schemes = dict(fs.SCHEMES)
    for f, name in enumerate(fs.FIELD_NAMES):
        schemes[f"field_{name}_m0"] = (fs.N_RING_DEFAULT, {f: [0]})
        schemes[f"field_{name}_m24"] = (fs.N_RING_DEFAULT, {f: fs.NO_M0})
    for name, sc in schemes.items():
        row = {}
        for k, c in libs.items():
            chi, nd = fs.scheme_chi2(cd, c, cov, valid, sc)
            row[k] = round(chi / nd, 3)
        row["N_data"] = nd
        res["chi2_per_N"][name] = row
    proposed = fs.SCHEMES["fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"]
    res["time_windows_proposed_scheme"] = {}
    for t_lo, t_hi in WINDOWS:
        ca = window_fit(traj_s, times, ok, t_lo, t_hi)
        cb = window_fit(traj_s, times, ok, t_lo, t_hi, sigma_dm)
        cc = window_fit(traj_t, times, ok, t_lo, t_hi, sigma_dm)
        vv = valid & np.all(np.isfinite(ca), axis=(1, 2)) & np.all(np.isfinite(cb), axis=(1, 2))
        vv &= np.all(np.isfinite(cc), axis=(1, 2))
        res["time_windows_proposed_scheme"][f"{t_lo:g}-{t_hi:g} Gyr"] = {
            k: round(float(np.divide(*fs.scheme_chi2(cd, c, cov, vv, proposed))), 3) for k, c in (("a", ca), ("b", cb), ("c", cc))}

    shape = (fs.NR, fs.NL, 6, NK)
    for layer, tag, label in LAYERS:
        fig, axes = plt.subplots(len(FIELDS), len(TERMS), figsize=(17.5, 12.5), gridspec_kw=dict(hspace=0.42, wspace=0.32))
        use = valid.reshape(fs.NR, fs.NL)[:, layer]
        r6 = use & (np.arange(fs.NR) < fs.N_RING_DEFAULT)
        for r, (name, f, kind) in enumerate(FIELDS):
            e = np.sqrt(np.diagonal(cov, axis1=1, axis2=2)).reshape(shape)[:, layer, f]
            dat, err = as_kms(cd.reshape(shape)[:, layer, f], e, kind)
            for k, term in enumerate(TERMS):
                ax = axes[r, k]
                ax.axvspan(6.0, 8.0, color=INK_MUTED, alpha=0.12, lw=0)
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                notes = []
                for key, (lab, colr, ls) in LINES.items():
                    mod, _ = as_kms(libs[key].reshape(shape)[:, layer, f], kind=kind)
                    ax.plot(RC[use], mod[use, k], color=colr, ls=ls, lw=1.7, label=f"({key}) {lab}" if (r, k) == (0, 0) else None)
                    notes.append(f"({key}) {np.mean(((dat[r6, k] - mod[r6, k]) / err[r6, k])**2):.2f}")
                ax.errorbar(RC[use], dat[use, k], yerr=err[use, k], fmt="o", ms=3.2, color=INK, ecolor=INK_2, elinewidth=0.9,
                            capsize=0, zorder=4, label="smeared data (σ_DM = 0.265 mag)" if (r, k) == (0, 0) else None)
                ax.text(0.02, 0.97, "χ²/N  " + "  ".join(notes), transform=ax.transAxes, fontsize=6.6, va="top",
                        color=INK_2, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
                ax.set_xlim(0, 8)
                ax.tick_params(labelsize=7.5)
                if r == 0:
                    ax.set_title(term, fontsize=10)
                if k == 0:
                    ax.set_ylabel(f"{name}\n" + ("km/s" if kind == "mean" else "σ [km/s]"), fontsize=9.5)
                if r == len(FIELDS) - 1:
                    ax.set_xlabel("R [kpc] (observed)", fontsize=9)
        fig.legend(loc="upper center", ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.965))
        fig.suptitle(f"Orbits from distance-smeared stars, true potential and pattern speed: Fourier coefficients per 0.3 kpc "
                     f"ring, {label}, bar frame at the true angle\nχ²/N per panel against the smeared data (diagonal "
                     f"errors, R < 6 kpc); σ terms: m = 0 is σ₀, others c_m / 2σ₀", fontsize=10.5, y=0.997)
        savefig(fig, f"test2b_smeared_initial_conditions_fourier_coefficients_{tag}")
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

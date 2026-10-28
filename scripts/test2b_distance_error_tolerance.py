"""Test 2b: how large a distance error can the current method take?

The current method is applied unchanged to data with distance errors: orbits are integrated from the observed (smeared)
stars, and the per-star Fourier coefficients of the observed stars are compared with the library's (proposed scheme,
joint bootstrap covariance, tempered logL = -chi2 / (2 sqrt(2 N_data)), N_data at the truth of each error level).
Error levels: fractional distance errors 0, 2.5, 5, 7.5, 10% and 12.22% (sigma_DM = FRAC / (ln 10 / 5) = 0.2654, the
DR3 Mira median), Gaussian in DM, with the same normal deviates for every level (those of
make_mock_distance_error.py, seed 2029), so the 12.22% level is that mock exactly.
Trial potentials: the truth and +-offsets in one parameter at a time (|Omega_b| +-1.2 km/s/kpc, phi_b +-5 deg,
S x exp(+-0.18), bar_frac x exp(+-0.30); about 2 sigma of the perfect-data conditional widths).
Per level and parameter, from L-, L0, L+ (in |Omega_b|, phi_b, ln S, ln bar_frac):
  width (conditional 1 sigma)   s^2 = delta^2 / (2 L0 - L+ - L-)
  shift of the maximum          b = s^2 (L+ - L-) / (2 delta)
Trajectories: ~/mira_dr4_runs/test2b/tolerance/traj_<level>_<point>.npy (re-used if present).
Figure: figures/test2b_distance_error_tolerance_bias_and_width.png
Numbers: results/test2b/distance_error_tolerance.json
"""
import json
import time
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import observe
from make_mock_distance_error import SEED as SMEAR_SEED
from test2b_smeared_initial_conditions import TRAJ as TRAJ_SMEARED

with h5py.File(MOCK_DIR / "nbody_mock_distance_error_sigmaDM0.265_R8_N50k.hdf5", "r") as _f:
    FRAC_MIRA = float(_f.attrs["frac_distance_error"])     # 0.1222: DR3 Miras, sigma_DM < 0.4 (median)
FRACS = (0.0, 0.025, 0.05, 0.075, 0.10, FRAC_MIRA)
TRUTH = dict(old.TRUTH)
OFFSETS = dict(Omega_abs=1.2, phi_b_deg=5.0, S=0.18, bar_frac=0.30)        # S, bar_frac: in ln
LABELS = {"Omega_abs": "|Ω_b|", "phi_b_deg": "φ_b", "S": "ln S", "bar_frac": "ln bar_frac"}
UNITS = {"Omega_abs": "km/s/kpc", "phi_b_deg": "deg", "S": "", "bar_frac": ""}
TRAJ_DIR = Path.home() / "mira_dr4_runs" / "test2b" / "tolerance"
RESULTS = ROOT / "results" / "test2b" / "distance_error_tolerance.json"
PROPOSED = fs.SCHEMES["fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"]


def points():
    pts = {"truth": TRUTH}
    for k, d in OFFSETS.items():
        for sgn, tag in ((-1, "minus"), (1, "plus")):
            p = dict(TRUTH)
            p[k] = TRUTH[k] * np.exp(sgn * d) if k in ("S", "bar_frac") else TRUTH[k] + sgn * d
            pts[f"{k}_{tag}"] = p
    return pts


def observed_stars(frac):
    """Galactocentric (MW frame) phase space of the mock observed with fractional distance error frac."""
    with h5py.File(old.MOCK, "r") as f:
        m = {k: f[k][:].astype(np.float64) for k in ("l", "b", "d", "pml", "pmb", "vlos")}
    eps = np.random.default_rng(SMEAR_SEED).standard_normal(len(m["d"]))      # the deviates of the 12.2% mock
    sigma_dm = frac / (np.log(10) / 5)
    d_obs = m["d"] * 10.0 ** (0.2 * sigma_dm * eps)
    return np.column_stack(observe.to_galactocentric(m["l"], m["b"], d_obs, m["pml"], m["pmb"], m["vlos"]))


def trajectory(frac, name, p, w_mw, log):
    path = TRAJ_DIR / f"traj_{frac:.3f}_{name}.npy"
    if frac == 0.0 and name == "truth":
        return np.load(old.ARR_DIR / "traj_000.npy")
    if frac == FRAC_MIRA and name == "truth":
        return np.load(TRAJ_SMEARED)
    if path.exists():
        return np.load(path)
    ssl.ACCURACY = old.ACCURACY
    t0 = time.time()
    pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], "nbody")
    _, traj = ssl.integrate(pot, ssl.rotate_z(w_mw, np.radians(p["phi_b_deg"])), -abs(p["Omega_abs"]))
    np.save(path, traj)
    log(f"  integrated {frac:.3f} {name} in {time.time() - t0:.0f} s")
    return traj


def main():
    set_style()
    TRAJ_DIR.mkdir(parents=True, exist_ok=True)
    log = lambda s: print(s, flush=True)
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(50000, np.full(50000, 1 / 50000), size=fs.N_BOOT).astype(float)
    pts = points()
    res = dict(fracs=list(FRACS), offsets=OFFSETS, levels={})
    for frac in FRACS:
        w_mw = observed_stars(frac)
        data = {}
        chi = {}
        for name, p in pts.items():
            phi = round(p["phi_b_deg"], 6)
            if phi not in data:
                data[phi] = fs.data_fit(w_mw, phi, boot_w)
            cd, cov, n, ok = data[phi]
            cl = fs.library_fit(trajectory(frac, name, p, w_mw, log), times, ok)
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            chi[name] = fs.scheme_chi2(cd, cl, cov, valid, PROPOSED)
        nd = chi["truth"][1]
        L = {k: -0.5 * c / np.sqrt(2 * nd) for k, (c, _) in chi.items()}
        lev = dict(sigma_DM=frac / (np.log(10) / 5), N_data=nd, truth_chi2_per_N=chi["truth"][0] / nd,
                   logL={k: round(v - L["truth"], 3) for k, v in L.items()}, params={})
        for k, d in OFFSETS.items():
            lm, l0, lp = L[f"{k}_minus"], L["truth"], L[f"{k}_plus"]
            curv = 2 * l0 - lp - lm
            s2 = d**2 / curv if curv > 0 else np.nan
            lev["params"][k] = dict(width=float(np.sqrt(s2)) if curv > 0 else None,
                                    shift=float(s2 * (lp - lm) / (2 * d)) if curv > 0 else None,
                                    curvature_ok=bool(curv > 0))
        res["levels"][f"{frac:.3f}"] = lev
        log(f"level {100 * frac:.1f}%: truth chi2/N {lev['truth_chi2_per_N']:.2f}; " + "; ".join(
            f"{LABELS[k]} shift {v['shift']:+.3f} width {v['width']:.3f}" if v["curvature_ok"] else f"{LABELS[k]} not curved"
            for k, v in lev["params"].items()))
        RESULTS.parent.mkdir(parents=True, exist_ok=True)
        with open(RESULTS, "w") as fh:
            json.dump(res, fh, indent=1)

    # shift in units of the perfect-data width
    w0 = {k: res["levels"]["0.000"]["params"][k]["width"] for k in OFFSETS}
    for lev in res["levels"].values():
        for k, v in lev["params"].items():
            v["shift_in_perfect_sigma"] = v["shift"] / w0[k] if v["shift"] is not None else None
    with open(RESULTS, "w") as fh:
        json.dump(res, fh, indent=1)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), gridspec_kw=dict(wspace=0.3))
    x = 100 * np.array(FRACS)
    colours = dict(Omega_abs=BLUE_RAMP[6], phi_b_deg=BLUE_RAMP[3], S=ACCENT, bar_frac=INK)
    marks = dict(Omega_abs="o", phi_b_deg="s", S="^", bar_frac="D")
    for k in OFFSETS:
        sh = [res["levels"][f"{f:.3f}"]["params"][k]["shift_in_perfect_sigma"] for f in FRACS]
        wd = [res["levels"][f"{f:.3f}"]["params"][k]["width"] for f in FRACS]
        axes[0].plot(x, [np.nan if s is None else s for s in sh], marker=marks[k], color=colours[k], lw=1.5, label=LABELS[k])
        axes[1].plot(x, [np.nan if v is None else v / w0[k] for v in wd], marker=marks[k], color=colours[k], lw=1.5,
                     label=LABELS[k])
    axes[0].axhspan(-0.5, 0.5, color=INK_MUTED, alpha=0.15, lw=0)
    axes[0].axhline(0, color=INK_MUTED, lw=0.8)
    axes[0].set_ylabel("shift of the best fit / perfect-data 1σ")
    axes[0].set_title("bias (conditional, other parameters at the truth)", fontsize=10)
    axes[1].axhline(1, color=INK_MUTED, lw=0.8)
    axes[1].set_ylabel("conditional 1σ width / perfect-data width")
    axes[1].set_title("width", fontsize=10)
    axes[2].plot(x, [res["levels"][f"{f:.3f}"]["truth_chi2_per_N"] for f in FRACS], marker="o", color=ACCENT, lw=1.5)
    axes[2].set_ylabel("χ²/N of the true potential")
    axes[2].set_title("fit quality at the truth", fontsize=10)
    for ax in axes:
        ax.set_xlabel("fractional distance error [%]")
        ax.axvline(100 * FRAC_MIRA, color=INK_2, lw=0.8, ls=":")
    axes[0].legend(fontsize=8.5, frameon=False)
    axes[0].text(12.0, axes[0].get_ylim()[1] * 0.9, "DR3 Miras\n(σ_DM < 0.4)", fontsize=8, color=INK_2, ha="right", va="top")
    fig.suptitle("Current method on data with distance errors (orbits from the observed stars, no correction); proposed "
                 "Fourier scheme, tempered logL; shifts and widths from ±offsets about the truth", fontsize=10.5, y=1.0)
    savefig(fig, "test2b_distance_error_tolerance_bias_and_width")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

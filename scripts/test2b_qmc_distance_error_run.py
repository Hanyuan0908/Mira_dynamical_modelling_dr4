"""Test 2b, brave test: the current method on data with a 7.5% distance error, over the focused Sobol grid.

Data: mock/nbody_mock_distance_error_sigmaDM0.163_R8_N50k.hdf5 (sigma_DM = 0.163 mag, 7.5%; the deviates of the 12.2%
mock, scaled). At every point of the focused global search (results/test1_nbody_qmc_focused/design.csv, 257 points)
the observed stars are rotated by the trial bar angle and integrated in the trial potential (accuracy 1e-5, 1 Gyr,
41 snapshots). Library: samples t >= 0.25 Gyr with 1/tau weights (tau: time fraction of the orbit at R < 8 kpc),
  observed  observed again with the same distance error (the user's choice; seed OBS_SEED for every point),
  asis      as it is (no extra smearing; costs nothing extra).
Per-star Fourier fit (test1_nbody_qmc_focused_fourier_stars.py) of the data at the trial angle (300 bootstraps) and of
the library; chi2 for its schemes; tempered logLs = -chi2 / (2 sqrt(2 N_data)), N_data at the truth.
Files: ~/mira_dr4_runs/test2b_qmc_7p5/traj_NNN.npy, fourier_NNN.npz; results/test2b_qmc_7p5/evals.csv. Resumable:
points with fourier_NNN.npz are re-used.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test1_nbody_qmc_focused_run import DESIGN_CSV
from plot_distance_error_effect_on_fourier_coefficients import stars, forward_smeared_library_fit

MOCK = MOCK_DIR / "nbody_mock_distance_error_sigmaDM0.163_R8_N50k.hdf5"
OUT_DIR = Path.home() / "mira_dr4_runs" / "test2b_qmc_7p5"
RES_DIR = ROOT / "results" / "test2b_qmc_7p5"
EVALS = RES_DIR / "evals.csv"
OBS_SEED = 11
VARIANTS = ("observed", "asis")


def write_evals(rows):
    df = pd.DataFrame(rows)
    t = df.kind == "truth"
    if t.any():
        for name in fs.SCHEMES:
            for v in VARIANTS:
                nd = int(df.loc[t, f"ndata_{name}_{v}"].iloc[0])
                df[f"logLs_{name}_{v}"] = -0.5 * df[f"chi2_{name}_{v}"] / np.sqrt(2 * nd)
    df.to_csv(EVALS, index=False)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)
    design = pd.read_csv(DESIGN_CSV).sort_values("idx").reset_index(drop=True)
    w_mw = stars(MOCK)
    sigma_dm = float(h5py.File(MOCK, "r").attrs["sigma_DM"])
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    boot_w = np.random.default_rng(fs.SEED).multinomial(len(w_mw), np.full(len(w_mw), 1 / len(w_mw)),
                                                        size=fs.N_BOOT).astype(float)
    ssl.ACCURACY = old.ACCURACY
    rows, t_start = [], time.time()
    for _, p in design.iterrows():
        i = int(p.idx)
        saved = OUT_DIR / f"fourier_{i:03d}.npz"
        row = {k: p[k] for k in ("idx", "kind", "batch", "Omega_abs", "phi_b_deg", "S", "bar_frac")}
        if saved.exists():
            s = np.load(saved)
            row.update({k: float(s[k]) for k in s.files if k.startswith(("chi2_", "ndata_"))})
            rows.append(row)
            continue
        t0 = time.time()
        path = OUT_DIR / f"traj_{i:03d}.npy"
        if path.exists():
            traj = np.load(path)
        else:
            pot = ssl.bar_aligned_potential(p.bar_frac, p.S, "nbody")
            _, traj = ssl.integrate(pot, ssl.rotate_z(w_mw, np.radians(p.phi_b_deg)), -abs(p.Omega_abs))
            np.save(path, traj)
        t_int = time.time() - t0
        cd, cov, n, ok = fs.data_fit(w_mw, float(p.phi_b_deg), boot_w)
        libs = dict(asis=fs.library_fit(traj, times, ok),
                    observed=forward_smeared_library_fit(traj, times, ok, sigma_dm, float(p.phi_b_deg), seed=OBS_SEED))
        out = {}
        for v, cl in libs.items():
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            for name, sc in fs.SCHEMES.items():
                out[f"chi2_{name}_{v}"], out[f"ndata_{name}_{v}"] = fs.scheme_chi2(cd, cl, cov, valid, sc)
        np.savez(saved, coef_data=cd, cov_data_joint=cov, n_stars=n, coef_lib_asis=libs["asis"],
                 coef_lib_observed=libs["observed"], **out)
        row.update(out)
        rows.append(row)
        write_evals(rows)
        key = "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"
        print(f"[{len(rows)}/{len(design)}] {p.kind:5s} Omega={p.Omega_abs:.2f} phi={p.phi_b_deg:.1f} S={p.S:.3f} "
              f"f={p.bar_frac:.3f}  chi2/N observed {out[f'chi2_{key}_observed'] / out[f'ndata_{key}_observed']:.2f}  "
              f"asis {out[f'chi2_{key}_asis'] / out[f'ndata_{key}_asis']:.2f}  integrate {t_int:.0f} s, total "
              f"{time.time() - t0:.0f} s (elapsed {(time.time() - t_start) / 60:.0f} min)", flush=True)
    write_evals(rows)


if __name__ == "__main__":
    main()

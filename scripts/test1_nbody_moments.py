"""Test 1b: moment-map likelihood with bootstrap errors for the N-body global search.

For every design point the observed stars are rotated by the trial bar angle and their means and
dispersions of (v_R, v_phi, v_z sign z) are computed in the folded 3D bar-frame cells of
steady_state_likelihood (R x phi mod pi x |z|). Errors on the observed moments come from N_BOOT
bootstrap resamplings of the stars (they include non-Gaussian tails, unlike sigma/sqrt(2N)). The
library moments saved by test1_nbody_qmc.py are compared with
    chi2 = sum over cells with N_obs >= N_MIN of (m_obs - m_lib)^2 / err_boot^2.
Columns added to results/test1_nbody/qmc_evals.csv: chi2b_<variant> for variant in means, means_R1,
meansdisp, meansdisp_R1, and the Schwarzschild-style tempered log-likelihood used from now on,
logLs_<variant> = -chi2b_<variant> / (2 sqrt(2 N_data)), N_data = (3 or 6 moments) x n_cells of that
variant at the truth, which treats the sqrt(2 N_data) scatter of chi2 itself as the 1-sigma threshold.
"""
import numpy as np
import pandas as pd
import h5py

import steady_state_likelihood as ssl
from test1_nbody_qmc import OUT_CSV, ARR_DIR, MOCK

N_BOOT = 200
N_MIN = 20
SEED = 0


def boot_errors(stars, rng):
    """Observed moments and their bootstrap errors in the folded cells."""
    N, m1, s = ssl.cell_moments(stars)
    bm, bs = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(stars), len(stars))
        _, a, b = ssl.cell_moments(stars[i])
        bm.append(a)
        bs.append(b)
    return N, m1, s, np.nanstd(bm, axis=0), np.nanstd(bs, axis=0)


def chi2(obs, lib_m1, lib_s, rmin, disp):
    N, m1, s, em, es = obs
    Rc = np.repeat(0.5 * (ssl.CELL_R[1:] + ssl.CELL_R[:-1]), (len(ssl.CELL_PHI) - 1) * (len(ssl.CELL_Z) - 1))
    good = (N >= N_MIN) & (Rc >= rmin) & np.all(np.isfinite(lib_m1), 1) & np.all(em > 0, 1) & np.all(es > 0, 1)
    c = np.sum(((m1 - lib_m1)[good] / em[good])**2)
    if disp:
        c += np.sum(((s - lib_s)[good] / es[good])**2)
    return float(c), int(good.sum())


def main():
    df = pd.read_csv(OUT_CSV).sort_values("idx").reset_index(drop=True)
    with h5py.File(MOCK, "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    cache = {}
    VARIANTS = (("chi2b_means", 0.0, False), ("chi2b_means_R1", 1.0, False),
                ("chi2b_meansdisp", 0.0, True), ("chi2b_meansdisp_R1", 1.0, True))
    out = {k: [] for k in [v[0] for v in VARIANTS] + [v[0] + "_ncells" for v in VARIANTS]}
    for _, row in df.iterrows():
        phi = round(float(row.phi_b_deg), 6)
        if phi not in cache:      # same bootstrap seed for every point: common random numbers
            cache[phi] = boot_errors(ssl.rotate_z(stars_mw, np.radians(phi)), np.random.default_rng(SEED))
        a = np.load(ARR_DIR / f"{int(row.idx):03d}.npz")
        for name, rmin, disp in VARIANTS:
            c, n = chi2(cache[phi], a["lib_m1"], a["lib_s"], rmin, disp)
            out[name].append(c)
            out[name + "_ncells"].append(n)
    for k, v in out.items():
        df[k] = v
    t = df.kind == "truth"
    for name, _, disp in VARIANTS:
        n_data = (6 if disp else 3) * int(df.loc[t, name + "_ncells"].iloc[0])
        df[name.replace("chi2b_", "logLs_")] = -0.5 * df[name] / np.sqrt(2 * n_data)
        print(f"tempering {name:20s}: N_data = {n_data:5d}, sqrt(2 N_data) = {np.sqrt(2 * n_data):5.1f}")
    df.to_csv(OUT_CSV, index=False)
    for k in ("logLs_means", "logLs_means_R1", "logLs_meansdisp", "logLs_meansdisp_R1"):
        dl = df[k] - df.loc[t, k].iloc[0]
        b = dl.idxmax()
        print(f"{k:20s} "
              f"truth rank {int((dl > 0).sum()) + 1}/65; best: Omega {df.Omega_abs[b]:.2f} phi {df.phi_b_deg[b]:.1f} "
              f"S {df.S[b]:.3f} f {df.bar_frac[b]:.3f} (dlogL {dl[b]:+.2f}); second best dlogL "
              f"{dl.drop(index=b).max():+.2f}")


if __name__ == "__main__":
    main()

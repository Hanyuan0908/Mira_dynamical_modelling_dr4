"""Rescore the 7.5% distance-error global search with the co-added Fourier error (user's recipe for a real observation).

Error at the true bar angle (27 deg): C = C_MC + C_boot, with C_MC the covariance over N_MC re-draws of every observed
star's distance (test2b_coadded_errors_truth_chi2.py) and C_boot the 300-bootstrap covariance of the unsmeared data
(Test 1). At another trial angle the data are the same stars rotated by d = phi_b - 27 deg: rings, layers, v_R and
v_phi are unchanged and the (cos, sin) pair of mode m rotates by m d, so the coefficients transform exactly as c' = T c
and the error as T C T^T (checked against the data fits saved at each angle). chi2 for the proposed scheme (its kept
terms are whole cos/sin pairs, so the transform is exact), Hartlap with 300; tempered logLs with N_data at the truth.
Writes results/test2b_qmc_7p5/evals_coadded.csv for the points finished so far.
"""
import json

import numpy as np
import pandas as pd
import h5py

from common import ROOT
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import observe
from test1_nbody_qmc_focused_run import OUT_DIR as FOCUSED_DIR, DESIGN_CSV
from test2b_qmc_distance_error_run import MOCK, OUT_DIR, RES_DIR

N_MC = 300
PHI0 = 27.0
SCHEME = fs.SCHEMES["fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"]
C_MC_FILE = OUT_DIR / "C_MC_27deg.npy"


def rotation(d_deg):
    """30 x 30 transform of the (field, term) coefficients when the bar-frame azimuth shifts by d."""
    T = np.eye(fs.NK)
    for m, (c, s) in ((2, (1, 2)), (4, (3, 4))):
        a = np.radians(m * d_deg)
        T[np.ix_([c, s], [c, s])] = [[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]
    return np.kron(np.eye(6), T)


def c_mc(ok):
    if C_MC_FILE.exists():
        return np.load(C_MC_FILE)
    sigma_dm = float(h5py.File(MOCK, "r").attrs["sigma_DM"])
    with h5py.File(MOCK, "r") as f:
        m = {k: f[k][:].astype(np.float64) for k in ("l", "b", "pml", "pmb", "vlos", "d")}
    rng = np.random.default_rng(21)
    mc = []
    for _ in range(N_MC):
        d = observe.scatter_distance_modulus(m["d"], sigma_dm, rng)
        obs = np.column_stack(observe.to_galactocentric(m["l"], m["b"], d, m["pml"], m["pmb"], m["vlos"]))
        b, X, v, _ = fs.prepare(ssl.rotate_z(obs, np.radians(PHI0)))
        mc.append(fs.fit(b, X, v, np.ones(len(b)), ok).reshape(fs.NB, -1))
    dev = np.array(mc) - np.mean(mc, 0)
    C = np.einsum("mbi,mbj->bij", dev, dev) / (N_MC - 1)
    np.save(C_MC_FILE, C)
    return C


def main():
    design = pd.read_csv(DESIGN_CSV).set_index("idx")
    s0 = np.load(OUT_DIR / "fourier_000.npz")
    ok = s0["n_stars"] >= fs.N_MIN_BLOCK
    p0 = np.load(FOCUSED_DIR / "fourier_stars_000.npz")
    ok &= p0["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK
    C0 = c_mc(s0["n_stars"] >= fs.N_MIN_BLOCK) + p0["cov_data_joint"].reshape(fs.NB, 30, 30)
    cd0 = s0["coef_data"].reshape(fs.NB, -1)
    rows, check = [], []
    for f in sorted(OUT_DIR.glob("fourier_[0-9][0-9][0-9].npz")):
        i = int(f.stem.split("_")[1])
        s = np.load(f)
        p = design.loc[i]
        T = rotation(float(p.phi_b_deg) - PHI0)
        cd = s["coef_data"].reshape(fs.NB, -1)
        check.append(float(np.nanmax(np.abs(cd[ok] - cd0[ok] @ T.T))))
        C = np.einsum("ij,bjk,lk->bil", T, C0, T)
        row = dict(idx=i, kind=p.kind, Omega_abs=p.Omega_abs, phi_b_deg=p.phi_b_deg, S=p.S, bar_frac=p.bar_frac)
        for lib in ("observed", "asis"):
            cl = s[f"coef_lib_{lib}"]
            valid = ok & np.all(np.isfinite(cl), axis=(1, 2))
            row[f"chi2_{lib}"], row[f"ndata_{lib}"] = fs.scheme_chi2(cd.reshape(fs.NB, 6, fs.NK), cl, C, valid, SCHEME)
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("idx")
    t = df.kind == "truth"
    for lib in ("observed", "asis"):
        nd = int(df.loc[t, f"ndata_{lib}"].iloc[0])
        df[f"chi2N_{lib}"] = df[f"chi2_{lib}"] / df[f"ndata_{lib}"]
        df[f"logLs_{lib}"] = -0.5 * df[f"chi2_{lib}"] / np.sqrt(2 * nd)
        df[f"dlogLs_{lib}"] = df[f"logLs_{lib}"] - df.loc[t, f"logLs_{lib}"].iloc[0]
    df.to_csv(RES_DIR / "evals_coadded.csv", index=False)
    print(f"rotation check: max |c_data(phi) - T c_data(27)| over {len(check)} points = {max(check):.2e} km/s or (km/s)^2")
    print(df[["idx", "kind", "Omega_abs", "phi_b_deg", "S", "bar_frac", "chi2N_observed", "dlogLs_observed",
              "chi2N_asis", "dlogLs_asis"]].sort_values("dlogLs_observed", ascending=False).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

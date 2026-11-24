"""Can the library noise floor be lowered, and does the pattern-speed signal then stand out?

Background (PROGRESS 2026-11-24): in the Milky Way family, a +-20-40% change of |Omega_b| changes the predicted
kinematics by D ~ 4-7 (tempered), only ~2-5 above the library noise floor D ~ 2.4 that separates two libraries of
near-identical potentials (+-0.8%). The floor comes from Monte Carlo sampling: each orbit is represented by 16 phases
(0.25-1 Gyr every 50 Myr) and each phase by K = 10 distance-error draws.

Here every |Omega_b| is integrated ONCE for 4 Gyr at 25 Myr cadence, and the libraries of several variants are built
from the same trajectories (columns always include t = 0):
  base      1 Gyr, 50 Myr, K = 10   (15 phases: the pipeline; the control)
  dt25      1 Gyr, 25 Myr, K = 10   (30 phases, same time span)
  T4_dt100  4 Gyr, 100 Myr, K = 10  (38 phases spread over 4x the time)
  T4_dt50   4 Gyr, 50 Myr, K = 10   (75 phases)
  K30       1 Gyr, 50 Myr, K = 30   (15 phases, 3x the distance draws)
For each variant: the noise floor (D at +-0.8%), the change of the prediction D at -30 ... +40%, the signal D - floor,
the data's Delta chi2, and chi2/N at the reference (Omega_ref = 32.5; for the steady-state mock that is the truth).
Star sets: OSARG and the steady-state OSARG-like mock (truth 32.5, phi_b 25, S = bar_frac = 1), both in Hunter+24 at
S = bar_frac = 1, phi_b 25 deg.

Usage: python library_noise_floor_test.py [--cells osarg:mw h24mock:mw]
Output: results/library_noise_floor_<stars>_<family>.npz and a printed table
"""
import argparse
import time

import numpy as np
import agama

from common import RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from potential_family import TIME_UNIT_GYR, ORBIT_KW
from rescore_empirical_selection import library_fit_empirical
from cross_response_stars_vs_potential import STARS, FAMILY, FRACS, family_potential, B

T_LONG, N_LONG = 4.0, 161                     # 4 Gyr at 25 Myr
VARIANTS = {  # name: (max time [Gyr], column step in units of 25 Myr, K)
    "base": (1.0, 2, 10), "dt25": (1.0, 1, 10), "T4_dt100": (4.0, 4, 10), "T4_dt50": (4.0, 2, 10), "K30": (1.0, 2, 30),
}


def load_stars(stars):
    import h5py
    from common import ROOT, MOCK_DIR
    from rescore_empirical_selection import data_lbdR
    st = STARS[stars]
    mock = MOCK_DIR / st["mock"]
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
        sigma_dm = float(f.attrs["sigma_DM"])
    dat = data_lbdR(mock)
    if st["keep"]:
        keep = np.load(ROOT / st["keep"])
        w, dat = w[keep], tuple(np.asarray(a)[keep] for a in dat)
    return w, dat, sigma_dm, st["phi"]


def run(stars, family):
    w, dat, sigma_dm, phi = load_stars(stars)
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    T = fs.rotation(phi - BAR_ANGLE_DEG)
    cd = (cd0.reshape(fs.NB, -1) @ T.T).reshape(fs.NB, 6, fs.NK)
    cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
    pot = family_potential(family)
    ic = ssl.rotate_z(w, np.radians(phi))
    om_ref = FAMILY[family]
    omegas = om_ref * (1 + FRACS)
    libs = {v: [] for v in VARIANTS}
    chi2 = {v: [] for v in VARIANTS}
    nd = None
    print(f"{stars} in {family}: {len(w):,} stars, phi_b {phi}, Omega_ref {om_ref}", flush=True)
    for om in omegas:
        t0 = time.time()
        orbs = agama.orbit(potential=pot, ic=ic, time=T_LONG / TIME_UNIT_GYR, trajsize=N_LONG, Omega=-abs(om),
                           accuracy=old.ACCURACY, **ORBIT_KW)
        # the pipeline's float32 time grid (ssl.snapshot_times): t = 0.25 Gyr rounds to 0.24999999 and is not kept,
        # so 'base' has the pipeline's 15 phases (0.30-1.0 Gyr)
        times = np.asarray(orbs[0, 0]).astype(np.float32) * TIME_UNIT_GYR
        traj = np.stack(orbs[:, 1]).astype(np.float32)
        del orbs
        t_int = time.time() - t0
        line = []
        for v, (tmax, step, K) in VARIANTS.items():
            cols = np.arange(0, N_LONG, step)
            cols = cols[times[cols] <= tmax + 1e-6]
            cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, cols]), phi, sigma_dm, dat, ok, n_draw=K,
                                          times=times[cols])
            cl = cl.reshape(fs.NB, 6, fs.NK)
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            x2, nd = fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[B])
            libs[v].append(cl)
            chi2[v].append(x2)
            line.append(f"{v} {x2 / nd:.3f}")
        print(f"   Omega {om:6.2f}: chi2/N " + ", ".join(line) + f"  (integrate {t_int:.0f} s, total {time.time() - t0:.0f} s)",
              flush=True)
        del traj
    iref = int(np.where(FRACS == 0)[0][0])
    temp = 2 * np.sqrt(2 * nd)
    out = dict(fracs=FRACS, omegas=omegas, ndata=nd, variants=list(VARIANTS))
    print(f"   {stars} in {family}: tempered by {temp:.1f};  frac " + " ".join(f"{x:+6.3f}" for x in FRACS))
    for v in VARIANTS:
        L = np.array(libs[v])
        D = []
        for cl in L:
            valid = ok & np.all(np.isfinite(cl), axis=(1, 2)) & np.all(np.isfinite(L[iref]), axis=(1, 2))
            D.append(fs.scheme_chi2(L[iref], cl, cov, valid, fs.SCHEMES[B])[0])
        D, c = np.array(D), np.array(chi2[v])
        out[f"D_{v}"], out[f"chi2_{v}"], out[f"coef_{v}"] = D, c, L
        floor = 0.5 * (D[np.argmin(np.abs(FRACS + 0.008))] + D[np.argmin(np.abs(FRACS - 0.008))]) / temp
        print(f"   {v:9s} floor {floor:5.2f} | D " + " ".join(f"{x / temp:6.2f}" for x in D)
              + f" | data " + " ".join(f"{(x - c[iref]) / temp:6.2f}" for x in c) + f" | chi2/N at ref {c[iref] / nd:.3f}",
              flush=True)
    np.savez(RESULTS_DIR / f"library_noise_floor_{stars}_{family}.npz", **out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", nargs="+", default=["h24mock:mw", "osarg:mw"])
    args = ap.parse_args()
    for c in args.cells:
        run(*c.split(":"))


if __name__ == "__main__":
    main()

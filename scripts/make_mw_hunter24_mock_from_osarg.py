"""An "OSARG sampled from a Milky Way simulation" mock: the OSARG stars phase-mixed in a known Hunter+24 bar.

The question it answers (user, 2026-11-23): at the same N, coverage, distance error and trial family, would a sample
drawn from a steady state of the Hunter+24 family constrain the pattern speed as well as the 25 simulated galaxies
do? Real OSARG data give 43% (16-84% width of |Omega_b|); the t = 7 mock with OSARG's coverage and N gives 8.8%.

Construction, keeping every ingredient of the real run except the data's origin:
  1. Initial conditions: the 26,630 OSARG stars exactly as observed (mock/osarg_zhang2024_observed_R8.hdf5), so the
     population has OSARG's orbits and kinematic temperature. No extra cloning or smearing.
  2. Truth: the Hunter+24 family at S = bar_frac = 1 (Hunter+24's own bar), |Omega_b| = OMEGA_TRUE,
     phi_b = PHI_TRUE; inside the OSARG design box (|Omega_b| 15-45, phi_b 20-30, S and bar_frac 0.5-1.5), away
     from its edges. The stars are rotated into the bar frame and integrated for T_MIX Gyr in the rotating frame
     (accuracy as the searches), then sampled every DT Myr from T_START on: N_SNAP_KEEP phases per orbit, a
     phase-mixed steady state of test particles in that potential.
  3. Observation: back to the Milky Way frame, (l, b, d, pm, v_los) computed exactly; the distance modulus is
     scattered by OSARG's sigma_DM (0.2084), the same scalar the pipeline assumes; velocities are recomputed at the
     observed distance. Angles, proper motions and v_los are exact, as in every mock so far.
  4. Selection: each instance is kept with probability proportional to rho_OSARG / rho_mock in OBSERVED (x, y, z)
     (32 nearest neighbours in each set), restricted to observed R < 8 kpc, normalised to OSARG's star count.
     So the mock has OSARG's observed spatial distribution and N.

Usage: python make_mw_hunter24_mock_from_osarg.py
Output: mock/mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5 (observed values under the usual names, truth as
<name>_true, orbit index, time), and a printed comparison with the OSARG data.
"""
import time

import numpy as np
import h5py
import agama
from scipy.spatial import cKDTree

from common import MOCK_DIR, rel
import observe
import steady_state_likelihood as ssl
from potential_family import TIME_UNIT_GYR, ORBIT_KW
import test1_nbody_qmc as old

OSARG = MOCK_DIR / "osarg_zhang2024_observed_R8.hdf5"
OUT = MOCK_DIR / "mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5"
OMEGA_TRUE, PHI_TRUE, S_TRUE, F_TRUE = 32.5, 25.0, 1.0, 1.0
T_MIX, T_START, DT = 6.0, 3.0, 0.05            # Gyr
K_NN = 32
SEED = 2311


def vphi_clockwise(x, y, vx, vy):
    return -(x * vy - y * vx) / np.hypot(x, y)


def summary(tag, x, y, z, vx, vy, vz):
    R = np.hypot(x, y)
    vR = (x * vx + y * vy) / R
    vp = vphi_clockwise(x, y, vx, vy)
    out = [f"{tag:6s} N {len(x):6d}; R<6 {np.sum(R < 6):6d}; x>0 {np.sum(x > 0):6d}; |z|>1 {np.sum(np.abs(z) > 1):5d}"]
    for lo in (0, 2, 4, 6):
        m = (R >= lo) & (R < lo + 2) & (np.abs(z) < 1)
        out.append(f"   R {lo}-{lo + 2}, |z|<1: n {m.sum():5d}, <v_phi> {vp[m].mean():6.1f}, sig_R {vR[m].std():5.1f}, "
                   f"sig_phi {vp[m].std():5.1f}, sig_z {vz[m].std():5.1f}")
    return "\n".join(out)


def main():
    rng = np.random.default_rng(SEED)
    with h5py.File(OSARG, "r") as f:
        ob = {k: f[k][:] for k in ("x", "y", "z", "vx", "vy", "vz")}
        sigma_dm = float(f.attrs["sigma_DM"])
    w_mw = np.column_stack([ob[k] for k in ("x", "y", "z", "vx", "vy", "vz")]).astype(np.float64)
    n_osarg = len(w_mw)

    # 1-2. phase-mix the OSARG stars in the truth potential (bar frame, rotating at -|Omega|, as the searches)
    pot = ssl.bar_aligned_potential(F_TRUE, S_TRUE, "mw")
    ic = ssl.rotate_z(w_mw, np.radians(PHI_TRUE))
    n_snap = int(round(T_MIX / DT)) + 1
    t0 = time.time()
    orbs = agama.orbit(potential=pot, ic=ic, time=T_MIX / TIME_UNIT_GYR, trajsize=n_snap, Omega=-OMEGA_TRUE,
                       accuracy=old.ACCURACY, **ORBIT_KW)
    times = orbs[0, 0] * TIME_UNIT_GYR
    traj = np.stack(orbs[:, 1])
    keep_t = times >= T_START - 1e-6
    print(f"integrated {n_osarg:,} orbits for {T_MIX} Gyr in {time.time() - t0:.0f} s; {keep_t.sum()} phases per orbit "
          f"from {times[keep_t][0]:.2f} to {times[keep_t][-1]:.2f} Gyr")
    # Jacobi integral E - Omega L_z in the rotating frame, as an integration check
    x0, xe = traj[:, 0], traj[:, -1]
    Om = -OMEGA_TRUE
    EJ = lambda w: pot.potential(w[:, :3]) + 0.5 * np.sum(w[:, 3:] ** 2, axis=1) - Om * (w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3])
    dE = np.abs(EJ(xe) / EJ(x0) - 1)
    print(f"Jacobi integral drift over {T_MIX} Gyr: median {np.median(dE):.1e}, 99% {np.percentile(dE, 99):.1e}")

    inst = traj[:, keep_t].reshape(-1, 6)
    orbit_id = np.repeat(np.arange(n_osarg), keep_t.sum())
    t_inst = np.tile(times[keep_t], n_osarg)
    w_true = ssl.rotate_z(inst, -np.radians(PHI_TRUE))         # back to the Milky Way frame

    # 3. observe
    g = observe.to_galactic(*w_true.T)
    d_obs = observe.scatter_distance_modulus(g["d"], sigma_dm, rng)
    obs = np.column_stack(observe.to_galactocentric(g["l"], g["b"], d_obs, g["pml"], g["pmb"], g["vlos"]))
    R_obs = np.hypot(obs[:, 0], obs[:, 1])

    # 4. select to OSARG's observed spatial distribution and N
    m = R_obs < 8.0
    pos_o = np.column_stack([ob["x"], ob["y"], ob["z"]])
    r_o = cKDTree(pos_o).query(obs[m, :3], k=K_NN, workers=-1)[0][:, -1]
    r_e = cKDTree(obs[m, :3]).query(obs[m, :3], k=K_NN + 1, workers=-1)[0][:, -1]   # +1: the point itself
    ratio = (r_e / r_o) ** 3                                    # rho_OSARG / rho_mock up to a constant
    p = np.zeros(len(obs))
    lo, hi = 0.0, 1.0 / ratio.min()
    for _ in range(60):                                         # scale so that sum(min(c * ratio, 1)) = N_OSARG
        c = 0.5 * (lo + hi)
        lo, hi = (c, hi) if np.minimum(c * ratio, 1).sum() < n_osarg else (lo, c)
    p[m] = np.minimum(c * ratio, 1)
    sel = rng.random(len(p)) < p
    print(f"instances {len(p):,}; observed R < 8: {m.sum():,}; keep probability mean {p[m].mean():.4f}, "
          f"clipped at 1: {np.mean(p[m] >= 1):.2%}; kept {sel.sum():,} (target {n_osarg:,}); distinct orbits "
          f"{len(np.unique(orbit_id[sel])):,}, orbits kept more than once {np.sum(np.bincount(orbit_id[sel]) > 1):,}")

    names = ("x", "y", "z", "vx", "vy", "vz")
    with h5py.File(OUT, "w") as f:
        for i, k in enumerate(names):
            f[k] = obs[sel, i]
            f[f"{k}_true"] = w_true[sel, i]
        for k in ("l", "b", "pml", "pmb", "vlos"):
            f[k] = g[k][sel]
        f["d"] = d_obs[sel]
        f["d_true"] = g["d"][sel]
        f["DM_true"] = 5 * np.log10(g["d"][sel] * 100.0)
        f["DM_obs"] = 5 * np.log10(d_obs[sel] * 100.0)
        f["DM_err"] = np.full(sel.sum(), sigma_dm)
        f["orbit_id"] = orbit_id[sel]
        f["source_id"] = orbit_id[sel]
        f["time_gyr"] = t_inst[sel]
        f.attrs.update(dict(
            description=__doc__, source=OSARG.name, galaxy="mw", sigma_DM=sigma_dm, R_max_obs=8.0,
            Omega_abs_true=OMEGA_TRUE, phi_b_deg_true=PHI_TRUE, S_true=S_TRUE, bar_frac_true=F_TRUE,
            T_mix_gyr=T_MIX, T_start_gyr=T_START, dt_gyr=DT, k_nn=K_NN, seed=SEED))
    print(f"wrote {rel(OUT)}")
    print(summary("OSARG", *(ob[k] for k in names)))
    print(summary("mock", *(obs[sel, i] for i in range(6))))
    # the same stars at t = 0 vs phase-mixed, both TRUE (no new distance error), all R: does H24 keep OSARG's kinematics?
    w0 = w_mw
    late = traj[:, keep_t][:, -1]
    late_mw = ssl.rotate_z(late, -np.radians(PHI_TRUE))
    print("same orbits, t = 0 (observed OSARG) vs t = 6 Gyr (mixed in the truth potential), no selection:")
    print(summary("t=0", *w0.T))
    print(summary("t=6", *late_mw.T))


if __name__ == "__main__":
    main()

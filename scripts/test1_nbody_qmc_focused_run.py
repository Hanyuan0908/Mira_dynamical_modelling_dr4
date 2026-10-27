"""Test 1, N-body: denser, focused global search for the Fourier moment-map likelihood (overnight run, 2026-10-27).

Design (results/test1_nbody_qmc_focused/design.csv):
  - the truth and the 64 Sobol points of test1_nbody_qmc.py, re-used from their saved trajectories (batch "old");
  - N_NEW new scrambled Sobol points (batch "new") in a box focused on the posterior in Omega_b and phi_b and wider in
    S and bar_frac, so that the S-bar_frac ridge (which ran into the old box corner S 0.75, f 1.6) is contained.
Integration as in test1_nbody_qmc.py (true N-body potential family, merged CylSpline, accuracy 1e-5, 1 Gyr, 41 snapshots).

Saved per point in ~/mira_dr4_runs/test1_nbody_qmc_focused/:
  lib_NNN.npz   library sums on a fine folded grid (R 0.1 kpc to 8 kpc x phi mod pi 7.5 deg x |z| 0.25 kpc to 4 kpc + one
                bin beyond): N (samples), W = sum w, W2 = sum w^2, S1 = sum w v for v = (v_R, v_phi, v_z sign z),
                S2 = sum w v_a v_b for (RR, pp, zz, Rp, Rz, pz); any coarser (R, |z|, phi) grid with edges on these is
                an exact sum. Also tau, the polar cell moments of the current likelihood (lib_m1, lib_s), params, timing.
  traj_NNN.npy  trajectories of the new points (bar frame, float32, (N, 41, 6)); old points: ../test1_nbody/traj_NNN.npy
Library weights 1/tau (R < 8 kpc selection), snapshots t >= 0.25 Gyr; folding (phi mod pi, |z|, v_z sign z) is equivalent
to the four-image library. Resumable: points with an existing lib_NNN.npz are skipped.
"""
import csv
import time
from pathlib import Path

import numpy as np
import h5py
from scipy.stats import qmc

from common import ROOT
import steady_state_likelihood as ssl
import test1_nbody_qmc as old

OUT_DIR = Path.home() / "mira_dr4_runs" / "test1_nbody_qmc_focused"
OLD_DIR = old.ARR_DIR
DESIGN_CSV = ROOT / "results" / "test1_nbody_qmc_focused" / "design.csv"
N_NEW = 192
SEED = 2028
ACCURACY = 1e-5
TRUTH = old.TRUTH
BOX = dict(Omega_abs=(20.0, 32.0), phi_b_deg=(15.0, 39.0), lnS=(np.log(0.6), np.log(1.65)),
           lnf=(np.log(0.45), np.log(2.0)))
FINE_R = np.round(np.arange(0.0, 8.0 + 1e-9, 0.1), 6)
FINE_PHI = np.radians(np.arange(0.0, 180.0 + 1e-9, 7.5))
FINE_Z = np.r_[np.round(np.arange(0.0, 4.0 + 1e-9, 0.25), 6), np.inf]
FINE_SHAPE = (len(FINE_R) - 1, len(FINE_PHI) - 1, len(FINE_Z) - 1)
PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
FIELDS = ["idx", "kind", "batch", "Omega_abs", "phi_b_deg", "S", "bar_frac", "traj_file", "tau_median", "t_integrate_s",
          "t_process_s"]


def design():
    """Old design (truth + 64 Sobol, trajectories exist) followed by N_NEW new Sobol points in BOX."""
    pts = [dict(kind=p["kind"], batch="old", traj_file=str(OLD_DIR / f"traj_{i:03d}.npy"),
                **{k: p[k] for k in ("Omega_abs", "phi_b_deg", "S", "bar_frac")}) for i, p in enumerate(old.design())]
    u = qmc.Sobol(d=4, scramble=True, seed=SEED).random(256)[:N_NEW]
    lo = np.array([v[0] for v in BOX.values()])
    hi = np.array([v[1] for v in BOX.values()])
    for a, b, c, d in lo + u * (hi - lo):
        pts.append(dict(kind="qmc", batch="new", Omega_abs=float(a), phi_b_deg=float(b), S=float(np.exp(c)),
                        bar_frac=float(np.exp(d)), traj_file=""))
    return pts


def fine_sums(w, weight):
    """Sums on the fine folded grid for phase-space points w (n, 6) with weights."""
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.mod(np.arctan2(w[:, 1], w[:, 0]), np.pi)
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5] * np.sign(w[:, 2])])
    i = np.digitize(R, FINE_R) - 1
    j = np.minimum(np.digitize(phi, FINE_PHI) - 1, FINE_SHAPE[1] - 1)
    k = np.digitize(np.abs(w[:, 2]), FINE_Z) - 1
    ok = (i >= 0) & (i < FINE_SHAPE[0]) & (j >= 0) & (k >= 0) & (k < FINE_SHAPE[2])
    flat = np.ravel_multi_index((i[ok], j[ok], k[ok]), FINE_SHAPE)
    n, wt, vv = int(np.prod(FINE_SHAPE)), weight[ok], v[ok]
    out = dict(N=np.bincount(flat, minlength=n), W=np.bincount(flat, weights=wt, minlength=n),
               W2=np.bincount(flat, weights=wt**2, minlength=n),
               S1=np.column_stack([np.bincount(flat, weights=wt * vv[:, a], minlength=n) for a in range(3)]),
               S2=np.column_stack([np.bincount(flat, weights=wt * vv[:, a] * vv[:, b], minlength=n) for a, b in PAIRS]))
    return {k: a.reshape(FINE_SHAPE + a.shape[1:]) for k, a in out.items()}


def library_samples(times, traj):
    """Unfolded library samples (t >= T_MIN) and their 1/tau weights, tau per orbit."""
    R = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])
    tau = np.clip(np.mean(R < ssl.R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
    keep = times >= ssl.T_MIN_GYR - 1e-9
    w = traj[:, keep].reshape(-1, 6).astype(np.float64)
    return w, np.repeat(1.0 / tau, keep.sum()), tau


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DESIGN_CSV.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(old.MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    times = np.load(OLD_DIR / "000.npz")["times"]
    ssl.ACCURACY = ACCURACY
    pts = design()
    new_csv = not DESIGN_CSV.exists()
    t_start = time.time()
    with open(DESIGN_CSV, "a", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        if new_csv:
            wr.writeheader()
        for i, p in enumerate(pts):
            out = OUT_DIR / f"lib_{i:03d}.npz"
            if out.exists():
                continue
            t0 = time.time()
            if p["batch"] == "old":
                traj = np.load(p["traj_file"])
                t_int = 0.0
            else:
                pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], "nbody")
                ic = ssl.rotate_z(stars, np.radians(p["phi_b_deg"]))
                tt, traj = ssl.integrate(pot, ic, -abs(p["Omega_abs"]))
                t_int = time.time() - t0
                np.save(OUT_DIR / f"traj_{i:03d}.npy", traj)
                p["traj_file"] = str(OUT_DIR / f"traj_{i:03d}.npy")
            t1 = time.time()
            w, wt, tau = library_samples(times, traj)
            sums = fine_sums(w, wt)
            _, lib_m1, lib_s = ssl.cell_moments(w, wt)
            np.savez(out, params=np.array([p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"]]), tau=tau.astype(np.float32),
                     lib_m1=lib_m1, lib_s=lib_s, fine_R=FINE_R, fine_phi=FINE_PHI, fine_z=FINE_Z,
                     **{f"fine_{k}": v for k, v in sums.items()})
            del traj, w, wt
            row = dict(idx=i, tau_median=float(np.median(tau)), t_integrate_s=round(t_int, 1),
                       t_process_s=round(time.time() - t1, 1), **p)
            wr.writerow(row)
            fh.flush()
            print(f"[{i + 1}/{len(pts)}] {p['batch']} {p['kind']:5s} Omega={p['Omega_abs']:.2f} phi={p['phi_b_deg']:.1f} "
                  f"S={p['S']:.3f} f={p['bar_frac']:.3f}  integrate {t_int:.0f} s  process {time.time() - t1:.0f} s  "
                  f"(elapsed {(time.time() - t_start) / 60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()

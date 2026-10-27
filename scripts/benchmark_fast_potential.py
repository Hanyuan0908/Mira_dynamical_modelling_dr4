"""Speed vs accuracy of fast potential set-ups for MCMC, integrating N mock orbits for T_GYR in the
bar frame (agama.orbit with Omega = Omega_b).

Set-ups: the exact make_potential (three CylSplines, mmax = 8) and make_static_fast re-expansions
with different mmax / symmetry / integrator accuracy. Accuracy measures:
  - force error vs the exact potential, relative to the bar (non-axisymmetric) force;
  - Jacobi-integral conservation;
  - after T_GYR, how far the integrated sample's bar-frame kinematics differ from those of the
    exact integration: Delta<v_R>, Delta<v_phi>, Delta sigma maps on 1 kpc pixels.
Usage: python benchmark_fast_potential.py [--n 5000]
"""
import argparse
import json
import resource
import time

import numpy as np
import h5py
import agama

from common import ROOT, MOCK_DIR
from potential_family import (make_potential, make_static_fast, load_axi, OMEGA_B_TRUE, TIME_UNIT_GYR)
from time_orbit_integration import jacobi

T_GYR = 1.0
MOCK = MOCK_DIR / "mock_perfect_R8_N50k.hdf5"
RESULTS_FILE = ROOT / "results" / "fast_potential_benchmark.json"
SETUPS = [  # name, builder, accuracy
    ("exact mmax8, acc1e-8", "exact", dict(), 1e-8),
    ("exact mmax8, acc1e-7 (chaos floor)", "exact", dict(), 1e-7),
    ("exact mmax8, acc1e-6 (chaos floor)", "exact", dict(), 1e-6),
    ("merged mmax8 none, acc1e-6", "fast", dict(mmax=8, symmetry="none"), 1e-6),
    ("merged mmax8 bisym, acc1e-6", "fast", dict(mmax=8, symmetry="bisymmetric"), 1e-6),
    ("merged mmax6 bisym, acc1e-6", "fast", dict(mmax=6, symmetry="bisymmetric"), 1e-6),
    ("merged mmax8 reflection, acc1e-6", "fast", dict(mmax=8, symmetry="reflection"), 1e-6),
    ("merged mmax4 none, acc1e-6", "fast", dict(mmax=4, symmetry="none"), 1e-6),
    ("merged mmax4 bisym, acc1e-6", "fast", dict(mmax=4, symmetry="bisymmetric"), 1e-6),
    ("merged mmax4 bisym, acc1e-5", "fast", dict(mmax=4, symmetry="bisymmetric"), 1e-5),
    ("merged mmax2 bisym, acc1e-5", "fast", dict(mmax=2, symmetry="bisymmetric"), 1e-5),
]


def kin_maps(w, edges):
    """Mass-free mean v_R, v_phi and their dispersions in face-on pixels (bar frame)."""
    x, y = w[:, 0], w[:, 1]
    R = np.hypot(x, y)
    vR = (x * w[:, 3] + y * w[:, 4]) / R
    vphi = -(x * w[:, 4] - y * w[:, 3]) / R
    N = np.histogram2d(x, y, bins=[edges, edges])[0]
    out = {}
    for k, q in (("vR", vR), ("vphi", vphi)):
        s1 = np.histogram2d(x, y, bins=[edges, edges], weights=q)[0]
        s2 = np.histogram2d(x, y, bins=[edges, edges], weights=q * q)[0]
        with np.errstate(invalid="ignore", divide="ignore"):
            m = s1 / N
            out[k] = np.where(N >= 30, m, np.nan)
            out["s" + k] = np.where(N >= 30, np.sqrt(s2 / N - m * m), np.nan)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=5000)
    args = p.parse_args()
    with h5py.File(MOCK, "r") as f:
        ic = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    ic = ic[np.sort(np.random.default_rng(1).choice(len(ic), args.n, replace=False))]   # random subset
    T = T_GYR / TIME_UNIT_GYR
    axi = load_axi()
    exact = make_potential(Omega_b=0.0, pot_axi=axi)

    rng = np.random.default_rng(0)
    R = np.sqrt(rng.uniform(0, 100, 20000))
    ph = rng.uniform(0, 2 * np.pi, R.size)
    pts = np.column_stack([R * np.cos(ph), R * np.sin(ph), rng.normal(0, 0.8, R.size)])
    Fe = exact.force(pts)
    Fb = np.linalg.norm(Fe - axi.force(pts), axis=1)

    edges = np.arange(-8, 8.001, 1.0)
    res, ref_maps = dict(N=len(ic), T_Gyr=T_GYR), None
    for name, kind, kw, acc in SETUPS:
        t0 = time.time()
        pot = exact if kind == "exact" else make_static_fast(**kw)
        t_build = time.time() - t0
        dF = np.linalg.norm(pot.force(pts) - Fe, axis=1)
        dF_R = {f"{a:g}-{b:g}": float(np.sqrt(np.mean(dF[(R >= a) & (R < b)]**2)) /
                                    np.sqrt(np.mean(Fb[(R >= a) & (R < b)]**2)))
                for a, b in ((0, 1), (1, 3), (3, 6), (6, 10))}
        t0 = time.time()
        orbs = agama.orbit(potential=pot, ic=ic, time=T, trajsize=2, Omega=OMEGA_B_TRUE, accuracy=acc,
                           verbose=False)
        dt = time.time() - t0
        end = np.array([o[-1] for o in orbs[:, 1]])
        dEJ = np.abs(jacobi(pot, end) / jacobi(pot, ic) - 1)
        maps = kin_maps(end, edges)
        if ref_maps is None:
            ref_maps = maps
        dmap = {k: float(np.sqrt(np.nanmean((maps[k] - ref_maps[k])**2))) for k in maps}
        res[name] = dict(build_s=round(t_build, 2), wall_s=round(dt, 2),
                         min_for_50k=round(dt * 50000 / len(ic) / 60, 2),
                         force_err_rms_over_bar=float(np.sqrt(np.mean(dF**2)) / np.sqrt(np.mean(Fb**2))),
                         force_err_max_over_bar=float(dF.max() / Fb.max()),
                         force_err_rms_over_bar_by_R=dF_R,
                         dEJ_median=float(np.median(dEJ)), dEJ_max=float(dEJ.max()),
                         rms_map_diff_vs_exact_kms=dmap)
        r = res[name]
        print(f"{name:36s} build {t_build:5.2f}s  {dt:6.1f}s/{len(ic)} -> 50k {r['min_for_50k']:5.2f} min | "
              f"dF/F_bar rms {r['force_err_rms_over_bar']:.1e} max {r['force_err_max_over_bar']:.1e} | "
              f"dE_J med {r['dEJ_median']:.0e} | map rms diff vR {dmap['vR']:.1f} vphi {dmap['vphi']:.1f} "
              f"svR {dmap['svR']:.1f} km/s | dF/F_bar by R " + " ".join(f"{k}:{v:.3f}" for k, v in dF_R.items()))
    # the 'chaos floor' row (exact potential, slightly different integrator accuracy) shows how much the
    # maps change from reshuffling chaotic orbits alone; fast set-ups are only worse if they exceed it
    res["peak_rss_GB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2)
    with open(RESULTS_FILE, "w") as f:
        json.dump(res, f, indent=2)


if __name__ == "__main__":
    main()

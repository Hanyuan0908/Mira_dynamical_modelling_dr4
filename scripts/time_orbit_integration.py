"""Wall-clock cost of integrating the Test 1 mock for T_GYR in the true potential, in two equivalent ways:
  inertial : make_potential() with the rotation modifier, integrated in the inertial frame
  rotating : static bar (Omega_b = 0) integrated in the bar frame with agama.orbit(Omega=Omega_b)
Accuracy check: the Jacobi integral E_J = v^2/2 + Phi - Omega_b L_z (bar frame) should be conserved.
Usage: python time_orbit_integration.py [--n N] [--methods inertial rotating]
"""
import argparse
import json
import resource
import time

import numpy as np
import h5py
import agama

from common import ROOT, MOCK_DIR
from potential_family import make_potential, load_axi, to_bar_frame, OMEGA_B_TRUE, TIME_UNIT_GYR

T_GYR = 1.0
TRAJSIZE = 101                 # snapshots every 10 Myr
MOCK = MOCK_DIR / "mock_perfect_R8_N50k.hdf5"
RESULTS_FILE = ROOT / "results" / "orbit_timing.json"


def jacobi(pot_static, w):
    """E_J in the bar frame (positions in bar axes, inertial velocities along bar axes)."""
    Lz = w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]
    return 0.5 * np.sum(w[:, 3:]**2, axis=1) + pot_static.potential(w[:, :3]) - OMEGA_B_TRUE * Lz


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=None)
    p.add_argument("--methods", nargs="+", default=["inertial", "rotating"])
    args = p.parse_args()

    with h5py.File(MOCK, "r") as f:
        ic = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    if args.n:
        ic = ic[: args.n]
    T = T_GYR / TIME_UNIT_GYR
    axi = load_axi()
    pots = {"inertial": make_potential(pot_axi=axi), "rotating": make_potential(Omega_b=0.0, pot_axi=axi)}
    out = dict(N=len(ic), T_Gyr=T_GYR, trajsize=TRAJSIZE, threads="all (OpenMP default)", accuracy=1e-8)
    for m in args.methods:
        kw = dict(Omega=OMEGA_B_TRUE) if m == "rotating" else {}
        t0 = time.time()
        orbs = agama.orbit(potential=pots[m], ic=ic, time=T, trajsize=TRAJSIZE, verbose=False, **kw)
        dt = time.time() - t0
        tgrid = orbs[0, 0]
        end = np.array([o[-1] for o in orbs[:, 1]])
        if m == "inertial":
            end = to_bar_frame(end, tgrid[-1])
        dEJ = np.abs(jacobi(pots["rotating"], end) / jacobi(pots["rotating"], ic) - 1)
        out[m] = dict(wall_s=round(dt, 1), orbits_per_s=round(len(ic) / dt, 1),
                      dEJ_rel_median=float(np.median(dEJ)), dEJ_rel_99pct=float(np.percentile(dEJ, 99)),
                      dEJ_rel_max=float(dEJ.max()))
        print(f"{m:9s}: {dt:7.1f} s for {len(ic):,} orbits x {T_GYR} Gyr ({len(ic) / dt:.0f} orbits/s); "
              f"|dE_J/E_J| median {np.median(dEJ):.1e}, 99% {np.percentile(dEJ, 99):.1e}, max {dEJ.max():.1e}")
        del orbs
    out["peak_rss_GB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2)
    print(f"peak RSS {out['peak_rss_GB']} GB")
    if not args.n:
        with open(RESULTS_FILE, "w") as f:
            json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()

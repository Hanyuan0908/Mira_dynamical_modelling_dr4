"""Independent check of the N-body bar's pattern speed: track the m = 2 phase between snapshots.

`measure_nbody_pattern_speed_evolution.py` uses the Dehnen, Semczuk & Schoenrich (2023) Fourier method, which
infers Omega_b from the velocities of a single snapshot. This script instead measures where the bar actually
*is* in consecutive snapshots and differences the phase - no kinematics, no method assumptions, just geometry.

The phase is defined modulo pi, and the snapshots are 0.125 model time units apart, so the measurement is
aliased: Omega and Omega + pi/dt = 25.13 km/s/kpc fit equally well. That is wide enough that a single solution
survives in any sensible prior range, so the aliasing is harmless here (it is only fatal if you try to
difference without a prior, which is why the Fourier method exists).

All snapshots are put in **one fixed frame** - the rotation is computed once, at the reference time, and
re-used - because recomputing it per snapshot lets the frame itself wobble in azimuth and that wobble is
indistinguishable from bar rotation.

Usage: python check_nbody_pattern_speed_independent.py [--time 7] [--half-window 2] [--comp all|disc|bulge]
"""
import argparse

import numpy as np

import auriga as au
import nbody
from common import rel
from prepare_au18_mwframe import shrinking_sphere, R_DISC

CENTRE_RADII = (20, 10, 5, 3, 2, 1, 0.5)
Z_SLAB = 2.0
ANNULI = ((1.0, 4.0), (0.5, 3.0), (1.0, 5.0), (2.0, 5.0))
OMEGA_RANGE = (12.0, 38.0)        # prior; alias spacing is pi/dt = 25.13 km/s/kpc for dt = 0.125


def load(t, rot=None, comp_sel="all"):
    """Positions and masses of the chosen component at model time t, in a frame fixed by `rot`."""
    d = nbody.read_nemo(nbody.snapshot(t), tags=("Mass", "Position", "Velocity"))
    comp = nbody.components(d["Mass"])
    st = {"all": comp < 2, "disc": comp == 0, "bulge": comp == 1}[comp_sel]
    pos = d["Position"].astype(np.float64)
    vel = d["Velocity"].astype(np.float64)
    m = d["Mass"].astype(np.float64) * nbody.MASS_UNIT_MSUN
    stars = comp < 2
    pos -= shrinking_sphere(pos[stars], m[stars], radii=CENTRE_RADII)
    if rot is None:                                     # define the frame once, from this snapshot
        s = stars & (np.linalg.norm(pos, axis=1) < R_DISC)
        vel -= np.average(vel[s], axis=0, weights=m[s])
        L = np.sum(m[s, None] * np.cross(pos[s], vel[s]), axis=0)
        rot = au.rotation_to(L, np.array([0.0, 0.0, -1.0]))
    return pos[st] @ rot.T, m[st], rot


def phase(x, m, r0, r1):
    """m = 2 phase [rad] of the stars in an annulus, in (-pi/2, pi/2]."""
    slab = np.abs(x[:, 2]) < Z_SLAB
    _, ph = au.fourier_m2(x[slab, 0], x[slab, 1], m[slab], np.array([r0, r1]))
    return float(ph[0])


def fit_omega(times, phases, lo=OMEGA_RANGE[0], hi=OMEGA_RANGE[1], n=260001):
    """Omega minimising the wrapped residuals of psi(t) = psi_0 - Omega t (clockwise frame).

    Returns Omega, the rms residual [deg], and the formal error from that scatter."""
    t = np.asarray(times) - np.mean(times)
    p = np.asarray(phases)
    grid = np.linspace(lo, hi, n)
    best, res = None, None
    for om in grid:
        pred = -om * t
        r = np.mod(p - pred + np.pi / 2, np.pi) - np.pi / 2
        r = r - np.mean(r)                                  # psi_0 is a free offset
        r = np.mod(r + np.pi / 2, np.pi) - np.pi / 2
        s = np.sum(r**2)
        if best is None or s < best:
            best, om_best, res = s, om, r
    rms = float(np.sqrt(np.mean(res**2)))
    # slope error for a least-squares line with this scatter
    err = rms / np.sqrt(np.sum(t**2)) if np.sum(t**2) > 0 else np.nan
    return om_best, np.degrees(rms), err


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--time", type=float, default=7.0)
    ap.add_argument("--half-window", type=int, default=2, help="snapshots either side of --time")
    ap.add_argument("--comp", default="all", choices=["all", "disc", "bulge"])
    args = ap.parse_args()

    step = 0.125
    times = [round(args.time + k * step, 6) for k in range(-args.half_window, args.half_window + 1)]
    times = [t for t in times if nbody.snapshot(t).exists()]
    print(f"frame fixed at t = {args.time:g}; snapshots {times[0]:g} to {times[-1]:g} "
          f"({len(times)}), component {args.comp}")
    print(f"alias spacing pi/dt = {np.pi / step:.2f} km/s/kpc, prior {OMEGA_RANGE}\n")

    rot = load(args.time, None, args.comp)[2]
    data = {}
    for t in times:
        x, m, _ = load(t, rot, args.comp)
        data[t] = {a: phase(x, m, *a) for a in ANNULI}
        del x, m

    print(f"{'annulus [kpc]':>14s} {'Omega_b':>9s} {'rms resid':>10s} {'formal err':>11s}")
    out = []
    for a in ANNULI:
        ph = [data[t][a] for t in times]
        om, rms, err = fit_omega(times, ph)
        out.append(om)
        print(f"{a[0]:6.1f}-{a[1]:<6.1f} {om:9.2f} {rms:9.1f} deg {err:10.2f}")
    print(f"\nmedian over annuli: {np.median(out):.2f} km/s/kpc  (spread {np.ptp(out):.2f})")


if __name__ == "__main__":
    main()

"""Orbital frequencies in the bar frame, for the dynamical bar length.

The saved orbit library cannot be re-used for a frequency-ratio classification: its 25/50 Myr
cadence puts Nyquist (20 / 10 per Gyr) BELOW the median Omega_R of the x1 orbits in every galaxy,
so 73-99% of the radial frequencies alias and the measured Omega_R/Omega_x is wrong by 0.4-1.2
against a tolerance of +-0.2 (`diagnose_x1_frequency_sampling.py`). The baseline is not the problem:
1 Gyr at 1 Myr recovers 83-93% of the reference x1 set at +-0.06-0.09. So this re-integrates the
same orbits over the same 1 Gyr, sampled 25x more finely, and keeps only the frequencies.

Cost, measured on 5000 orbits: trajsize 41 -> 1001 is +2% wall clock (the force evaluations are
unchanged; AGAMA's `trajsize` only sets how often the interpolator is sampled), and the five
frequencies add 0.40 ms per orbit, i.e. ~10% on top. Nothing dense is stored - the output is one
small table per potential.

Frequencies follow Zhang et al. (2024, arXiv:2406.06678) / `/data/hz420/OSARGs/x1_orbits_obs.py`:
the highest peak of the FFT periodogram, taken of the COMPLEX signal (coordinate + i velocity) so
that the frequency is signed, and of the bare coordinate as well for x, y, z and R. Omega_phi has
no bare counterpart because phi wraps; its complex signal is the unit vector (x + iy)/R, whose
periodogram peak is the azimuthal frequency in the bar frame (circulation for a passing orbit,
libration for a trapped one). The `_withv` set is the one Zhang et al. use for the ratios.

Two modes.
  --run RUN     the FITTED galaxy. Potentials are drawn from the posterior re-weighted by the
                kinematic phi_b prior (`plot_corner_with_kinematic_prior.weighted_posterior`), so
                every draw is a (Omega_b, phi_b, S, bar_frac) the data actually supports. The
                orbits are the run's own mock stars, OBSERVED - distance errors and selection
                function included, exactly the sample the model was fitted to and the analogue of
                what Zhang et al. integrate for the real Miras.
  --truth       the TRUE galaxy. Potential is the family at (bar_frac, S) = (1, 1) with the
                tabulated Omega_b, which is identically the galaxy's own potential, and the orbits
                are the original snapshot particles: no distance error, no selection function.

The star pool is R < 8 kpc in both modes (`--rmax`), so that a difference between truth and model
is a difference in the potential and not in which stars were available.

Usage:
  python integrate_orbit_frequencies.py --run test2_au26_sf_12p5pct [--n-draw 30] [--n-star 10000]
  python integrate_orbit_frequencies.py --truth --galaxy au26 [--n-star 50000]
  ... --stride 4 --offset 0        one worker of four (draws are split; each writes its own file)
Output: results/bar_dynamics/<run>/draw_NNN.npz,  results/bar_dynamics/truth_<galaxy>.npz
"""
import argparse
import json
import time

import numpy as np
import h5py
import agama
from scipy.signal import find_peaks

from common import RESULTS_DIR, MOCK_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy as galaxy_entry
import steady_state_likelihood as ssl
import potential_family as pf
import test1_gp_emulator as ge
from plot_corner_with_kinematic_prior import (SCHEME_B, kinematic_prior, measure_from_mock,
                                              weighted_posterior)
# the pipeline's own truth, NOT potential_family.omega_true: the N-body snapshots have no
# <gal>_ground_truth.json and their pattern speed lives in nbody_pattern_speed_evolution.csv.
# Using omega_true here killed 11 of 25 truth integrations on 2026-11-07 - the same trap recorded
# on 2026-11-05. Always go through the function the pipeline scores against.
from test2b_qmc_distance_error_run import truth_omega

T_GYR = 1.0           # same baseline as the orbit library
N_SNAP = 1001         # dt = 1 Myr
CHUNK = 2000          # orbits per agama.orbit call, to bound memory
OUT_DIR = RESULTS_DIR / "bar_dynamics"

COLS = ["Omega_x", "Omega_y", "Omega_z", "Omega_R",
        "Omega_x_withv", "Omega_y_withv", "Omega_z_withv", "Omega_R_withv", "Omega_phi_withv",
        "r_apo", "R_apo", "R_apo_mean", "n_apo", "x_max", "y_max", "z_max", "R_ini"]

# r_apo / R_apo are the largest 3D and cylindrical radii reached; R_apo_mean averages the radii of
# the individual apocentre passages (find_peaks on R(t), as in OSARGs x1_orbits_obs.py), which is
# the steadier extent for an orbit that is not closed. n_apo is how many passages that average is
# over - a frequency measured from fewer than a handful of radial periods should not be trusted.


def _peak(t, y):
    """Highest peak of the FFT periodogram (OSARGs `extract_orbital_period`, n_peaks = 1)."""
    n = len(t)
    f = np.fft.fftfreq(n, t[1] - t[0])
    p = (1.0 / n) * abs(np.fft.fft(y)) ** 2
    m = f > 0
    f, p = f[m], p[m]
    pk, _ = find_peaks(p)
    return float(f[pk][np.argmax(p[pk])]) if len(pk) else np.nan


def frequencies(tr, t):
    """(n, len(COLS)) table from a (n, n_snap, 6) bar-frame trajectory block."""
    x, y, z, vx, vy, vz = (tr[:, :, k] for k in range(6))
    R = np.hypot(x, y)
    Rs = np.where(R > 0, R, 1.0)
    vR = (vx * x + vy * y) / Rs
    out = np.empty((len(tr), len(COLS)), dtype=np.float32)
    for i in range(len(tr)):
        out[i, :4] = [_peak(t, s[i]) for s in (x, y, z, R)]
        out[i, 4:9] = [_peak(t, a[i] + 1j * b[i]) for a, b in
                       ((x, vx), (y, vy), (z, vz), (R, vR))] + \
                      [_peak(t, (x[i] + 1j * y[i]) / Rs[i])]
    out[:, 9] = np.sqrt(x ** 2 + y ** 2 + z ** 2).max(axis=1)
    out[:, 10] = R.max(axis=1)
    for i in range(len(tr)):
        ap, _ = find_peaks(R[i])
        out[i, 11] = R[i][ap].mean() if len(ap) else R[i].max()
        out[i, 12] = len(ap)
    out[:, 13] = np.abs(x).max(axis=1)
    out[:, 14] = np.abs(y).max(axis=1)
    out[:, 15] = np.abs(z).max(axis=1)
    out[:, 16] = R[:, 0]
    return out


def dense_pass(pot, ic, Omega_abs, n_snap=N_SNAP, thin_to=0, chunk=CHUNK):
    """Integrate once at `n_snap` output points and return (thinned trajectory, frequency table).

    Sampling densely and thinning is BITWISE identical to integrating at the coarse trajesize -
    verified on 800 orbits, max |diff| exactly 0 - because AGAMA's `trajsize` only sets how often
    the interpolator is sampled, not the integration. So a pipeline run can record frequencies
    without changing any fitted number, for +2% integration and ~10% FFT.

    thin_to = 0 returns no trajectory (the frequency-only use); otherwise n_snap - 1 must be a
    multiple of thin_to - 1, and the returned trajectory is the usual coarse grid.
    """
    step = 0
    if thin_to:
        if (n_snap - 1) % (thin_to - 1):
            raise SystemExit(f"n_snap - 1 = {n_snap - 1} is not a multiple of thin_to - 1 = "
                             f"{thin_to - 1}; the coarse grid would not be a subset")
        step = (n_snap - 1) // (thin_to - 1)
    t = np.linspace(0.0, T_GYR, n_snap)
    tab = np.empty((len(ic), len(COLS)), dtype=np.float32)
    traj = np.empty((len(ic), thin_to, 6), dtype=np.float32) if thin_to else None
    for a in range(0, len(ic), chunk):
        b = min(a + chunk, len(ic))
        orbs = agama.orbit(potential=pot, ic=ic[a:b], time=T_GYR / ssl.TIME_UNIT_GYR,
                           trajsize=n_snap, Omega=-abs(Omega_abs), accuracy=ssl.ACCURACY,
                           **ssl.ORBIT_KW)
        dense = np.stack(orbs[:, 1])
        tab[a:b] = frequencies(dense, t)
        if thin_to:
            traj[a:b] = dense[:, ::step].astype(np.float32)
        del dense, orbs
    return traj, tab


def run_one(galaxy, ic, bar_frac, S, Omega_abs, label=""):
    """Integrate `ic` (inertial frame, bar on x) for T_GYR and return the frequency table."""
    t0 = time.time()
    _, tab = dense_pass(ssl.bar_aligned_potential(bar_frac, S, galaxy), ic, Omega_abs)
    print(f"    {label} bar_frac {bar_frac:.3f} S {S:.3f} Omega {Omega_abs:.2f}: "
          f"{len(ic):,} orbits in {time.time() - t0:.0f} s", flush=True)
    return tab


def mock_stars(mock_name, rmax, n_star, rng):
    """Observed phase space of the run's own mock: distance errors and selection included."""
    with h5py.File(MOCK_DIR / mock_name, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    w = w[np.hypot(w[:, 0], w[:, 1]) < rmax]
    if n_star and n_star < len(w):
        w = w[np.sort(rng.choice(len(w), n_star, replace=False))]
    return w


def snapshot_stars(gal, rmax, n_star, rng):
    """Original snapshot particles: no distance error, no selection function."""
    with h5py.File(galaxy_entry(gal)["stars"], "r") as f:
        pool = np.flatnonzero(np.hypot(f["x"][:], f["y"][:]) < rmax)
        if n_star and n_star < len(pool):
            pool = np.sort(rng.choice(pool, n_star, replace=False))
        return np.column_stack([f[k][:][pool].astype(np.float64)
                                for k in ("x", "y", "z", "vx", "vy", "vz")])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run")
    ap.add_argument("--galaxy")
    ap.add_argument("--truth", action="store_true", help="true potential, snapshot particles")
    ap.add_argument("--column", default=f"logLs_{SCHEME_B}_sfratio")
    ap.add_argument("--n-draw", type=int, default=30)
    ap.add_argument("--n-star", type=int, default=10000)
    ap.add_argument("--rmax", type=float, default=8.0)
    ap.add_argument("--rule", default="self:0.3")
    ap.add_argument("--est", default="absvR")
    ap.add_argument("--flat-posterior", action="store_true",
                    help="draw from the unweighted posterior instead of the phi_b-prior one")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)

    if args.truth:
        gal = args.galaxy or (json.load(open(RESULTS_DIR / args.run / "run_meta.json"))["galaxy"])
        out = OUT_DIR / f"truth_{gal}.npz"
        if out.exists():
            print(f"  {rel(out)} exists; nothing to do")
            return
        Om = abs(truth_omega(gal))
        w = snapshot_stars(gal, args.rmax, args.n_star, rng)
        print(f"truth {gal}: Omega_b {Om:.2f}, {len(w):,} snapshot particles with R < {args.rmax:g} "
              f"kpc, no distance error, no selection", flush=True)
        ic = ssl.rotate_z(w, np.deg2rad(BAR_ANGLE_DEG))
        tab = run_one(gal, ic, 1.0, 1.0, Om, label="truth")
        np.savez_compressed(out, table=tab, cols=np.array(COLS), galaxy=gal, Omega_abs=Om,
                            bar_frac=1.0, S=1.0, phi_b_deg=BAR_ANGLE_DEG, n_star=len(w),
                            rmax=args.rmax, kind="truth", T_Gyr=T_GYR, n_snap=N_SNAP)
        print(f"wrote {rel(out)}")
        return

    if not args.run:
        raise SystemExit("give --run RUN (fitted galaxy) or --truth --galaxy GAL")
    meta = json.load(open(RESULTS_DIR / args.run / "run_meta.json"))
    gal, mock = meta["galaxy"], meta["mock"]
    ge.use_run(args.run)
    meas = measure_from_mock(mock, args.rule, args.est)
    phi_kin, sig_kin, info = kinematic_prior(gal, args.rule, args.est, measured=meas)
    post, w8, ess = weighted_posterior(args.column, phi_kin, sig_kin)
    if args.flat_posterior:
        w8 = np.ones(len(post))
    idx = rng.choice(len(post), args.n_draw, replace=True, p=w8 / w8.sum())
    draws = post[idx]                                   # columns: Omega_abs, phi_b_deg, S, bar_frac

    d = OUT_DIR / args.run
    d.mkdir(parents=True, exist_ok=True)
    print(f"{args.run} ({gal}): kinematic phi_b = {phi_kin:.2f} +- {sig_kin:.2f} deg, "
          f"ESS {ess:.0f}; {args.n_draw} draws, {args.n_star:,} stars each", flush=True)
    stars = mock_stars(mock, args.rmax, args.n_star, rng)
    for k in range(args.offset, args.n_draw, args.stride):
        out = d / f"draw_{k:03d}.npz"
        if out.exists():
            continue
        Om, phi, S, bf = draws[k]
        ic = ssl.rotate_z(stars, np.deg2rad(phi))
        tab = run_one(gal, ic, bf, S, Om, label=f"[{k + 1:3d}/{args.n_draw}]")
        np.savez_compressed(out, table=tab, cols=np.array(COLS), galaxy=gal, run=args.run,
                            Omega_abs=Om, phi_b_deg=phi, S=S, bar_frac=bf, n_star=len(stars),
                            rmax=args.rmax, kind="posterior", phi_kin=phi_kin, sigma_kin=sig_kin,
                            flat_posterior=bool(args.flat_posterior), T_Gyr=T_GYR, n_snap=N_SNAP)
    print(f"wrote {len(list(d.glob('draw_*.npz')))}/{args.n_draw} draws to {rel(d)}")


if __name__ == "__main__":
    main()

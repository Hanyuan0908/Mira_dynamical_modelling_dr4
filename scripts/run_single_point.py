"""Run the pipeline at one parameter point on one mock: integrate, fit the Fourier modes, report chi2/N.

This is the per-point work of the global searches (`test1_nbody_qmc_focused_run.py`,
`test2b_qmc_distance_error_run.py`) for a single trial potential, so that one data set can be scored without
launching a whole Sobol design. The default point is the galaxy's ground truth, which is the reference every
search is measured against: chi2/N there is the floor the method can reach on that data.

Steps (all reusing the pipeline modules):
  1. `ssl.bar_aligned_potential(bar_frac, S, galaxy)`  - the trial potential (ground truth at bar_frac = S = 1)
  2. `ssl.integrate(...)`                              - 50k orbits, 1 Gyr, 41 snapshots, in the bar frame
  3. `fs.data_fit(...)`                                - per-star Fourier coefficients of the data + N_BOOT bootstraps
  4. `fs.library_fit(...)`                             - the same fit to the time-averaged orbit library
  5. `fs.scheme_chi2(...)`                             - chi2 and N_data for each scheme in fs.SCHEMES

Two library variants are scored whenever the mock carries a distance error (as in
`test2b_qmc_distance_error_run.py`):
  asis      the library as integrated - the data are smeared but the model is not, so part of the chi2 is
            just that mismatch;
  observed  the library observed again with the same distance error, which is the like-for-like comparison
            and the recipe chosen for the 7.5% search.

Saved under $MIRA_DR4_RUNS/<run>/ (default ROOT/runs):
  traj_<tag>.npy      trajectories, bar frame, float32 (N, 41, 6)
  point_<tag>.npz     coefficients (data and library), the joint bootstrap covariance, star counts,
                      chi2/ndata for every scheme, the parameters and the timings
  results/<run>/single_point.csv   one row per point run, appended

Usage:
    python run_single_point.py --mock nbody_mock_distance_error_sigmaDM0.217_R8_N50k.hdf5
    python run_single_point.py --mock <file> --omega 28 --phi 24 --S 0.9 --bar-frac 1.2
"""
import argparse
import json
import time

import numpy as np
import h5py

from common import MOCK_DIR, RESULTS_DIR, RUNS_DIR, rel
from galaxies import galaxy
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
from plot_distance_error_effect_on_fourier_coefficients import forward_smeared_library_fit

ACCURACY = 1e-5          # as in the global searches (ssl's own default 1e-6 is for the pilot runs)
MAIN_SCHEME = "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"


def truth_params(name):
    """Ground truth (|Omega_b|, phi_b, S, bar_frac) of a galaxy, from its measured pattern speed."""
    g = galaxy(name)
    if g["truth"].exists():
        omega = float(json.load(open(g["truth"]))["Omega_p_kms_kpc"])
    else:                                    # fall back on the snapshot evolution table
        import pandas as pd
        csv_path = RESULTS_DIR / "nbody_pattern_speed_evolution.csv"
        if not csv_path.exists():
            raise SystemExit(f"no truth for {name}: neither {rel(g['truth'])} nor {rel(csv_path)} exists")
        df = pd.read_csv(csv_path)
        hit = df[np.isclose(df.t_model, g.get("time", np.nan))]
        if not len(hit):
            raise SystemExit(f"{rel(csv_path)} has no entry for snapshot time {g.get('time')}")
        omega = float(hit.Omega_p.iloc[0])
    from common import BAR_ANGLE_DEG
    return dict(Omega_abs=omega, phi_b_deg=BAR_ANGLE_DEG, S=1.0, bar_frac=1.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", required=True, help="mock file name in mock/")
    ap.add_argument("--galaxy", default="nbody", help="galaxy whose potential is used (default nbody)")
    ap.add_argument("--omega", type=float, default=None, help="|Omega_b| [km/s/kpc] (default: the truth)")
    ap.add_argument("--phi", type=float, default=None, help="trial bar angle [deg] (default: the truth, 27)")
    ap.add_argument("--S", type=float, default=None, help="bar length scale (default 1)")
    ap.add_argument("--bar-frac", type=float, default=None, help="bar strength (default 1)")
    ap.add_argument("--run", default="single_point", help="output sub-directory name")
    ap.add_argument("--tag", default=None, help="file tag (default: from the mock name)")
    ap.add_argument("--accuracy", type=float, default=ACCURACY)
    ap.add_argument("--n-boot", type=int, default=fs.N_BOOT)
    ap.add_argument("--obs-seed", type=int, default=11, help="seed for re-observing the library")
    ap.add_argument("--fresh", action="store_true", help="re-integrate even if the trajectory is saved")
    args = ap.parse_args()

    mock = MOCK_DIR / args.mock
    if not mock.exists():
        raise SystemExit(f"no such mock: {rel(mock)}")
    p = truth_params(args.galaxy)
    for key, val in (("Omega_abs", args.omega), ("phi_b_deg", args.phi), ("S", args.S), ("bar_frac", args.bar_frac)):
        if val is not None:
            p[key] = val
    tag = args.tag or mock.stem
    out_dir = RUNS_DIR / args.run
    res_dir = RESULTS_DIR / args.run
    out_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(mock, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    print(f"mock     {rel(mock)}  ({len(stars):,} stars)")
    print(f"galaxy   {args.galaxy}  ({rel(galaxy(args.galaxy)['pot_dir'])})")
    print(f"point    |Omega_b| = {p['Omega_abs']:.2f}  phi_b = {p['phi_b_deg']:.1f} deg  S = {p['S']:.3f}  "
          f"bar_frac = {p['bar_frac']:.3f}   (accuracy {args.accuracy:g})", flush=True)

    ssl.ACCURACY = args.accuracy
    t0 = time.time()
    pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], args.galaxy)
    t_pot = time.time() - t0
    print(f"potential built in {t_pot:.1f} s", flush=True)

    traj_file = out_dir / f"traj_{tag}.npy"
    t1 = time.time()
    if traj_file.exists() and not args.fresh:
        traj, times = np.load(traj_file), ssl.snapshot_times()
        t_int = time.time() - t1
        print(f"re-used {rel(traj_file)} ({traj.shape}) in {t_int:.1f} s "
              f"(--fresh to re-integrate)", flush=True)
    else:
        ic = ssl.rotate_z(stars, np.radians(p["phi_b_deg"]))
        times, traj = ssl.integrate(pot, ic, -abs(p["Omega_abs"]))
        t_int = time.time() - t1
        np.save(traj_file, traj)
        print(f"integrated {len(stars):,} orbits x {ssl.T_GYR} Gyr, {ssl.N_SNAP} snapshots in {t_int:.1f} s "
              f"({len(stars) / t_int:.0f} orbits/s); wrote {rel(traj_file)} "
              f"({traj_file.stat().st_size / 1e6:.0f} MB)", flush=True)

    t2 = time.time()
    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(len(stars), np.full(len(stars), 1 / len(stars)), size=args.n_boot).astype(float)
    cd, cov, n_stars, ok = fs.data_fit(stars, float(p["phi_b_deg"]), boot_w)
    t_data = time.time() - t2
    print(f"data Fourier fit + {args.n_boot} bootstraps in {t_data:.1f} s "
          f"({int(ok.sum())} of {fs.NB} blocks have >= {fs.N_MIN_BLOCK} stars)", flush=True)

    t3 = time.time()
    with h5py.File(mock, "r") as f:
        sigma_dm = float(f.attrs.get("sigma_DM", 0.0))
    libs = {"asis": fs.library_fit(traj, times, ok)}
    if sigma_dm > 0:        # the like-for-like comparison: smear the model the way the data are smeared
        libs["observed"] = forward_smeared_library_fit(traj, times, ok, sigma_dm, float(p["phi_b_deg"]),
                                                       seed=args.obs_seed)
    out = {}
    for v, cl in libs.items():
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        for name, sc in fs.SCHEMES.items():
            out[f"chi2_{name}_{v}"], out[f"ndata_{name}_{v}"] = fs.scheme_chi2(cd, cl, cov, valid, sc)
    t_lib = time.time() - t3
    t_total = time.time() - t0

    npz = out_dir / f"point_{tag}.npz"
    np.savez(npz, coef_data=cd, cov_data_joint=cov, n_stars=n_stars, valid=valid, sigma_DM=sigma_dm,
             **{f"coef_lib_{v}": c for v, c in libs.items()},
             ring_edges=fs.RING_EDGES, layer_edges=fs.LAYER_EDGES, terms=fs.TERMS, fields=fs.FIELD_NAMES,
             times=times, n_boot=args.n_boot, accuracy=args.accuracy,
             params=np.array([p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"]]), **out)
    print(f"library fit + chi2 in {t_lib:.1f} s; wrote {rel(npz)}", flush=True)

    print(f"\n{'scheme':52s} {'library':>9s} {'chi2':>11s} {'N_data':>7s} {'chi2/N':>8s}")
    for name in fs.SCHEMES:
        for v in libs:
            c, nd = out[f"chi2_{name}_{v}"], out[f"ndata_{name}_{v}"]
            star = "  <-- main" if name == MAIN_SCHEME else ""
            print(f"{name:52s} {v:>9s} {c:11.1f} {nd:7d} {c / nd:8.3f}{star}")
    print(f"\ntimings [s]: potential {t_pot:.1f}, integrate {t_int:.1f}, data fit {t_data:.1f}, "
          f"library + chi2 {t_lib:.1f}, total {t_total:.1f} ({t_total / 60:.1f} min)")

    row = dict(mock=args.mock, galaxy=args.galaxy, n_stars=len(stars), **p, accuracy=args.accuracy,
               n_boot=args.n_boot, n_blocks=int(valid.sum()), sigma_DM=sigma_dm,
               **{k: round(v, 4) if isinstance(v, float) else v for k, v in out.items()},
               **{f"chi2N_main_{v}": round(out[f"chi2_{MAIN_SCHEME}_{v}"] / out[f"ndata_{MAIN_SCHEME}_{v}"], 4)
                  for v in libs},
               t_potential_s=round(t_pot, 1), t_integrate_s=round(t_int, 1), t_data_fit_s=round(t_data, 1),
               t_library_chi2_s=round(t_lib, 1), t_total_s=round(t_total, 1), traj_file=str(rel(traj_file)))
    # append, tolerating a changing set of columns: a perfect mock has no "observed" library, so its
    # row is narrower than a distance-error one and a plain DictWriter would write a short header
    csv_path = res_dir / "single_point.csv"
    import pandas as pd
    df = pd.concat([pd.read_csv(csv_path), pd.DataFrame([row])], ignore_index=True) \
        if csv_path.exists() else pd.DataFrame([row])
    df.to_csv(csv_path, index=False)
    print(f"appended {rel(csv_path)} ({len(df)} rows)")


if __name__ == "__main__":
    main()

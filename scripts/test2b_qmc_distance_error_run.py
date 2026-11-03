"""Test 2b: the current method on data with distance errors (and optionally a selection function), over the
focused Sobol grid.

At every point of the focused global search (results/test1_nbody_qmc_focused/design.csv, 257 points) the observed
stars are rotated by the trial bar angle and integrated in the trial potential (accuracy 1e-5, 1 Gyr, 41 snapshots).
Library: samples t >= 0.25 Gyr, weighted
  asis      1/tau, tau = time fraction of the orbit at R < 8 kpc; library not smeared
  observed  the same, but every library sample observed again with the mock's distance error (seed OBS_SEED)
  sf        selection-function weighting S(x_j)/tau_i, tau_i = <S> along the orbit, S evaluated on the observed
            library (selection.py). Added automatically for a mock that carries selection tables, and the only
            one of the three that is right for such a mock; the other two ignore the selection and are kept for
            comparison.
Per-star Fourier fit (test1_nbody_qmc_focused_fourier_stars.py) of the data and of the library; chi2 for its
schemes; tempered logLs = -chi2 / (2 sqrt(2 N_data)), N_data at the truth.

The data fit is done **once** at the true bar angle and rotated to every trial angle instead of being refitted:
rings, layers, v_R and v_phi are unchanged by a rotation about z and the (cos, sin) pair of mode m turns by m*dphi,
so c' = T c and C' = T C T^T exactly (`test2b_qmc_7p5_rescore_coadded.rotation`, checked there to 2e-10 and
re-checked at start-up here). That removes a 300-bootstrap refit from every point (PIPELINE section 5.4).

Mock, galaxy and output directory are arguments (PIPELINE section 7.8), and --stride/--offset give each worker a
subset of the design so job arrays - or two machines sharing one output directory - can run it (section 7.7).
Points are claimed with an atomic lock file, so workers never duplicate work. Each worker writes only per-point
files; `--merge` (or a final run with nothing left to do) builds evals.csv from them.

Files: $MIRA_DR4_RUNS/<run>/{traj_NNN.npy, fourier_NNN.npz}; results/<run>/{evals.csv, run_meta.json}.
Resumable: points with fourier_NNN.npz are re-used.

Usage:
    python test2b_qmc_distance_error_run.py                                  # the original 7.5% run
    python test2b_qmc_distance_error_run.py --mock <file> --run <name> --galaxy nbody_t4
    python test2b_qmc_distance_error_run.py --run <name> --stride 12 --offset 3   # one worker of twelve
    python test2b_qmc_distance_error_run.py --run <name> --merge                  # build evals.csv
"""
import argparse
import json
import os
import shutil
import time

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR, RUNS_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import selection
from test1_nbody_qmc_focused_run import DESIGN_CSV
from plot_distance_error_effect_on_fourier_coefficients import stars, forward_smeared_library_fit
from test1_nbody_qmc_focused_fourier_stars import rotation

MOCK = MOCK_DIR / "nbody_mock_distance_error_sigmaDM0.163_R8_N50k.hdf5"
OUT_DIR = RUNS_DIR / "test2b_qmc_7p5"
RES_DIR = ROOT / "results" / "test2b_qmc_7p5"
EVALS = RES_DIR / "evals.csv"
OBS_SEED = 11
VARIANTS = ("observed", "asis")
OMEGA_REF = 26.02          # the pattern speed the focused design was built around (galaxy "nbody")
MIN_FREE_GB = 25.0         # below this, stop saving trajectories rather than fill the disc


def truth_omega(name):
    """Measured pattern speed of a galaxy: its ground-truth file, else the snapshot evolution table."""
    g = galaxy(name)
    if g["truth"].exists():
        return float(json.load(open(g["truth"]))["Omega_p_kms_kpc"])
    df = pd.read_csv(ROOT / "results" / "nbody_pattern_speed_evolution.csv")
    hit = df[np.isclose(df.t_model, g.get("time", np.nan))]
    if not len(hit):
        raise SystemExit(f"no truth pattern speed for {name}; run measure_nbody_pattern_speed_evolution.py")
    return float(hit.Omega_p.iloc[0])


def scaled_design(galaxy_name):
    """The focused design, with the Omega axis scaled to this galaxy's pattern speed.

    The 257 Sobol points stay exactly where they are in the unit cube, so every mock is searched on the same
    design; only the Omega box moves, from (20, 32) around 26.02 to the same relative range around the galaxy's
    own Omega. Returns the table and the scale factor."""
    d = pd.read_csv(DESIGN_CSV).sort_values("idx").reset_index(drop=True)
    scale = truth_omega(galaxy_name) / OMEGA_REF
    d["Omega_abs"] = d["Omega_abs"] * scale
    return d, scale


STALE_LOCK_H = 3.0        # a lock older than this, with no result, belonged to a worker that died


def claim(path):
    """Atomically claim a design point; False if another live worker holds it.

    Several workers - job-array tasks, or two machines sharing one output directory - compete for the same
    design, so the claim has to be atomic (O_EXCL). A lock older than STALE_LOCK_H with no .npz beside it is
    taken to be a crashed worker's and is re-claimed, otherwise a single crash would strand that point."""
    try:
        os.close(os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        try:
            if (time.time() - os.path.getmtime(path)) / 3600 > STALE_LOCK_H:
                os.utime(path, None)              # re-stamp first, so two reclaimers cannot both win
                print(f"reclaimed stale lock {os.path.basename(path)}", flush=True)
                return True
        except FileNotFoundError:
            return claim(path)
        return False


def merge_evals(design, out_dir, evals_path, variants):
    """Build evals.csv from whatever per-point fourier_NNN.npz exist (PIPELINE section 7.7)."""
    rows = []
    for f in sorted(out_dir.glob("fourier_[0-9][0-9][0-9].npz")):
        i = int(f.stem.split("_")[1])
        s = np.load(f)
        p = design.set_index("idx").loc[i]
        row = dict(idx=i, kind=p.kind, batch=p.batch, Omega_abs=p.Omega_abs, phi_b_deg=p.phi_b_deg,
                   S=p.S, bar_frac=p.bar_frac)
        row.update({k: float(s[k]) for k in s.files if k.startswith(("chi2_", "ndata_"))})
        rows.append(row)
    if not rows:
        print("nothing to merge")
        return None
    df = pd.DataFrame(rows).sort_values("idx")
    t = df.kind == "truth"
    if t.any():
        for name in fs.SCHEMES:
            for v in variants:
                col = f"ndata_{name}_{v}"
                if col not in df:
                    continue
                nd = int(df.loc[t, col].iloc[0])
                df[f"chi2N_{name}_{v}"] = df[f"chi2_{name}_{v}"] / df[col]
                df[f"logLs_{name}_{v}"] = -0.5 * df[f"chi2_{name}_{v}"] / np.sqrt(2 * nd)
    df.to_csv(evals_path, index=False)
    print(f"merged {len(df)} points -> {rel(evals_path)}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", default=MOCK.name, help="mock file name in mock/")
    ap.add_argument("--run", default=OUT_DIR.name, help="output directory name under $MIRA_DR4_RUNS and results/")
    ap.add_argument("--galaxy", default="nbody", help="galaxy whose potential and truth Omega are used")
    ap.add_argument("--stride", type=int, default=1, help="run every stride-th design point")
    ap.add_argument("--offset", type=int, default=0, help="this worker's residue, 0 .. stride-1")
    ap.add_argument("--no-traj", action="store_true", help="do not save trajectories")
    ap.add_argument("--time-budget-h", type=float, default=None, help="stop cleanly after this many hours")
    ap.add_argument("--merge", action="store_true", help="only rebuild evals.csv from the per-point files")
    ap.add_argument("--cadence-step", type=int, default=1, metavar="K",
                    help="keep every K-th snapshot of the trajectory (2 = 50 Myr instead of 25). The "
                         "library is scored on the thinned trajectory AND only the thinned trajectory is "
                         "saved, so the stored file and every later rescoring are on the same grid. "
                         "Halves the trajectory files; measured cost on the posterior is < 0.2 km/s/kpc "
                         "(scripts/test_time_cadence.py).")
    args = ap.parse_args()

    mock = MOCK_DIR / args.mock
    out_dir, res_dir = RUNS_DIR / args.run, ROOT / "results" / args.run
    evals = res_dir / "evals.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)
    design, scale = scaled_design(args.galaxy)

    if args.merge:
        tab0 = selection.load_tables(mock)
        with h5py.File(mock, "r") as f:
            sig = float(f.attrs.get("sigma_DM", 0.0))
        merge_evals(design, out_dir, evals, list(VARIANTS) + (["sf"] if tab0 else []))
        return

    if not mock.exists():
        raise SystemExit(f"no such mock: {rel(mock)}")
    w_mw = stars(mock)
    with h5py.File(mock, "r") as f:
        sigma_dm = float(f.attrs.get("sigma_DM", 0.0))
    tab = selection.load_tables(mock)
    variants = list(VARIANTS) + (["sf"] if tab else [])
    omega_truth = truth_omega(args.galaxy)
    print(f"mock     {rel(mock)}  ({len(w_mw):,} stars, sigma_DM {sigma_dm:.4f}"
          f"{', with selection function' if tab else ''})")
    print(f"galaxy   {args.galaxy}: truth |Omega_b| = {omega_truth:.2f} (design Omega scaled by {scale:.4f}, "
          f"box {design.Omega_abs.min():.1f}-{design.Omega_abs.max():.1f})")
    print(f"library  variants: {', '.join(variants)}")
    print(f"worker   stride {args.stride} offset {args.offset}; output {rel(out_dir)}", flush=True)
    with open(res_dir / "run_meta.json", "w") as fh:
        json.dump(dict(run=args.run, mock=args.mock, galaxy=args.galaxy, omega_truth=omega_truth,
                       omega_scale=scale, variants=variants, sigma_DM=sigma_dm, accuracy=old.ACCURACY,
                       design=str(rel(DESIGN_CSV)), n_points=len(design),
                       cadence_step=args.cadence_step,
                       snapshot_dt_myr=round(float(ssl.T_GYR / (ssl.N_SNAP - 1) * args.cadence_step * 1000), 3)),
                  fh, indent=2)

    ssl.ACCURACY = old.ACCURACY
    times_full = ssl.snapshot_times()
    times = times_full[::args.cadence_step]
    if args.cadence_step != 1:
        print(f"cadence  keeping every {args.cadence_step} snapshots: {len(times)} of {len(times_full)} "
              f"({(times[1] - times[0]) * 1000:.0f} Myr), library keeps "
              f"{int((times >= ssl.T_MIN_GYR - 1e-9).sum())} of "
              f"{int((times_full >= ssl.T_MIN_GYR - 1e-9).sum())}", flush=True)
    boot_w = np.random.default_rng(fs.SEED).multinomial(len(w_mw), np.full(len(w_mw), 1 / len(w_mw)),
                                                        size=fs.N_BOOT).astype(float)
    t0 = time.time()
    cd0, cov0, n_stars, ok = fs.data_fit(w_mw, BAR_ANGLE_DEG, boot_w)
    cd0f = cd0.reshape(fs.NB, -1)
    print(f"data fit + {fs.N_BOOT} bootstraps at phi_b = {BAR_ANGLE_DEG}deg in {time.time() - t0:.1f} s "
          f"({int(ok.sum())} usable blocks)", flush=True)
    chk = 33.0                                     # verify the rotation identity before relying on it
    cd_direct = fs.data_fit(w_mw, chk, boot_w[:2])[0].reshape(fs.NB, -1)
    err = float(np.nanmax(np.abs(cd_direct[ok] - cd0f[ok] @ rotation(chk - BAR_ANGLE_DEG).T)))
    print(f"rotation check at {chk}deg: max |direct - rotated| = {err:.2e}", flush=True)
    if err > 1e-6:
        raise SystemExit("rotation identity failed; not safe to re-use the data fit")

    todo = [int(i) for i in design.idx if int(i) % args.stride == args.offset % args.stride]
    t_start, done = time.time(), 0
    for i in todo:
        saved = out_dir / f"fourier_{i:03d}.npz"
        if saved.exists():
            continue
        lock = out_dir / f"fourier_{i:03d}.lock"
        if not claim(str(lock)):
            continue
        if args.time_budget_h and (time.time() - t_start) / 3600 > args.time_budget_h:
            os.unlink(lock)
            print(f"time budget reached after {done} points", flush=True)
            break
        p = design.set_index("idx").loc[i]
        t1 = time.time()
        path = out_dir / f"traj_{i:03d}.npy"
        if path.exists():
            traj, t_int = np.load(path), 0.0
            if traj.shape[1] == len(times_full) and args.cadence_step != 1:
                traj = np.ascontiguousarray(traj[::1, ::args.cadence_step])   # full file, thin on read
            elif traj.shape[1] != len(times):
                raise SystemExit(f"{rel(path)} has {traj.shape[1]} snapshots; expected "
                                 f"{len(times)} (step {args.cadence_step}) or {len(times_full)}")
        else:
            pot = ssl.bar_aligned_potential(p.bar_frac, p.S, args.galaxy)
            _, traj = ssl.integrate(pot, ssl.rotate_z(w_mw, np.radians(p.phi_b_deg)), -abs(p.Omega_abs))
            t_int = time.time() - t1
            if args.cadence_step != 1:              # thin before scoring AND before saving, so the
                traj = np.ascontiguousarray(traj[:, ::args.cadence_step])   # stored file matches the fit
            if not args.no_traj and shutil.disk_usage(out_dir).free / 2**30 > MIN_FREE_GB:
                np.save(path, traj)
        t2 = time.time()
        T = rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        cd = (cd0f @ T.T).reshape(fs.NB, 6, fs.NK)
        cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
        libs, taus = {}, {}
        for v in variants:
            if v == "asis":
                libs[v] = fs.library_fit(traj, times, ok)
            elif v == "observed":
                libs[v] = forward_smeared_library_fit(traj, times, ok, sigma_dm, float(p.phi_b_deg),
                                                      seed=OBS_SEED)
            else:
                libs[v], taus[v] = selection.library_fit_sf(traj, times, ok, tab, float(p.phi_b_deg),
                                                            seed=OBS_SEED, fs=fs)
        out = {}
        for v, cl in libs.items():
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            for name, sc in fs.SCHEMES.items():
                out[f"chi2_{name}_{v}"], out[f"ndata_{name}_{v}"] = fs.scheme_chi2(cd, cl, cov, valid, sc)
        t_score = time.time() - t2
        np.savez(saved, coef_data=cd, cov_data_joint=cov, n_stars=n_stars,
                 params=np.array([p.Omega_abs, p.phi_b_deg, p.S, p.bar_frac]), kind=str(p.kind),
                 t_integrate_s=t_int, t_score_s=t_score,
                 **{f"coef_lib_{v}": c for v, c in libs.items()},
                 **{f"tau_{v}": a.astype(np.float32) for v, a in taus.items()}, **out)
        del traj
        done += 1
        key = f"fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024_{variants[-1]}"
        free = shutil.disk_usage(out_dir).free / 2**30
        print(f"[{done}/{len(todo)}] idx {i:3d} {p.kind:5s} Om={p.Omega_abs:.2f} phi={p.phi_b_deg:.1f} "
              f"S={p.S:.3f} f={p.bar_frac:.3f}  chi2/N({variants[-1]}) "
              f"{out['chi2_' + key] / out['ndata_' + key]:.2f}  int {t_int:.0f}s score {t_score:.0f}s  "
              f"(elapsed {(time.time() - t_start) / 60:.0f} min, {free:.0f} GB free)", flush=True)
    merge_evals(design, out_dir, evals, variants)
    print(f"worker done: {done} points in {(time.time() - t_start) / 3600:.2f} h", flush=True)


if __name__ == "__main__":
    main()

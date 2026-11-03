"""Does a 50 Myr trajectory cadence cost anything? Rescores a finished search on every other snapshot.

The runs save 41 snapshots over 1 Gyr (25 Myr); the library keeps the 30 with t >= 0.25 Gyr. Storing every
*other* snapshot would halve the trajectory files, which are 98% of a run's 12.3 GB. The question is whether
the posterior notices.

Halving the cadence halves the number of library samples per orbit (30 -> 16), so the library's own sampling
noise grows by about sqrt(2). PROGRESS 2026-11-01 measured that noise at C_lib/C_data ~ 0.09 and showed it
accounts for ~37% of the excess chi2/N, so doubling it is not obviously harmless - hence this test rather
than an argument.

Nothing is re-integrated: the saved trajectories are re-scored with traj[:, ::step] and times[::step], which
is exactly what a coarser run would have produced (agama samples the same grid; the kept subset is a subset).
Writes `cad<step>_NNN.npz` beside the run's files and the usual chi2_/ndata_/logLs_ columns under the variant
name `cad<step>`, so test1_gp_emulator.py reads it directly.

Usage: python test_time_cadence.py --run test2_t7_sf_12p5pct [--step 2] [--stride N --offset K] [--merge]
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import h5py

from common import ROOT, MOCK_DIR, RUNS_DIR, BAR_ANGLE_DEG, rel
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import selection
from test1_nbody_qmc_focused_fourier_stars import rotation
from test2b_qmc_distance_error_run import scaled_design, claim
from plot_distance_error_effect_on_fourier_coefficients import forward_smeared_library_fit

OBS_SEED = 11


def merge(design, out_dir, path, variant, res_dir):
    rows = []
    for f in sorted(out_dir.glob(f"{variant}_[0-9][0-9][0-9].npz")):
        i = int(f.stem.split("_")[1])
        s = np.load(f)
        p = design.set_index("idx").loc[i]
        row = dict(idx=i, kind=p.kind, batch=p.batch, Omega_abs=p.Omega_abs, phi_b_deg=p.phi_b_deg,
                   S=p.S, bar_frac=p.bar_frac)
        row.update({k: float(s[k]) for k in s.files if k.startswith(("chi2_", "ndata_"))})
        rows.append(row)
    if not rows:
        print("nothing to merge")
        return
    df = pd.DataFrame(rows).sort_values("idx")
    t = df.kind == "truth"
    if t.any():                       # one set of columns per (scheme, library variant)
        for c in [c for c in df.columns if c.startswith("chi2_")]:
            tail = c[len("chi2_"):]
            nd = int(df.loc[t, f"ndata_{tail}"].iloc[0])
            df[f"chi2N_{tail}"] = df[c] / df[f"ndata_{tail}"]
            df[f"logLs_{tail}"] = -0.5 * df[c] / np.sqrt(2 * nd)
    df.to_csv(path, index=False)
    print(f"merged {len(df)} points -> {rel(path)}")
    ev_p = res_dir / "evals.csv"
    if ev_p.exists():                       # fold the new columns in, which is what the emulator reads
        ev = pd.read_csv(ev_p)
        new = [c for c in df.columns if f"_{variant}_" in c]
        ev = ev.drop(columns=[c for c in new if c in ev], errors="ignore")
        ev = ev.merge(df[["idx"] + new], on="idx", how="left")
        ev.to_csv(ev_p, index=False)
        print(f"  folded {len(new)} {variant} columns into {rel(ev_p)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--step", type=int, default=2, help="keep every STEP-th snapshot (2 = 50 Myr)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--merge", action="store_true")
    args = ap.parse_args()
    variant = f"cad{args.step}"

    res_dir, out_dir = ROOT / "results" / args.run, RUNS_DIR / args.run
    meta = json.load(open(res_dir / "run_meta.json"))
    design, _ = scaled_design(meta["galaxy"])
    evals = res_dir / f"evals_{variant}.csv"
    if args.merge:
        merge(design, out_dir, evals, variant, res_dir)
        return

    ssl.ACCURACY = old.ACCURACY
    mock = MOCK_DIR / meta["mock"]
    tab = selection.load_tables(mock)
    sigma_dm = float(meta["sigma_DM"])
    stars = np.column_stack([h5py.File(mock, "r")[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    boot_w = np.random.default_rng(fs.SEED).multinomial(len(stars), np.full(len(stars), 1 / len(stars)),
                                                        size=fs.N_BOOT).astype(float)
    cd0, cov0, n_stars, ok = fs.data_fit(stars, BAR_ANGLE_DEG, boot_w)
    cd0f = cd0.reshape(fs.NB, -1)
    times_full = ssl.snapshot_times()
    times = times_full[::args.step]
    keep_full = int((times_full >= ssl.T_MIN_GYR - 1e-9).sum())
    keep_cut = int((times >= ssl.T_MIN_GYR - 1e-9).sum())
    print(f"{args.run}: step {args.step} -> {len(times)} of {len(times_full)} snapshots "
          f"({(times[1] - times[0]) * 1000:.0f} Myr), library keeps {keep_cut} of {keep_full}", flush=True)

    todo = [int(i) for i in design.idx if int(i) % args.stride == args.offset % args.stride]
    t0, done = time.time(), 0
    for i in todo:
        out = out_dir / f"{variant}_{i:03d}.npz"
        traj_f = out_dir / f"traj_{i:03d}.npy"
        if out.exists() or not traj_f.exists():
            continue
        if not claim(str(out_dir / f"{variant}_{i:03d}.lock")):
            continue
        p = design.set_index("idx").loc[i]
        traj = np.ascontiguousarray(np.load(traj_f)[:, ::args.step])
        libs = {}
        if tab:
            libs["sf"] = selection.library_fit_sf(traj, times, ok, tab, float(p.phi_b_deg),
                                                  seed=OBS_SEED, fs=fs)[0]
        libs["asis"] = fs.library_fit(traj, times, ok)
        libs["observed"] = forward_smeared_library_fit(traj, times, ok, sigma_dm, float(p.phi_b_deg),
                                                       seed=OBS_SEED)
        T = rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        cd = (cd0f @ T.T).reshape(fs.NB, 6, fs.NK)
        cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
        res = {}
        for v, cl in libs.items():
            valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
            for name, sc in fs.SCHEMES.items():
                res[f"chi2_{name}_{variant}_{v}"], res[f"ndata_{name}_{variant}_{v}"] = \
                    fs.scheme_chi2(cd, cl, cov, valid, sc)
        np.savez(out, step=args.step, n_snap=len(times), n_keep=keep_cut,
                 params=np.array([p.Omega_abs, p.phi_b_deg, p.S, p.bar_frac]), kind=str(p.kind),
                 **{f"coef_lib_{v}": c for v, c in libs.items()}, **res)
        del traj
        done += 1
        if done % 20 == 0 or done == 1:
            print(f"  [{done}/{len(todo)}] idx {i:3d}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
    merge(design, out_dir, evals, variant, res_dir)
    print(f"done: {done} points in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()

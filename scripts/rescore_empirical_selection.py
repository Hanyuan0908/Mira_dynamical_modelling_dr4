"""Rescore a finished search with an *empirical* selection function, estimated from the data itself.

The `sf` weighting of test2b_qmc_distance_error_run.py uses the selection tables stored in the mock, i.e.
perfect knowledge of S, which no real survey has. Here S is instead read off the ratio of the observed sample
to the (re-observed) orbit library in (l, b, d, R): if the trial potential is right, that is the only reason
the two differ. It is fitted in the same separable form the mock uses, by weighted iterative proportional
fitting against the data's three margins, with the library as the base population (`selection.empirical_tables`).

This costs no discriminating power in principle - the likelihood fits <v>(phi) by least squares and so does
not use the density, and reweighting by a function of position cannot change the library's mean velocity at a
given position. The risk is that S, refitted for every trial potential, absorbs genuine model error.

Settings come from a scan at the truth point (PROGRESS): HEALPix nside 32 and a *single* IPF pass. Iterating
to self-consistency in the 1/tau feedback makes S_emp worse, not better (correlation with the true S falls
from 0.73 to 0.55), and coarser pixels are worse too.

Needs only the saved trajectories - no new integrations. Writes `sfemp_NNN.npz` beside the run's own files and
`results/<run>/evals_sfemp.csv`, with the same chi2_/ndata_/logLs_ columns under the variant name `sfemp`,
so test1_gp_emulator.py can read it directly.

Usage: python rescore_empirical_selection.py --run test2_t7_sf_10pct [--stride N --offset K] [--merge]
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
import selection
from test2b_qmc_distance_error_run import load_design, reference_ndata, claim
from test1_nbody_qmc_focused_fourier_stars import rotation

NSIDE = 32
VARIANT = "sfemp"
ESTIMATORS = {"ipf": "sfemp", "ratio": "sfratio"}


def data_lbdR(mock):
    with h5py.File(mock, "r") as f:
        return (f["l"][:], f["b"][:], f["d"][:], np.hypot(f["x"][:], f["y"][:]))


def merge(design, out_dir, path):
    rows = []
    for f in sorted(out_dir.glob(f"{VARIANT}_[0-9][0-9][0-9].npz")):
        i = int(f.stem.rsplit("_", 1)[1])
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
    for name in fs.SCHEMES:
        nd = reference_ndata(df, f"ndata_{name}_{VARIANT}")
        df[f"chi2N_{name}_{VARIANT}"] = df[f"chi2_{name}_{VARIANT}"] / df[f"ndata_{name}_{VARIANT}"]
        df[f"logLs_{name}_{VARIANT}"] = -0.5 * df[f"chi2_{name}_{VARIANT}"] / np.sqrt(2 * nd)
    df.to_csv(path, index=False)
    print(f"merged {len(df)} points -> {rel(path)}")
    # also fold the new columns into the run's own evals.csv, which is what test1_gp_emulator reads
    main = path.parent / "evals.csv"
    if not main.exists():          # a --keep sub-sample's own results directory: the scored table is the evals
        df.to_csv(main, index=False)
        print(f"  wrote {rel(main)}")
    else:
        base = pd.read_csv(main)
        new_cols = [c for c in df.columns if c.endswith(f"_{VARIANT}")]
        base = base.drop(columns=[c for c in new_cols if c in base.columns])
        base = base.merge(df[["idx"] + new_cols], on="idx", how="left")
        base.to_csv(main, index=False)
        print(f"  folded {len(new_cols)} {VARIANT} columns into {rel(main)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--nside", type=int, default=NSIDE)
    ap.add_argument("--dbin", type=float, default=1.0, help="distance bin [kpc] for --estimator ratio")
    ap.add_argument("--estimator", default="ipf", choices=list(ESTIMATORS),
                    help="ipf: separable S fitted to the data margins; ratio: direct rho_obs/rho_lib in (l,b,d)")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--cadence-step", type=int, default=1, metavar="K",
                    help="thin the saved trajectory to every K-th snapshot before scoring (2 = 50 Myr "
                         "from a 25 Myr file). Runs made after 2026-11-03 are already saved thinned, "
                         "so this is only for re-scoring older ones on the same grid.")
    ap.add_argument("--n-draw", type=int, default=10, metavar="K",
                    help="distance-error draws per library sample, pooled into one fit. DEFAULT 10 "
                         "since 2026-11-03: a one-sample library (K=1) biases chi2 high by a "
                         "parameter-dependent amount, and K=10 cuts the emulator's leave-one-out "
                         "error 1.47 -> 0.65. Pass --n-draw 1 to reproduce results from before that "
                         "date. See selection.observed_library.")
    ap.add_argument("--ring-width", type=float, default=0.3,
                    help="Fourier ring width inside 6 kpc (default 0.3, the validated grid)")
    ap.add_argument("--z-max", type=float, default=3.0,
                    help="top of the |z| layers, 1 kpc each (default 3)")
    ap.add_argument("--keep", default=None, metavar="NPY",
                    help="boolean mask over the mock's stars: score only these, dropping the same rows from the "
                         "data and from every saved trajectory. Orbits are test particles, so this is exactly a "
                         "search on the sub-sample, with no new integration. Needs --keep-tag.")
    ap.add_argument("--keep-tag", default=None,
                    help="name of the sub-sample; results go to results/<run>_<tag>/ under variant sfratio_<tag>")
    ap.add_argument("--sigma-dm", type=float, default=None, metavar="MAG",
                    help="smear the library with this distance-modulus error instead of the run's own (e.g. to "
                         "ask what the data imply if their distances are worse than assumed); variant sfratio_sdm<MAG>")
    args = ap.parse_args()
    if bool(args.keep) != bool(args.keep_tag):
        raise SystemExit("--keep and --keep-tag go together")
    global VARIANT
    # a non-default grid gets its own variant name (e.g. sfratio_r0.5_z2), so default results are never overwritten
    VARIANT = ESTIMATORS[args.estimator] + fs.set_grid(args.ring_width, args.z_max) + \
        (f"_{args.keep_tag}" if args.keep else "") + (f"_sdm{args.sigma_dm:g}" if args.sigma_dm else "")

    res_dir = ROOT / "results" / args.run
    out_dir = RUNS_DIR / args.run
    meta = json.load(open(res_dir / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    keep = None
    if args.keep:
        keep = np.load(ROOT / args.keep)
        sub_dir = ROOT / "results" / f"{args.run}_{args.keep_tag}"
        sub_dir.mkdir(parents=True, exist_ok=True)
        json.dump(dict(meta, parent_run=args.run, keep_mask=args.keep, n_keep=int(keep.sum()), n_mock=len(keep)),
                  open(sub_dir / "run_meta.json", "w"), indent=2)
        res_dir = sub_dir
    evals = res_dir / f"evals_{VARIANT}.csv"
    if args.merge:
        merge(design, out_dir, evals)
        return

    mock = MOCK_DIR / meta["mock"]
    dat = data_lbdR(mock)
    sigma_dm = float(args.sigma_dm or meta["sigma_DM"])
    stars = np.column_stack([h5py.File(mock, "r")[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    if keep is not None:
        if len(keep) != len(stars):
            raise SystemExit(f"--keep has {len(keep)} entries for {len(stars)} mock stars")
        dat = tuple(np.asarray(a)[keep] for a in dat)
        stars = stars[keep]
        print(f"keeping {int(keep.sum()):,} of {len(keep):,} stars ({args.keep})", flush=True)
    boot_w = np.random.default_rng(fs.SEED).multinomial(len(stars), np.full(len(stars), 1 / len(stars)),
                                                        size=fs.N_BOOT).astype(float)
    cd0, cov0, n_stars, ok = fs.data_fit(stars, BAR_ANGLE_DEG, boot_w)
    cd0f = cd0.reshape(fs.NB, -1)
    print(f"{args.run}: {len(stars):,} stars, sigma_DM {sigma_dm:.4f}, nside {args.nside}", flush=True)

    todo = [int(i) for i in design.idx if int(i) % args.stride == args.offset % args.stride]
    t0, done = time.time(), 0
    for i in todo:
        out = out_dir / f"{VARIANT}_{i:03d}.npz"
        traj_f = out_dir / f"traj_{i:03d}.npy"
        if out.exists() or not traj_f.exists():
            continue
        if not claim(str(out_dir / f"{VARIANT}_{i:03d}.lock")):
            continue
        p = design.set_index("idx").loc[i]
        traj = np.load(traj_f)
        if keep is not None:
            traj = traj[keep]                   # row i of a trajectory file is star i of the mock
        if args.cadence_step > 1:
            traj = np.ascontiguousarray(traj[:, ::args.cadence_step])
        times = ssl.times_for_traj(traj)        # the file may already be thinned
        K = max(1, args.n_draw)
        _obs_bar, lbdR, nk = selection.observed_library(traj, times, float(p.phi_b_deg), sigma_dm,
                                                        seed=11, n_draw=K)
        R0 = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])
        tau0 = np.clip(np.mean(R0 < ssl.R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
        w0 = np.tile(np.repeat(1.0 / tau0, nk), K) / K
        if args.estimator == "ratio":
            tab = selection.ratio_tables(dat, lbdR, w0, nside=args.nside, d_bin=args.dbin)
            Se = selection.evaluate_ratio(tab, *lbdR)
        else:
            tab = selection.empirical_tables(dat, lbdR, w0, nside=args.nside)
            Se = selection.evaluate_raw(tab, *lbdR)
        obs_bar = _obs_bar
        tau = np.clip(Se.reshape(K, traj.shape[0], nk).mean(axis=(0, 2)), selection.TAU_FLOOR, None)
        wt = Se / np.tile(np.repeat(tau, nk), K) / K
        b, X, v, m = fs.prepare(obs_bar)
        cl = fs.fit(b, X, v, wt[m], ok, leverage=False)
        T = rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        cd = (cd0f @ T.T).reshape(fs.NB, 6, fs.NK)
        cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        res = {}
        for name, sc in fs.SCHEMES.items():
            res[f"chi2_{name}_{VARIANT}"], res[f"ndata_{name}_{VARIANT}"] = fs.scheme_chi2(cd, cl, cov, valid, sc)
        np.savez(out, coef_lib=cl, tau=tau.astype(np.float32), nside=args.nside, n_draw=K,
                 estimator=args.estimator, dbin=args.dbin,
                 params=np.array([p.Omega_abs, p.phi_b_deg, p.S, p.bar_frac]), kind=str(p.kind), **res)
        del traj
        done += 1
        if done % 10 == 0 or done == 1:
            print(f"  [{done}/{len(todo)}] idx {i:3d}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
    merge(design, out_dir, evals)
    print(f"done: {done} points in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()

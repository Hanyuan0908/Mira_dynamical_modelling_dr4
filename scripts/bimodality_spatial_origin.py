"""Is the low-|Omega_b| / weak-bar part of a Hunter+24 posterior real, and which part of the sample produces it?
(user, 2026-11-25)

Stage "points" (emulator, no integration). The run's chi2 on the default grid, at that grid's T1 (logL = -dchi2/2T1),
emulated as the pipeline does. Its maximum is found separately in the window of the second mode (|Omega_b| < --split
with --second low, > --split with --second high) and in the rest (the main mode). A second-window maximum that sits ON
the window edge is the slope of the main mode (a tail); one inside the window is a separate mode. Printed with the
chain fraction in each window. (The first run, h24mock, labelled the second mode "low".)

Stage "direct". At each point (the two window maxima, and the truth for a mock), N_REP runs are integrated: the point
itself plus N_REP - 1 tiny random perturbations (the widths of chi2_jitter_small_perturbations.py). Each is scored
exactly as the pipeline (data fitted at 27 deg on that grid with the 300-bootstrap joint covariance; library with the
sfratio selection, K = 10, rotated to 27 deg; scheme B; Hartlap) on every grid in GRIDS. Per block, chi2 is split
EXACTLY into terms, chi2_block = sum_t h d_t (C^-1 d)_t (individual terms can be negative), and summed into the
four fields. Control: the per-block sums reproduce fs.scheme_chi2.

Delta(block) = mean over runs of chi2_block(point B) - chi2_block(point A), with its standard error from the run
scatter, on the blocks valid in every run. Negative Delta: that block's data prefer B.

Usage: python bimodality_spatial_origin.py --case h24mock --split 28 [--n-rep 12]
       python bimodality_spatial_origin.py --case osarg_s22 --split 34.5 --second high
Output: results/bimodality_origin_<case>.npz and .json (points, totals, per-block Delta), printed tables.
        Trajectories: $MIRA_DR4_RUNS (default ~/mira_dr4_runs)/bimodality_traj_<case>/.
"""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
from scipy.optimize import minimize

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import test1_gp_emulator as ge
from compare_grid_strategies import posterior
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from scan_omega_at_fixed_bar import CASES
from chi2_jitter_small_perturbations import RANDOM_WIDTH

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
GRIDS = [((0.3, 3.0, 1.0, None), "default: 0.3 kpc rings, 3 layers"),
         ((0.3, 2.0, 2.0, "geo20"), "geometric 20 rings, |z| < 2"),
         ((0.3, 2.0, 1.0, "geo20"), "geometric 20 rings, |z| 0-1, 1-2")]
FIELDS = {0: "<v_R> m2,4", 1: "<v_phi> m0,2,4", 3: "sigma_R^2 m2,4", 4: "sigma_phi^2 m2,4"}
TRUTH = {"h24mock": (32.5, 25.0, 1.0, 1.0)}
R_BANDS = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 6)]


def t1_default(case):
    """T1 of the default grid: from chi2_jitter_per_grid_<case>.csv if it exists, else the 100-run random file."""
    p = RESULTS_DIR / f"chi2_jitter_per_grid_{case}.csv"
    if p.exists():
        j = pd.read_csv(p)
        j = j[j.grid == "default"]
        if len(j) > 2:
            return float(np.std(j.chi2.values[1:] - j.chi2.iloc[0], ddof=1)), p.name
    j = pd.read_csv(RESULTS_DIR / f"chi2_jitter_{case}_random.csv")
    x = j.loc[j.param == "random", "C_data"].values - j.loc[j.param == "repeat", "C_data"].iloc[0]
    return float(np.std(x, ddof=1)), f"chi2_jitter_{case}_random.csv"


def find_points(run, case, split, second="low"):
    ge.use_run(run)
    ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    chi2 = tr[f"chi2_{B}_sfratio"].values
    T1, src = t1_default(case)
    u, y = ge.to_unit(tr), -(chi2 - chi2.min()) / (2 * T1)
    emu = ge.Emulator(u, y)
    us = (split - ge.LO[0]) / (ge.HI[0] - ge.LO[0])
    out = {}
    wins = (("second", 0.0, us), ("main", us, 1.0)) if second == "low" else (("main", 0.0, us), ("second", us, 1.0))
    for lab, lo, hi in wins:
        best = None
        inside = (u[:, 0] >= lo) & (u[:, 0] <= hi)
        starts = np.vstack([u[inside][np.argsort(-y[inside])[:10]],
                            np.column_stack([np.random.default_rng(3).uniform(lo, hi, 30), np.random.default_rng(4).random((30, 3))])])
        for x0 in starts:
            r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(lo, hi)] + [(0, 1)] * 3, method="L-BFGS-B")
            if best is None or r.fun < best.fun:
                best = r
        p = ge.from_unit(best.x[None])[0]
        edge = min(abs(best.x[0] - lo), abs(best.x[0] - hi)) < 0.01
        out[lab] = dict(Omega=float(p[0]), phi=float(p[1]), S=float(p[2]), f=float(p[3]), logL=float(-best.fun),
                        on_window_edge=bool(edge))
    P = posterior(u, y)
    m = P[:, 0] < split if second == "low" else P[:, 0] > split
    info = dict(T1=T1, T1_source=src, split=split, second_side=second, chain_frac_second=float(m.mean()),
                chain_second=dict(S=float(np.median(P[m, 2])), f=float(np.median(P[m, 3]))),
                chain_main=dict(S=float(np.median(P[~m, 2])), f=float(np.median(P[~m, 3]))))
    return out, info


def block_terms(cd, cl, cov, valid, spec):
    """Per block and field: the exact additive split of the joint-covariance scheme chi2 (NaN for invalid blocks)."""
    idx_all = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    field_of = np.array([f for f, terms in spec.items() for _ in terms])
    out = np.full((fs.NB, len(spec)), np.nan)
    h = (fs.N_BOOT - len(idx_all) - 2) / (fs.N_BOOT - 1)
    for r in range(spec_rings(spec)):
        for z in range(fs.NL):
            i = r * fs.NL + z
            if not valid[i]:
                continue
            d = (cd[i] - cl[i]).reshape(-1)[idx_all]
            t = h * d * np.linalg.solve(cov[i][np.ix_(idx_all, idx_all)], d)
            out[i] = [t[field_of == f].sum() for f in spec]
    return out


def spec_rings(spec):
    return fs.SCHEMES[B][0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, choices=list(CASES))
    ap.add_argument("--split", type=float, required=True, help="|Omega_b| separating the low window from the rest")
    ap.add_argument("--n-rep", type=int, default=12)
    ap.add_argument("--second", choices=["low", "high"], default="low", help="side of --split the second mode lies on")
    args = ap.parse_args()
    c = CASES[args.case]
    run = c["run"]
    pts, info = find_points(run, args.case, args.split, args.second)
    side = "below" if args.second == "low" else "above"
    print(f"[{args.case}] default grid, T1 = {info['T1']:.1f} ({info['T1_source']}); chain {side} {args.split}: "
          f"{info['chain_frac_second']:.2f} (S {info['chain_second']['S']:.2f}, f {info['chain_second']['f']:.2f}; rest S "
          f"{info['chain_main']['S']:.2f}, f {info['chain_main']['f']:.2f})")
    for k, p in pts.items():
        print(f"   emulator max, {k:6s} window: Omega {p['Omega']:.2f} phi {p['phi']:.2f} S {p['S']:.3f} f {p['f']:.3f}  "
              f"logL {p['logL']:+.2f}  {'ON THE WINDOW EDGE (a tail, not a separate maximum)' if p['on_window_edge'] else 'interior (a separate maximum)'}",
              flush=True)
    points = {"main": pts["main"], "second": pts["second"]}
    if args.case in TRUTH:
        o, ph, s, f = TRUTH[args.case]
        points["truth"] = dict(Omega=o, phi=ph, S=s, f=f)

    meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as fh:
        w = np.column_stack([fh[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    tdir = Path(os.environ.get("MIRA_DR4_RUNS", Path.home() / "mira_dr4_runs")) / f"bimodality_traj_{args.case}"
    tdir.mkdir(parents=True, exist_ok=True)
    ssl.ACCURACY = old.ACCURACY
    runs = {}
    for lab, p in points.items():
        rng = np.random.default_rng({"main": 101, "second": 202, "low": 202, "truth": 303}[lab])
        runs[lab] = []
        for r in range(args.n_rep):
            q = dict(Omega_abs=p["Omega"], phi_b_deg=p["phi"], S=p["S"], bar_frac=p["f"])
            if r:
                for k, wdt in RANDOM_WIDTH.items():
                    q[k] = q[k] + wdt * rng.standard_normal()
            path = tdir / f"{lab}_r{r:02d}.npy"
            if not path.exists():
                t0 = time.time()
                pot = ssl.bar_aligned_potential(float(q["bar_frac"]), float(q["S"]), c["galaxy"])
                _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(float(q["phi_b_deg"]))), -abs(float(q["Omega_abs"])))
                np.save(path, np.ascontiguousarray(traj[:, ::2]))
                print(f"integrated {lab} run {r} ({time.time() - t0:.0f} s)", flush=True)
            runs[lab].append((float(q["phi_b_deg"]), path))

    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    spec = fs.SCHEMES[B][1]
    res = dict(points=points, emulator_info=info, grids={})
    arrays = {}
    for gargs, gdesc in GRIDS:
        tag = fs.set_grid(*gargs) or "default"
        cd0, cov0, nstar, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
        cd0 = cd0.reshape(fs.NB, 6, fs.NK)
        per = {}
        for lab, lst in runs.items():
            per[lab] = []
            for phi, path in lst:
                cl, _ = library_fit_empirical(np.load(path), phi, sigma_dm, dat, ok)
                cl = (cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - phi).T).reshape(fs.NB, 6, fs.NK)
                valid = ok & np.all(np.isfinite(cd0), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
                bt = block_terms(cd0, cl, cov0, valid, spec)
                tot, _ = fs.scheme_chi2(cd0, cl, cov0, valid, fs.SCHEMES[B])
                assert abs(np.nansum(bt) - tot) < 1e-6 * max(1.0, tot), (np.nansum(bt), tot)
                per[lab].append(bt)
            per[lab] = np.array(per[lab])                   # (n_rep, NB, 4)
        common = np.all([np.all(np.isfinite(a), axis=(0, 2)) for a in per.values()], axis=0)
        edges = fs.RING_EDGES
        rc = 0.5 * (edges[:-1] + edges[1:])
        ring_of = np.arange(fs.NB) // fs.NL
        layer_of = np.arange(fs.NB) % fs.NL
        g = dict(desc=gdesc, n_blocks_common=int(common.sum()), totals={}, layers=fs.LAYER_EDGES.tolist())
        print(f"\n== grid {tag} ({gdesc}); blocks valid in every run: {common.sum()}; data stars per block median "
              f"{int(np.median(nstar[common]))}")
        for lab, a in per.items():
            t = np.nansum(a[:, common], axis=(1, 2))
            g["totals"][lab] = dict(mean=float(t.mean()), sd=float(t.std(ddof=1)), se=float(t.std(ddof=1) / np.sqrt(len(t))))
            print(f"   total chi2 over the common blocks, {lab:5s}: mean {t.mean():8.1f}  sd {t.std(ddof=1):5.1f}  "
                  f"s.e. {t.std(ddof=1) / np.sqrt(len(t)):4.1f}")
        T1g = None
        pj = RESULTS_DIR / f"chi2_jitter_per_grid_{args.case}.csv"
        if pj.exists():
            jj = pd.read_csv(pj)
            jj = jj[jj.grid == tag]
            if len(jj) > 2:
                T1g = float(np.std(jj.chi2.values[1:] - jj.chi2.iloc[0], ddof=1))
        g["T1"] = T1g
        pairs = [("second", "main")] + ([("main", "truth"), ("second", "truth")] if "truth" in per else [])
        g["pairs"] = {}
        for b_, a_ in pairs:
            D = per[b_][:, common] .mean(0) - per[a_][:, common].mean(0)          # (n_common, 4)
            V = (per[b_][:, common].var(0, ddof=1) + per[a_][:, common].var(0, ddof=1)) / args.n_rep
            Dtot, Stot = D.sum(), np.sqrt(np.sum(per[b_][:, common].sum(axis=(1, 2)).var(ddof=1) / args.n_rep
                                                + per[a_][:, common].sum(axis=(1, 2)).var(ddof=1) / args.n_rep))
            unit = f" = {Dtot / T1g:+.2f} T1" if T1g else ""
            print(f"   Delta chi2 ({b_} - {a_}) = {Dtot:+.1f} +- {Stot:.1f}{unit}  [negative: the data prefer '{b_}']")
            by = {}
            print(f"      {'R band':10s} " + " ".join(f"|z| {lo:g}-{hi:g}".rjust(16) for lo, hi in zip(fs.LAYER_EDGES[:-1], fs.LAYER_EDGES[1:])))
            for lo, hi in R_BANDS:
                cells = []
                for z in range(fs.NL):
                    sel = (rc[ring_of[common]] >= lo) & (rc[ring_of[common]] < hi) & (layer_of[common] == z)
                    if sel.any():
                        dv = D[sel].sum()
                        sv = np.sqrt(np.sum((per[b_][:, common][:, sel].sum(axis=(1, 2)).var(ddof=1) +
                                             per[a_][:, common][:, sel].sum(axis=(1, 2)).var(ddof=1)) / args.n_rep))
                        cells.append(f"{dv:+7.1f} +-{sv:5.1f}")
                        by[f"R{lo}-{hi}_z{z}"] = (float(dv), float(sv))
                    else:
                        cells.append(" " * 16)
                print(f"      {f'{lo}-{hi} kpc':10s} " + " ".join(cells))
            fld = {FIELDS[f]: float(D[:, k].sum()) for k, f in enumerate(spec)}
            print("      by field: " + ", ".join(f"{k} {v:+.1f}" for k, v in fld.items()))
            g["pairs"][f"{b_}-{a_}"] = dict(delta=float(Dtot), se=float(Stot), by_band_layer=by, by_field=fld)
            full = np.full((fs.NB, len(spec)), np.nan)
            full[common] = D
            arrays[f"{tag}__{b_}-{a_}__delta_block_field"] = full
        arrays[f"{tag}__ring_edges"] = edges
        arrays[f"{tag}__layer_edges"] = fs.LAYER_EDGES
        arrays[f"{tag}__n_data_stars"] = nstar
        for lab, a in per.items():
            arrays[f"{tag}__{lab}__chi2_rep_block_field"] = a
        res["grids"][tag] = g
    fs.set_grid()
    np.savez(RESULTS_DIR / f"bimodality_origin_{args.case}.npz", **arrays)
    json.dump(res, open(RESULTS_DIR / f"bimodality_origin_{args.case}.json", "w"), indent=1)


if __name__ == "__main__":
    main()

"""Statistical and systematic errors on (Omega_b, phi_b, S, bar_frac) for the Test 2 searches.

Two errors, kept separate (PROGRESS, 2026-11-01: the tempering conflates them):

  statistical - the Fisher error of the likelihood itself, sigma = sqrt(diag(F^-1)) with
                F = h (dE/dtheta)(dE/dtheta)^T and E the library predictions whitened by the data
                covariance. dE/dtheta is a local linear fit over the N_NEAR design points nearest the truth.
                This is what the error would be with a correctly specified model and T = 1.

  systematic  - the offset of the recovered posterior median from the truth, across the five snapshots.
                The mean is the systematic shift; the scatter divided by sqrt(n) is its uncertainty.

Truth values are the simulation directory's finite-difference column (user's choice).
Parameters are handled in the units the search uses: Omega_b [km/s/kpc], phi_b [deg], ln S, ln bar_frac, so
the S and bar_frac numbers are fractional.

Usage: python quantify_test2_errors.py [--variant sf] [--scheme B]
"""
import argparse, json
import numpy as np, pandas as pd
from common import ROOT, RUNS_DIR, RESULTS_DIR, BAR_ANGLE_DEG
from nbody import SNAPSHOT_DIR
import test1_nbody_qmc_focused_fourier_stars as fs
from test1_nbody_qmc_focused_fourier_stars import rotation
from test2b_qmc_distance_error_run import scaled_design

SCHEMES = {"A": "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024",
           "B": "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"}
RUNS = [("test2_t2_sf_10pct", 2.0), ("test2_t3_sf_10pct", 3.0), ("test2_t4_sf_10pct", 4.0),
        ("test2_t5.5_sf_10pct", 5.5), ("test2_t7_sf_10pct", 7.0)]
TABLE = SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv"   # ships with the snapshots, in mock/Nbody
KEYS = ["Omega_b", "phi_b", "ln S", "ln bar_frac"]
N_NEAR = 60


def fisher(run, variant, scheme):
    """Marginal 1-sigma on each parameter from the likelihood's own curvature."""
    meta = json.load(open(ROOT / "results" / run / "run_meta.json"))
    design, _ = scaled_design(meta["galaxy"])
    out = RUNS_DIR / run
    s0 = np.load(out / "fourier_000.npz")
    cov0 = s0["cov_data_joint"].reshape(fs.NB, 6 * fs.NK, 6 * fs.NK)
    ok = s0["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK
    n_ring, spec = fs.SCHEMES[scheme]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    h = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)
    blocks = [r * fs.NL + z for r in range(n_ring) for z in range(fs.NL) if ok[r * fs.NL + z]]
    Li = [np.linalg.inv(np.linalg.cholesky(cov0[b][np.ix_(idx, idx)])) for b in blocks]
    E, P, kinds = [], [], []
    for i in design.idx:
        i = int(i)
        f = out / f"fourier_{i:03d}.npz"
        if not f.exists():
            continue
        cl = np.load(f)[f"coef_lib_{variant}"].reshape(fs.NB, -1)
        p = design.set_index("idx").loc[i]
        clr = cl @ rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        E.append(np.concatenate([Li[j] @ clr[b][idx] for j, b in enumerate(blocks)]))
        P.append([p.Omega_abs, p.phi_b_deg, np.log(p.S), np.log(p.bar_frac)])
        kinds.append(p.kind)
    E, P, kinds = np.array(E), np.array(P), np.array(kinds)
    t = np.flatnonzero(kinds == "truth")[0]
    q = np.flatnonzero(kinds == "qmc")
    dn = np.linalg.norm((P[q] - P[t]) / P[q].std(0), axis=1)
    near = q[np.argsort(dn)[:N_NEAR]]
    A = np.column_stack([np.ones(len(near)), P[near] - P[t]])
    J = np.linalg.lstsq(A, E[near], rcond=None)[0][1:]
    return np.sqrt(np.diag(np.linalg.inv(h * (J @ J.T))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", default="sf")
    ap.add_argument("--scheme", default="B", choices=list(SCHEMES))
    args = ap.parse_args()
    sch = SCHEMES[args.scheme]
    tab = pd.read_csv(TABLE)

    rows = []
    for run, tm in RUNS:
        om = float(tab[np.isclose(tab.snap_time, tm)].pattern_speed_fd.iloc[0])
        truth = np.array([om, BAR_ANGLE_DEG, 0.0, 0.0])
        d = json.load(open(RESULTS_DIR / run / f"emulator_results_logLs_{sch}_{args.variant}.json"))
        p = d["posterior"]
        med = np.array([p["Omega_abs"]["p50"], p["phi_b_deg"]["p50"],
                        np.log(p["S"]["p50"]), np.log(p["bar_frac"]["p50"])])
        wid = np.array([0.5 * (p[k]["p84"] - p[k]["p16"]) for k in ("Omega_abs", "phi_b_deg")] +
                       [0.5 * (np.log(p[k]["p84"]) - np.log(p[k]["p16"])) for k in ("S", "bar_frac")])
        rows.append(dict(run=run, stat=fisher(run, args.variant, sch), bias=med - truth,
                         tempered=wid, truth=truth))

    print(f"scheme {args.scheme}, {args.variant} weighting, 10% distance error + selection function")
    print(f"truth: Omega_b = tabulated finite-difference column, phi_b = {BAR_ANGLE_DEG}, S = bar_frac = 1\n")
    print(f"{'parameter':<12s} {'statistical':>12s} {'tempered':>10s} {'ratio':>7s} "
          f"{'systematic':>20s} {'sys/stat':>9s}")
    res = {}
    for j, k in enumerate(KEYS):
        st = np.mean([r["stat"][j] for r in rows])
        tp = np.mean([r["tempered"][j] for r in rows])
        b = np.array([r["bias"][j] for r in rows])
        sysm, syse = b.mean(), b.std(ddof=1) / np.sqrt(len(b))
        unit = "" if j < 2 else " (frac)"
        print(f"{k:<12s} {st:12.3f} {tp:10.3f} {tp/st:7.1f} "
              f"{sysm:+9.3f} +/- {syse:5.3f} {abs(sysm)/st:9.1f}")
        res[k] = dict(statistical=st, tempered=tp, systematic=sysm, systematic_err=syse,
                      per_run_bias={r["run"]: float(r["bias"][j]) for r in rows})
    print(f"\n  statistical = Fisher 1-sigma (correctly specified model, T = 1), mean over the 5 runs")
    print(f"  tempered    = the width actually quoted, from logLs = -chi2/(2 sqrt(2 N))")
    print(f"  systematic  = mean(recovered - truth) over the 5 runs, +/- scatter/sqrt(5)")
    print(f"\n  fractional: Omega_b {100*res['Omega_b']['systematic']/np.mean([r['truth'][0] for r in rows]):+.1f}% "
          f"systematic, {100*res['Omega_b']['statistical']/np.mean([r['truth'][0] for r in rows]):.1f}% statistical")
    print(f"              S {100*res['ln S']['systematic']:+.1f}%, bar_frac {100*res['ln bar_frac']['systematic']:+.1f}%")
    with open(RESULTS_DIR / f"test2_error_budget_{args.variant}_scheme{args.scheme}.json", "w") as fh:
        json.dump({k: {kk: (vv if not isinstance(vv, np.floating) else float(vv))
                       for kk, vv in v.items()} for k, v in res.items()}, fh, indent=2, default=float)
    print(f"\n  wrote results/test2_error_budget_{args.variant}_scheme{args.scheme}.json")


if __name__ == "__main__":
    main()

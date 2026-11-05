"""Leave-one-out validation of the logL emulator, with calibration of its error bars.

The emulator already reports a LOO rms, but only for the 32 best training points and only as a single
number. This does the full test and keeps the per-point results:

  for each training point i, refit the whole model (quadratic mean + GP hyperparameters) on the other
  256 points and predict point i. Nothing about point i enters its own prediction.

Three things come out, and they answer different questions:

  residual = predicted - actual        how big is the error, in tempered logL (1 sigma of the
                                       posterior is dlogL ~ 0.5)
  pull     = residual / predicted sd   is the emulator's *own* error bar honest? A calibrated
                                       emulator gives pulls with mean 0 and sd 1. sd < 1 means it
                                       over-states its uncertainty, > 1 means it under-states it.
  residual against logL                is the error uniform over the design, or worse near the peak
                                       where the posterior is actually determined?

Writes results/emulator_loo_<tag>.csv (one row per training point) and
figures/emulator_validation/<tag>_loo.png, plus a summary table.

Usage: python validate_emulator_loo.py [--runs ...] [--variant sfratio] [--restarts 3]
"""
import argparse
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from common import ROOT, RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from test1_nbody_qmc_focused_run import BOX
from test1_gp_emulator import quad_features, KEYS

SCHEME_B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
DEFAULT_RUNS = ["test2_t2_sf_12p5pct", "test2_t3_sf_12p5pct", "test2_t4_sf_12p5pct",
                "test2_t5.5_sf_12p5pct", "test2_t7_sf_12p5pct", "test2_t2.375_sf_12p5pct"]
LO = np.array([BOX["Omega_abs"][0], BOX["phi_b_deg"][0], BOX["lnS"][0], BOX["lnf"][0]])
HI = np.array([BOX["Omega_abs"][1], BOX["phi_b_deg"][1], BOX["lnS"][1], BOX["lnf"][1]])


def to_unit(df, scale):
    p = np.column_stack([df["Omega_abs"] / scale, df["phi_b_deg"], np.log(df["S"]), np.log(df["bar_frac"])])
    return (p - LO) / (HI - LO)


def fit(u, y, restarts):
    """The emulator of test1_gp_emulator.Emulator, returning (predict_fn)."""
    ymax = y.max(); y = y - ymax
    beta, *_ = np.linalg.lstsq(quad_features(u), y, rcond=None)
    r = y - quad_features(u) @ beta
    k = (ConstantKernel(np.var(r), (1e-3 * np.var(r) + 1e-6, 1e3 * np.var(r) + 1)) *
         Matern(length_scale=[0.3] * 4, length_scale_bounds=(0.02, 10), nu=2.5) +
         WhiteKernel(100.0, (1e-2, 1e5)))
    gp = GaussianProcessRegressor(k, normalize_y=False, n_restarts_optimizer=restarts,
                                  random_state=0).fit(u, r)

    def predict(q):
        m, s = gp.predict(np.atleast_2d(q), return_std=True)
        return quad_features(np.atleast_2d(q)) @ beta + m + ymax, s
    return predict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=DEFAULT_RUNS)
    ap.add_argument("--variant", default="sfratio")
    ap.add_argument("--scheme", default=SCHEME_B)
    ap.add_argument("--restarts", type=int, default=3,
                    help="hyperparameter restarts per fold (the production emulator uses 15; 3 is "
                         "enough for a LOO sweep and 5x faster)")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()
    col = f"logLs_{args.scheme}_{args.variant}"
    rows = []
    for run in args.runs:
        f = RESULTS_DIR / run / "evals.csv"
        if not f.exists():
            print(f"  skipping {run}: no evals.csv")
            continue
        ev = pd.read_csv(f)
        if col not in ev:
            print(f"  skipping {run}: no {col}")
            continue
        meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
        tr = ev[(ev.kind == "qmc") & np.isfinite(ev[col])]
        u, y = to_unit(tr, meta.get("omega_scale", 1.0)), tr[col].values
        print(f"{run}: {len(u)} training points, LOO with {args.restarts} restarts ...", flush=True)
        pred = np.empty(len(u)); sd = np.empty(len(u))
        for i in range(len(u)):
            p = fit(np.delete(u, i, 0), np.delete(y, i), args.restarts)
            m, s = p(u[i:i + 1]); pred[i], sd[i] = m[0], s[0]
        for j in range(len(u)):
            rows.append(dict(run=run, idx=int(tr.idx.iloc[j]), logL=y[j], loo_pred=pred[j],
                             loo_sd=sd[j], residual=pred[j] - y[j],
                             pull=(pred[j] - y[j]) / sd[j] if sd[j] > 0 else np.nan,
                             rank=int((y > y[j]).sum())))
    if not rows:
        raise SystemExit("nothing to do")
    df = pd.DataFrame(rows)
    out = RESULTS_DIR / f"emulator_loo_{args.variant}{args.tag}.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {rel(out)}  ({len(df)} rows)")

    # --- summary table ---------------------------------------------------------------------------
    print(f"\n{'run':24s} {'n':>4s} {'rms all':>8s} {'rms top32':>10s} {'pull mean':>10s} {'pull sd':>8s}")
    print("-" * 70)
    for run, g in df.groupby("run", sort=False):
        t = g.nsmallest(32, "rank")
        print(f"{run:24s} {len(g):4d} {np.sqrt((g.residual**2).mean()):8.2f} "
              f"{np.sqrt((t.residual**2).mean()):10.2f} {g.pull.mean():+10.2f} {g.pull.std():8.2f}")
    print("-" * 70)
    t = df.nsmallest(32 * df.run.nunique(), "rank")
    print(f"{'ALL':24s} {len(df):4d} {np.sqrt((df.residual**2).mean()):8.2f} "
          f"{np.sqrt((t.residual**2).mean()):10.2f} {df.pull.mean():+10.2f} {df.pull.std():8.2f}")

    # --- figure ----------------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    top = df["rank"] < 32
    ax = axes[0]
    lim = [df.logL.min() - 2, df.logL.max() + 2]
    ax.plot(lim, lim, color=INK_MUTED, ls="--", lw=1)
    ax.plot(df.logL[~top], df.loo_pred[~top], "o", ms=3, color=BLUE_RAMP[2], label="rest of the design")
    ax.plot(df.logL[top], df.loo_pred[top], "o", ms=5, color=ACCENT, label="32 best per run")
    ax.set_xlabel("actual logLs"); ax.set_ylabel("left-out prediction")
    ax.set_title("predicted against actual", color=INK); ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.axhline(0, color=INK_MUTED, ls="--", lw=1)
    ax.plot(df.logL[~top], df.residual[~top], "o", ms=3, color=BLUE_RAMP[2])
    ax.plot(df.logL[top], df.residual[top], "o", ms=5, color=ACCENT)
    # running rms, to show whether the error grows towards the peak
    o = np.argsort(df.logL.values); xs = df.logL.values[o]; rs = df.residual.values[o]
    w = max(20, len(xs) // 12)
    mid = np.array([xs[i:i + w].mean() for i in range(0, len(xs) - w, w // 2)])
    rms = np.array([np.sqrt((rs[i:i + w]**2).mean()) for i in range(0, len(xs) - w, w // 2)])
    ax.plot(mid, rms, "-", color=INK, lw=1.6, label="running rms")
    ax.plot(mid, -rms, "-", color=INK, lw=1.6)
    ax.set_xlabel("actual logLs"); ax.set_ylabel("predicted - actual")
    ax.set_title("error against position in the design", color=INK); ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    b = np.linspace(-4, 4, 33)
    ax.hist(df.pull.clip(-4, 4), bins=b, color=BLUE_RAMP[3], alpha=0.85, density=True,
            label=f"all ({df.pull.std():.2f})")
    ax.hist(df.pull[top].clip(-4, 4), bins=b, histtype="step", color=ACCENT, lw=1.8, density=True,
            label=f"32 best ({df.pull[top].std():.2f})")
    g = np.linspace(-4, 4, 200)
    ax.plot(g, np.exp(-g**2 / 2) / np.sqrt(2 * np.pi), color=INK, lw=1.4, ls="--",
            label="N(0, 1) if calibrated")
    ax.set_xlabel("pull = (predicted - actual) / predicted sd"); ax.set_ylabel("density")
    ax.set_title("are the emulator's error bars honest?", color=INK); ax.legend(frameon=False, fontsize=9)

    fig.suptitle(f"Emulator leave-one-out validation — {args.variant}, scheme B, 12.5% distance error "
                 f"({df.run.nunique()} runs, {len(df)} points)", color=INK, y=1.0)
    fig.tight_layout()
    savefig(fig, f"emulator_validation/loo_{args.variant}{args.tag}")


if __name__ == "__main__":
    main()

"""Recovered bar strength profile A2/A0(R) against the truth, for every modelled galaxy.

Rather than trying to read a bar strength off the bar_frac parameter - which is an amplitude on a
difference of two potentials and has no clean physical meaning on its own - this builds the trial
potential at each posterior draw and measures the thing directly: the m = 2 Fourier amplitude of
the midplane density, A2/A0, as a function of radius.

The trial potential's shape depends only on (S, bar_frac): phi_b rotates the stars, not the
potential, and Omega_b is a pattern speed, so neither changes A2/A0(R). The truth is the family at
(bar_frac, S) = (1, 1), which is identically the galaxy's own potential with its bar along x.

Two bands per galaxy, from the same chain:
  flat      the model's own posterior
  + prior   the same posterior with the kinematic prior on phi_b (plot_corner_with_kinematic_prior)
drawn separately rather than by reweighting one set, because the prior keeps only ~25% of the
chain's effective weight and reweighting a small sample would thin the band too far.

Sample size: 300 draws per band. For a roughly Gaussian posterior the 16th/84th percentile has a
standard error of about 1.3 sigma / sqrt(n), so n = 300 pins the band edges to ~7% of the posterior
width - well below the band width itself. Each draw costs ~0.35 s (one CylSpline re-expansion).

Usage: python plot_bar_strength_profile_recovery.py [--n-draw 300] [--workers 6]
Figure: figures/test2_bar_strength_profile_recovery.png
Table:  results/bar_strength_profile_recovery.npz
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import emcee
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import (RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT,
                    BLUE_RAMP)
import test1_gp_emulator as ge
import steady_state_likelihood as ssl
from plot_corner_with_kinematic_prior import kinematic_prior, SCHEME_B

R_GRID = np.arange(0.3, 8.0 + 1e-9, 0.2)
N_PHI = 64
_RR, _PP = np.meshgrid(R_GRID, np.linspace(0, 2 * np.pi, N_PHI, endpoint=False), indexing="ij")
_PTS = np.column_stack([(_RR * np.cos(_PP)).ravel(), (_RR * np.sin(_PP)).ravel(),
                        np.zeros(_RR.size)])


def a2_over_a0(galaxy, bar_frac, S):
    """m = 2 amplitude of the midplane density, normalised by the azimuthal mean, vs radius."""
    p = ssl.bar_aligned_potential(bar_frac, S, galaxy)
    rho = p.density(_PTS).reshape(_RR.shape)
    a0 = rho.mean(axis=1)
    a2 = 2 * np.abs((rho * np.exp(-2j * _PP)).mean(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(a0 > 0, a2 / a0, np.nan)


def one_galaxy(job):
    run, column, n_draw = job
    meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
    gal = meta["galaxy"]
    phi_kin, sig_kin, info = kinematic_prior(gal)

    ge.use_run(run)
    df, _ = ge.load_table(column, None, None)
    tr = df[df.kind == "qmc"]
    u, y = ge.to_unit(tr), tr["logL"].values
    emu = ge.Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    sam = emcee.EnsembleSampler(
        32, 4, lambda q: -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0])
    sam.random_state = np.random.RandomState(2).get_state()
    sam.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1),
                 3000, progress=False)
    post = ge.from_unit(sam.get_chain(discard=1000, thin=5, flat=True))

    rng = np.random.default_rng(7)
    w = np.exp(-0.5 * ((post[:, 1] - phi_kin) / sig_kin) ** 2)
    idx_flat = rng.choice(len(post), n_draw, replace=True)
    idx_prio = rng.choice(len(post), n_draw, replace=True, p=w / w.sum())

    out = {}
    for tag, idx in (("flat", idx_flat), ("prior", idx_prio)):
        curves = np.array([a2_over_a0(gal, post[i, 3], post[i, 2]) for i in idx])
        out[tag] = np.nanpercentile(curves, [16, 50, 84], axis=0)
    out["truth"] = a2_over_a0(gal, 1.0, 1.0)
    out["galaxy"] = gal
    out["ess"] = float(w.sum() ** 2 / np.sum(w ** 2))
    return run, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draw", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--column", default=f"logLs_{SCHEME_B}_sfratio")
    ap.add_argument("--runs", nargs="+", default=None)
    args = ap.parse_args()
    set_style()

    runs = args.runs or [f"test2_{t}_sf_12p5pct" for t in
                         ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5",
                          "t5", "t5.5", "t6.25", "t7",
                          "au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24", "au26", "au27", "au28"]]
    runs = [r for r in runs
            if (RESULTS_DIR / r / f"emulator_results_{args.column}.json").exists()]
    print(f"{len(runs)} galaxies, {args.n_draw} draws per band "
          f"(~{2 * args.n_draw * 0.35 / 60:.0f} min each serial)", flush=True)

    res = {}
    with ProcessPoolExecutor(args.workers) as ex:
        for k, (run, o) in enumerate(ex.map(one_galaxy, [(r, args.column, args.n_draw)
                                                         for r in runs]), 1):
            res[run] = o
            pk = np.nanmax(o["truth"])
            pm = np.nanmax(o["prior"][1])
            print(f"  [{k:2d}/{len(runs)}] {o['galaxy']:10s} peak A2/A0 truth {pk:.3f}, "
                  f"model {pm:.3f} ({100*(pm-pk)/pk:+.0f}%)", flush=True)

    np.savez_compressed(RESULTS_DIR / "bar_strength_profile_recovery.npz", R=R_GRID,
                        **{f"{r}_{k}": v for r, o in res.items() for k, v in o.items()
                           if k in ("truth", "flat", "prior")})

    n = len(runs)
    nc = 5
    nr = int(np.ceil(n / nc))
    fig, AX = plt.subplots(nr, nc, figsize=(3.2 * nc, 2.7 * nr), sharex=True, sharey=True)
    AX = np.atleast_2d(AX)
    for i, run in enumerate(runs):
        ax = AX[i // nc][i % nc]
        o = res[run]
        ax.fill_between(R_GRID, o["flat"][0], o["flat"][2], color=INK_MUTED, alpha=0.30, lw=0,
                        label="model, flat prior (16–84%)")
        ax.fill_between(R_GRID, o["prior"][0], o["prior"][2], color=BLUE_RAMP[3], alpha=0.60,
                        lw=0, label="model, $\\phi_b$ prior (16–84%)")
        ax.plot(R_GRID, o["prior"][1], "-", color=BLUE_RAMP[6], lw=1.6, label="model median")
        ax.plot(R_GRID, o["truth"], "-", color=INK, lw=2.0, label="truth")
        ax.set_title(o["galaxy"], color=INK, fontsize=9)
        ax.set_ylim(0, 0.8)
        ax.set_xlim(R_GRID[0], R_GRID[-1])
        if i == 0:
            ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    for j in range(n, nr * nc):
        AX[j // nc][j % nc].axis("off")
    for j in range(nc):
        AX[-1][j].set_xlabel("R [kpc]")
    for j in range(nr):
        AX[j][0].set_ylabel(r"$A_2/A_0$ (midplane density)")
    fig.suptitle("Recovered bar strength profile against the truth — trial potential rebuilt at "
                 f"{args.n_draw} posterior draws per galaxy", color=INK, y=0.995, fontsize=12)
    fig.tight_layout()
    savefig(fig, "test2_bar_strength_profile_recovery")
    print(f"wrote {rel(RESULTS_DIR / 'bar_strength_profile_recovery.npz')}")


if __name__ == "__main__":
    main()

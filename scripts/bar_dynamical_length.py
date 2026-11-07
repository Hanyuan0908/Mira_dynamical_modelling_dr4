"""Dynamical bar length from trapped x1 orbits, following Petersen et al. (2021).

The photometric length (the radius where the density m = 2 amplitude falls off) is not what a
dynamicist means by the bar. Petersen, Weinberg & Katz (2021) argue for the DYNAMICAL length: the
radial extent of the orbits actually trapped by the bar, which is a property of the potential and
pattern speed rather than of the light.

That is directly computable from what this project already fits. The trial potential's shape is
set by (S, bar_frac) and the frame by Omega_b, so a posterior sample gives a potential and a
pattern speed, and the dynamical length follows from integrating orbits in it - no new mock, no
re-scoring, no orbit library. The truth is the family at (bar_frac, S) = (1, 1) with the galaxy's
true Omega_b, which is identically the galaxy's own potential.

Method. Orbits are launched from the bar major axis (the trial potential always puts the bar along
x) over a grid of radii and a range of speeds about circular, and integrated in the frame rotating
at Omega_b. An orbit is bar-supporting if its azimuth LIBRATES rather than circulates: 2*phi stays
within 2*pi instead of winding. The dynamical length is the largest apocentre among the trapped
orbits.

Two AGAMA conventions matter and both were verified numerically before use (a circular orbit in an
axisymmetric potential must keep r fixed, and a star whose angular speed equals the pattern speed
must sit at a fixed azimuth in the rotating frame):
  * agama.orbit(..., Omega=) takes the initial conditions in the INERTIAL frame;
  * it returns the trajectory in the ROTATING frame.
Passing rotating-frame velocities instead sends orbits plunging to the centre and reports no
trapped orbits at all.

Usage: python bar_dynamical_length.py [--n-draw 120] [--workers 8]
Figure: figures/test2_bar_dynamical_length.png
Table:  results/bar_dynamical_length.csv
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import agama
import emcee
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import potential_family as pf
import test1_gp_emulator as ge
from plot_corner_with_kinematic_prior import SCHEME_B

N_R, N_V, T_BAR, N_SAMP = 44, 15, 6, 400
V_FAC = (0.35, 1.45)


def dynamical_length(galaxy, bar_frac, S, Omega):
    """Largest apocentre among orbits trapped by the bar. Omega in AGAMA's sign convention."""
    p = ssl.bar_aligned_potential(bar_frac, S, galaxy)
    R = np.linspace(0.3, 14.0, N_R)
    vc = np.sqrt(np.maximum(-p.force(np.column_stack([R, 0 * R, 0 * R]))[:, 0] * R, 0))
    fac = np.linspace(*V_FAC, N_V)
    # INERTIAL ICs, co-rotating sense (sign of Omega); the integrator handles the frame
    ic = np.array([[r, 0, 0, 0, np.sign(Omega) * f * v, 0] for r, v in zip(R, vc) for f in fac])
    orb = agama.orbit(potential=p, ic=ic, time=abs(2 * np.pi / Omega) * T_BAR,
                      trajsize=N_SAMP, Omega=Omega)
    best = 0.0
    for i in range(len(ic)):
        xy = np.asarray(orb[i][1])[:, :2]
        ph = np.unwrap(2 * np.arctan2(xy[:, 1], xy[:, 0]))
        if np.ptp(ph) < 2 * np.pi:                       # librates -> trapped by the bar
            best = max(best, float(np.hypot(*xy.T).max()))
    return best


def one(job):
    run, galaxy, column, n_draw = job
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
    idx = np.random.default_rng(5).choice(len(post), n_draw, replace=True)
    # Omega_abs is the magnitude; the family's true pattern speed is negative (clockwise)
    L = np.array([dynamical_length(galaxy, post[i, 3], post[i, 2], -post[i, 0]) for i in idx])
    Ltrue = dynamical_length(galaxy, 1.0, 1.0, pf.omega_true(galaxy))
    q = np.percentile(L[L > 0], [16, 50, 84]) if (L > 0).any() else [np.nan] * 3
    return dict(run=run, galaxy=galaxy, L_true=Ltrue, L16=q[0], L50=q[1], L84=q[2],
                n_valid=int((L > 0).sum()), n_draw=n_draw)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-draw", type=int, default=120)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--column", default=f"logLs_{SCHEME_B}_sfratio")
    ap.add_argument("--tags", nargs="+", default=None)
    args = ap.parse_args()
    set_style()
    tags = args.tags or ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5",
                         "t5.5", "t6.25", "t7", "au5", "au9", "au10", "au12", "au14", "au17",
                         "au18", "au20", "au22", "au24", "au26", "au27", "au28"]
    jobs = []
    for t in tags:
        run = f"test2_{t}_sf_12p5pct"
        meta = RESULTS_DIR / run / "run_meta.json"
        if not (RESULTS_DIR / run / f"emulator_results_{args.column}.json").exists():
            print(f"  skip {t}: no posterior")
            continue
        jobs.append((run, json.load(open(meta))["galaxy"], args.column, args.n_draw))
    print(f"{len(jobs)} galaxies, {args.n_draw} draws each", flush=True)

    import pandas as pd
    rows = []
    with ProcessPoolExecutor(args.workers) as ex:
        for k, r in enumerate(ex.map(one, jobs), 1):
            rows.append(r)
            print(f"  [{k:2d}/{len(jobs)}] {r['galaxy']:10s} L_dyn truth {r['L_true']:5.2f}  "
                  f"model {r['L50']:5.2f} (+{r['L84']-r['L50']:.2f}/-{r['L50']-r['L16']:.2f})",
                  flush=True)
    D = pd.DataFrame(rows)
    D["bias"] = 100 * (D.L50 / D.L_true - 1)
    D.to_csv(RESULTS_DIR / "bar_dynamical_length.csv", index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 5.5), gridspec_kw=dict(wspace=0.26))
    au = D.galaxy.str.startswith("au")
    lim = (0.85 * D.L_true.min(), 1.15 * D.L_true.max())
    a1.plot(lim, lim, "--", color=INK, lw=1.1, zorder=1, label="perfect recovery")
    for m, c, mk, lab in ((au, ACCENT, "s", "Auriga"), (~au, BLUE_RAMP[5], "o", "N-body")):
        a1.errorbar(D.L_true[m], D.L50[m], yerr=[(D.L50 - D.L16)[m], (D.L84 - D.L50)[m]],
                    fmt=mk, ms=7, color=c, mec=INK, mew=0.8, capsize=3, lw=1.2, zorder=4, label=lab)
    a1.set_xlim(*lim); a1.set_ylim(*lim); a1.set_aspect("equal")
    a1.set_xlabel("true dynamical bar length [kpc]")
    a1.set_ylabel("recovered dynamical bar length [kpc]")
    a1.legend(loc="upper left", fontsize=9)
    a1.set_title(f"{len(D)} galaxies, trapped-orbit length", color=INK_2, fontsize=10)

    a2.axhline(0, color=INK, lw=1.1, ls="--", zorder=1)
    for m, c, mk in ((au, ACCENT, "s"), (~au, BLUE_RAMP[5], "o")):
        a2.plot(D.L_true[m], D.bias[m], mk, ms=8, color=c, mec=INK, mew=0.8, zorder=4)
    for _, r in D.iterrows():
        a2.annotate(r.galaxy.replace("nbody_", "").replace("nbody", "t7"),
                    (r.L_true, r.bias), xytext=(0, 9), textcoords="offset points",
                    ha="center", fontsize=7, color=INK_2)
    a2.set_xlabel("true dynamical bar length [kpc]")
    a2.set_ylabel("bias in dynamical length [%]")
    a2.set_title(f"mean bias {D.bias.mean():+.1f}%, mean |bias| {D.bias.abs().mean():.1f}%",
                 color=INK_2, fontsize=10)
    fig.suptitle("Dynamical bar length from trapped x1 orbits (Petersen et al. 2021), "
                 "recovered against truth", color=INK, y=0.98, fontsize=12)
    fig.tight_layout()
    savefig(fig, "test2_bar_dynamical_length")
    print(f"\n  mean bias {D.bias.mean():+.1f}%, mean |bias| {D.bias.abs().mean():.1f}%")
    print(f"wrote {rel(RESULTS_DIR / 'bar_dynamical_length.csv')}")


if __name__ == "__main__":
    main()

"""Is the second pattern-speed peak of the OSARG posterior in Sormani+22 (|Omega_b| ~ 37-40, next to the main one at ~30)
real, or an emulator artefact? (user, 2026-11-24)

Four checks, on the pipeline's tempered column (logLs, scheme B, sfratio):
  1. The evaluated design points alone, no emulator: the best logLs per 2.5 km/s/kpc bin of |Omega_b| (other
     parameters free) and where the 20 best points lie. Single-run jitter is sigma_1 / (2 sqrt(2N)) = 0.42 tempered.
  2. The emulator's two modes: the chain split at the histogram valley between them, and the emulator maximum inside
     each mode's |Omega_b| range (the parameters used in check 4), with the predicted logL difference.
  3. Robustness: 10 refits, each on a random 90% of the design points (jackknife). Fraction of the chain in the
     second mode, and its peak position.
  4. Direct: integrate and score both mode maxima exactly as the pipeline, 12 runs each (the point itself and 11 tiny
     random perturbations, as chi2_jitter_small_perturbations.py), so the jitter averages down. Compare the mean chi2
     difference, with its standard error, against the emulator's prediction.

Usage: python check_bimodality_sormani22.py
Output: results/osarg_zhang2024_sormani22/bimodality_check.json and a printed report
"""
import json
import time

import numpy as np
import pandas as pd
import h5py
import emcee
from scipy.optimize import minimize

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import test1_gp_emulator as ge
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from chi2_jitter_small_perturbations import RANDOM_WIDTH

RUN = "osarg_zhang2024_sormani22"
MOCK = "osarg_zhang2024_observed_R8.hdf5"
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
N_REP = 12


def sample(emu, u, y, seed=2):
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    s = emcee.EnsembleSampler(32, 4, lambda q: -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0])
    s.random_state = np.random.RandomState(seed).get_state()
    s.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1), 3000, progress=False)
    return s.get_chain(discard=1000, thin=5, flat=True)


def max_in_range(emu, u, y, lo_u, hi_u):
    best = None
    bounds = [(lo_u, hi_u), (0, 1), (0, 1), (0, 1)]
    inside = (u[:, 0] >= lo_u) & (u[:, 0] <= hi_u)
    starts = np.vstack([u[inside][np.argsort(-y[inside])[:8]], np.random.default_rng(3).random((16, 4))])
    starts[:, 0] = np.clip(starts[:, 0], lo_u, hi_u)
    for x0 in starts:
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=bounds, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun


def main():
    ge.use_run(RUN)
    ev = pd.read_csv(RESULTS_DIR / RUN / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    y = tr[f"logLs_{B}_sfratio"].values
    y = y - y.max()
    u = ge.to_unit(tr)
    om = tr.Omega_abs.values
    nd = int(tr[f"ndata_{B}_sfratio"].mode().iloc[0])
    temp = 2 * np.sqrt(2 * nd)
    report = {}

    # 1. evaluated points
    edges = np.arange(15, 45.01, 2.5)
    prof = [float(y[(om >= a) & (om < b)].max()) for a, b in zip(edges[:-1], edges[1:])]
    top = np.argsort(-y)[:20]
    print("1. evaluated points: best logLs per 2.5 km/s/kpc bin, 15 -> 45:", " ".join(f"{v:5.2f}" for v in prof))
    print(f"   20 best points: |Omega_b| " + " ".join(f"{om[i]:.1f}" for i in top))
    print(f"   their logLs: " + " ".join(f"{y[i]:.2f}" for i in top))
    print(f"   single-run jitter sigma_1/T = 31.7/{temp:.1f} = {31.7 / temp:.2f} tempered; difference of two runs "
          f"{44.8 / temp:.2f}")
    report["profile_bins"], report["profile"] = edges.tolist(), prof

    # 2. emulator modes
    emu = ge.Emulator(u, y)
    ch = sample(emu, u, y)
    o = ge.from_unit(ch)[:, 0]
    h, e = np.histogram(o, bins=np.arange(15, 45.01, 1.0))
    c = 0.5 * (e[1:] + e[:-1])
    pkA = float(c[(c >= 25) & (c <= 33.5)][np.argmax(h[(c >= 25) & (c <= 33.5)])])
    pkB = float(c[(c >= 35.5) & (c <= 42.5)][np.argmax(h[(c >= 35.5) & (c <= 42.5)])])
    sel = (c > pkA) & (c < pkB)
    valley = float(c[sel][np.argmin(h[sel])])
    fracB = float(np.mean(o > valley))
    print(f"   chain histogram peaks at {pkA:.1f} and {pkB:.1f} (counts {h[c == pkA][0]}, {h[c == pkB][0]}); valley {valley:.1f} "
          f"(count {h[c == valley][0]})")
    tu = lambda v: (v - ge.LO[0]) / (ge.HI[0] - ge.LO[0])
    xa, la = max_in_range(emu, u, y, tu(25.0), tu(33.5))     # mode A window
    xb, lb = max_in_range(emu, u, y, tu(35.5), tu(42.5))     # mode B window
    xg, lg = max_in_range(emu, u, y, 0.0, 1.0)
    print(f"   emulator global maximum: Omega {ge.from_unit(xg[None])[0][0]:.2f}, logL {lg:.2f}")
    pa, pb = ge.from_unit(xa[None])[0], ge.from_unit(xb[None])[0]
    print(f"2. emulator: chain valley at |Omega_b| {valley:.1f}; fraction above it {fracB:.2f}")
    print(f"   mode A maximum: Omega {pa[0]:.2f}, phi {pa[1]:.2f}, S {pa[2]:.3f}, f {pa[3]:.3f}; predicted logL {la:.2f}")
    print(f"   mode B maximum: Omega {pb[0]:.2f}, phi {pb[1]:.2f}, S {pb[2]:.3f}, f {pb[3]:.3f}; predicted logL {lb:.2f}")
    print(f"   predicted difference B - A: {lb - la:+.2f} tempered ({(lb - la) * temp:+.0f} in chi2, sign flipped)")
    report.update(valley=valley, frac_mode_B=fracB, mode_A=pa.tolist(), mode_B=pb.tolist(), pred_logL_A=la, pred_logL_B=lb)

    # 3. jackknife
    rng = np.random.default_rng(11)
    jk = []
    for k in range(10 if __import__("os").environ.get("SKIP_JK") is None else 0):
        keep = rng.random(len(y)) < 0.9
        ek = ge.Emulator(u[keep], y[keep] - y[keep].max())
        ok_ = ge.from_unit(sample(ek, u[keep], y[keep] - y[keep].max(), seed=2 + k))[:, 0]
        hk, _ = np.histogram(ok_, bins=e)
        above = c > valley
        jk.append(dict(frac_B=float(np.mean(ok_ > valley)), peak_B=float(c[above][np.argmax(hk[above])]),
                       peak_A=float(c[~above][np.argmax(hk[~above])]), median=float(np.median(ok_)),
                       p16=float(np.percentile(ok_, 16)), p84=float(np.percentile(ok_, 84))))
        print(f"3. jackknife {k}: fraction above {valley:.1f}: {jk[-1]['frac_B']:.2f}; histogram peaks {jk[-1]['peak_A']:.1f} | "
              f"{jk[-1]['peak_B']:.1f}; Omega {jk[-1]['median']:.2f} ({jk[-1]['p16']:.2f}-{jk[-1]['p84']:.2f})", flush=True)
    report["jackknife"] = jk

    # 4. direct evaluation of both maxima
    meta = json.load(open(RESULTS_DIR / RUN / "run_meta.json"))
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / MOCK
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    cd0, Cd, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    cd0 = cd0.reshape(fs.NB, -1)
    nr, spec = fs.SCHEMES[B]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    hart = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)
    W = {b: np.linalg.inv(np.linalg.cholesky(Cd[b][np.ix_(idx, idx)]))
         for b in (r * fs.NL + z for r in range(nr) for z in range(fs.NL)) if ok[b]}
    ssl.ACCURACY = old.ACCURACY
    rng = np.random.default_rng(2124)
    direct = {}
    for tag, p in (("A", pa), ("B", pb)):
        vals = []
        for r in range(N_REP):
            q = dict(Omega_abs=p[0], phi_b_deg=p[1], S=p[2], bar_frac=p[3])
            if r:
                q = {k: v + RANDOM_WIDTH[k] * rng.standard_normal() for k, v in q.items()}
            t0 = time.time()
            pot = ssl.bar_aligned_potential(q["bar_frac"], q["S"], "sormani22")
            _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(q["phi_b_deg"])), -abs(q["Omega_abs"]))
            cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, ::2]), q["phi_b_deg"], sigma_dm, dat, ok)
            cl = cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - q["phi_b_deg"]).T
            blocks = [b for b in W if np.all(np.isfinite(cl[b][idx]))]
            vals.append(hart * sum(float(np.sum((W[b] @ (cd0[b][idx] - cl[b][idx])) ** 2)) for b in blocks))
            print(f"4. mode {tag} run {r:2d}: chi2 {vals[-1]:.1f} ({time.time() - t0:.0f} s)", flush=True)
        direct[tag] = vals
    a, b = np.array(direct["A"]), np.array(direct["B"])
    d, se = b.mean() - a.mean(), np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    print(f"4. direct: mode A chi2 {a.mean():.1f} +- {a.std(ddof=1):.1f} (sem {a.std(ddof=1) / np.sqrt(len(a)):.1f}); "
          f"mode B {b.mean():.1f} +- {b.std(ddof=1):.1f}; B - A = {d:+.1f} +- {se:.1f} in chi2 = {-d / temp:+.2f} +- "
          f"{se / temp:.2f} tempered logL (emulator predicted {lb - la:+.2f})")
    report.update(direct_A=a.tolist(), direct_B=b.tolist(), direct_dchi2=d, direct_se=se)
    json.dump(report, open(RESULTS_DIR / RUN / "bimodality_check.json", "w"), indent=2)


if __name__ == "__main__":
    main()

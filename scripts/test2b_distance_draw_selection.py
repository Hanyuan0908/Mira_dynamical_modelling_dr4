"""Test 2b: several candidate distances per star, then choose among them by fitting the observed kinematic field.

Each smeared star (sigma_DM = 0.265 mag) gets K candidate distance moduli at fixed quantiles of its measurement error,
DM_obs + sigma_DM z_k with z_k = Phi^-1((k - 1/2) / K) (equal prior probability each). With its exact l, b, pm and v_los
every candidate is a full initial condition; all N K orbits are integrated in a trial potential (accuracy 1e-5, 1 Gyr,
41 snapshots), for three trial potentials of the focused global search: #0 (truth), #17 (|Omega_b| 33.4, rejected on
perfect data) and #15 (S 0.82, bar_frac 1.50, on the S - bar_frac ridge).

Statistic (per 0.3 kpc ring and |z| layer, R < 6 kpc, blocks with >= 50 data stars): the least-squares Fourier
coefficients c = G^-1 H of <v_R>, <v_phi>, <v_R^2>, <v_phi^2> on (1, cos 2phi, sin 2phi, cos 4phi, sin 4phi); kept terms
<v_R> m = 2, 4 and the other three m = 0, 2, 4 (19 per block); chi2 with the joint bootstrap covariance of the data's
coefficients (300 bootstraps, Hartlap). <v^2> rather than sigma^2 so that the model sums are linear in the weights.
Model: every library sample (t >= 0.25 Gyr, true R < 8 kpc, weight 1 / tau) is observed with the distance errors
(fixed seed); per orbit and block the 35 sums (15 of G, 20 of H) are stored, so the model's block sums for any orbit
weights are a weighted total of these.

Fits (no regularisation, to measure the flexibility of each):
  pick1   one candidate per star carries the star's whole weight; greedy coordinate descent from a random start: stars
          in random order, the exact chi2 change of each of the K candidates (only the blocks its orbits touch change),
          switch to the best if it lowers chi2; sweeps until < 0.1% of stars switch (max MAX_SWEEPS).
  simplex continuous weights, softmax per star (sum 1 per star); L-BFGS on all N K logits with the analytic gradient.
  free    free positive weights w = exp(theta) (Schwarzschild-like); L-BFGS.
Reference: the truth library from the true initial conditions (one orbit per star), observed the same way.
Usage: python test2b_distance_draw_selection.py --stage integrate | fit
Files: ~/mira_dr4_runs/test2b/draws_K10_traj_NNN.npy, draws_K10_features_NNN.npz
Numbers: results/test2b/distance_draw_selection.json
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import h5py
from scipy.stats import norm
from scipy.optimize import minimize

from common import ROOT
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import observe
from test1_nbody_qmc_focused_run import DESIGN_CSV
from plot_distance_error_effect_on_fourier_coefficients import SMEARED, stars, observe_library_samples
from test2b_smeared_initial_conditions import OUT_DIR

K = 10
POINTS = (0, 17, 15)
MAX_SWEEPS = 30
NB, NK = fs.NB, fs.NK
PAIRS = fs.PAIRS                                    # 15 (k, l) pairs of G
NF = len(PAIRS) + 4 * NK                            # 35 sums per (orbit, block)
KEEP = [0 * NK + k for k in fs.NO_M0] + [f * NK + k for f in (1, 2, 3) for k in fs.ALL]   # 19 of the 20 coefficients
RESULTS = ROOT / "results" / "test2b" / "distance_draw_selection.json"


def params(idx):
    p = pd.read_csv(DESIGN_CSV).set_index("idx").loc[idx]
    return dict(Omega_abs=float(p.Omega_abs), phi_b_deg=float(p.phi_b_deg), S=float(p.S), bar_frac=float(p.bar_frac))


def candidate_ics():
    """(N K, 6) MW-frame initial conditions, star-major (orbit o = star * K + k), and the candidate distances."""
    with h5py.File(SMEARED, "r") as f:
        m = {k: f[k][:].astype(np.float64) for k in ("l", "b", "pml", "pmb", "vlos", "DM_obs")}
        sigma_dm = float(f.attrs["sigma_DM"])
    z = norm.ppf((np.arange(K) + 0.5) / K)
    d = 10.0 ** (0.2 * (m["DM_obs"][:, None] + sigma_dm * z[None, :]) - 2.0)
    rep = lambda a: np.repeat(a, K)
    ic = np.column_stack(observe.to_galactocentric(rep(m["l"]), rep(m["b"]), d.ravel(), rep(m["pml"]), rep(m["pmb"]),
                                                   rep(m["vlos"])))
    return ic, d, sigma_dm


def integrate_all():
    ic, _, _ = candidate_ics()
    ssl.ACCURACY = old.ACCURACY
    for idx in POINTS:
        path = OUT_DIR / f"draws_K{K}_traj_{idx:03d}.npy"
        if path.exists():
            continue
        p = params(idx)
        t0 = time.time()
        pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], "nbody")
        _, traj = ssl.integrate(pot, ssl.rotate_z(ic, np.radians(p["phi_b_deg"])), -abs(p["Omega_abs"]))
        np.save(path, traj)
        print(f"#{idx}: {len(ic)} orbits integrated in {time.time() - t0:.0f} s", flush=True)
        del traj


# ---------------------------------------------------------------------------------------------------------- statistics
def sums_per_row(w, wt):
    """Per sample: block index, the 35 sums' summands (G pairs, then <v_R>, <v_phi>, <v_R^2>, <v_phi^2> times f_k), mask."""
    b, X, v, m = fs.prepare(w)
    wt = wt[m]
    Y = np.column_stack([v[:, 0], v[:, 1], v[:, 0]**2, v[:, 1]**2])
    cols = [wt * X[:, k] * X[:, l] for k, l in PAIRS] + [wt * Y[:, a] * X[:, k] for a in range(4) for k in range(NK)]
    return b, np.column_stack(cols), m


def coefficients(S):
    """Block sums (..., 35) -> coefficients (..., 20) [field-major: <v_R>, <v_phi>, <v_R^2>, <v_phi^2>] and G^-1."""
    G = np.zeros(S.shape[:-1] + (NK, NK))
    for j, (k, l) in enumerate(PAIRS):
        G[..., k, l] = G[..., l, k] = S[..., j]
    H = S[..., len(PAIRS):].reshape(S.shape[:-1] + (4, NK))
    Ginv = np.linalg.inv(G)
    return np.einsum("...kl,...al->...ak", Ginv, H).reshape(S.shape[:-1] + (4 * NK,)), Ginv


def data_statistics(phi_deg, boot_w):
    """Data coefficients (NB, 20), Cinv of the kept terms with the Hartlap factor (NB, 19, 19), usable blocks."""
    w = ssl.rotate_z(stars(SMEARED), np.radians(phi_deg))
    b, F, m = sums_per_row(w, np.ones(len(w)))
    n = np.bincount(b, minlength=NB)
    ring = np.arange(NB) // fs.NL
    ok = (n >= fs.N_MIN_BLOCK) & (ring < fs.N_RING_DEFAULT)
    S = np.stack([np.bincount(b, weights=F[:, j], minlength=NB) for j in range(NF)], -1)
    c = np.full((NB, 4 * NK), np.nan)
    c[ok] = coefficients(S[ok])[0]
    boots = []
    for bw in boot_w:
        bwm = bw[m]
        Sb = np.stack([np.bincount(b, weights=F[:, j] * bwm, minlength=NB) for j in range(NF)], -1)
        cb = np.full_like(c, np.nan)
        cb[ok] = coefficients(Sb[ok])[0]
        boots.append(cb)
    boots = np.array(boots)[..., KEEP]
    dev = boots - boots.mean(0)
    cov = np.einsum("sbi,sbj->bij", dev, dev) / (len(boot_w) - 1)
    hartlap = (len(boot_w) - len(KEEP) - 2) / (len(boot_w) - 1)
    Cinv = np.zeros((NB, len(KEEP), len(KEEP)))
    Cinv[ok] = hartlap * np.linalg.inv(cov[ok])
    return c, Cinv, ok


def library_features(traj, times, phi_deg, sigma_dm, n_orbits):
    """Sparse per-(orbit, block) sums: orbit id, block id, F (rows, 35), for the observed library."""
    R = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])
    tau = np.clip(np.mean(R < ssl.R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
    keep = times >= ssl.T_MIN_GYR - 1e-9
    nt = int(keep.sum())
    rows_o, rows_b, rows_F = [], [], []
    for s in range(0, n_orbits, 50000):
        e = min(s + 50000, n_orbits)
        w = traj[s:e][:, keep].reshape(-1, 6).astype(np.float64)
        oid = np.repeat(np.arange(s, e), nt)
        wt = np.repeat(1.0 / tau[s:e], nt)
        sel = np.hypot(w[:, 0], w[:, 1]) < ssl.R_SEL
        wo, wto = observe_library_samples(w, wt, sigma_dm, phi_deg, seed=s)       # applies the same R < 8 cut
        oid = oid[sel]
        b, F, m = sums_per_row(wo, wto)
        key = oid[m] * NB + b
        u, inv = np.unique(key, return_inverse=True)
        rows_o.append(u // NB)
        rows_b.append(u % NB)
        rows_F.append(np.stack([np.bincount(inv, weights=F[:, j], minlength=len(u)) for j in range(NF)], -1)
                      .astype(np.float32))
    return np.concatenate(rows_o), np.concatenate(rows_b), np.concatenate(rows_F)


# ---------------------------------------------------------------------------------------------------------- chi2
class Objective:
    def __init__(self, ro, rb, F, c_data, Cinv, ok, n_orbits):
        order = np.argsort(ro, kind="stable")
        self.ro, self.rb, self.F = ro[order], rb[order], F[order].astype(np.float64)
        self.ptr = np.searchsorted(self.ro, np.arange(n_orbits + 1))
        self.cd, self.Cinv, self.ok = c_data[:, KEEP], Cinv, ok
        self.n_orbits = n_orbits
        self.N = int(ok.sum()) * len(KEEP)

    def block_sums(self, w):
        return np.stack([np.bincount(self.rb, weights=self.F[:, j] * w[self.ro], minlength=NB) for j in range(NF)], -1)

    def block_chi2(self, S, blocks):
        """chi2 of the given blocks for block sums S (..., len(blocks), 35)."""
        c, _ = coefficients(S)
        d = self.cd[blocks] - c[..., KEEP]
        return np.einsum("...i,...ij,...j->...", d, self.Cinv[blocks], d)

    def chi2(self, w):
        S = self.block_sums(w)
        bl = np.flatnonzero(self.ok)
        return float(self.block_chi2(S[bl], bl).sum())

    def chi2_and_grad_w(self, w):
        """chi2 and d chi2 / d w_o for every orbit (analytic, through c = G^-1 H)."""
        S = self.block_sums(w)
        bl = np.flatnonzero(self.ok)
        c, Ginv = coefficients(S[bl])
        d = self.cd[bl] - c[:, KEEP]
        chi = float(np.einsum("bi,bij,bj->", d, self.Cinv[bl], d))
        gk = np.zeros((len(bl), 4 * NK))
        gk[:, KEEP] = -2 * np.einsum("bij,bj->bi", self.Cinv[bl], d)         # d chi2 / d c
        g = gk.reshape(len(bl), 4, NK)
        eta = np.einsum("bkl,bal->bak", Ginv, g)                              # d chi2 / d H_a = G^-1 g_a
        cc = c.reshape(len(bl), 4, NK)
        D = -np.einsum("bak,bal->bkl", eta, cc)                               # d chi2 / d G (as a full matrix)
        dS = np.zeros((NB, NF))
        dS[bl, :len(PAIRS)] = np.stack([D[:, k, l] + (D[:, l, k] if k != l else 0) for k, l in PAIRS], -1)
        dS[bl, len(PAIRS):] = eta.reshape(len(bl), -1)
        grow = np.einsum("rj,rj->r", self.F, dS[self.rb])
        return chi, np.bincount(self.ro, weights=grow, minlength=self.n_orbits)

    # one candidate per star -----------------------------------------------------------------------------------------
    def pick_one(self, n_stars, rng, log):
        choice = rng.integers(0, K, n_stars)
        w = np.zeros(self.n_orbits)
        w[np.arange(n_stars) * K + choice] = 1.0
        S = self.block_sums(w)
        bchi = np.zeros(NB)
        bl = np.flatnonzero(self.ok)
        bchi[bl] = self.block_chi2(S[bl], bl)
        history = [float(bchi.sum())]
        for sweep in range(MAX_SWEEPS):
            t0, n_switch = time.time(), 0
            for i in rng.permutation(n_stars):
                lo, hi = self.ptr[i * K], self.ptr[(i + 1) * K]
                rb, F, ro = self.rb[lo:hi], self.F[lo:hi], self.ro[lo:hi] - i * K
                U = np.unique(rb)
                U = U[self.ok[U]]
                if len(U) == 0:
                    continue
                pos = np.searchsorted(U, rb)
                inU = (pos < len(U)) & (U[np.minimum(pos, len(U) - 1)] == rb)
                contrib = np.zeros((K, len(U), NF))
                np.add.at(contrib, (ro[inU], pos[inU]), F[inU])
                trial = S[U][None] - contrib[choice[i]][None] + contrib           # (K, |U|, 35)
                new = self.block_chi2(trial, U)                                   # (K, |U|)
                delta = new.sum(1) - bchi[U].sum()
                k = int(np.argmin(delta))
                if k != choice[i] and delta[k] < -1e-9:
                    S[U] = trial[k]
                    bchi[U] = new[k]
                    choice[i] = k
                    n_switch += 1
            history.append(float(bchi.sum()))
            log(f"    pick1 sweep {sweep + 1}: chi2/N {history[-1] / self.N:.4f}, switched {n_switch} "
                f"({100 * n_switch / n_stars:.2f}%), {time.time() - t0:.0f} s")
            if n_switch < 0.001 * n_stars:
                break
        return choice, history

    # continuous weights -----------------------------------------------------------------------------------------------
    def fit_simplex(self, n_stars, maxiter):
        def f(theta):
            th = theta.reshape(n_stars, K)
            e = np.exp(th - th.max(1, keepdims=True))
            w = e / e.sum(1, keepdims=True)
            chi, g = self.chi2_and_grad_w(w.ravel())
            g = g.reshape(n_stars, K)
            return chi, (w * (g - np.sum(w * g, 1, keepdims=True))).ravel()
        r = minimize(f, np.zeros(n_stars * K), jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter, maxfun=2 * maxiter))
        th = r.x.reshape(n_stars, K)
        e = np.exp(th - th.max(1, keepdims=True))
        return e / e.sum(1, keepdims=True), r

    def fit_free(self, maxiter):
        def f(theta):
            w = np.exp(np.clip(theta, -30, 30))
            chi, g = self.chi2_and_grad_w(w)
            return chi, g * w
        r = minimize(f, np.zeros(self.n_orbits), jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter, maxfun=2 * maxiter))
        return np.exp(np.clip(r.x, -30, 30)), r


def fit_all(maxiter):
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    ic, d_cand, sigma_dm = candidate_ics()
    n_stars = len(d_cand)
    with h5py.File(SMEARED, "r") as f:
        d_true, d_obs = f["d_true"][:], f["d"][:]
    rng_b = np.random.default_rng(fs.SEED)
    boot_w = rng_b.multinomial(n_stars, np.full(n_stars, 1 / n_stars), size=fs.N_BOOT).astype(float)
    res = dict(K=K, z_nodes=[round(float(z), 4) for z in norm.ppf((np.arange(K) + 0.5) / K)], points={})
    # how far the observed distances and the nearest candidate are from the truth
    res["observed_frac_error_median"] = round(float(np.median(np.abs(d_obs / d_true - 1))), 4)
    res["nearest_candidate_frac_error_median"] = round(float(np.median(np.min(np.abs(d_cand / d_true[:, None] - 1), axis=1))), 4)
    log = lambda s: print(s, flush=True)
    for idx in POINTS:
        p = params(idx)
        log(f"#{idx} {p}")
        t0 = time.time()
        feat = OUT_DIR / f"draws_K{K}_features_{idx:03d}.npz"
        if feat.exists():
            z = np.load(feat)
            ro, rb, F = z["ro"], z["rb"], z["F"]
        else:
            traj = np.load(OUT_DIR / f"draws_K{K}_traj_{idx:03d}.npy", mmap_mode="r")
            ro, rb, F = library_features(traj, times, p["phi_b_deg"], sigma_dm, n_stars * K)
            np.savez(feat, ro=ro, rb=rb, F=F)
            del traj
        c_data, Cinv, ok = data_statistics(p["phi_b_deg"], boot_w)
        obj = Objective(ro, rb, F, c_data, Cinv, ok, n_stars * K)
        log(f"  features: {len(ro):,} (orbit, block) rows, data: {int(ok.sum())} blocks, N_data {obj.N}; {time.time() - t0:.0f} s")
        out = dict(params=p, N_data=obj.N)
        out["uniform_weights"] = round(obj.chi2(np.ones(n_stars * K)) / obj.N, 4)
        # reference: the library from the true initial conditions (true positions) of this potential, for #0 only
        if idx == 0:
            tr = np.load(old.ARR_DIR / "traj_000.npy")
            ro0, rb0, F0 = library_features(tr, times, p["phi_b_deg"], sigma_dm, n_stars)
            out["true_initial_conditions"] = round(Objective(ro0, rb0, F0, c_data, Cinv, ok, n_stars).chi2(np.ones(n_stars)) / obj.N, 4)
        t0 = time.time()
        choice, hist = obj.pick_one(n_stars, np.random.default_rng(1), log)
        dsel = d_cand[np.arange(n_stars), choice]
        out["pick1"] = dict(start=round(hist[0] / obj.N, 4), final=round(hist[-1] / obj.N, 4), sweeps=len(hist) - 1,
                            seconds=round(time.time() - t0), history=[round(h / obj.N, 4) for h in hist],
                            chosen_frac_error_median=round(float(np.median(np.abs(dsel / d_true - 1))), 4),
                            chosen_node_histogram=np.bincount(choice, minlength=K).tolist())
        log(f"  pick1: {out['pick1']['start']} -> {out['pick1']['final']}")
        t0 = time.time()
        wts, r = obj.fit_simplex(n_stars, maxiter)
        dmean = np.sum(wts * d_cand, 1)
        out["simplex"] = dict(final=round(r.fun / obj.N, 4), iterations=int(r.nit), seconds=round(time.time() - t0),
                              message=str(r.message), weighted_mean_frac_error_median=round(float(np.median(np.abs(dmean / d_true - 1))), 4),
                              effective_candidates_per_star_median=round(float(np.median(1 / np.sum(wts**2, 1))), 3))
        log(f"  simplex: {out['simplex']['final']} after {r.nit} iterations")
        t0 = time.time()
        wf, r = obj.fit_free(maxiter)
        out["free"] = dict(final=round(r.fun / obj.N, 4), iterations=int(r.nit), seconds=round(time.time() - t0),
                           message=str(r.message), weight_spread_16_50_84=[round(float(x), 4) for x in
                                                                           np.percentile(wf / wf.mean(), [16, 50, 84])])
        log(f"  free: {out['free']['final']} after {r.nit} iterations")
        res["points"][str(idx)] = out
        RESULTS.parent.mkdir(parents=True, exist_ok=True)
        with open(RESULTS, "w") as fh:
            json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("integrate", "fit"), required=True)
    ap.add_argument("--maxiter", type=int, default=300)
    args = ap.parse_args()
    if args.stage == "integrate":
        integrate_all()
    else:
        fit_all(args.maxiter)


if __name__ == "__main__":
    main()

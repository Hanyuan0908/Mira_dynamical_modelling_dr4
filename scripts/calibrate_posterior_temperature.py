"""How wide should the posterior be? Measure the Delta-chi2 fluctuation and calibrate by coverage.

The likelihood is tempered as logLs = -chi2 / (2 sqrt(2 N_data)). The justification offered for sqrt(2N) is
that it is the standard deviation of chi2 for N degrees of freedom. That is the scatter of the *absolute*
chi2 over data realisations; it is not the scatter of the *difference* between two models on the same data,
which is what a posterior width depends on. For two models with predictions m1, m2 and data covariance C,

    Delta = sqrt((m1-m2)^T C^-1 (m1-m2)),   E[dchi2] = +/- Delta^2,   std[dchi2] = 2 Delta

so the noise on dchi2 is set by the model separation, not by N. Part 1 measures 2*Delta for models near the
peak and compares it with sqrt(2N).

Part 2 calibrates the temperature by coverage, which is the thing that actually matters: generate synthetic
data from the truth model, build the posterior at each temperature from the 256 Sobol points by importance
weighting, and count how often the truth falls inside the 68% interval. The temperature that gives 68%
coverage is the faithful one.

Everything is done at the reference bar angle: c(phi) = T c(27) with T orthogonal, so rotating the library
back by T^T and using the fixed C(27) is exact (and much cheaper than rotating per point).

This uses only saved coefficients - no integrations.
Usage: python calibrate_posterior_temperature.py --run test2_t7_sf_10pct --variant sf
"""
import argparse, json
import numpy as np, pandas as pd
from common import ROOT, RUNS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test1_nbody_qmc_focused_fourier_stars import rotation
from test2b_qmc_distance_error_run import scaled_design

SCHEME = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="test2_t7_sf_10pct")
    ap.add_argument("--variant", default="sf")
    ap.add_argument("--n-real", type=int, default=400)
    args = ap.parse_args()

    meta = json.load(open(ROOT / "results" / args.run / "run_meta.json"))
    design, _ = scaled_design(meta["galaxy"])
    out = RUNS_DIR / args.run
    s0 = np.load(out / "fourier_000.npz")
    cov0 = s0["cov_data_joint"].reshape(fs.NB, 6 * fs.NK, 6 * fs.NK)
    ok = s0["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK
    n_ring, spec = fs.SCHEMES[SCHEME]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    h = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)

    # whitening per block at the reference angle
    blocks, Linv = [], []
    for r in range(n_ring):
        for z in range(fs.NL):
            i = r * fs.NL + z
            if ok[i]:
                blocks.append(i)
                Linv.append(np.linalg.inv(np.linalg.cholesky(cov0[i][np.ix_(idx, idx)])))

    E, P, kinds = [], [], []
    for i in design.idx:
        i = int(i)
        f = out / f"fourier_{i:03d}.npz"
        if not f.exists():
            continue
        cl = np.load(f)[f"coef_lib_{args.variant}"].reshape(fs.NB, -1)
        p = design.set_index("idx").loc[i]
        T = rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
        clr = cl @ T                                        # T^T cl, back to the reference angle
        E.append(np.concatenate([Linv[b] @ clr[j][idx] for b, j in zip(range(len(blocks)), blocks)]))
        P.append([p.Omega_abs, p.phi_b_deg, p.S, p.bar_frac]); kinds.append(p.kind)
    E = np.array(E); P = np.array(P); kinds = np.array(kinds)
    nd = len(blocks) * len(idx)
    truth = np.flatnonzero(kinds == "truth")[0]
    qmc = np.flatnonzero(kinds == "qmc")
    print(f"{args.run} / {args.variant} / scheme B: {len(E)} models, N_data {nd}, sqrt(2N) = {np.sqrt(2*nd):.1f}")

    # --- Part 1: how much does Delta-chi2 actually fluctuate? ---
    d = np.linalg.norm(E[qmc] - E[truth], axis=1)           # Delta for each model vs the truth
    near = qmc[np.argsort(d)[:20]]
    print(f"\n  20 models closest to the truth: Delta = {d[np.argsort(d)[:20]].min():.1f} to "
          f"{d[np.argsort(d)[:20]].max():.1f}")
    print(f"  so std(Delta chi2) = 2*Delta = {2*d[np.argsort(d)[:20]].min():.1f} to "
          f"{2*d[np.argsort(d)[:20]].max():.1f}")
    print(f"  the tempering instead assumes the scale is sqrt(2N) = {np.sqrt(2*nd):.1f}")

    # --- Part 2: coverage calibration ---
    rng = np.random.default_rng(0)
    temps = [1.0, 2.0, 4.0, 8.0, 16.0, np.sqrt(2*nd), 64.0]
    cover = {T: 0 for T in temps}
    width = {T: [] for T in temps}
    om = P[qmc, 0]; om_true = P[truth, 0]
    for _ in range(args.n_real):
        dsyn = E[truth] + rng.standard_normal(E.shape[1])    # synthetic data, whitened units
        chi2 = h * np.sum((dsyn - E[qmc])**2, axis=1)
        for T in temps:
            w = np.exp(-(chi2 - chi2.min()) / (2 * T))
            w /= w.sum()
            o = np.argsort(om); c = np.cumsum(w[o])
            lo, hi = np.interp([0.16, 0.84], c, om[o])
            cover[T] += (lo <= om_true <= hi); width[T].append(0.5 * (hi - lo))
    print(f"\n  {'temperature':>12s} {'coverage of 68% interval':>26s} {'median width':>14s}")
    for T in temps:
        tag = "  <- current sqrt(2N)" if abs(T - np.sqrt(2*nd)) < 1e-6 else ""
        print(f"  {T:12.1f} {100*cover[T]/args.n_real:25.0f}% {np.median(width[T]):14.2f}{tag}")
    print(f"\n  (target 68%; {args.n_real} synthetic data realisations drawn from the truth model)")


if __name__ == "__main__":
    main()

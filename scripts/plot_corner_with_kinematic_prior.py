"""Corner plot with a second posterior added: the same likelihood under a kinematic prior on phi_b.

The dynamical model constrains the bar angle weakly, and that weakness propagates: phi_b is
degenerate with bar_frac, so a poorly determined angle drags the bar strength with it. But phi_b is
also measurable independently, from the m = 2 phase of <|v_R|/v_tot> - no orbit library, no
potential, just the kinematics - and that measurement is far tighter than the model's own.

So this takes the finished posterior and multiplies in a Gaussian prior

    phi_b ~ N(phi_kin, sigma_kin)

where phi_kin is the kinematic angle inverted through the Auriga calibration band and sigma_kin is
that band's half-width divided by its slope (+-1.88 deg for absvR at self:0.3).

Done by importance weighting the existing chain rather than re-sampling: the target is the same
likelihood times a bounded prior, so w_i = exp(-0.5 ((phi_i - phi_kin)/sigma_kin)^2) is exact. The
effective sample size is reported, since a prior much narrower than the posterior would thin the
chain too far to trust.

CIRCULARITY: the calibration is built FROM the 13 Auriga galaxies. For an Auriga run the prior is
therefore derived leave-one-out, from the other twelve. For an N-body run the full Auriga
calibration is used and the prior is genuinely blind.

Usage: python plot_corner_with_kinematic_prior.py --run test2_t5_sf_12p5pct
Figure: figures/<run>_corner_kinprior_<column>.png
Table:  results/<run>/kinematic_prior_<column>.json
"""
import argparse
import json

import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge
from scan_bar_angle_windows import phi_from, CACHE

SCHEME_B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
NATIVE_ANGLE = 27.0                 # every mock is built with its bar at BAR_ANGLE_DEG


def measure_from_mock(mock_name, rule="self:0.3", est="absvR"):
    """Measured bar angle from a mock's own observed stars - one dataset, as for real data.

    The ring cache holds angles measured on the FIDUCIAL uniform-error mocks. A run whose data is
    a different mock (e.g. the per-star distance-error test) must be measured on that mock, or the
    prior describes data that was never modelled. The calibration it is inverted through stays the
    external Auriga relation either way.
    """
    import h5py
    from common import MOCK_DIR
    from bar_angle_vecsum_verification import phase
    with h5py.File(MOCK_DIR / mock_name, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    win = f"kin{rule.split(':')[1]}" if rule.startswith("self:") else rule
    # phase() returns the bar axis, which is MINUS the angle from the Sun-GC line
    return -float(phase(w, est, win))


def kinematic_prior(run_galaxy, rule="self:0.3", est="absvR", measured=None):
    """(phi_kin, sigma_kin) for one galaxy, from the Auriga calibration band.

    Auriga galaxies get a leave-one-out calibration so the prior is not derived from themselves.
    """
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    ang = np.array(angles)
    ei = ["vR", "absvR"].index(est)
    if run_galaxy not in snaps:
        raise SystemExit(f"{run_galaxy} is not in the ring cache")
    P = np.array([[np.nanmean([phi_from(Z[g, a, s, ei], rule, Z[g, a, s, 0], est)
                               for s in range(Z.shape[2])])
                   for a in range(len(ang))] for g in range(len(snaps))])
    iau = [i for i, g in enumerate(snaps) if g.startswith("au")]
    cal = [i for i in iau if snaps[i] != run_galaxy]          # leave-one-out for Auriga
    loo = len(cal) < len(iau)
    A = P[cal]
    med = np.nanmedian(A, 0)
    width = float(np.mean(np.nanpercentile(A, 84, 0) - np.nanpercentile(A, 16, 0)))
    slope = float(np.polyfit(ang, med, 1)[0])
    meas = float(P[snaps.index(run_galaxy)][list(ang).index(NATIVE_ANGLE)]) if measured is None \
        else float(measured)
    inv = (np.interp(meas, med[::-1], ang[::-1]) if med[0] > med[-1]
           else np.interp(meas, med, ang))
    return float(inv), float(0.5 * width / slope), dict(
        measured=meas, calibrators=len(cal), leave_one_out=loo, band=width, slope=slope,
        rule=rule, est=est)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--column", default=f"logLs_{SCHEME_B}_sfratio")
    ap.add_argument("--rule", default="self:0.3")
    ap.add_argument("--est", default="absvR")
    ap.add_argument("--loo-top", type=int, default=32)
    ap.add_argument("--measure-from-mock", action="store_true",
                    help="measure the kinematic angle on this run's own mock instead of taking it "
                         "from the ring cache (required for any run whose data is not the "
                         "fiducial uniform-error mock)")
    args = ap.parse_args()
    set_style()

    run = ge.use_run(args.run)
    TRUTH, pre = run["truth"], run["prefix"]
    meta = json.load(open(RESULTS_DIR / args.run / "run_meta.json"))
    meas = measure_from_mock(meta["mock"], args.rule, args.est) if args.measure_from_mock else None
    phi_kin, sig_kin, info = kinematic_prior(meta["galaxy"], args.rule, args.est, measured=meas)
    if meas is not None:
        info["measured_from"] = meta["mock"]
    print(f"  {args.run}: galaxy {meta['galaxy']}")
    print(f"    kinematic angle {info['measured']:.2f} deg measured -> phi_b = {phi_kin:.2f} "
          f"+- {sig_kin:.2f} deg  ({'leave-one-out, ' if info['leave_one_out'] else ''}"
          f"{info['calibrators']} calibrators)")

    df, tag = ge.load_table(args.column, None, None)
    tr = df[df.kind == "qmc"]
    u, y = ge.to_unit(tr), tr["logL"].values
    emu = ge.Emulator(u, y)

    def logpost(q):
        return -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0]
    from scipy.optimize import minimize
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    nw = 32
    p0 = np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((nw, 4)), 0, 1)
    s = emcee.EnsembleSampler(nw, 4, logpost)
    s.random_state = np.random.RandomState(2).get_state()
    s.run_mcmc(p0, 3000, progress=False)
    post = ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True))

    # the prior enters as an importance weight; the likelihood is untouched
    w = np.exp(-0.5 * ((post[:, 1] - phi_kin) / sig_kin) ** 2)
    ess = float(w.sum() ** 2 / np.sum(w ** 2))
    print(f"    chain {len(post)} samples -> effective {ess:.0f} under the prior")

    def pct(a, weights=None):
        if weights is None:
            return np.percentile(a, [16, 50, 84], axis=0)
        o = np.argsort(a, axis=0)
        out = []
        for k in range(a.shape[1]):
            v, ww = a[o[:, k], k], weights[o[:, k]]
            c = np.cumsum(ww) / ww.sum()
            out.append(np.interp([0.16, 0.50, 0.84], c, v))
        return np.array(out).T

    q0, q1 = pct(post), pct(post, w)
    res = dict(run=args.run, column=args.column, galaxy=meta["galaxy"],
               phi_kin=phi_kin, sigma_kin=sig_kin, prior_info=info, ess=ess,
               truth=TRUTH,
               flat={k: dict(p16=float(a), p50=float(b), p84=float(c))
                     for k, a, b, c in zip(ge.KEYS, *q0)},
               with_prior={k: dict(p16=float(a), p50=float(b), p84=float(c))
                           for k, a, b, c in zip(ge.KEYS, *q1)})
    out = RESULTS_DIR / args.run / f"kinematic_prior_{args.column}.json"
    json.dump(res, open(out, "w"), indent=2)

    lab = [r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", r"$\phi_b$ [deg]", "$S$", "bar_frac"]
    rng = [(min(post[:, i].min(), TRUTH[k]), max(post[:, i].max(), TRUTH[k]))
           for i, k in enumerate(ge.KEYS)]
    fig = corner.corner(post, labels=lab, range=rng, color=BLUE_RAMP[5], plot_datapoints=False,
                        levels=(0.39, 0.68, 0.95), truths=[TRUTH[k] for k in ge.KEYS],
                        truth_color=INK_2, hist_kwargs=dict(density=True))
    corner.corner(post, weights=w, labels=lab, range=rng, color=ACCENT, fig=fig,
                  plot_datapoints=False, levels=(0.39, 0.68, 0.95),
                  hist_kwargs=dict(density=True))
    for i, k in enumerate(ge.KEYS):
        a, b, c = q1[:, i]
        fig.axes[i * 5].set_title(f"{k}: {b:.2f}$^{{+{c-b:.2f}}}_{{-{b-a:.2f}}}$",
                                  color=ACCENT, fontsize=9)
    fig.suptitle(f"{pre}: flat prior (blue) vs kinematic prior on $\\phi_b$ (orange)\n"
                 f"$\\phi_b$ = {phi_kin:.1f} $\\pm$ {sig_kin:.1f}$^\\circ$ from "
                 f"{args.est} / {args.rule}"
                 f"{' (leave-one-out)' if info['leave_one_out'] else ' (blind)'}; "
                 f"truth in grey; effective samples {ess:.0f}",
                 color=INK, fontsize=11, y=1.02)
    savefig(fig, f"{pre}_corner_kinprior_{args.column}")

    print(f"    {'param':10s} {'truth':>8s} | {'flat':>22s} | {'with prior':>22s}")
    for i, k in enumerate(ge.KEYS):
        a0, b0, c0 = q0[:, i]
        a1, b1, c1 = q1[:, i]
        print(f"    {k:10s} {TRUTH[k]:8.2f} | {b0:8.2f} +{c0-b0:5.2f} -{b0-a0:5.2f} | "
              f"{b1:8.2f} +{c1-b1:5.2f} -{b1-a1:5.2f}")
    print(f"  wrote {rel(out)}")


if __name__ == "__main__":
    main()

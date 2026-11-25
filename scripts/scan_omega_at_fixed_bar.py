"""Do the predicted kinematics react to the pattern speed at a FIXED bar? A direct test, no emulator.

For each case, (phi_b, S, bar_frac) are held at the best evaluated point of its search and only |Omega_b| is varied.
At every |Omega_b| the stars are integrated (1 Gyr, accuracy 1e-5, thinned to 50 Myr) and scored exactly as the
pipeline does (rescore_empirical_selection.library_fit_empirical: sfratio, K = 10; scheme B; data fitted once at
27 deg and rotated). Recorded per |Omega_b|:
  chi2 against the data                        -> does the DATA prefer one pattern speed?
  D = chi2(library(Omega) vs library(Omega_ref) | data covariance)
                                               -> do the MODEL moments change, in units of the data errors?
and the library <v_R> m = 2 amplitude per ring (|z| < 1) for the plot. Omega_ref is the scan value nearest the
search's best point.
Cases:
  osarg_h24  OSARG data, Hunter+24 family, best point of osarg_zhang2024_mw_hunter24
  t7_own     t = 7 mock (12.5%), its own family, best point of test2_t7_sf_12p5pct
  t7_h24     t = 7 mock (12.5%), Hunter+24 family, best point of test2_t7_sf_12p5pct_mwH24_osargdesign

Usage: python scan_omega_at_fixed_bar.py [--cases osarg_h24 t7_h24 t7_own] [--plot-only]
Output: results/omega_scan_fixed_bar_<case>.npz; figures/test2_omega_scan_at_fixed_bar_osarg_vs_mock.png
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from test2b_qmc_distance_error_run import load_design
from rescore_empirical_selection import data_lbdR, library_fit_empirical

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
CASES = {
    "osarg_h24": dict(label="OSARG data, Hunter+24 family", run="osarg_zhang2024_mw_hunter24", galaxy="mw",
                      mock="osarg_zhang2024_observed_R8.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "osarg_h24w": dict(label="OSARG data, Hunter+24 family, wide bar box (S 0.4-1.7, bar_frac 0.1-2)", run="osarg_zhang2024_mw_hunter24_widebar", galaxy="mw",
                      mock="osarg_zhang2024_observed_R8.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "t7_h24": dict(label="t = 7 mock, Hunter+24 family", run="test2_t7_sf_12p5pct_mwH24_osargdesign", galaxy="mw",
                   mock="nbody_mock_selection_sigmaDM0.271_N50k.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "osarg_s22": dict(label="OSARG data, Sormani+22 family", run="osarg_zhang2024_sormani22", galaxy="sormani22",
                      mock="osarg_zhang2024_observed_R8.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "osarg_s22w": dict(label="OSARG data, Sormani+22 family, wide bar box (S 0.4-1.7, bar_frac 0.1-2)", run="osarg_zhang2024_sormani22_widebar", galaxy="sormani22",
                      mock="osarg_zhang2024_observed_R8.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "h24mock": dict(label="steady-state OSARG-like mock, Hunter+24 family (truth 32.5)", run="mw_hunter24_mock_from_osarg_ics",
                    galaxy="mw", mock="mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5", omegas=np.arange(15.0, 45.01, 2.5)),
    "t7_own": dict(label="t = 7 mock, own family (truth 25.20)", run="test2_t7_sf_12p5pct", galaxy="nbody",
                   mock="nbody_mock_selection_sigmaDM0.271_N50k.hdf5", omegas=np.arange(15.0, 40.01, 2.5)),
}
RINGS = [r for r in range(fs.N_RING_DEFAULT) if 0.5 * (fs.RING_EDGES[r] + fs.RING_EDGES[r + 1]) <= 5.0]


def run_case(key, s_override=None):
    c = CASES[key]
    meta = json.load(open(RESULTS_DIR / c["run"] / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    ev = pd.read_csv(RESULTS_DIR / c["run"] / "evals.csv").set_index("idx")
    ibest = int(ev[f"logLs_{B}_sfratio"].idxmax())
    p = design.loc[ibest].copy()
    out_key = key
    if s_override is not None:            # same phi_b and bar_frac, a different bar length
        p["S"] = s_override
        out_key = f"{key}_S{s_override:g}"
    sigma_dm = float(meta["sigma_DM"])
    mock = MOCK_DIR / c["mock"]
    dat = data_lbdR(mock)
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    T = fs.rotation(float(p.phi_b_deg) - BAR_ANGLE_DEG)
    cd = (cd0.reshape(fs.NB, -1) @ T.T).reshape(fs.NB, 6, fs.NK)
    cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
    print(f"{key}: best point idx {ibest}: Omega {p.Omega_abs:.2f}, phi {p.phi_b_deg:.2f}, S {p.S:.3f}, "
          f"f {p.bar_frac:.3f}; {len(w):,} stars, sigma_DM {sigma_dm}", flush=True)
    ssl.ACCURACY = old.ACCURACY
    pot = ssl.bar_aligned_potential(float(p.bar_frac), float(p.S), c["galaxy"])
    ic = ssl.rotate_z(w, np.radians(float(p.phi_b_deg)))
    libs, chi2, nd = [], [], None
    for om in c["omegas"]:
        t0 = time.time()
        _, traj = ssl.integrate(pot, ic, -abs(om))
        traj = np.ascontiguousarray(traj[:, ::2])            # 50 Myr, as the searches
        cl, _ = library_fit_empirical(traj, float(p.phi_b_deg), sigma_dm, dat, ok)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        x2, nd = fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[B])
        libs.append(cl)
        chi2.append(x2)
        print(f"  Omega {om:5.1f}: chi2/N {x2 / nd:.3f}  ({time.time() - t0:.0f} s)", flush=True)
    libs = np.array(libs)
    iref = int(np.argmin(np.abs(c["omegas"] - p.Omega_abs)))
    D = []
    for cl in libs:
        valid = ok & np.all(np.isfinite(cl), axis=(1, 2)) & np.all(np.isfinite(libs[iref]), axis=(1, 2))
        D.append(fs.scheme_chi2(libs[iref], cl, cov, valid, fs.SCHEMES[B])[0])
    np.savez(RESULTS_DIR / f"omega_scan_fixed_bar_{out_key}.npz", omegas=c["omegas"], chi2=np.array(chi2), D=np.array(D),
             ndata=nd, coef_lib=libs, coef_data=cd, cov_data=cov, best=np.array([p.Omega_abs, p.phi_b_deg, p.S,
                                                                                 p.bar_frac]), iref=iref)


NAMES = {0: "<v_R> m2,4", 1: "<v_phi> m0,2,4", 3: "sig_R^2 m2,4", 4: "sig_phi^2 m2,4"}


def field_breakdown(key, step=10.0):
    """Tempered model change and data change per moment, R < 3 | R >= 3 kpc, for a +step change of |Omega_b| from the
    reference (field-by-field covariance: an apportioning approximation of the joint chi2)."""
    s = np.load(RESULTS_DIR / f"omega_scan_fixed_bar_{key}.npz")
    om, iref, cov, cd, L = s["omegas"], int(s["iref"]), s["cov_data"], s["coef_data"], s["coef_lib"]
    temp = 2 * np.sqrt(2 * int(s["ndata"]))
    j = int(np.argmin(np.abs(om - (om[iref] + step))))
    nr, spec = fs.SCHEMES[B]
    inner = [r for r in range(nr) if fs.RING_EDGES[r + 1] <= 3.0 + 1e-9]
    outer = [r for r in range(nr) if fs.RING_EDGES[r] >= 3.0 - 1e-9]

    def chi2_field(a, b, f, rings):
        tot, idx = 0.0, [f * fs.NK + k for k in spec[f]]
        for r in rings:
            for z in range(fs.NL):
                i = r * fs.NL + z
                if not (np.all(np.isfinite(a[i])) and np.all(np.isfinite(b[i]))):
                    continue
                d, C = (a[i] - b[i]).reshape(-1)[idx], cov[i][np.ix_(idx, idx)]
                if np.all(np.isfinite(C)):
                    tot += (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1) * float(d @ np.linalg.solve(C, d))
        return tot

    print(f"{key}: Omega {om[iref]:.1f} -> {om[j]:.1f}, tempered, R<3 | R>=3 kpc")
    for f in spec:
        m = [chi2_field(L[j], L[iref], f, rr) / temp for rr in (inner, outer)]
        d = [(chi2_field(cd, L[j], f, rr) - chi2_field(cd, L[iref], f, rr)) / temp for rr in (inner, outer)]
        print(f"   {NAMES[f]:16s} model change {m[0]:6.2f} | {m[1]:6.2f}    data change {d[0]:6.2f} | {d[1]:6.2f}")


def plot(keys, tag=""):
    set_style()
    fig, axes = plt.subplots(2, len(keys), figsize=(5.4 * len(keys), 9.0), squeeze=False,
                             gridspec_kw=dict(hspace=0.32, wspace=0.28))
    R = np.array([0.5 * (fs.RING_EDGES[r] + fs.RING_EDGES[r + 1]) for r in RINGS])
    for j, key in enumerate(keys):
        s = np.load(RESULTS_DIR / f"omega_scan_fixed_bar_{key}.npz")
        om, chi2, D, nd, best = s["omegas"], s["chi2"], s["D"], int(s["ndata"]), s["best"]
        temp = 2 * np.sqrt(2 * nd)
        a = axes[0, j]
        cmap = plt.get_cmap("viridis")
        for k, o in enumerate(om):
            amp = [np.hypot(*s["coef_lib"][k, r * fs.NL, 0, 1:3]) for r in RINGS]
            a.plot(R, amp, "-", color=cmap(k / (len(om) - 1)), lw=1.4, label=f"{o:g}" if k % 2 == 0 else None)
        cd, cov = s["coef_data"], s["cov_data"]
        damp = [np.hypot(*cd[r * fs.NL, 0, 1:3]) for r in RINGS]
        derr = [np.sqrt(0.5 * (cov[r * fs.NL, 1, 1] + cov[r * fs.NL, 2, 2])) for r in RINGS]
        a.errorbar(R, damp, yerr=derr, fmt="o", color=INK, ms=5, capsize=2, label="data", zorder=5)
        a.set(xlabel="R [kpc]  (|z| < 1)", ylabel=r"$\langle v_R\rangle$ quadrupole amplitude [km s$^{-1}$]",
              ylim=(0, 45))
        a.set_title(f"{CASES[key.split('_S')[0]]['label']}\nbar fixed: φ_b {best[1]:.1f}, S {best[2]:.2f}, bar_frac {best[3]:.2f}",
                    color=INK, fontsize=10)
        a.legend(title=r"$|\Omega_b|$", fontsize=7.5, title_fontsize=8, ncol=2, frameon=False, loc="upper right")
        b = axes[1, j]
        b.plot(om, -(chi2 - chi2.min()) / temp, "o-", color=INK, lw=1.8, label="data: χ² against the data")
        b.plot(om, -D / temp, "s--", color=ACCENT, lw=1.6, label="models: change of the prediction")
        b.axvline(best[0], color=BLUE_RAMP[3], lw=1, ls=":", label=f"search best {best[0]:.1f}")
        if key.startswith("t7_own"):
            b.axvline(25.20, color=BLUE_RAMP[5], lw=1.2, ls="--", label="truth 25.20")
        for lev in (-0.5, -2, -4.5):
            b.axhline(lev, color=INK_2, lw=0.5, ls=":")
        b.set(xlabel=r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", ylabel="tempered Δ log L  (χ² / 2√(2N))",
              ylim=(-40, 1))
        b.legend(fontsize=8, frameon=False, loc="lower left")
        print(f"{key}: data dlogL at each Omega: " + " ".join(f"{x:.1f}" for x in -(chi2 - chi2.min()) / temp))
        print(f"{' ' * len(key)}  model change D:      " + " ".join(f"{x:.1f}" for x in -D / temp))
    fig.suptitle("Pattern speed varied with the bar held fixed: do the predicted kinematics change, and does the data "
                 "notice?", color=INK, y=0.995, fontsize=12)
    savefig(fig, "test2_omega_scan_at_fixed_bar_osarg_vs_mock" + tag)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", nargs="+", default=["osarg_h24", "t7_h24", "t7_own"])
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--S", type=float, default=None, help="override the bar length S (phi_b and bar_frac stay at the "
                    "best point); outputs are tagged _S<value>")
    ap.add_argument("--figure-tag", default="", help="suffix for the figure name")
    args = ap.parse_args()
    keys = [k if args.S is None else f"{k}_S{args.S:g}" for k in args.cases]
    if not args.plot_only:
        for k in args.cases:
            run_case(k, args.S)
    plot(keys, args.figure_tag)
    for k in keys:
        field_breakdown(k)


if __name__ == "__main__":
    main()

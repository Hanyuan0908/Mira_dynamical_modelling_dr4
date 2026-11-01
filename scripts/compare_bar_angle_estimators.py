"""Two ways to read the bar angle off the kinematics, compared against distance error.

  vR quadrupole   the m = 2 node of <v_R>(phi), which for a bar lies along the major axis. This is what
                  plot_bar_angle_inversion.py uses; a 10% distance error rotates it by about +9 deg.

  |vR/vtot|       the m = 2 phase of <|v_R| / |v|>(phi). Zhang et al. (2024) propose this as a bar-orientation
                  proxy and argue it is more robust to heliocentric distance uncertainty than a <v_R> map;
                  the GaiaNIR follow-up (arXiv:2604.18329) quotes 1-7 deg accuracy with a slight over-estimate.

Both are fitted the same way - per-star least squares of the quantity against 1, cos2phi, sin2phi, cos4phi,
sin4phi in each (ring, |z| layer) block, so neither is biased by how the stars are spread in azimuth - and both
get the same zero-distance-error calibration, so the comparison is like for like. What is being compared is how
far each drifts as the distance error grows, and how much it scatters.

Usage: python compare_bar_angle_estimators.py [--galaxy nbody] [--n-real 15]
Figure: figures/test2_bar_angle_estimators_<prefix>.png
"""
import argparse
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy
import observe, selection, steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs

FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
R_FIT = (1.0, 4.0)


def fourier_phase(w_mw, which, layer=0):
    """m = 2 amplitude and phase [deg, MW frame] of the chosen quantity, median over rings in R_FIT.

    which = "vR"   : <v_R>(phi); the quadrupole NODE is 45 deg behind the m=2 maximum (bar major axis).
            "absvR": <|v_R|/|v|>(phi); the phase of the m=2 maximum is returned as is and calibrated later.
    """
    b, X, v, _ = fs.prepare(ssl.rotate_z(w_mw, np.radians(BAR_ANGLE_DEG)))
    n = np.bincount(b, minlength=fs.NB)
    ok = n >= fs.N_MIN_BLOCK
    if which == "vR":
        c = fs.fit(b, X, v, np.ones(len(b)), ok)[:, 0, :]
    else:
        q = np.abs(v[:, 0]) / np.maximum(np.linalg.norm(v, axis=1), 1e-9)
        G = np.zeros((fs.NB, fs.NK, fs.NK))
        for k in range(fs.NK):
            for l in range(fs.NK):
                G[:, k, l] = np.bincount(b, weights=X[:, k] * X[:, l], minlength=fs.NB)
        H = np.stack([np.bincount(b, weights=X[:, k] * q, minlength=fs.NB) for k in range(fs.NK)], -1)
        Gi = np.full_like(G, np.nan)
        Gi[ok] = np.linalg.inv(G[ok])
        c = np.einsum("bkl,bl->bk", Gi, H)
    c = c.reshape(fs.NR, fs.NL, fs.NK)[:, layer, :]
    amp = np.hypot(c[:, 1], c[:, 2])
    psi = 0.5 * np.degrees(np.arctan2(c[:, 2], c[:, 1]))          # m=2 maximum, bar-aligned frame
    ang = (psi - (45.0 if which == "vR" else 0.0) - BAR_ANGLE_DEG + 90.0) % 180.0 - 90.0
    Rc = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
    m = (Rc > R_FIT[0]) & (Rc < R_FIT[1]) & np.isfinite(amp)
    return float(np.median(ang[m])), float(np.median(amp[m]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", default="nbody")
    ap.add_argument("--fracs", nargs="+", type=float, default=FRACS)
    ap.add_argument("--n-real", type=int, default=15)
    ap.add_argument("--selection", action="store_true", help="also apply the selection function")
    args = ap.parse_args()
    set_style()
    g = galaxy(args.galaxy)
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    pool = np.flatnonzero(np.hypot(W[:, 0], W[:, 1]) < 8.0)
    sub = np.sort(np.random.default_rng(2026).choice(pool, 50_000, replace=False))
    o = observe.to_galactic(*W[sub].T)
    tab = None
    if args.selection:
        s10 = 0.10 * 5 / np.log(10)
        tab = selection.load_tables(MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{s10:.3f}_N50k.hdf5")
        oa = observe.to_galactic(*W.T)

    rows = []
    for fr in args.fracs:
        sig = fr * 5 / np.log(10)
        got = {k: [] for k in ("vR", "absvR")}
        for s in range(args.n_real if fr > 0 else 1):
            d = (observe.scatter_distance_modulus(o["d"], sig, np.random.default_rng(700 + s))
                 if fr > 0 else o["d"])
            w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
            if tab is not None:
                da = (observe.scatter_distance_modulus(oa["d"], sig, np.random.default_rng(800 + s))
                      if fr > 0 else oa["d"])
                wa = np.column_stack(observe.to_galactocentric(oa["l"], oa["b"], da, oa["pml"], oa["pmb"],
                                                               oa["vlos"]))
                S = selection.evaluate(tab, oa["l"], oa["b"], da, np.hypot(wa[:, 0], wa[:, 1]))
                w = wa[np.random.default_rng(900 + s).random(len(S)) < np.minimum(S, 1.0)]
            for k in got:
                got[k].append(fourier_phase(w, k)[0])
        r = dict(frac=fr)
        for k, vv in got.items():
            r[f"{k}_mean"] = float(np.mean(vv))
            r[f"{k}_sd"] = float(np.std(vv, ddof=1)) if len(vv) > 1 else 0.0
        rows.append(r)
        print(f"  frac {fr:5.3f}: vR node {r['vR_mean']:+7.2f} +/- {r['vR_sd']:.2f}   "
              f"|vR/vtot| phase {r['absvR_mean']:+7.2f} +/- {r['absvR_sd']:.2f}", flush=True)
    df = pd.DataFrame(rows)
    # both estimators get the same zero-error calibration, so only the DRIFT is compared
    for k in ("vR", "absvR"):
        df[f"{k}_bias"] = df[f"{k}_mean"] - df[f"{k}_mean"].iloc[0]
    tag = "_sf" if args.selection else ""
    out = RESULTS_DIR / f"{g['prefix']}_bar_angle_estimators{tag}.csv"
    df.to_csv(out, index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.2, 5.0), gridspec_kw=dict(wspace=0.26))
    x = 100 * df.frac
    for k, c, mk, lab in (("vR", BLUE_RAMP[5], "o", r"$\langle v_R\rangle$ quadrupole node"),
                          ("absvR", ACCENT, "s", r"$\langle |v_R|/v_{\rm tot}\rangle$  $m=2$ phase")):
        a1.errorbar(x, df[f"{k}_mean"], yerr=df[f"{k}_sd"], fmt=mk + "-", color=c, lw=1.6, ms=6,
                    capsize=3, label=lab)
        a2.errorbar(x, df[f"{k}_bias"], yerr=df[f"{k}_sd"], fmt=mk + "-", color=c, lw=1.6, ms=6, capsize=3,
                    label=f"{lab}   ({df[f'{k}_bias'].iloc[-1]:+.1f}° at 20%)")
    a1.axhline(-BAR_ANGLE_DEG, color=INK, lw=1.2, ls="--", label=f"true bar axis {-BAR_ANGLE_DEG:g}°")
    a1.set_xlabel("fractional distance error [%]")
    a1.set_ylabel("measured angle [deg]")
    a1.legend(fontsize=8.5)
    a1.set_title("raw angle from each estimator", color=INK_2, fontsize=10)
    a2.axhline(0, color=INK, lw=1.2, ls="--")
    a2.set_xlabel("fractional distance error [%]")
    a2.set_ylabel("drift from the zero-error value [deg]")
    a2.legend(fontsize=8.5)
    a2.set_title("drift caused by the distance error (both zero-calibrated)", color=INK_2, fontsize=10)
    fig.suptitle(f"{g['label']}: which bar-angle estimator resists distance error?"
                 + (" (with selection)" if args.selection else " (no selection)"), fontsize=11.5, y=0.99)
    savefig(fig, f"test2_bar_angle_estimators_{g['prefix']}{tag}")
    print(f"\n{'frac':>6s} {'vR drift':>10s} {'|vR/vtot| drift':>17s}")
    for _, r in df.iterrows():
        print(f"{r.frac:6.3f} {r.vR_bias:+10.2f} {r.absvR_bias:+17.2f}")
    print(f"\nmedian |drift| over 5-20%: vR {df.vR_bias[df.frac>=0.05].abs().median():.2f}, "
          f"|vR/vtot| {df.absvR_bias[df.frac>=0.05].abs().median():.2f} deg")
    print(f"median scatter: vR {df.vR_sd[df.frac>0].median():.2f}, "
          f"|vR/vtot| {df.absvR_sd[df.frac>0].median():.2f} deg")
    print(f"wrote {rel(out)}")


if __name__ == "__main__":
    main()

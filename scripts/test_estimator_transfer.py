"""Do both bar-angle estimators transfer between snapshots, or only one?

The selection function is what breaks the transferability of the <v_R> quadrupole calibration: switching it
off reduces the t = 7 -> t = 4 offset from -1.86 to -0.35 deg (PROGRESS, 2026-11-01). The <|v_R|/v_tot> phase
of Zhang et al. (2024) is far less sensitive to distance errors; the question here is whether it is also less
sensitive to the galaxy, which is what a simulation-calibrated method needs.

For each snapshot the bar is rotated through a range of angles, observed with a 10% distance error and the
DR3-Mira-like selection, and both estimators are measured. The relation fitted on the reference snapshot is
then applied to the other one.

Usage: python test_estimator_transfer.py [--ref nbody] [--test nbody_t4] [--n-real 10]
"""
import argparse
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy
import observe, selection, steady_state_likelihood as ssl
from compare_bar_angle_estimators import fourier_phase

ANGLES = [15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35.]


def scan(gname, angles, frac, n_real):
    g = galaxy(gname)
    sig = frac * 5 / np.log(10)
    tab = selection.load_tables(MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{sig:.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    rows = []
    for phi in angles:
        o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - phi)).T)
        got = {"vR": [], "absvR": []}
        for s in range(n_real):
            d = observe.scatter_distance_modulus(o["d"], sig, np.random.default_rng(1000 + s))
            w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
            S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
            sel = w[np.random.default_rng(1100 + s).random(len(S)) < np.minimum(S, 1.0)]
            for k in got:
                got[k].append(fourier_phase(sel, k)[0])
        rows.append(dict(phi_true=phi, **{f"{k}_mean": np.mean(v) for k, v in got.items()},
                         **{f"{k}_sd": np.std(v, ddof=1) for k, v in got.items()}))
        print(f"  {gname} phi {phi:5.1f}: vR {rows[-1]['vR_mean']:+7.2f}  "
              f"|vR/vtot| {rows[-1]['absvR_mean']:+7.2f}", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="nbody")
    ap.add_argument("--test", default="nbody_t4")
    ap.add_argument("--frac", type=float, default=0.10)
    ap.add_argument("--n-real", type=int, default=10)
    args = ap.parse_args()
    set_style()
    R = scan(args.ref, ANGLES, args.frac, args.n_real)
    T = scan(args.test, ANGLES, args.frac, args.n_real)
    R.to_csv(RESULTS_DIR / f"{galaxy(args.ref)['prefix']}_estimator_scan.csv", index=False)
    T.to_csv(RESULTS_DIR / f"{galaxy(args.test)['prefix']}_estimator_scan.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), gridspec_kw=dict(wspace=0.26))
    print(f"\n{'estimator':<14s} {'ref slope':>10s} {'test slope':>11s} {'transfer mean':>14s} {'rms':>7s}")
    for ax, (k, c, lab) in zip(axes, (("vR", BLUE_RAMP[5], r"$\langle v_R\rangle$ node"),
                                      ("absvR", ACCENT, r"$\langle |v_R|/v_{\rm tot}\rangle$ phase"))):
        a, b = np.polyfit(R.phi_true, R[f"{k}_mean"], 1)
        at, bt = np.polyfit(T.phi_true, T[f"{k}_mean"], 1)
        inf = (T[f"{k}_mean"] - b) / a
        err = inf - T.phi_true
        print(f"{lab[7:30]:<14s} {a:10.3f} {at:11.3f} {err.mean():+14.2f} {np.sqrt((err**2).mean()):7.2f}")
        lim = (12, 38)
        ax.plot(lim, lim, color=INK, lw=1.2, ls="--", label="perfect recovery")
        ax.errorbar(T.phi_true, inf, yerr=T[f"{k}_sd"] / abs(a), fmt="o", color=c, ms=7, capsize=3, lw=1.4,
                    label=f"mean {err.mean():+.2f}°, rms {np.sqrt((err**2).mean()):.2f}°")
        ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
        ax.set_xlabel("true bar angle [deg]"); ax.set_ylabel("inferred bar angle [deg]")
        ax.legend(fontsize=8.5, loc="upper left")
        ax.set_title(lab, color=INK_2, fontsize=10)
    fig.suptitle(f"{galaxy(args.ref)['label']} calibration applied to {galaxy(args.test)['label']}: "
                 f"which estimator transfers? ({args.frac:.0%} error + selection)", fontsize=11.5, y=0.99)
    savefig(fig, f"test2_estimator_transfer_{galaxy(args.ref)['prefix']}_to_{galaxy(args.test)['prefix']}")


if __name__ == "__main__":
    main()

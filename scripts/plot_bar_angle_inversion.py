"""Observed v_R quadrupole angle against the true bar angle, at a fixed distance error.

This is the curve you would actually invert: measure the node in the data, read off the bar angle. It is not
the 1:1 line, because the line-of-sight distance error rotates the node by an amount that itself depends on
the bar's orientation to the Sun (PROGRESS, 2026-11-01), so the relation is flatter than 1:1 and a given node
error maps to a larger angle error.

The galaxy is rotated about z so the bar sits at the chosen angle, then re-observed from the fixed Sun. The
selection function is applied as a property of the *survey*: the same S(l, b, d) tables are used at every bar
angle, which is physically right - a survey's footprint and depth do not know where the bar is.

Usage: python plot_bar_angle_inversion.py [--frac 0.10] [--n-real 5]
Figure: figures/test2_bar_angle_inversion_<prefix>_<frac>pct.png
"""
import argparse
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy
import observe, selection, steady_state_likelihood as ssl
from plot_bar_angle_vs_distance_error import node_and_amp

ANGLES = [15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35.]
N_STARS, R_MAX, SEED = 50_000, 8.0, 2026


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", default="nbody")
    ap.add_argument("--frac", type=float, default=0.10)
    ap.add_argument("--angles", nargs="+", type=float, default=ANGLES)
    ap.add_argument("--n-real", type=int, default=5)
    args = ap.parse_args()
    set_style()
    g = galaxy(args.galaxy)
    sig = args.frac * 5.0 / np.log(10)
    tab = selection.load_tables(MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{sig:.3f}_N50k.hdf5")

    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    pool = np.flatnonzero(np.hypot(W[:, 0], W[:, 1]) < R_MAX)
    sub = np.sort(np.random.default_rng(SEED).choice(pool, N_STARS, replace=False))
    print(f"{g['label']}: {len(W):,} particles; {args.frac:.1%} distance error; "
          f"{'selection tables loaded' if tab else 'NO selection tables'}")

    rows = []
    for phi in args.angles:
        rot = np.radians(BAR_ANGLE_DEG - phi)                    # bar -27 -> -phi
        o_sub = observe.to_galactic(*ssl.rotate_z(W[sub], rot).T)
        o_all = observe.to_galactic(*ssl.rotate_z(W, rot).T) if tab else None
        n_no, n_sf = [], []
        for s in range(args.n_real):
            rng = np.random.default_rng(200 + s)
            d = observe.scatter_distance_modulus(o_sub["d"], sig, rng)
            w = np.column_stack(observe.to_galactocentric(o_sub["l"], o_sub["b"], d, o_sub["pml"],
                                                          o_sub["pmb"], o_sub["vlos"]))
            n_no.append(node_and_amp(w)[0])
            if tab:
                rng2 = np.random.default_rng(300 + s)
                da = observe.scatter_distance_modulus(o_all["d"], sig, rng2)
                wa = np.column_stack(observe.to_galactocentric(o_all["l"], o_all["b"], da, o_all["pml"],
                                                               o_all["pmb"], o_all["vlos"]))
                S = selection.evaluate(tab, o_all["l"], o_all["b"], da, np.hypot(wa[:, 0], wa[:, 1]))
                keep = np.random.default_rng(400 + s).random(len(S)) < np.minimum(S, 1.0)
                n_sf.append(node_and_amp(wa[keep])[0])
        r = dict(phi_true=phi, node_nosf=np.mean(n_no), node_nosf_sd=np.std(n_no, ddof=1))
        if tab:
            r.update(node_sf=np.mean(n_sf), node_sf_sd=np.std(n_sf, ddof=1), n_sel=int(keep.sum()))
        rows.append(r)
        print(f"  phi_true {phi:5.1f}: node no-SF {r['node_nosf']:+7.2f} +/- {r['node_nosf_sd']:.2f}"
              + (f"   SF {r['node_sf']:+7.2f} +/- {r['node_sf_sd']:.2f} ({r['n_sel']:,} stars)" if tab else ""),
              flush=True)
    df = pd.DataFrame(rows)
    out = RESULTS_DIR / f"{g['prefix']}_bar_angle_inversion_{100*args.frac:g}pct.csv"
    df.to_csv(out, index=False)

    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    lim = (min(args.angles) - 3, max(args.angles) + 3)
    ax.plot(lim, [-l for l in lim], color=INK, lw=1.3, ls="--", label="node = true bar axis (1:1)")
    for col, sdc, c, mk, lab in (("node_nosf", "node_nosf_sd", BLUE_RAMP[5], "o", "no selection function"),
                                 ("node_sf", "node_sf_sd", ACCENT, "s", "with selection function")):
        if col not in df:
            continue
        ax.errorbar(df.phi_true, df[col], yerr=df[sdc], fmt=mk, color=c, ms=7, capsize=3, lw=1.4, label=lab)
        p = np.polyfit(df.phi_true, df[col], 1)
        xx = np.array(lim)
        ax.plot(xx, np.polyval(p, xx), "-", color=c, lw=1.5, alpha=0.75,
                label=f"   slope {p[0]:.2f}  (1:1 would be −1)")
    ax.set_xlim(*lim)
    ax.set_xlabel("true bar angle [deg]")
    ax.set_ylabel(r"observed $m=2$ node of $\langle v_R \rangle$ [deg]")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title(f"{g['label']}, {args.frac:.0%} distance error: the node under-reports the bar angle,\n"
                 f"and the relation is flatter than 1:1 so it must be inverted with this curve",
                 fontsize=10, color=INK_2)
    savefig(fig, f"test2_bar_angle_inversion_{g['prefix']}_{100*args.frac:g}pct")

    for col, sdc, lab in (("node_nosf", "node_nosf_sd", "no selection"), ("node_sf", "node_sf_sd", "selection")):
        if col in df:
            p = np.polyfit(df.phi_true, df[col], 1)
            sd = df[sdc].median()
            print(f"\n{lab}: node = {p[0]:.3f} x phi_true + {p[1]:+.2f};  node scatter {sd:.2f} deg")
            print(f"   -> phi_true recoverable to {sd/abs(p[0]):.2f} deg")
    print(f"wrote {rel(out)}")


if __name__ == "__main__":
    main()

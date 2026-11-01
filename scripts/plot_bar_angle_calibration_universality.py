"""Is the distance-error -> bar-angle bias the same for every true bar angle?

`plot_bar_angle_vs_distance_error.py` calibrates the rotation of the v_R quadrupole node against the distance
error, at the one bar angle the mocks are built with (-27 deg). That calibration is only usable if it does not
itself depend on the bar angle - otherwise you would need to know the answer before you could apply it.

The test: rotate the galaxy about z so the bar sits at a chosen angle to the Sun-GC line, re-observe from the
(fixed) Sun, smear the distances, and measure the node. The same 50k stars are used at every angle - the
R < 8 kpc pool is rotation-invariant - so the only thing that changes is the geometry relative to the Sun,
which is exactly what the distance error couples to.

Usage: python plot_bar_angle_calibration_universality.py [--angles 15 20 ...] [--n-real 5]
Figure: figures/test2_bar_angle_calibration_universality_<prefix>.png
"""
import argparse
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, CMAP_SEQ
from galaxies import galaxy
import observe, steady_state_likelihood as ssl
from plot_bar_angle_vs_distance_error import node_and_amp

ANGLES = [15.0, 20.0, 25.0, 27.0, 30.0, 35.0]
FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
N_STARS, R_MAX, SEED = 50_000, 8.0, 2026


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", default="nbody")
    ap.add_argument("--angles", nargs="+", type=float, default=ANGLES)
    ap.add_argument("--fracs", nargs="+", type=float, default=FRACS)
    ap.add_argument("--n-real", type=int, default=5)
    args = ap.parse_args()
    set_style()
    g = galaxy(args.galaxy)

    with h5py.File(g["stars"], "r") as f:
        x, y = f["x"][:], f["y"][:]
        pool = np.flatnonzero(np.hypot(x, y) < R_MAX)
        idx = np.sort(np.random.default_rng(SEED).choice(pool, N_STARS, replace=False))
        w0 = np.column_stack([f[k][:][idx].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    print(f"{g['label']}: {len(w0):,} stars, bar at {-BAR_ANGLE_DEG:g} deg; "
          f"rotating to {args.angles}")

    rows = []
    for phi in args.angles:
        w = ssl.rotate_z(w0, np.radians(BAR_ANGLE_DEG - phi))     # bar -27 -> -phi
        o = observe.to_galactic(*w.T)
        for fr in args.fracs:
            sig = fr * 5.0 / np.log(10)
            nodes = []
            for s in range(args.n_real if fr > 0 else 1):
                d = (observe.scatter_distance_modulus(o["d"], sig, np.random.default_rng(100 + s))
                     if fr > 0 else o["d"])
                wo = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
                nodes.append(node_and_amp(wo)[0])
            n = np.array(nodes)
            shift = np.mod(n.mean() - (-phi) + 90, 180) - 90       # node minus the true axis
            rows.append(dict(phi_true=phi, frac=fr, node=n.mean(),
                             node_sd=float(n.std(ddof=1)) if len(n) > 1 else 0.0, shift_deg=shift))
        print(f"  phi_true {phi:5.1f}: shift {rows[-len(args.fracs)]['shift_deg']:+5.2f} deg at 0%, "
              f"{[round(r['shift_deg'],1) for r in rows[-len(args.fracs):]][-4]:+5.1f} at 10%", flush=True)
    df = pd.DataFrame(rows)
    out = RESULTS_DIR / f"{g['prefix']}_bar_angle_calibration_universality.csv"
    df.to_csv(out, index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.2, 5.0), gridspec_kw=dict(wspace=0.26))
    cols = CMAP_SEQ(np.linspace(0.25, 0.95, len(args.angles)))
    for c, phi in zip(cols, args.angles):
        d = df[df.phi_true == phi]
        a1.errorbar(100*d.frac, d['shift_deg'], yerr=d.node_sd, fmt="o-", color=c, lw=1.5, ms=5, capsize=2,
                    label=f"{phi:g}°")
    a1.axhline(0, color=INK, lw=1.2, ls="--")
    a1.set_xlabel("fractional distance error [%]")
    a1.set_ylabel("node − true bar angle [deg]")
    a1.legend(fontsize=8, title="true bar angle", title_fontsize=8, ncol=2)
    a1.set_title("calibration curve at each true bar angle", color=INK_2, fontsize=10)

    for c, phi in zip(cols, args.angles):
        d = df[df.phi_true == phi]
        ref = df[df.phi_true == BAR_ANGLE_DEG].set_index("frac")["shift_deg"]
        a2.plot(100*d.frac, d["shift_deg"].values - ref.reindex(d.frac).values, "o-", color=c, lw=1.5, ms=5)
    a2.axhline(0, color=INK, lw=1.2, ls="--")
    a2.set_xlabel("fractional distance error [%]")
    a2.set_ylabel(f"difference from the {BAR_ANGLE_DEG:g}° curve [deg]")
    a2.set_title("is the calibration universal? (flat at zero = yes)", color=INK_2, fontsize=10)
    fig.suptitle(f"{g['label']}: does the distance-error bar-angle bias depend on the true bar angle?",
                 fontsize=11.5, y=0.99)
    savefig(fig, f"test2_bar_angle_calibration_universality_{g['prefix']}")

    piv = df.pivot_table(index="frac", columns="phi_true", values="shift_deg")
    print(f"\nshift [deg] vs true bar angle:\n{piv.round(2).to_string()}")
    sp = piv.max(axis=1) - piv.min(axis=1)
    print(f"\nspread across bar angles, per error level: median {sp.median():.2f} deg, max {sp.max():.2f} deg")
    print(f"typical single-realisation scatter: {df[df.node_sd>0].node_sd.median():.2f} deg")
    print(f"wrote {rel(out)}")


if __name__ == "__main__":
    main()

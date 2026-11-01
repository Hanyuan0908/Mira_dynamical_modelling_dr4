"""Bar angle from the v_R quadrupole node, as a function of the distance error.

The orbital-library posterior recovers phi_b poorly - not for lack of precision (Fisher 1.2 deg) but because
of a -4.8 deg systematic. The bar angle, though, is directly visible: a bar drives a quadrupole in <v_R> whose
node lies along the bar's major axis. The distance error rotates that node, so the node alone is biased - but
the rotation is a smooth, monotonic function of the error, so if sigma_DM is known the bias can be calibrated
away. This builds that calibration curve.

For each fractional distance error a mock is generated at the chosen snapshot, with and without the
DR3-Mira-like selection function, and the node is fitted per ring with
`test1_nbody_qmc_focused_fourier_stars.fit` - the per-star least-squares fit of <v_R>(phi), which is insensitive
to how the stars are spread in azimuth and so is not biased by the selection.

Usage: python plot_bar_angle_vs_distance_error.py [--galaxy nbody] [--fracs 0 0.025 ...] [--no-regenerate]
Figure: figures/test2_bar_angle_vs_distance_error_<prefix>.png;  table: results/<prefix>_bar_angle_vs_distance_error.csv
"""
import argparse, subprocess, sys
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs

FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
R_FIT = (1.0, 4.0)       # kpc, rings used for the node (inside the bar)
PY = sys.executable


def node_and_amp(w_mw, layer=0):
    """Median m = 2 node angle [deg, MW frame] and amplitude [km/s] of <v_R> over R_FIT."""
    b, X, v, _ = fs.prepare(ssl.rotate_z(w_mw, np.radians(BAR_ANGLE_DEG)))
    n = np.bincount(b, minlength=fs.NB)
    c = fs.fit(b, X, v, np.ones(len(b)), n >= fs.N_MIN_BLOCK)[:, 0, :].reshape(fs.NR, fs.NL, fs.NK)[:, layer, :]
    amp = np.hypot(c[:, 1], c[:, 2])
    node = (0.5 * np.degrees(np.arctan2(c[:, 2], c[:, 1])) - 45.0 - BAR_ANGLE_DEG + 90.0) % 180.0 - 90.0
    Rc = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
    m = (Rc > R_FIT[0]) & (Rc < R_FIT[1]) & np.isfinite(amp)
    return float(np.median(node[m])), float(np.median(amp[m]))


def load(path, true=False):
    sfx = "_true" if true else ""
    with h5py.File(path, "r") as f:
        return np.column_stack([f[f"{k}{sfx}"][:].astype(np.float64) for k in
                                ("x", "y", "z", "vx", "vy", "vz")])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", default="nbody")
    ap.add_argument("--fracs", nargs="+", type=float, default=FRACS)
    ap.add_argument("--no-regenerate", action="store_true")
    ap.add_argument("--n-real", type=int, default=10,
                    help="independent distance-error draws per level, for the no-selection error bar")
    args = ap.parse_args()
    set_style()
    g = galaxy(args.galaxy)
    pre = g["prefix"]
    rows = []
    for fr in args.fracs:
        sig = fr * 5.0 / np.log(10)
        nosf = MOCK_DIR / f"{pre}_mock_distance_error_sigmaDM{sig:.3f}_R8_N50k.hdf5"
        sfm = MOCK_DIR / f"{pre}_mock_selection_sigmaDM{sig:.3f}_N50k.hdf5"
        if not args.no_regenerate:
            if fr > 0 and not nosf.exists():
                subprocess.run([PY, "-W", "ignore", "make_mock_distance_error.py", "--galaxy", args.galaxy,
                                "--frac", str(fr)], check=True, capture_output=True)
            if not sfm.exists():
                subprocess.run([PY, "-W", "ignore", "make_mock_selection.py", "--galaxy", args.galaxy,
                                "--frac", str(max(fr, 1e-9))], check=True, capture_output=True)
        row = dict(frac=fr, sigma_DM=sig)
        if fr == 0:
            perfect = MOCK_DIR / f"{pre}_mock_perfect_R8_N50k.hdf5"
            row["node_nosf"], row["amp_nosf"] = node_and_amp(load(perfect))
        elif nosf.exists():
            row["node_nosf"], row["amp_nosf"] = node_and_amp(load(nosf))
            row["node_nosf_true"], _ = node_and_amp(load(nosf, true=True))
        if sfm.exists():
            row["node_sf"], row["amp_sf"] = node_and_amp(load(sfm))
            row["node_sf_true"], _ = node_and_amp(load(sfm, true=True))
        # repeat selection-mock realisations: each needs the whole mock rebuilt (new error draw, refit
        # selection, new Bernoulli draw), so they are made, measured and deleted one at a time
        if args.n_real > 1 and not args.no_regenerate:
            ns = []
            for s in range(1, args.n_real):
                se, ss = 3000 + s, 4000 + s
                tmp = MOCK_DIR / (f"{pre}_mock_selection_sigmaDM{sig:.3f}_N50k_s{se}.hdf5")
                try:
                    subprocess.run([PY, "-W", "ignore", "make_mock_selection.py", "--galaxy", args.galaxy,
                                    "--frac", str(max(fr, 1e-9)), "--seed-err", str(se),
                                    "--seed-sel", str(ss)], check=True, capture_output=True)
                    ns.append(node_and_amp(load(tmp))[0])
                finally:
                    tmp.unlink(missing_ok=True)
            if ns:
                allv = ns + [row["node_sf"]]
                row["node_sf_mean"], row["node_sf_sd"] = float(np.mean(allv)), float(np.std(allv, ddof=1))
        # scatter of the node over independent error draws (no-selection case, done in memory)
        if args.n_real > 1 and fr > 0:
            import observe
            with h5py.File(MOCK_DIR / f"{pre}_mock_perfect_R8_N50k.hdf5", "r") as f:
                mm = {k: f[k][:].astype(np.float64) for k in ("l", "b", "d", "pml", "pmb", "vlos")}
            nn = []
            for s in range(args.n_real):
                dd = observe.scatter_distance_modulus(mm["d"], sig, np.random.default_rng(100 + s))
                w = np.column_stack(observe.to_galactocentric(mm["l"], mm["b"], dd, mm["pml"],
                                                              mm["pmb"], mm["vlos"]))
                nn.append(node_and_amp(w)[0])
            row["node_nosf_mean"], row["node_nosf_sd"] = float(np.mean(nn)), float(np.std(nn, ddof=1))
        rows.append(row)
        print(f"  frac {fr:5.3f} (sigma_DM {sig:.3f}): node no-SF {row.get('node_nosf', np.nan):+7.2f}, "
              f"SF {row.get('node_sf', np.nan):+7.2f} deg", flush=True)
    df = pd.DataFrame(rows)
    out = RESULTS_DIR / f"{pre}_bar_angle_vs_distance_error.csv"
    df.to_csv(out, index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw=dict(wspace=0.26))
    x = 100 * df.frac
    a1.axhline(-BAR_ANGLE_DEG, color=INK, lw=1.3, ls="--", label=f"true bar axis, {-BAR_ANGLE_DEG:g}°")
    if "node_nosf_mean" in df:
        a1.errorbar(x, df.node_nosf_mean, yerr=df.node_nosf_sd, fmt="o-", color=BLUE_RAMP[5], lw=1.7,
                    ms=6, capsize=3, label=f"no selection ({args.n_real} error draws)")
    elif "node_nosf" in df:
        a1.plot(x, df.node_nosf, "o-", color=BLUE_RAMP[5], lw=1.7, ms=6, label="no selection function")
    if "node_sf_mean" in df:
        a1.errorbar(x, df.node_sf_mean, yerr=df.node_sf_sd, fmt="s-", color=ACCENT, lw=1.7, ms=6,
                    capsize=3, label=f"with selection ({args.n_real} realisations)")
    elif "node_sf" in df:
        a1.plot(x, df.node_sf, "s-", color=ACCENT, lw=1.7, ms=6, label="with selection function")
    a1.set_xlabel("fractional distance error [%]")
    a1.set_ylabel(r"measured $m=2$ node of $\langle v_R\rangle$ [deg]")
    a1.legend(fontsize=8.5, loc="upper left")
    a1.set_title("the distance error rotates the apparent bar angle", color=INK_2, fontsize=10)

    col = "node_nosf_mean" if "node_nosf_mean" in df else "node_nosf"
    sdc = df.node_nosf_sd if "node_nosf_sd" in df else None
    sh = df[col] - (-BAR_ANGLE_DEG)
    a2.errorbar(x, sh, yerr=sdc, fmt="o", color=BLUE_RAMP[5], ms=6, capsize=3, label="no selection")
    k = np.isfinite(sh) & (df.frac <= 0.125)
    pf = np.polyfit(df.frac[k], sh[k], 1)
    xx = np.linspace(0, 0.125, 50)
    a2.plot(100*xx, np.polyval(pf, xx), "-", color=BLUE_RAMP[5], lw=1.4,
            label=f"{pf[0]/100:.2f}°  per 1% error")
    if "node_sf_mean" in df:
        a2.errorbar(x, df.node_sf_mean - (-BAR_ANGLE_DEG), yerr=df.node_sf_sd, fmt="s", color=ACCENT,
                    ms=6, capsize=3, label="with selection")
    elif "node_sf" in df:
        a2.plot(x, df.node_sf - (-BAR_ANGLE_DEG), "s", color=ACCENT, ms=6, label="with selection")
    a2.axhline(0, color=INK, lw=1.3, ls="--")
    a2.set_xlabel("fractional distance error [%]")
    a2.set_ylabel("node − true bar angle [deg]")
    a2.legend(fontsize=8.5, loc="upper left")
    a2.set_title("the calibration: subtract this to recover the bar angle", color=INK_2, fontsize=10)
    fig.suptitle(f"{g['label']}: bar angle from the $v_R$ quadrupole node against distance error",
                 fontsize=11.5, y=0.99)
    savefig(fig, f"test2_bar_angle_vs_distance_error_{pre}")
    print(f"\nwrote {rel(out)}")
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

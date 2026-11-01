"""Bar-angle estimators with amplitude-weighted ring combination: bias vs distance error, and transfer.

Earlier versions combined the per-ring m = 2 phases with a plain median over a fixed 1-4 kpc window, which
gives a ring with no bar signal the same vote as one with a strong signal. Here the complex m = 2
coefficients are summed instead, so each ring weights itself by its own amplitude. That halves the
snapshot-to-snapshot spread of the <|v_R|/v_tot> drift, and - more usefully - makes the <v_R> node's drift
nearly universal (PROGRESS, 2026-11-01).

Seeds are independent per (angle, error level, realisation); sharing them across points was found to inflate
the scatter of fitted slopes by ~2x and was the reason two earlier transfer runs disagreed.

Two figures:
  ..._drift   measured angle and bias against distance error, both snapshots, both estimators
  ..._transfer  the relation fitted on one snapshot applied to the other, in both directions

Usage: python bar_angle_vecsum_verification.py [--n-real 12]
"""
import argparse
import numpy as np, pandas as pd, h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy, nbody_name
import observe, selection, steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs

Rc = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
FRACS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20]
ANGLES = [15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35.]
EST = {"vR": r"$\langle v_R\rangle$ node", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}


def ring_coeffs(w, which, layer=0):
    b, X, v, _ = fs.prepare(ssl.rotate_z(w, np.radians(BAR_ANGLE_DEG)))
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
    return c[:, 1], c[:, 2]


def phase(w, which):
    """Amplitude-weighted m=2 phase over 1-4 kpc: the argument of the summed complex coefficients."""
    a2, b2 = ring_coeffs(w, which)
    m = (Rc > 1) & (Rc < 4) & np.isfinite(a2) & np.isfinite(b2)
    z = np.nansum(a2[m] + 1j * b2[m])
    off = 45.0 if which == "vR" else 0.0
    return float((0.5 * np.degrees(np.angle(z)) - off - BAR_ANGLE_DEG + 90) % 180 - 90)


def measure(gname, phi_true, frac, n_real, seed0):
    """Mean and SEM of both estimators at one (bar angle, distance error)."""
    g = galaxy(gname)
    s10 = 0.10 * 5 / np.log(10)
    tab = selection.load_tables(MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{s10:.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - phi_true)).T)
    sig = frac * 5 / np.log(10)
    got = {k: [] for k in EST}
    for s in range(n_real):
        rng = np.random.default_rng(seed0 + s)
        d = observe.scatter_distance_modulus(o["d"], sig, rng) if frac > 0 else o["d"]
        w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
        S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
        sel = w[np.random.default_rng(seed0 + 500000 + s).random(len(S)) < np.minimum(S, 1.0)]
        for k in got:
            got[k].append(phase(sel, k))
    return {k: (float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(n_real))) for k, v in got.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-real", type=int, default=12)
    args = ap.parse_args()
    set_style()
    snaps = [("nbody", 7.0), ("nbody_t4", 4.0)]

    # --- 1. bias against distance error, at the native bar angle ---
    rows = []
    for gn, tm in snaps:
        for i, fr in enumerate(FRACS):
            r = measure(gn, BAR_ANGLE_DEG, fr, args.n_real, seed0=10_000 + 7919 * i + 104729 * int(tm * 2))
            rows.append(dict(snap=f"t={tm:g}", frac=fr,
                             **{f"{k}_{s}": r[k][j] for k in EST for j, s in enumerate(("mean", "err"))}))
            print(f"  {gn} frac {fr:5.3f}: " + "  ".join(
                f"{k} {r[k][0]:+7.2f}+/-{r[k][1]:.2f}" for k in EST), flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS_DIR / "bar_angle_vecsum_drift.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.0), gridspec_kw=dict(wspace=0.26))
    for ax, k in zip(axes, EST):
        for (gn, tm), c, mk in zip(snaps, (BLUE_RAMP[5], ACCENT), ("o", "s")):
            d = D[D.snap == f"t={tm:g}"]
            ax.errorbar(100 * d.frac, d[f"{k}_mean"] + BAR_ANGLE_DEG, yerr=d[f"{k}_err"], fmt=mk + "-",
                        color=c, lw=1.6, ms=6, capsize=3, label=f"{galaxy(gn)['label']}")
        ax.axhline(0, color=INK, lw=1.2, ls="--")
        ax.set_xlabel("fractional distance error [%]")
        ax.set_ylabel("measured angle − true bar angle [deg]")
        ax.set_title(EST[k], color=INK_2, fontsize=10)
        ax.legend(fontsize=8.5)
    fig.suptitle("Bar-angle bias against distance error, amplitude-weighted ring combination",
                 fontsize=11.5, y=0.99)
    savefig(fig, "test2_bar_angle_vecsum_drift")

    # --- 2. transfer, both directions, at 10% ---
    scans = {}
    for gn, tm in snaps:
        out = []
        for i, phi in enumerate(ANGLES):
            r = measure(gn, phi, 0.10, args.n_real, seed0=900_000 + 31337 * i + 7717 * int(tm * 2))
            out.append(dict(phi_true=phi,
                            **{f"{k}_{s}": r[k][j] for k in EST for j, s in enumerate(("mean", "err"))}))
            print(f"  scan {gn} phi {phi:5.1f}: " + "  ".join(f"{k} {r[k][0]:+7.2f}" for k in EST), flush=True)
        scans[gn] = pd.DataFrame(out)
        scans[gn].to_csv(RESULTS_DIR / f"{galaxy(gn)['prefix']}_vecsum_scan.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4), gridspec_kw=dict(wspace=0.26))
    print(f"\n{'estimator':<12s} {'direction':<18s} {'mean err':>9s} {'rms':>7s}")
    for ax, k in zip(axes, EST):
        lim = (12, 38)
        ax.plot(lim, lim, color=INK, lw=1.2, ls="--", label="perfect recovery")
        for (src, dst), c, mk in ((("nbody", "nbody_t4"), BLUE_RAMP[5], "o"),
                                  (("nbody_t4", "nbody"), ACCENT, "s")):
            R, T = scans[src], scans[dst]
            a, b = np.polyfit(R.phi_true, R[f"{k}_mean"], 1)
            inf = (T[f"{k}_mean"] - b) / a
            e = inf - T.phi_true
            lab = f"{galaxy(src)['label'][-6:-1]} → {galaxy(dst)['label'][-6:-1]}"
            print(f"{k:<12s} {lab:<18s} {e.mean():+9.2f} {np.sqrt((e**2).mean()):7.2f}")
            ax.errorbar(T.phi_true, inf, yerr=T[f"{k}_err"] / abs(a), fmt=mk, color=c, ms=6, capsize=3,
                        lw=1.3, label=f"{lab}: {e.mean():+.2f}° ± {e.std(ddof=1):.2f}")
        ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
        ax.set_xlabel("true bar angle [deg]"); ax.set_ylabel("inferred bar angle [deg]")
        ax.set_title(EST[k], color=INK_2, fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("Calibration transferred between snapshots, 10% distance error + selection",
                 fontsize=11.5, y=0.99)
    savefig(fig, "test2_bar_angle_vecsum_transfer")


if __name__ == "__main__":
    main()

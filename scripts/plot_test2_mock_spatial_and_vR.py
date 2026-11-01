"""Spatial distribution (x-y, x-z, l-b) and v_R quadrupole of the Test 2 mocks.

The Test 2 mock set (PROGRESS, 2026-10-29) is one galaxy seen at five snapshot times, i.e. five pattern speeds,
each observed with a 10% distance error and a DR3-Mira-like selection function, plus a no-selection benchmark at
the fiducial snapshot. This script puts them side by side.

Figure 1 (`test2_mock_spatial`): counts in x-y, x-z and l-b, one column per data set, with the DR3 Miras as the
reference column. Log10 of the fraction of the sample per pixel, with a shared colour range per projection row so
that the columns are comparable.

Figure 2 (`test2_mock_vR_quadrupole`): the bar's signature in <v_R>.
  - Top: face-on <v_R> maps (|z| < 1 kpc). A bar drives a quadrupole: gas and stars stream inwards on one pair of
    quadrants and outwards on the other, with nodes along the bar's major and minor axes. The dashed line is the
    true bar axis (-27 deg); the measured node is printed on each panel.
  - Bottom: the m = 2 amplitude of <v_R> against R, and the node angle it implies, for every set. These come from
    `test1_nbody_qmc_focused_fourier_stars` (`fs.prepare` / `fs.fit`), the same per-star least-squares fit the
    likelihood uses: it fits <v_R>(phi) rather than summing v_R cos 2phi over stars, so it is insensitive to how the
    stars are spread in phi and therefore unbiased by the selection function. Solid lines are the observed (smeared)
    data, dashed the same stars at their true positions, so the gap is what the 10% distance error costs.

Usage:
    python plot_test2_mock_spatial_and_vR.py                  # the standard 10% Test 2 set
    python plot_test2_mock_spatial_and_vR.py --frac 0.075     # another error level
    python plot_test2_mock_spatial_and_vR.py --mocks a.hdf5 b.hdf5 --labels "A" "B"
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import (MOCK_DIR, MIRA_DR3_FILE, RESULTS_DIR, R0_KPC, Z_SUN_KPC, BAR_ANGLE_DEG, set_style, savefig,
                    rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP, CMAP_SEQ_W)
from galaxies import galaxy, nbody_name
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs

TIMES = (2, 3, 4, 5.5, 7)            # snapshot model times of the selection mocks
SIGMA_DM_CUT = 0.4                   # DR3 Mira cut, as in make_mock_selection.py
R_MAX = 8.0
PIX = 0.2                            # kpc, spatial pixels
SKY_PIX = 1.0                        # deg
PIX_V = 0.4                          # kpc, velocity-map pixels
MIN_STARS = 12                       # minimum stars per velocity pixel
AMP_MIN = 5.0                        # km/s, below this the m=2 amplitude is noise, not a quadrupole
V_R_LIM = 40.0                       # km/s, colour range of the v_R maps
Z_CUT = 1.0                          # kpc, slab for the face-on v_R maps
EVOLUTION_CSV = RESULTS_DIR / "nbody_pattern_speed_evolution.csv"


def default_sets(frac):
    """(file, label) of the standard Test 2 set: the selection mocks, then the no-selection benchmark."""
    sigma = frac * 5.0 / np.log(10)
    out = []
    for t in TIMES:
        g = galaxy(nbody_name(t))
        out.append((MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{sigma:.3f}_N50k.hdf5", f"t = {t:g}"))
    g = galaxy(nbody_name(TIMES[-1]))
    out.append((MOCK_DIR / f"{g['prefix']}_mock_distance_error_sigmaDM{sigma:.3f}_R8_N50k.hdf5",
                f"t = {TIMES[-1]:g}, no selection"))
    return out


def omega_of(t):
    """Measured Omega_b of a snapshot, if the evolution table has been made."""
    if not EVOLUTION_CSV.exists():
        return None
    df = pd.read_csv(EVOLUTION_CSV)
    hit = df[np.isclose(df.t_model, float(t))]
    return float(hit.Omega_p.iloc[0]) if len(hit) else None


def load(path):
    """Observed and true phase space of a mock, plus (l, b)."""
    with h5py.File(path, "r") as f:
        obs = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
        true = np.column_stack([f[f"{k}_true"][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
        l, b = f["l"][:].astype(np.float64), f["b"][:].astype(np.float64)
        attrs = dict(f.attrs)
    return obs, true, ((l + 180.0) % 360.0) - 180.0, b, attrs


def mira_reference():
    df = pd.read_csv(MIRA_DR3_FILE, usecols=["l", "b", "DM_err", "x", "y", "z"]).dropna()
    df = df[(df.DM_err < SIGMA_DM_CUT) & (np.hypot(df.x, df.y) < R_MAX)]
    return (df.x.values, df.y.values, df.z.values, ((df.l.values + 180.0) % 360.0) - 180.0, df.b.values), len(df)


def projections():
    xy = np.arange(-8.6, 8.6 + 1e-9, PIX)
    zz = np.arange(-4.0, 4.0 + 1e-9, PIX)
    ll = np.arange(-90.0, 90.0 + 1e-9, SKY_PIX)       # everything with R < 8 kpc lies within |l| < 90 deg
    bb = np.arange(-20.0, 20.0 + 1e-9, SKY_PIX)
    return [(0, 1, xy, xy, "x [kpc]", "y [kpc]", "face-on (x–y)"),
            (0, 2, xy, zz, "x [kpc]", "z [kpc]", "edge-on (x–z)"),
            (3, 4, ll, bb, "l [deg]", "b [deg]", "sky (l, b)")]


def vR_map(w, edges):
    """Face-on <v_R> in a |z| < Z_CUT slab, and the star count per pixel."""
    s = np.abs(w[:, 2]) < Z_CUT
    x, y = w[s, 0], w[s, 1]
    R = np.hypot(x, y)
    vR = (x * w[s, 3] + y * w[s, 4]) / R
    N = np.histogram2d(x, y, bins=[edges, edges])[0]
    V = np.histogram2d(x, y, bins=[edges, edges], weights=vR)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = V / N
    mean[N < MIN_STARS] = np.nan
    return mean.T


def vR_m2(w_mw, layer=0):
    """m = 2 amplitude [km/s] and node angle [deg] of <v_R> per ring, from the likelihood's per-star fit.

    `fs.fit` least-squares fits <v_R>(phi) = a0 + a2 cos 2phi + b2 sin 2phi + (m = 4), so it measures the mean
    field and not the density-weighted one: the selection function changes how stars are spread in phi, which
    would bias a plain sum but does not bias this. Returns (R centres, amplitude, node angle in the MW frame)."""
    b, X, v, _ = fs.prepare(ssl.rotate_z(w_mw, np.radians(BAR_ANGLE_DEG)))     # bar on the x-axis
    n = np.bincount(b, minlength=fs.NB)
    ok = n >= fs.N_MIN_BLOCK
    c = fs.fit(b, X, v, np.ones(len(b)), ok)[:, 0, :]          # field 0 = <v_R>
    c = c.reshape(fs.NR, fs.NL, fs.NK)[:, layer, :]
    amp = np.hypot(c[:, 1], c[:, 2])
    # the quadrupole node sits 45 deg behind the maximum of the m = 2 term; -27 deg back to the MW frame
    node = (0.5 * np.degrees(np.arctan2(c[:, 2], c[:, 1])) - 45.0 - BAR_ANGLE_DEG + 90.0) % 180.0 - 90.0
    Rc = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
    good = ok.reshape(fs.NR, fs.NL)[:, layer] & np.isfinite(amp)
    return Rc[good], amp[good], node[good]


def figure_spatial(sets, data, tag, sigma):
    """Counts in three projections; one column per data set, DR3 Miras first."""
    mira, n_mira = mira_reference()
    cols = [(f"DR3 Miras\n{n_mira:,}", mira)] + \
           [(f"{lab}\n{len(d[0]):,}", (d[0][:, 0], d[0][:, 1], d[0][:, 2], d[2], d[3]))
            for (_, lab), d in zip(sets, data)]
    projs = projections()
    # row heights follow the drawn aspect of each projection, so the panels fill their rows:
    # x-y is square, x-z is 8/17.2 of that, and l-b is free-aspect.
    pw = 2.45
    fig, axes = plt.subplots(3, len(cols), figsize=(pw * len(cols) + 1.5, pw * 2.35 + 1.5),
                             gridspec_kw=dict(height_ratios=[1.0, 0.47, 0.6], hspace=0.42, wspace=0.28))
    for r, (iu, iv, eu, ev, xl, yl, title) in enumerate(projs):
        hs = [np.histogram2d(c[iu], c[iv], bins=[eu, ev])[0].T / len(c[iu]) for _, c in cols]
        pos = np.concatenate([h[h > 0] for h in hs])
        lo, hi = np.log10(np.percentile(pos, 1)), np.log10(pos.max())
        for k, ((label, _), h) in enumerate(zip(cols, hs)):
            ax = axes[r, k]
            im = ax.pcolormesh(eu, ev, np.log10(np.where(h > 0, h, np.nan)), cmap=CMAP_SEQ_W,
                               vmin=lo, vmax=hi, rasterized=True)
            if r < 2:
                ax.plot(-R0_KPC, Z_SUN_KPC if r == 1 else 0.0, marker="*", color=ACCENT, ms=9, ls="none", zorder=5)
                ax.set_aspect("equal")
            else:
                ax.set_xlim(eu[-1], eu[0])            # l increases to the left
            ax.set_xlabel(xl, fontsize=8.5)
            ax.set_ylabel(yl if k == 0 else "", fontsize=8.5)
            ax.tick_params(labelsize=7.5)
            if r == 0:
                ax.set_title(label, fontsize=9.5, color=INK)
            if k == len(cols) - 1:
                cb = fig.colorbar(im, ax=axes[r, :].tolist(), shrink=0.85, pad=0.012, aspect=12)
                cb.set_label(f"log10 fraction per pixel\n{title}", fontsize=7.5)
                cb.ax.tick_params(labelsize=7)
    fig.suptitle(f"Test 2 mocks, spatial distribution ({100 * np.log(10) / 5 * sigma:.3g}% distance error): "
                 f"selection-function mocks at five snapshot times, the no-selection benchmark, and the DR3 Miras. "
                 f"{PIX} kpc and {SKY_PIX:g}° pixels; star = Sun.", fontsize=10.5, y=0.985)
    return savefig(fig, f"test2_mock_spatial{tag}")


def figure_vR(sets, data, tag, sigma):
    """Face-on <v_R> maps and the m = 2 (quadrupole) profiles."""
    edges = np.arange(-8.0, 8.0 + 1e-9, PIX_V)
    ncol = 3
    nrow = int(np.ceil(len(sets) / ncol))
    fig = plt.figure(figsize=(4.4 * ncol, 4.0 * nrow + 4.2))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[1.0] * nrow + [0.92], hspace=0.38, wspace=0.30)

    for k, ((path, lab), (obs, true, _, _, attrs)) in enumerate(zip(sets, data)):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        img = vR_map(obs, edges)
        im = ax.pcolormesh(edges, edges, img, cmap=CMAP_DIV, vmin=-V_R_LIM, vmax=V_R_LIM, rasterized=True)
        phi = np.radians(-BAR_ANGLE_DEG)                       # bar at -27 deg in the MW frame
        Rb = 5.0
        ax.plot([-Rb * np.cos(phi), Rb * np.cos(phi)], [-Rb * np.sin(phi), Rb * np.sin(phi)],
                color=INK, lw=1.3, ls=(0, (5, 3)), zorder=6)
        ax.plot(-R0_KPC, 0.0, marker="*", color=ACCENT, ms=10, ls="none", zorder=6)
        Rc, amp, node = vR_m2(obs)
        inner = (Rc > 1.0) & (Rc < 4.0)
        ax.set_aspect("equal")
        ax.set_xlim(-8, 8)
        ax.set_ylim(-8, 8)
        ax.set_xlabel("x [kpc]", fontsize=8.5)
        ax.set_ylabel("y [kpc]", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        om = omega_of(attrs.get("snapshot_time", np.nan))
        head = lab + (f",  $\\Omega_b$ = {om:.1f}" if om else "")
        ax.set_title(f"{head}    N = {len(obs):,}\nm=2 node {np.median(node[inner]):+.1f}° (bar −27°), "
                     f"amp {np.median(amp[inner]):.1f} km/s", fontsize=9, color=INK_2)
        fig.colorbar(im, ax=ax, shrink=0.86, pad=0.02).set_label(r"$\langle v_R \rangle$ [km s$^{-1}$]", fontsize=8)

    colours = [BLUE_RAMP[6], BLUE_RAMP[5], BLUE_RAMP[4], BLUE_RAMP[3], BLUE_RAMP[1], ACCENT]
    a_amp = fig.add_subplot(gs[nrow, 0])
    a_damp = fig.add_subplot(gs[nrow, 1])
    a_node = fig.add_subplot(gs[nrow, 2])
    for k, ((_, lab), (obs, true, _, _, _)) in enumerate(zip(sets, data)):
        c = colours[k % len(colours)]
        Rc, amp, node = vR_m2(obs)
        Rt, amp_t, _ = vR_m2(true)
        a_amp.plot(Rc, amp, "-", color=c, lw=1.7, label=f"{lab}  (N = {len(obs):,})")
        if len(Rt) == len(Rc) and np.allclose(Rt, Rc):
            # only where the quadrupole is actually measurable: a ratio of two small, noisy
            # amplitudes carries no information
            m = amp_t > AMP_MIN
            a_damp.plot(Rc[m], (amp / amp_t)[m], "-", color=c, lw=1.7)
        a_node.plot(Rc, node, "-", color=c, lw=1.5)
    a_amp.set_xlabel("R [kpc]")
    a_amp.set_ylabel(r"$m=2$ amplitude of $\langle v_R \rangle$ [km s$^{-1}$]")
    a_amp.set_xlim(0, 8)
    a_amp.set_ylim(0, None)
    a_amp.legend(fontsize=7.5, loc="upper right", framealpha=0.92, facecolor="white", edgecolor="none")
    a_amp.set_title("quadrupole strength (observed, |z| < 1 kpc)", fontsize=9, color=INK_2)

    a_damp.axhline(1.0, color=INK_MUTED, lw=1.0, ls="--")
    a_damp.set_xlabel("R [kpc]")
    a_damp.set_ylabel("observed / true amplitude")
    a_damp.set_xlim(0, 8)
    a_damp.set_ylim(0.5, 1.6)
    a_damp.set_title(f"what the {100 * np.log(10) / 5 * sigma:.3g}% distance error does to it\n"
                     f"(same stars, true positions; rings with amplitude > {AMP_MIN:g} km/s)",
                     fontsize=9, color=INK_2)

    a_node.axhline(-BAR_ANGLE_DEG, color=INK, lw=1.2, ls="--", zorder=1)
    a_node.text(0.2, -BAR_ANGLE_DEG - 7, "bar axis, −27°", fontsize=7.5, color=INK, ha="left")
    a_node.set_xlabel("R [kpc]")
    a_node.set_ylabel("node angle [deg]")
    a_node.set_xlim(0, 8)
    a_node.set_ylim(-60, 10)
    a_node.set_title("quadrupole node vs the bar axis\n(observed; the offset is a distance-error bias)",
                     fontsize=9, color=INK_2)
    fig.suptitle(r"Test 2 mocks: the bar's $\langle v_R \rangle$ quadrupole "
                 f"({100 * np.log(10) / 5 * sigma:.3g}% distance error, observed positions)", fontsize=11.5,
                 y=0.995)
    return savefig(fig, f"test2_mock_vR_quadrupole{tag}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frac", type=float, default=0.10, help="fractional distance error of the set (default 0.10)")
    ap.add_argument("--mocks", nargs="+", default=None, help="mock file names in mock/ (overrides --frac)")
    ap.add_argument("--labels", nargs="+", default=None, help="labels for --mocks")
    args = ap.parse_args()
    set_style()

    if args.mocks:
        sets = [(MOCK_DIR / m, (args.labels[i] if args.labels else m)) for i, m in enumerate(args.mocks)]
        tag = ""
    else:
        sets = default_sets(args.frac)
        tag = f"_{100 * args.frac:g}pct"
    missing = [p for p, _ in sets if not p.exists()]
    if missing:
        raise SystemExit("missing mocks:\n  " + "\n  ".join(str(rel(m)) for m in missing))

    data = []
    for path, lab in sets:
        data.append(load(path))
        print(f"loaded {rel(path)}: {len(data[-1][0]):,} stars", flush=True)
    sigma = float(data[0][4].get("sigma_DM", np.nan))

    figure_spatial(sets, data, tag, sigma)
    figure_vR(sets, data, tag, sigma)

    mira, n_mira = mira_reference()
    print("\nm = 2 amplitude of <v_R>, median over 1 < R < 4 kpc (|z| < 1 kpc):")
    print(f"  {'set':26s} {'N':>8s} {'N(|z|<1)':>9s} {'obs':>8s} {'true':>7s} {'change':>7s} "
          f"{'node obs':>10s} {'node true':>9s} {'shift':>8s}")
    for (_, lab), (obs, true, _, _, _) in zip(sets, data):
        Rc, a_o, n_o = vR_m2(obs)
        _, a_t, n_t = vR_m2(true)
        s = (Rc > 1) & (Rc < 4)
        mo, mt = np.median(a_o[s]), np.median(a_t[s])
        no, nt = np.median(n_o[s]), np.median(n_t[s])
        nz = int(np.sum(np.abs(obs[:, 2]) < Z_CUT))
        print(f"  {lab:26s} {len(obs):8,d} {nz:9,d} {mo:8.2f} {mt:7.2f} {100 * (mo / mt - 1):+6.1f}% "
              f"{no:+9.1f}° {nt:+8.1f}° {no - nt:+7.1f}°")
    print(f"  {'DR3 Miras (reference)':26s} {n_mira:8,d}")


if __name__ == "__main__":
    main()

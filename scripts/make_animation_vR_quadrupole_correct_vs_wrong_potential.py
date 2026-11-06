"""Animation of the project's central idea: the bar-frame <v_R> quadrupole of a stellar sample stays put only when the
orbits are integrated in the correct potential.

Stars, rotated by the true bar angle (27 deg), are integrated for 1 Gyr (41 snapshots, every 25 Myr) in two potentials
of the trial family (steady_state_likelihood.bar_aligned_potential, accuracy 1e-5):
  left   the true potential (|Omega_b| 26.02 km/s/kpc, S 1, bar_frac 1);
  right  a wrong one: --wrong pattern_speed (|Omega_b| 33, else true) or --wrong bar_strength (bar_frac 0.5).
Star sets (--stars):
  mock   the 50k stars of the perfect N-body mock (mock/nbody_mock_perfect_R8_N50k.hdf5);
  pool   every N-body star particle with R < 8 kpc (1,493,610; the pool the mock was drawn from), in random chunks of
         CHUNK stars integrated for all three potentials in turn (--stage integrate), so that any number of finished
         chunks gives the same random subset for every potential; --stage render uses the chunks finished so far.
Both panels start from the same stars, so the first frame is identical. Each frame (one snapshot): mean v_R in the bar
frame (bar along x), |z| < 1 kpc, folded with the bar's point symmetry ((x, y) -> (-x, -y); v_R unchanged), binned on a
--pix grid with the velocity sums and the counts each smoothed by a Gaussian of --smooth before dividing; grey where
the smoothed count is below N_MIN stars, and (default) beyond R = --rmax 6 kpc, the radial range the likelihood uses.
Beyond it the libraries grow an outer m = 2 pattern (3-7 km/s in 6-8 kpc) that turns by 40-60 deg over 1 Gyr even in
the correct potential (non-stationary there, likely trial-potential mismatch near corotation); the ring-averaged
<v_R> stays within +-1.6 km/s (correct potential). Contours guide the eye: solid at +-CONTOUR_LEVELS km/s for the
current frame, dashed at the same levels for t = 0 (both on maps smoothed to at least CONTOUR_SMOOTH kpc), so a
weaker pattern shows as solid contours shrinking inside the dashed ones. A strip below the maps plots the quadrupole
amplitude relative to t = 0 against time for both potentials.
The box in each panel gives the m = 2 term of <v_R> in 1 < R < 4 kpc
(least-squares fit over the stars on 1, cos 2phi, sin 2phi, cos 4phi, sin 4phi), as its amplitude relative to t = 0
and its turn since t = 0.
Trajectories: ~/mira_dr4_runs/animations/traj_<case>.npy (mock), traj_pool_<case>_cNN.npy (pool).
Output: figures/animations/vR_quadrupole_correct_vs_wrong_<pattern_speed|bar_strength>[_N<n>k].{mp4,gif}
"""
import argparse
import time
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.ndimage import gaussian_filter

from common import FIG_DIR, set_style, INK, INK_2, INK_MUTED, ACCENT
from galaxies import galaxy
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
import test1_nbody_qmc as old

OUT_DIR = Path.home() / "mira_dr4_runs" / "animations"
ANIM_DIR = FIG_DIR / "animations"
TRUTH = dict(Omega_abs=26.02, phi_b_deg=27.0, S=1.0, bar_frac=1.0)
CASES = {"truth": TRUTH,
         "pattern_speed": dict(TRUTH, Omega_abs=33.0),
         "bar_strength": dict(TRUTH, bar_frac=0.5)}
LABEL = {"truth": "correct potential\n|Ω_b| = 26.0 km/s/kpc, bar_frac = 1",
         "pattern_speed": "wrong pattern speed\n|Ω_b| = 33.0 km/s/kpc (+27%), bar_frac = 1",
         "bar_strength": "wrong bar strength\n|Ω_b| = 26.0 km/s/kpc, bar_frac = 0.5"}
CHUNK = 100_000
POOL_SEED = 7
Z_CUT = 1.0
N_MIN = 20
VLIM = 40.0
CONTOUR_LEVELS = (10.0, 20.0)
CONTOUR_SMOOTH = 0.4      # kpc, minimum smoothing of the maps the contours are drawn on
FPS = 6
HOLD = 6                 # extra frames shown at the start and at the end


def integrate(stars_mw, case):
    p = CASES[case]
    ssl.ACCURACY = old.ACCURACY
    pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], "nbody")
    return ssl.integrate(pot, ssl.rotate_z(stars_mw, np.radians(p["phi_b_deg"])), -abs(p["Omega_abs"]))[1]


def load_stars(path, R_max=None):
    with h5py.File(path, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    return w if R_max is None else w[np.hypot(w[:, 0], w[:, 1]) < R_max]


def mock_trajectory(case):
    path = OUT_DIR / f"traj_{case}.npy"
    if not path.exists():
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        np.save(path, integrate(load_stars(old.MOCK), case))
        print(f"integrated {case} in {time.time() - t0:.0f} s", flush=True)
    return np.load(path)


def integrate_pool():
    """All star particles with R < 8 kpc, in random chunks, each chunk for all three potentials before the next."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pool = load_stars(galaxy("nbody")["stars"], R_max=8.0)
    order = np.random.default_rng(POOL_SEED).permutation(len(pool))
    n_chunks = int(np.ceil(len(pool) / CHUNK))
    print(f"pool: {len(pool):,} stars, {n_chunks} chunks of {CHUNK:,}", flush=True)
    for k in range(n_chunks):
        for case in CASES:
            path = OUT_DIR / f"traj_pool_{case}_c{k:02d}.npy"
            if path.exists():
                continue
            t0 = time.time()
            np.save(path, integrate(pool[order[k * CHUNK:(k + 1) * CHUNK]], case))
            print(f"chunk {k + 1}/{n_chunks} {case}: {time.time() - t0:.0f} s", flush=True)


def pool_trajectories(cases):
    """Concatenated chunks finished for every one of the cases (the same stars for each)."""
    k = 0
    while all((OUT_DIR / f"traj_pool_{c}_c{k:02d}.npy").exists() for c in cases):
        k += 1
    if k == 0:
        raise SystemExit("no finished pool chunks yet")
    return {c: np.concatenate([np.load(OUT_DIR / f"traj_pool_{c}_c{j:02d}.npy") for j in range(k)]) for c in cases}


def vR_map(w, edges, smooth, rmax=None):
    """Folded, smoothed mean v_R map (NaN where < N_MIN stars or beyond rmax) for bar-frame points w (n, 6)."""
    w = w[np.abs(w[:, 2]) < Z_CUT]
    R = np.hypot(w[:, 0], w[:, 1])
    vR = (w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R
    x, y, v = np.r_[w[:, 0], -w[:, 0]], np.r_[w[:, 1], -w[:, 1]], np.r_[vR, vR]
    sig = smooth / (edges[1] - edges[0])
    n = gaussian_filter(np.histogram2d(x, y, bins=[edges, edges])[0], sig, mode="constant")
    s = gaussian_filter(np.histogram2d(x, y, bins=[edges, edges], weights=v)[0], sig, mode="constant")
    n_eff = n * 2 * np.pi * sig**2                     # stars within the kernel's effective area
    c = 0.5 * (edges[1:] + edges[:-1])
    outside = np.hypot(*np.meshgrid(c, c, indexing="ij")) > rmax if rmax else np.zeros(n.shape, bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((n_eff >= N_MIN) & ~outside, s / n, np.nan).T


def quadrupole(w):
    """m = 2 term of <v_R>(phi) in 1 < R < 4 kpc, |z| < 1: complex amplitude a_cos - i a_sin (least squares)."""
    R = np.hypot(w[:, 0], w[:, 1])
    s = (np.abs(w[:, 2]) < Z_CUT) & (R > 1) & (R < 4)
    phi = np.arctan2(w[s, 1], w[s, 0])
    vR = (w[s, 0] * w[s, 3] + w[s, 1] * w[s, 4]) / R[s]
    X = np.column_stack([np.ones_like(phi), np.cos(2 * phi), np.sin(2 * phi), np.cos(4 * phi), np.sin(4 * phi)])
    a = np.linalg.lstsq(X, vR, rcond=None)[0]
    return complex(a[1], -a[2])


def render(wrong, stars, pix, smooth, rmax):
    set_style()
    times = np.linspace(0.0, ssl.T_GYR, ssl.N_SNAP)
    cases = ("truth", wrong)
    trajs = {c: mock_trajectory(c) for c in cases} if stars == "mock" else pool_trajectories(cases)
    n_stars = len(trajs["truth"])
    edges = np.arange(-7.0, 7.0 + 1e-9, pix)
    maps = {c: [vR_map(t[:, k].astype(np.float64), edges, smooth, rmax) for k in range(len(times))] for c, t in trajs.items()}
    cs = max(smooth, CONTOUR_SMOOTH)
    cmaps = {c: [vR_map(t[:, k].astype(np.float64), edges, cs, rmax) for k in range(len(times))] for c, t in trajs.items()}
    centres = 0.5 * (edges[1:] + edges[:-1])
    quad = {c: np.array([quadrupole(t[:, k].astype(np.float64)) for k in range(len(times))]) for c, t in trajs.items()}

    fig = plt.figure(figsize=(12.4, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[4.2, 1.25], hspace=0.32, wspace=0.08, left=0.07, right=0.88,
                          top=0.86, bottom=0.12)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    strip = fig.add_subplot(gs[1, :])
    cmap = CMAP_DIV.copy()
    cmap.set_bad("#cfcfcc")
    levels = [-v for v in CONTOUR_LEVELS[::-1]] + list(CONTOUR_LEVELS)
    widths = [1.7, 1.0, 1.0, 1.7]                      # the +-20 contours drawn thicker
    meshes, boxes, lines = {}, {}, {}
    for ax, c in zip(axes, cases):
        meshes[c] = ax.pcolormesh(edges, edges, maps[c][0], cmap=cmap, vmin=-VLIM, vmax=VLIM, rasterized=True)
        ax.axhline(0, color=INK, lw=0.9, ls=":", alpha=0.7)
        ax.text(0.97, 0.505, "bar axis", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=INK)
        ax.contour(centres, centres, cmaps[c][0], levels=levels, colors=[INK_2], linewidths=widths, linestyles="--")
        ax.set_aspect("equal")
        lim = min(edges[-1], 1.05 * rmax) if rmax else edges[-1]
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("x [kpc] (bar frame)")
        ax.set_title(LABEL[c], fontsize=11, color=INK if c == "truth" else ACCENT)
        boxes[c] = ax.text(0.03, 0.03, "", transform=ax.transAxes, fontsize=9, color=INK_2, va="bottom",
                           bbox=dict(fc="white", ec="none", alpha=0.85, pad=3))
    axes[0].set_ylabel("y [kpc]")
    axes[1].set_yticklabels([])
    cb = fig.colorbar(meshes["truth"], cax=axes[1].inset_axes([1.04, 0.0, 0.045, 1.0]))
    cb.set_label("mean radial velocity ⟨v_R⟩ [km/s]")

    amp = {c: np.abs(quad[c]) / np.abs(quad[c][0]) for c in cases}
    tmyr = 1000 * times
    strip.axhline(1.0, color=INK_MUTED, lw=0.8, ls=":")
    strip.text(1005, 1.0, " start", va="center", fontsize=8, color=INK_MUTED)
    curves, dots = {}, {}
    for c, colr in zip(cases, (INK, ACCENT)):
        curves[c], = strip.plot([], [], color=colr, lw=2.0, label=LABEL[c].split("\n")[0])
        dots[c], = strip.plot([], [], "o", color=colr, ms=6)
    strip.set_xlim(0, 1000)
    strip.set_ylim(min(0.4, min(a.min() for a in amp.values()) - 0.05), 1.12)
    strip.set_xlabel("time [Myr]")
    strip.set_ylabel("quadrupole\namplitude / start", fontsize=9)
    strip.set_title("strength of the ⟨v_R⟩ quadrupole (m = 2, 1 < R < 4 kpc) relative to t = 0", fontsize=9.5,
                    color=INK_2)
    strip.legend(loc="lower left", fontsize=8.5, frameon=False, ncol=2)

    stamp = fig.text(0.475, 0.955, "", ha="center", fontsize=14, color=INK)
    fig.text(0.475, 0.915, f"same {n_stars:,} stars, integrated forward; the pattern stays put only in the correct potential",
             ha="center", fontsize=10, color=INK_2)
    fig.text(0.475, 0.025, f"|z| < 1 kpc, folded (x, y) → (−x, −y), {pix:g} kpc pixels, Gaussian-smoothed σ = {smooth:g} kpc"
             + (f", R < {rmax:g} kpc" if rmax else "") + "; one snapshot per frame, every 25 Myr\ncontours at "
             + " and ".join(f"±{v:g}" for v in CONTOUR_LEVELS) + " km/s: solid = this frame, dashed = t = 0 (thick = ±"
             + f"{CONTOUR_LEVELS[-1]:g})", ha="center", fontsize=8.5, color=INK_2)

    order = [0] * HOLD + list(range(len(times))) + [len(times) - 1] * HOLD

    def draw(frame):
        k = order[frame]
        for c in meshes:
            meshes[c].set_array(maps[c][k].ravel())
            if c in lines:
                lines[c].remove()
            lines[c] = meshes[c].axes.contour(centres, centres, cmaps[c][k], levels=levels,
                                              colors=["#1f4e96"] * len(CONTOUR_LEVELS) + ["#b2182b"] * len(CONTOUR_LEVELS),
                                              linewidths=widths)
            q0, q = quad[c][0], quad[c][k]
            turn = np.degrees(np.angle(q / q0)) / 2          # rotation of the m = 2 pattern since t = 0
            boxes[c].set_text(f"quadrupole (1 < R < 4 kpc):\namplitude {abs(q) / abs(q0):.2f} × start\n"
                              f"turned {turn:+.1f}° since t = 0")
            curves[c].set_data(tmyr[:k + 1], amp[c][:k + 1])
            dots[c].set_data([tmyr[k]], [amp[c][k]])
        stamp.set_text(f"t = {1000 * times[k]:.0f} Myr")
        return list(meshes.values()) + list(boxes.values()) + [stamp]

    draw(0)

    ANIM_DIR.mkdir(parents=True, exist_ok=True)
    anim = animation.FuncAnimation(fig, draw, frames=len(order), interval=1000 / FPS, blit=False)
    suffix = "" if stars == "mock" else f"_N{n_stars // 1000}k"
    base = ANIM_DIR / f"vR_quadrupole_correct_vs_wrong_{wrong}{suffix}"
    anim.save(f"{base}.mp4", writer=animation.FFMpegWriter(fps=FPS, bitrate=6000), dpi=200)
    anim.save(f"{base}.gif", writer=animation.PillowWriter(fps=FPS), dpi=90)
    for k in (0, len(times) // 2, len(times) - 1):            # stills for checking
        draw(order.index(k) if k else 0)
        fig.savefig(f"{base}_t{1000 * times[k]:04.0f}Myr.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    for c in quad:
        q = quad[c]
        print(f"{c:14s} amplitude at 0.5 / 1 Gyr: {abs(q[20]) / abs(q[0]):.2f} / {abs(q[-1]) / abs(q[0]):.2f}; "
              f"turn at 0.5 / 1 Gyr: {np.degrees(np.angle(q[20] / q[0])) / 2:+.1f} / {np.degrees(np.angle(q[-1] / q[0])) / 2:+.1f} deg")
    print(f"saved {base}.mp4 and .gif ({n_stars:,} stars)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", choices=("render", "integrate"), default="render", help="integrate: the pool chunks")
    ap.add_argument("--wrong", choices=("pattern_speed", "bar_strength"))
    ap.add_argument("--stars", choices=("mock", "pool"), default="mock")
    ap.add_argument("--pix", type=float, default=0.25, help="pixel size [kpc]")
    ap.add_argument("--smooth", type=float, default=0.5, help="Gaussian smoothing sigma [kpc]")
    ap.add_argument("--rmax", type=float, default=6.0, help="show R < rmax only (0: no limit) [kpc]")
    args = ap.parse_args()
    if args.stage == "integrate":
        integrate_pool()
    else:
        render(args.wrong, args.stars, args.pix, args.smooth, args.rmax or None)


if __name__ == "__main__":
    main()

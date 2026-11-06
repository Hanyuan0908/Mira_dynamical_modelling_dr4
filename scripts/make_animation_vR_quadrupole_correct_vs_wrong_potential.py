"""Animation of the project's central idea: the bar-frame <v_R> quadrupole of a stellar sample stays put only when the
orbits are integrated in the correct potential.

The 50k stars of the perfect N-body mock (mock/nbody_mock_perfect_R8_N50k.hdf5) are rotated by the true bar angle
(27 deg) and integrated for 1 Gyr (41 snapshots, every 25 Myr) in two potentials of the trial family
(steady_state_likelihood.bar_aligned_potential, accuracy 1e-5):
  left   the true potential (|Omega_b| 26.02 km/s/kpc, S 1, bar_frac 1);
  right  a wrong one: --wrong pattern_speed (|Omega_b| 33, else true) or --wrong bar_strength (bar_frac 0.5).
Both start from the same stars, so the first frame is identical. Each frame (one snapshot): mean v_R in the bar
frame (bar along x), |z| < 1 kpc, folded with the bar's point symmetry ((x, y) -> (-x, -y); v_R unchanged), binned on
a 0.25 kpc grid with the velocity sums and the counts each smoothed by a Gaussian of SMOOTH_KPC before dividing;
grey where the smoothed count is below N_MIN stars. The box in each panel gives the m = 2 term of <v_R> in 1 < R < 4 kpc (least-squares fit
over the stars on 1, cos 2phi, sin 2phi, cos 4phi, sin 4phi), as its amplitude relative to t = 0 and its turn since
t = 0.
Trajectories: ~/mira_dr4_runs/animations/traj_<case>.npy (re-used if present).
Output: figures/animations/vR_quadrupole_correct_vs_wrong_<pattern_speed|bar_strength>.{mp4,gif}
"""
import argparse
import time
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.ndimage import gaussian_filter

from common import FIG_DIR, set_style, INK, INK_2, ACCENT
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
EDGES = np.arange(-7.0, 7.0 + 1e-9, 0.25)
SMOOTH_KPC = 0.5
Z_CUT = 1.0
N_MIN = 20
VLIM = 40.0
FPS = 6
HOLD = 6                 # extra frames shown at the start and at the end


def trajectory(case):
    path = OUT_DIR / f"traj_{case}.npy"
    if path.exists():
        return np.load(path)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p = CASES[case]
    with h5py.File(old.MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    ssl.ACCURACY = old.ACCURACY
    t0 = time.time()
    pot = ssl.bar_aligned_potential(p["bar_frac"], p["S"], "nbody")
    _, traj = ssl.integrate(pot, ssl.rotate_z(stars, np.radians(p["phi_b_deg"])), -abs(p["Omega_abs"]))
    np.save(path, traj)
    print(f"integrated {case} in {time.time() - t0:.0f} s", flush=True)
    return traj


def vR_map(w):
    """Folded mean v_R map (NaN where < N_MIN stars) for bar-frame phase-space points w (n, 6)."""
    w = w[np.abs(w[:, 2]) < Z_CUT]
    R = np.hypot(w[:, 0], w[:, 1])
    vR = (w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R
    x, y, v = np.r_[w[:, 0], -w[:, 0]], np.r_[w[:, 1], -w[:, 1]], np.r_[vR, vR]
    sig = SMOOTH_KPC / (EDGES[1] - EDGES[0])
    n = gaussian_filter(np.histogram2d(x, y, bins=[EDGES, EDGES])[0], sig, mode="constant")
    s = gaussian_filter(np.histogram2d(x, y, bins=[EDGES, EDGES], weights=v)[0], sig, mode="constant")
    n_eff = n * 2 * np.pi * sig**2                     # stars within the kernel's effective area
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n_eff >= N_MIN, s / n, np.nan).T


def quadrupole(w):
    """m = 2 term of <v_R>(phi) in 1 < R < 4 kpc, |z| < 1: complex amplitude a_cos - i a_sin (least squares)."""
    R = np.hypot(w[:, 0], w[:, 1])
    s = (np.abs(w[:, 2]) < Z_CUT) & (R > 1) & (R < 4)
    phi = np.arctan2(w[s, 1], w[s, 0])
    vR = (w[s, 0] * w[s, 3] + w[s, 1] * w[s, 4]) / R[s]
    X = np.column_stack([np.ones_like(phi), np.cos(2 * phi), np.sin(2 * phi), np.cos(4 * phi), np.sin(4 * phi)])
    a = np.linalg.lstsq(X, vR, rcond=None)[0]
    return complex(a[1], -a[2])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wrong", choices=("pattern_speed", "bar_strength"), required=True)
    args = ap.parse_args()
    set_style()
    times = np.linspace(0.0, ssl.T_GYR, ssl.N_SNAP)
    trajs = {c: trajectory(c) for c in ("truth", args.wrong)}
    maps = {c: [vR_map(t[:, k].astype(np.float64)) for k in range(len(times))] for c, t in trajs.items()}
    quad = {c: np.array([quadrupole(t[:, k].astype(np.float64)) for k in range(len(times))]) for c, t in trajs.items()}

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.6), gridspec_kw=dict(wspace=0.08))
    cmap = CMAP_DIV.copy()
    cmap.set_bad("#cfcfcc")
    meshes, boxes = {}, {}
    for ax, c in zip(axes, ("truth", args.wrong)):
        meshes[c] = ax.pcolormesh(EDGES, EDGES, maps[c][0], cmap=cmap, vmin=-VLIM, vmax=VLIM, rasterized=True)
        ax.axhline(0, color=INK, lw=0.9, ls="--", alpha=0.7)
        ax.text(6.7, 0.25, "bar axis of the potential", ha="right", va="bottom", fontsize=8.5, color=INK)
        ax.set_aspect("equal")
        ax.set_xlim(EDGES[0], EDGES[-1])
        ax.set_ylim(EDGES[0], EDGES[-1])
        ax.set_xlabel("x [kpc] (bar frame)")
        ax.set_title(LABEL[c], fontsize=11, color=INK if c == "truth" else ACCENT)
        boxes[c] = ax.text(0.03, 0.03, "", transform=ax.transAxes, fontsize=9, color=INK_2, va="bottom",
                           bbox=dict(fc="white", ec="none", alpha=0.85, pad=3))
    axes[0].set_ylabel("y [kpc]")
    axes[1].set_yticklabels([])
    cb = fig.colorbar(meshes["truth"], ax=axes, shrink=0.82, pad=0.02)
    cb.set_label("mean radial velocity ⟨v_R⟩ [km/s]")
    stamp = fig.text(0.44, 0.955, "", ha="center", fontsize=14, color=INK)
    fig.text(0.44, 0.905, "same 50,000 stars, integrated forward; the pattern stays put only in the correct potential",
             ha="center", fontsize=10, color=INK_2)
    fig.text(0.44, 0.03, f"|z| < 1 kpc, folded (x, y) → (−x, −y), Gaussian-smoothed σ = {SMOOTH_KPC} kpc; one snapshot per "
             f"frame, every 25 Myr", ha="center", fontsize=8.5, color=INK_2)

    order = [0] * HOLD + list(range(len(times))) + [len(times) - 1] * HOLD

    def draw(frame):
        k = order[frame]
        for c in meshes:
            meshes[c].set_array(maps[c][k].ravel())
            q0, q = quad[c][0], quad[c][k]
            turn = np.degrees(np.angle(q / q0)) / 2          # rotation of the m = 2 pattern since t = 0
            boxes[c].set_text(f"quadrupole (1 < R < 4 kpc):\namplitude {abs(q) / abs(q0):.2f} × start\n"
                              f"turned {turn:+.1f}° since t = 0")
        stamp.set_text(f"t = {1000 * times[k]:.0f} Myr")
        return list(meshes.values()) + list(boxes.values()) + [stamp]

    ANIM_DIR.mkdir(parents=True, exist_ok=True)
    anim = animation.FuncAnimation(fig, draw, frames=len(order), interval=1000 / FPS, blit=False)
    base = ANIM_DIR / f"vR_quadrupole_correct_vs_wrong_{args.wrong}"
    anim.save(f"{base}.mp4", writer=animation.FFMpegWriter(fps=FPS, bitrate=3000), dpi=150)
    anim.save(f"{base}.gif", writer=animation.PillowWriter(fps=FPS), dpi=80)
    for k in (0, len(times) // 2, len(times) - 1):            # stills for checking
        draw(order.index(k) if k else 0)
        fig.savefig(f"{base}_t{1000 * times[k]:04.0f}Myr.png", dpi=110, bbox_inches="tight")
    for c in quad:
        q = quad[c]
        print(f"{c:14s} amplitude at 0.5 / 1 Gyr: {abs(q[20]) / abs(q[0]):.2f} / {abs(q[-1]) / abs(q[0]):.2f}; "
              f"turn at 0.5 / 1 Gyr: {np.degrees(np.angle(q[20] / q[0])) / 2:+.1f} / {np.degrees(np.angle(q[-1] / q[0])) / 2:+.1f} deg")
    print(f"saved {base}.mp4 and .gif")


if __name__ == "__main__":
    main()

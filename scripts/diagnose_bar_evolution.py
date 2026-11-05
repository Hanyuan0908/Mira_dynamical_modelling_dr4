"""Bar strength, growth and buckling against time for the N-body model, and maps of one snapshot.

t = 2 is an outlier in every test we have run — the largest pattern-speed bias, the largest bar-angle
bias and the largest zero-error offset of the kinematic estimators. This asks what is different about
it, by tracking the bar through the whole run and mapping the suspect snapshot against a control.

Diagnostics, all on the disc particles, measured in each snapshot's own bar frame:

  A2(R)        |sum e^(2 i phi)| / N in radial bins; A2_max is the bar strength, R_bar the radius
               where A2 falls to 0.5 A2_max beyond the peak
  dA2/dt       how fast the bar is growing or weakening - a bar caught mid-growth is not in the
               steady state the modelling assumes
  A1_z         |sum z e^(i phi)| / sum |z| inside R_bar. The buckling instability bends one half of
               the bar up and the other down, which is exactly an m = 1 pattern in z(phi). This is
               the standard signature; it spikes during buckling and decays after.
  A2_z         the m = 2 equivalent, which responds to a boxy/peanut bulge rather than to bending
  <|z|>        disc thickness inside R_bar: it jumps when buckling converts a flat bar into a B/P

Usage: python diagnose_bar_evolution.py [--maps 2.0 5.0] [--tmax 8] [--no-series]
Figures: figures/nbody_bar_evolution.png, figures/nbody_bar_maps_t<T>.png
Table:   results/nbody_bar_evolution.csv
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP, \
    CMAP_SEQ_W
import nbody

R_EDGES = np.arange(0.0, 10.01, 0.25)
Rc = 0.5 * (R_EDGES[1:] + R_EDGES[:-1])


def load_disc(t):
    """Disc particles of one snapshot, centred and rotated so the bar lies along x."""
    d = nbody.read_nemo(nbody.snapshot(t), tags=("Mass", "Position", "Velocity"))
    comp = nbody.components(d["Mass"])
    x, v = d["Position"][comp == 0], d["Velocity"][comp == 0]
    c = np.median(x[np.linalg.norm(x, axis=1) < 5], axis=0)
    x = x - c
    v = v - np.median(v[np.linalg.norm(x, axis=1) < 5], axis=0)
    R = np.hypot(x[:, 0], x[:, 1])
    m = (R > 0.5) & (R < 5) & (np.abs(x[:, 2]) < 2)          # the bar region, for the phase
    phi = np.arctan2(x[m, 1], x[m, 0])
    ang = 0.5 * np.angle(np.sum(np.exp(2j * phi)))
    c_, s_ = np.cos(-ang), np.sin(-ang)
    xr = np.column_stack([c_ * x[:, 0] - s_ * x[:, 1], s_ * x[:, 0] + c_ * x[:, 1], x[:, 2]])
    vr = np.column_stack([c_ * v[:, 0] - s_ * v[:, 1], s_ * v[:, 0] + c_ * v[:, 1], v[:, 2]])
    return xr, vr


def bar_metrics(x):
    R = np.hypot(x[:, 0], x[:, 1])
    phi = np.arctan2(x[:, 1], x[:, 0])
    thin = np.abs(x[:, 2]) < 2
    A2 = np.full(len(Rc), np.nan)
    for i in range(len(Rc)):
        s = thin & (R >= R_EDGES[i]) & (R < R_EDGES[i + 1])
        if s.sum() > 200:
            A2[i] = np.abs(np.sum(np.exp(2j * phi[s]))) / s.sum()
    if np.all(np.isnan(A2)):
        return dict(A2_max=np.nan, R_bar=np.nan, A1_z=np.nan, A2_z=np.nan, absz=np.nan)
    k = int(np.nanargmax(A2))
    A2max = float(A2[k])
    half = np.where(np.arange(len(Rc)) > k, A2 < 0.5 * A2max, False)
    Rbar = float(Rc[np.argmax(half)]) if half.any() else float(Rc[-1])
    inb = (R < Rbar) & (R > 0.5)
    z, p = x[inb, 2], phi[inb]
    sz = np.sum(np.abs(z))
    return dict(A2_max=A2max, R_bar=Rbar,
                A1_z=float(np.abs(np.sum(z * np.exp(1j * p))) / sz),
                A2_z=float(np.abs(np.sum(z * np.exp(2j * p))) / sz),
                absz=float(np.mean(np.abs(z))), A2_prof=A2)


def maps(t, tag=""):
    x, v = load_disc(t)
    R = np.hypot(x[:, 0], x[:, 1])
    vR = (x[:, 0] * v[:, 0] + x[:, 1] * v[:, 1]) / np.maximum(R, 1e-9)
    vtot = np.linalg.norm(v, axis=1)
    q = np.abs(vR) / np.maximum(vtot, 1e-9)
    thin = np.abs(x[:, 2]) < 1.5
    b = np.linspace(-8, 8, 121)
    bz = np.linspace(-3, 3, 61)
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.5))

    H, *_ = np.histogram2d(x[thin, 0], x[thin, 1], bins=[b, b])
    ax[0].imshow(np.log10(H.T + 1), origin="lower", extent=[-8, 8, -8, 8], cmap=CMAP_SEQ_W)
    ax[0].set_title("face-on density (log)", color=INK); ax[0].set_xlabel("x [kpc]"); ax[0].set_ylabel("y [kpc]")

    He, *_ = np.histogram2d(x[np.abs(x[:, 1]) < 1.5, 0], x[np.abs(x[:, 1]) < 1.5, 2], bins=[b, bz])
    ax[1].imshow(np.log10(He.T + 1), origin="lower", extent=[-8, 8, -3, 3], cmap=CMAP_SEQ_W,
                 aspect="auto")
    ax[1].set_title("edge-on density, |y| < 1.5 (log)", color=INK)
    ax[1].set_xlabel("x [kpc]"); ax[1].set_ylabel("z [kpc]")

    for j, (val, lab, cm, lim) in enumerate(
            ((vR, r"$\langle v_R\rangle$ [km s$^{-1}$]", "RdBu_r", 40),
             (q, r"$\langle |v_R|/v_{\rm tot}\rangle$", "viridis", None))):
        n, *_ = np.histogram2d(x[thin, 0], x[thin, 1], bins=[b, b])
        s, *_ = np.histogram2d(x[thin, 0], x[thin, 1], bins=[b, b], weights=val[thin])
        M = np.where(n > 20, s / np.maximum(n, 1), np.nan)
        kw = dict(norm=TwoSlopeNorm(0, -lim, lim)) if lim else {}
        im = ax[2 + j].imshow(M.T, origin="lower", extent=[-8, 8, -8, 8], cmap=cm, **kw)
        plt.colorbar(im, ax=ax[2 + j], fraction=0.046)
        ax[2 + j].set_title(lab, color=INK); ax[2 + j].set_xlabel("x [kpc]")
    for a in ax:
        a.set_xticks([-8, -4, 0, 4, 8])
    fig.suptitle(f"N-body model at t = {t:g} ({t * nbody.TIME_UNIT_GYR:.2f} Gyr), "
                 f"disc particles, bar along x", color=INK, y=1.0)
    fig.tight_layout()
    savefig(fig, f"nbody_bar_maps_t{t:g}{tag}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maps", nargs="*", type=float, default=[2.0, 5.0])
    ap.add_argument("--tmax", type=float, default=8.0)
    ap.add_argument("--no-series", action="store_true")
    args = ap.parse_args()
    set_style()

    if not args.no_series:
        ts = [t for t in nbody.available_times() if 0.25 <= t <= args.tmax]
        rows = []
        for t in ts:
            x, _ = load_disc(t)
            m = bar_metrics(x)
            m.pop("A2_prof", None)
            rows.append(dict(t=t, t_Gyr=t * nbody.TIME_UNIT_GYR, **m))
            print(f"  t={t:5.3g}  A2 {m['A2_max']:.3f}  R_bar {m['R_bar']:.2f}  "
                  f"A1_z {m['A1_z']:.4f}  A2_z {m['A2_z']:.4f}  <|z|> {m['absz']:.3f}", flush=True)
        D = pd.DataFrame(rows)
        D["dA2dt"] = np.gradient(D.A2_max, D.t_Gyr)
        D.to_csv(RESULTS_DIR / "nbody_bar_evolution.csv", index=False)
        print(f"wrote {rel(RESULTS_DIR / 'nbody_bar_evolution.csv')}")

        fig, ax = plt.subplots(4, 1, figsize=(9.5, 10), sharex=True,
                               gridspec_kw=dict(hspace=0.12))
        for a, (c, lab, col) in zip(ax, (("A2_max", "bar strength $A_2^{\\rm max}$", BLUE_RAMP[5]),
                                         ("dA2dt", "growth rate $dA_2/dt$ [Gyr$^{-1}$]", "#1baf7a"),
                                         ("A1_z", "buckling  $A_{1,z}$  (bending mode)", ACCENT),
                                         ("absz", r"$\langle |z|\rangle$ inside $R_{\rm bar}$ [kpc]",
                                          "#d4a017"))):
            a.plot(D.t, D[c], "-", color=col, lw=1.8)
            a.set_ylabel(lab, fontsize=9.5)
            for tv in args.maps:
                a.axvline(tv, color=INK_MUTED, lw=1, ls="--")
            if c == "dA2dt":
                a.axhline(0, color=INK_MUTED, lw=1, ls=":")
        ax[-1].set_xlabel("snapshot model time  (1 unit = 0.978 Gyr)")
        for tv in args.maps:
            ax[0].annotate(f"t = {tv:g}", (tv, ax[0].get_ylim()[1]), xytext=(0, 3),
                           textcoords="offset points", ha="center", fontsize=9, color=INK_2)
        fig.suptitle("Bar evolution of the N-body model", color=INK, y=0.995)
        savefig(fig, "nbody_bar_evolution")

    for t in args.maps:
        maps(t)


if __name__ == "__main__":
    main()

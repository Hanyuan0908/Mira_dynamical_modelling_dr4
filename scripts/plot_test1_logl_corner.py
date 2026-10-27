"""Corner plot of the Test 1 log-likelihood surface over (|Omega_b|, phi_b, S, bar_frac).

Upper triangle: the evaluated Sobol points, coloured by logL - max. Lower triangle and diagonal: the
profile log-likelihood of the GP emulator (test1_gp_emulator.Emulator, fitted to the Sobol points only),
i.e. the maximum over the parameters not shown, found on a grid. Truth, best evaluated point and
emulator maximum are marked.
Usage: python plot_test1_logl_corner.py [--run au18|nbody] [--column logL|logL_gauss|chi2_...]
       [--rmin R] [--rmax R]   (radial cuts need per-star files)
"""
import argparse

import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
from scipy.optimize import minimize

from common import set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge
from test1_gp_emulator import Emulator, to_unit, from_unit, KEYS

LABELS = [r"$|\Omega_b|$", r"$\phi_b$ [deg]", r"$S$", r"bar_frac"]
LEVELS = [-5000, -2000, -1000, -500, -200, -100, -50, -20, 0.1]
LOG_TICKS = {2: [0.8, 1.0, 1.2, 1.4, 1.6], 3: [0.5, 0.7, 1.0, 1.4]}


def log_axis(axis, j):
    axis.set_major_locator(FixedLocator(LOG_TICKS[j]))
    axis.set_major_formatter(FormatStrFormatter("%g"))
    axis.set_minor_formatter(NullFormatter())
N2, NOTHER, N1 = 25, 12, 60


def phys(j, u):
    """Unit-cube coordinate -> physical value of parameter j."""
    v = ge.LO[j] + np.asarray(u) * (ge.HI[j] - ge.LO[j])
    return np.exp(v) if j >= 2 else v


def profile(emu, dims):
    """Max of the emulator over the dimensions not in `dims`, on a grid."""
    g = np.linspace(0, 1, N2 if len(dims) == 2 else N1)
    go = np.linspace(0, 1, NOTHER)
    others = [d for d in range(4) if d not in dims]
    mesh = np.meshgrid(*([g] * len(dims) + [go] * len(others)), indexing="ij")
    U = np.zeros((mesh[0].size, 4))
    for d, m in zip(list(dims) + others, mesh):
        U[:, d] = m.ravel()
    pred = emu.predict(U).reshape(mesh[0].shape)
    return g, pred.max(axis=tuple(range(len(dims), 4)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="au18", choices=list(ge.RUNS))
    ap.add_argument("--column", default="logL")
    ap.add_argument("--rmin", type=float, default=None)
    ap.add_argument("--rmax", type=float, default=None)
    args = ap.parse_args()
    set_style()
    run = ge.use_run(args.run)
    TRUTH = run["truth"]
    df, tag = ge.load_table(args.column, args.rmin, args.rmax)
    tr = df[df.kind == "qmc"]
    u, y = to_unit(tr), tr["logL"].values
    emu = Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        best = r if best is None or r.fun < best.fun else best
    ymax = max(-best.fun, y.max())
    u_truth = to_unit(pd.DataFrame([TRUTH]))[0]
    u_bestpt = u[np.argmax(y)]
    t_row = df[df.kind == "truth"]
    print(f"{len(tr)} Sobol points; best evaluated logL - truth logL = "
          f"{y.max() - t_row.logL.iloc[0]:+.1f}" if len(t_row) else "")
    print("emulator max:", dict(zip(KEYS, np.round(from_unit(best.x[None])[0], 3))),
          f"(logL_max - logL_truth,emulated = {-best.fun - emu.predict(u_truth[None])[0]:+.1f})")

    cmap = LinearSegmentedColormap.from_list("lr", ["#f7fafe"] + BLUE_RAMP)
    norm = BoundaryNorm(LEVELS, cmap.N)
    fig, axes = plt.subplots(4, 4, figsize=(12.5, 12.5), gridspec_kw=dict(hspace=0.08, wspace=0.08))
    marks = [(u_truth, ACCENT, "*", 15, "truth"), (u_bestpt, INK, "D", 6, "best evaluated"),
             (best.x, "#1baf7a", "X", 9, "emulator max")]
    for i in range(4):
        for j in range(4):
            ax = axes[i, j]
            if j > i:            # upper: raw evaluations
                sc = ax.scatter(phys(j, u[:, j]), phys(i, u[:, i]), c=np.clip(y - ymax, LEVELS[0], 0), cmap=cmap,
                                norm=norm, s=28, edgecolor=INK_2, lw=0.3)
                for um, c, mk, ms, _ in marks:
                    ax.plot(phys(j, um[j]), phys(i, um[i]), mk, color=c, ms=ms, mec="white", mew=0.6)
            elif j < i:          # lower: emulator profile
                g, P = profile(emu, (j, i))
                ax.contourf(phys(j, g), phys(i, g), (P - ymax).T, levels=LEVELS, cmap=cmap, norm=norm, extend="min")
                ax.contour(phys(j, g), phys(i, g), (P - ymax).T, levels=[-50, -20], colors=INK_2, linewidths=0.6)
                ax.plot(phys(j, u[:, j]), phys(i, u[:, i]), ".", color=INK_2, ms=2)
                for um, c, mk, ms, _ in marks:
                    ax.plot(phys(j, um[j]), phys(i, um[i]), mk, color=c, ms=ms, mec="white", mew=0.6)
            else:                # diagonal: 1D profile + evaluated points
                g, P = profile(emu, (i,))
                ax.plot(phys(i, g), P - ymax, color=BLUE_RAMP[5], lw=1.8)
                ax.plot(phys(i, u[:, i]), y - ymax, ".", color=INK_2, ms=4)
                ax.axvline(phys(i, u_truth[i]), color=ACCENT, lw=1.2, ls="--")
                ax.axvline(phys(i, best.x[i]), color="#1baf7a", lw=1.0, ls=":")
                ax.set_ylim(-1500, 100)
                ax.yaxis.tick_right()
                ax.tick_params(labelsize=7)
            if i == 3:
                ax.set_xlabel(LABELS[j])
            else:
                ax.tick_params(labelbottom=False)
            if j == 0 and i > 0:
                ax.set_ylabel(LABELS[i])
            elif j > 0 and i != j:
                ax.tick_params(labelleft=False)
            if i != j:
                ax.set_xlim(phys(j, 0), phys(j, 1))
                ax.set_ylim(phys(i, 0), phys(i, 1))
                if j >= 2:
                    ax.set_xscale("log")
                    log_axis(ax.xaxis, j)
                if i >= 2:
                    ax.set_yscale("log")
                    log_axis(ax.yaxis, i)
            elif i >= 2:
                ax.set_xscale("log")
                log_axis(ax.xaxis, i)
            if i == j and i < 3:
                ax.tick_params(labelbottom=False)
            if j > i or (j < i and j > 0):
                ax.tick_params(labelleft=False)
    axes[0, 0].set_ylabel(r"profile $\log L - \max$", fontsize=8)
    cb = fig.colorbar(sc, ax=axes, shrink=0.35, location="right", pad=0.06)
    cb.set_label(r"$\log L - \max$")
    handles = [plt.Line2D([], [], marker=mk, color=c, ls="", ms=ms * 0.8, label=lab) for _, c, mk, ms, lab in marks]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.86, 0.9), fontsize=9)
    fig.suptitle(f"{run['prefix']} log-likelihood ({args.column}{', ' + tag[1:] if tag else ''}): "
                 f"upper = {len(tr)} Sobol evaluations, lower/diagonal = GP profile (max over hidden parameters); "
                 "contours at −20, −50", y=0.92, fontsize=10.5)
    savefig(fig, f"{run['prefix']}_logl_corner" + tag)


if __name__ == "__main__":
    main()

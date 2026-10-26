"""Column-normalised period distribution of DR3 Miras as a function of Galactocentric R,
i.e. P(period | R), with the number of stars per R column shown above."""
import numpy as np
import matplotlib.pyplot as plt

from common import (load_dr3_mira, set_style, savefig, CMAP_SEQ, ACCENT, INK_2,
                    INK_MUTED, R0_KPC)

R_EDGES = np.arange(0.0, 16.0 + 1e-9, 0.25)       # kpc
P_EDGES = np.arange(80.0, 820.0 + 1e-9, 10.0)      # d; covers the full period range
P_PLOT = (80.0, 650.0)
MIN_PER_COLUMN = 30                                 # columns with fewer stars are masked


def main():
    set_style()
    df = load_dr3_mira()

    H, _, _ = np.histogram2d(df["R"], df["period"], bins=[R_EDGES, P_EDGES])
    n_col = H.sum(axis=1)
    dP = np.diff(P_EDGES)
    with np.errstate(invalid="ignore", divide="ignore"):
        cond = H / (n_col[:, None] * dP[None, :])   # P(period | R) in d^-1
    cond[n_col < MIN_PER_COLUMN] = np.nan
    cond[cond == 0] = np.nan                          # empty cells drawn as white

    # running percentiles of period per R column
    R_mid = 0.5 * (R_EDGES[1:] + R_EDGES[:-1])
    idx = np.digitize(df["R"], R_EDGES) - 1
    pct = np.full((len(R_mid), 3), np.nan)
    for i in range(len(R_mid)):
        if n_col[i] >= MIN_PER_COLUMN:
            pct[i] = np.percentile(df["period"].values[idx == i], [16, 50, 84])

    fig, (ax_n, ax) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 3.2], hspace=0.05))

    ax_n.stairs(n_col, R_EDGES, color=INK_2, lw=1.2)
    ax_n.axhline(MIN_PER_COLUMN, color=INK_MUTED, lw=0.8, ls=":")
    ax_n.text(0.3, MIN_PER_COLUMN * 1.15, f"columns with N < {MIN_PER_COLUMN} masked",
              ha="left", va="bottom", fontsize=8, color=INK_MUTED)
    ax_n.set_yscale("log")
    ax_n.set_ylabel(r"$N$ per $\Delta R$ = 0.25 kpc")
    ax_n.set_title(f"Gaia DR3 Miras (N = {len(df):,}): P(period | R)")

    vmax = np.nanpercentile(cond[n_col >= 10 * MIN_PER_COLUMN], 99.5)   # not set by noisy columns
    m = ax.pcolormesh(R_EDGES, P_EDGES, cond.T, cmap=CMAP_SEQ, vmin=0, vmax=vmax,
                      rasterized=True)
    ax.plot(R_mid, pct[:, 1], color=ACCENT, lw=2, label="median")
    ax.plot(R_mid, pct[:, 0], color=ACCENT, lw=1, ls="--", label="16th / 84th pct")
    ax.plot(R_mid, pct[:, 2], color=ACCENT, lw=1, ls="--")
    ax.legend(loc="upper right", fontsize=8, labelcolor=INK_2)

    for a in (ax_n, ax):
        a.axvline(R0_KPC, color=INK_MUTED, lw=0.8, ls="-.")
    ax.text(R0_KPC + 0.15, P_PLOT[0] + 12, r"$R_\odot$", color=INK_2, fontsize=9)

    ax.set_xlim(R_EDGES[0], R_EDGES[-1])
    ax.set_ylim(*P_PLOT)
    ax.set_xlabel(r"Galactocentric $R$ [kpc]")
    ax.set_ylabel("Period [d]")

    cax = ax.inset_axes([1.02, 0.0, 0.03, 1.0])
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r"$p(P\,|\,R)$ [d$^{-1}$]")
    cb.outline.set_linewidth(0.5)

    savefig(fig, "dr3_mira_period_vs_R")


if __name__ == "__main__":
    main()

"""An equilibrium test: does the <v_R> quadrupole align with the density bar?

For a steady bar, stars on closed x1 orbits sit at apocentre on the major axis and pericentre on the
minor axis, so v_R = 0 on both: the <v_R> quadrupole nodes MUST lie along the bar axes. A nonzero
offset means the orbits are not closed in the rotating frame - there is a net radial mass flux -
which is exactly the steady state the orbit-library modelling assumes.

Measured on ALL stars with no selection and no distance error, so it is a property of the galaxy,
not of the observation. It is also measurable in real data: the density bar angle comes from
photometry and the <v_R> node from kinematics, so a survey can apply this test to its own galaxy.

Run on all twelve campaign snapshots, so the criterion is defined independently of which snapshot
turned out to be an outlier.

Usage: python diagnose_equilibrium_offset.py [--kin 0.7]
Table: results/nbody_equilibrium_offset.csv;  figure: figures/test2_equilibrium_offset.png
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from scipy import stats

from common import RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, rel, INK, INK_2, INK_MUTED, \
    ACCENT, BLUE_RAMP
from galaxies import galaxy, nbody_name
import test1_nbody_qmc_focused_fourier_stars as fs
from bar_angle_vecsum_verification import ring_coeffs
from nbody import SNAPSHOT_DIR
import steady_state_likelihood as ssl

TIMES = [2.0, 2.375, 2.625, 3.0, 3.375, 3.625, 4.0, 4.5, 5.0, 5.5, 6.25, 7.0]
BIAS = {2.0: -4.1, 2.375: +2.0, 2.625: +0.3, 3.0: +1.1, 3.375: +2.5, 3.625: +3.4,
        4.0: +5.7, 4.5: +2.8, 5.0: -4.5, 5.5: -5.1, 6.25: -0.9, 7.0: -1.5}
RING_C = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])


def ring_sums(gname):
    """Complex m = 2 coefficient of <v_R> and of the density, per ring, in the nominal bar frame.

    The velocity coefficients come from the LEAST-SQUARES fit of the field to
    [1, cos 2phi, sin 2phi, cos 4phi, sin 4phi] (ring_coeffs -> fs.fit), not from a discrete sum
    sum_i v_i exp(2 i phi_i). The two differ: the discrete sum is the m = 2 of the product n*v, so it
    picks up (m=0 of v) x (m=2 of the density) and therefore mixes the bar's own density quadrupole -
    and the selection function - into what is meant to be a property of the velocity field. For
    <v_R> the m = 0 term is ~0 by symmetry and the two agree to ~1 deg, but the fit is the right
    quantity and is what the production estimator uses. The density coefficient is a discrete sum,
    which for the density is exactly right."""
    g = galaxy(gname)
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    a2, b2 = ring_coeffs(W, "vR")                  # ring_coeffs rotates to the bar frame itself
    zv_fs = (a2 + 1j * b2) * np.exp(-1j * np.pi / 2)     # v_R extrema sit 45 deg from the node
    w = ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG))
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.arctan2(w[:, 1], w[:, 0])
    thin = np.abs(w[:, 2]) < 1.5
    zv = np.full(len(RING_C), np.nan, complex)
    zd = np.full(len(RING_C), np.nan, complex)
    nn = np.zeros(len(RING_C), int)
    for i in range(len(RING_C)):
        s = thin & (R >= fs.RING_EDGES[i]) & (R < fs.RING_EDGES[i + 1])
        nn[i] = s.sum()
        if nn[i] < 400 or not np.isfinite(zv_fs[i]):
            continue
        zv[i] = zv_fs[i]
        zd[i] = np.sum(np.exp(2j * phi[s])) / nn[i]
    return zv, zd, nn


def combine(zv, zd, nn, kin):
    """Offset over the kinematically-defined bar region, and the region itself."""
    ok = np.isfinite(zv) & np.isfinite(zd)
    amp = np.where(ok, np.abs(zv), np.nan)
    m = ok & (amp > kin * np.nanmax(amp))
    off = np.degrees(np.angle(np.nansum(zv[m] * nn[m]) * np.conj(np.nansum(zd[m] * nn[m])))) / 2
    return float(off), float(RING_C[m].min()), float(RING_C[m].max())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kin", type=float, default=0.7,
                    help="keep rings whose kinematic m=2 amplitude exceeds this fraction of its peak")
    args = ap.parse_args()
    set_style()
    tb = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
    rows = []
    for t in TIMES:
        gn = nbody_name(t)
        zv, zd, nn = ring_sums(gn)
        ok = np.isfinite(zv)
        fixed = np.degrees(np.angle(
            np.nansum(zv[ok & (RING_C > 1) & (RING_C < 4)] * nn[ok & (RING_C > 1) & (RING_C < 4)]) *
            np.conj(np.nansum(zd[ok & (RING_C > 1) & (RING_C < 4)] * nn[ok & (RING_C > 1) & (RING_C < 4)])))) / 2
        kino, rlo, rhi = combine(zv, zd, nn, args.kin)
        rows.append(dict(t=t, galaxy=gn,
                         Omega_b=float(tb[np.isclose(tb.snap_time, t)]["pattern_speed_fd"].iloc[0]),
                         A2=float(tb[np.isclose(tb.snap_time, t)]["bar_strength"].iloc[0]),
                         off_fixed=float(fixed), off_kin=kino, R_lo=rlo, R_hi=rhi, bias=BIAS[t]))
        print(f"  t = {t:<6g} offset: fixed 1-4 kpc {fixed:+6.2f}   "
              f"kinematic ({rlo:.2f}-{rhi:.2f} kpc) {kino:+6.2f}   bias {BIAS[t]:+5.1f}%", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS_DIR / "nbody_equilibrium_offset.csv", index=False)
    print(f"\nwrote {rel(RESULTS_DIR / 'nbody_equilibrium_offset.csv')}")
    for c in ("off_fixed", "off_kin"):
        print(f"  {c:9s}: median |offset| {np.median(np.abs(D[c])):5.2f} deg, "
              f"max {np.abs(D[c]).max():5.2f} (t = {D.t[np.abs(D[c]).idxmax()]:g})")
    r1 = stats.pearsonr(np.abs(D.off_kin), np.abs(D.bias))
    rs = stats.spearmanr(np.abs(D.off_kin), np.abs(D.bias))
    print(f"\n  |kinematic offset| vs |Omega_b bias|:  r = {r1.statistic:+.2f} (p = {r1.pvalue:.3f}), "
          f"Spearman {rs.statistic:+.2f} (p = {rs.pvalue:.3f})")

    fig, (a, b) = plt.subplots(1, 2, figsize=(12.4, 4.8))
    a.axhline(0, color=INK_MUTED, lw=1, ls="--")
    a.plot(D.Omega_b, D.off_fixed, "o-", ms=5, lw=1.2, color=BLUE_RAMP[3], label="fixed 1–4 kpc")
    a.plot(D.Omega_b, D.off_kin, "D-", ms=6, lw=1.6, color=ACCENT,
           label=f"kinematic, >{args.kin:.0%} of peak")
    for _, r in D.iterrows():
        a.annotate(f"{r.t:g}", (r.Omega_b, r.off_kin), textcoords="offset points", xytext=(0, 8),
                   ha="center", fontsize=8, color=INK_2)
    a.set_xlabel(r"true $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    a.set_ylabel(r"$\langle v_R\rangle$ node $-$ density bar axis [deg]")
    a.set_title("equilibrium test (no selection, no distance error)", color=INK, fontsize=10.5)
    a.legend(frameon=False, fontsize=9)

    b.axhline(0, color=INK_MUTED, lw=1, ls="--")
    b.axvline(0, color=INK_MUTED, lw=1, ls="--")
    b.plot(D.off_kin, D.bias, "o", ms=7, color=ACCENT)
    for _, r in D.iterrows():
        b.annotate(f"{r.t:g}", (r.off_kin, r.bias), textcoords="offset points", xytext=(0, 8),
                   ha="center", fontsize=8, color=INK_2)
    b.set_xlabel(r"$\langle v_R\rangle$ node $-$ density bar axis [deg]")
    b.set_ylabel(r"$\Omega_b$ bias [%]")
    b.set_title(f"does disequilibrium predict the bias?   "
                f"r(|x|,|y|) = {r1.statistic:+.2f}, p = {r1.pvalue:.2f}", color=INK, fontsize=10.5)
    fig.suptitle("Is the kinematic bar aligned with the density bar? A steady bar requires zero",
                 color=INK, y=0.99)
    fig.tight_layout()
    savefig(fig, "test2_equilibrium_offset")


if __name__ == "__main__":
    main()

"""Bar pattern speed and strength of the isolated N-body run against time (mock/Nbody/t_t0_*).

Reuses the mock-data machinery: `nbody.read_nemo` / `nbody.components` read each NEMO snapshot,
`prepare_au18_mwframe.shrinking_sphere` and `measure_bar` give the centre and the m = 2 bar, and
Omega_b comes from the Dehnen, Semczuk & Schoenrich (2023) Fourier method exactly as in
`measure_au18_pattern_speed.py` (disc mirrored to rotate anticlockwise, masses as weights).

Each snapshot is put in the same MW-like frame as the mocks: centred on the stars, bulk velocity of
the stars within R_DISC removed, and the disc spun so that L_z < 0 (clockwise seen from +z). The bar
is *not* rotated to -27 deg here; the measured phase is kept as it is, so that it can be compared
between snapshots.

Omega(R), and so the corotation radius, comes from the snapshot's own accelerations (azimuthal
average of a_R in a thin midplane slab), not from a spherical enclosed-mass estimate: at t = 7 this
gives R_CR = 8.61 kpc, matching the 8.61 kpc of the AGAMA ground-truth potential in
results/nbody_ground_truth.json, while the spherical estimate is 2 kpc too small.

Why a single-snapshot method. The snapshots are 0.125 model time units apart (0.122 Gyr), over which
a 26 km/s/kpc bar turns by 3.2 rad, i.e. more than the pi period of the m = 2 phase. Differencing the
phase between snapshots is therefore aliased and cannot give Omega_b. The Fourier method instead uses
the velocities and needs only one snapshot. `phase_check` turns the aliasing around and uses it as a
weak consistency test (see its docstring).

Usage:
    python measure_nbody_pattern_speed_evolution.py                 # 2 to 8 Gyr
    python measure_nbody_pattern_speed_evolution.py --tmin 0 --tmax 9.5 --step 0.25
    python measure_nbody_pattern_speed_evolution.py --model-time    # --tmin/--tmax in model units
Times are in Gyr unless --model-time is given; snapshot file names are in model time units
(1 unit = 0.9778 Gyr).

Outputs (both repo-relative):
    results/nbody_pattern_speed_evolution.csv
    figures/nbody_pattern_speed_evolution.png
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import auriga as au
import nbody
from common import RESULTS_DIR, EXTERNAL_DIR, ROOT, rel, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from prepare_au18_mwframe import shrinking_sphere, measure_bar, R_DISC

PATTERN_SPEED_DIR = os.environ.get("PATTERN_SPEED_DIR") or str(EXTERNAL_DIR / "patternSpeed-main")
if not hasattr(np, "asfarray"):        # removed in NumPy 2.0; patternSpeed.py still calls it
    np.asfarray = lambda a: np.asarray(a, dtype=np.float64)
sys.path.insert(0, PATTERN_SPEED_DIR)
import patternSpeed  # noqa: E402

OUT_CSV = RESULTS_DIR / "nbody_pattern_speed_evolution.csv"
CENTRE_RADII = (20, 10, 5, 3, 2, 1, 0.5)      # shrinking-sphere radii, as in prepare_nbody_mwframe.py
TRUTH_TIME = nbody.DEFAULT_TIME                # the snapshot the mock tests use (t = 7)
RAD_PER_GYR = 1.0 / nbody.TIME_UNIT_GYR        # 1 km/s/kpc in rad/Gyr


def mw_like_frame(d, comp):
    """Centre, de-boost and spin one snapshot into the mocks' MW-like frame (L_z < 0).
    Returns the stars' positions, velocities and masses, and the accelerations of all particles
    (the full potential gradient, used for Omega(R))."""
    st = comp < 2
    pos = d["Position"].astype(np.float64)
    vel = d["Velocity"].astype(np.float64)
    m = d["Mass"].astype(np.float64) * nbody.MASS_UNIT_MSUN
    pos -= shrinking_sphere(pos[st], m[st], radii=CENTRE_RADII)
    s = st & (np.linalg.norm(pos, axis=1) < R_DISC)
    vel -= np.average(vel[s], axis=0, weights=m[s])
    L = np.sum(m[s, None] * np.cross(pos[s], vel[s]), axis=0)
    R = au.rotation_to(L, np.array([0.0, 0.0, -1.0]))
    return pos[st] @ R.T, vel[st] @ R.T, m[st], pos @ R.T, d["Acceleration"].astype(np.float64) @ R.T


CR_EDGES = np.arange(0.2, 25.0 + 1e-9, 0.2)      # radial bins for Omega(R)
CR_SLAB_KPC = 0.25                               # |z| slab for the midplane average


def circular_frequency(x, a):
    """Omega(R) [km/s/kpc] from the snapshot's own accelerations: the azimuthal average of the
    radial acceleration in a thin midplane slab, Omega = sqrt(-<a_R> / R). This uses the real
    (flattened) potential of the disc, unlike a spherical enclosed-mass estimate, and includes the
    gravitational softening exactly as the simulation does. Returns (R_centres, Omega)."""
    R = np.hypot(x[:, 0], x[:, 1])
    aR = (x[:, 0] * a[:, 0] + x[:, 1] * a[:, 1]) / R
    Rc = 0.5 * (CR_EDGES[1:] + CR_EDGES[:-1])
    i = np.digitize(R, CR_EDGES) - 1
    ok = (np.abs(x[:, 2]) < CR_SLAB_KPC) & (i >= 0) & (i < len(Rc))
    n = np.bincount(i[ok], minlength=len(Rc))
    s = np.bincount(i[ok], weights=aR[ok], minlength=len(Rc))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_aR = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    return Rc, np.sqrt(np.maximum(-mean_aR, 0.0) / Rc)


def corotation_radius(Rc, Omega, Omega_p):
    """Innermost radius where Omega(R) falls through Omega_p."""
    f = Omega - Omega_p
    i = np.flatnonzero(np.diff(np.sign(f)) != 0)
    if not len(i):
        return np.nan
    j = i[0]
    return float(Rc[j] - f[j] * (Rc[j + 1] - Rc[j]) / (f[j + 1] - f[j]))


def analyse(t_model):
    """Bar strength, length, phase and pattern speed of the snapshot at model time t_model."""
    d = nbody.read_nemo(nbody.snapshot(t_model), tags=("Mass", "Position", "Velocity", "Acceleration"))
    comp = nbody.components(d["Mass"])
    x, v, m, x_all, a_all = mw_like_frame(d, comp)
    del d

    bar = measure_bar(x[:, 0], x[:, 1], x[:, 2], m)
    # MW-like frame rotates clockwise: mirror y -> -y to give the method an anticlockwise disc
    pat, _ = patternSpeed.FourierMethod.patternSpeed(x[:, 0], -x[:, 1], v[:, 0], -v[:, 1], mu=m)
    Om_p = float(pat["Om"])
    Rc, Omega = circular_frequency(x_all, a_all)
    return dict(
        t_model=float(t_model), t_Gyr=float(t_model) * nbody.TIME_UNIT_GYR,
        Omega_p=Om_p, Omega_p_err=float(pat["Om_e"]),
        bar_R0_kpc=float(pat["R0"]), bar_R1_kpc=float(pat["R1"]), bar_Rm_kpc=float(pat["Rm"]),
        A2_max=float(bar["A2_max"]), R_peak_kpc=float(bar["R_peak"]), R_bar_A2_kpc=float(bar["R_bar"]),
        phase_bar_rad=float(bar["phase_bar"]),
        R_CR_kpc=corotation_radius(Rc, Omega, Om_p),
        v_c_8kpc=float(np.interp(8.0, Rc, Omega * Rc)),
        n_stars=int(len(m)), M_stars=float(m.sum()),
    )


def phase_check(df):
    """Per-step consistency of Omega_b with the measured m = 2 phase: the phase the bar actually turned
    between consecutive snapshots minus Omega_b dt, wrapped to (-90, 90] deg.

    The frame rotates clockwise, so the phase decreases: dphi = -Omega_b dt (trapezoidal in Omega).
    The m = 2 phase is defined modulo pi, hence the wrap, and that sets what this can test: one step is
    0.122 Gyr, over which the ambiguity limit is dOmega = 12.6 km/s/kpc, so the check catches a gross
    error in Omega_b but not a small one. It must be done per step: integrating across 2-8 Gyr
    accumulates 60 pi of phase, where even a 1% error in Omega_b wraps.

    Returns the per-step residuals (length N-1) at the step mid-times. Their median tests Omega_b for
    bias; their scatter is dominated by the bar's own phase wobble (the m = 2 phase of a live N-body
    bar is not a rigid-body angle), not by the method's quoted error."""
    t, Om, ph = df.t_Gyr.values, df.Omega_p.values, df.phase_bar_rad.values
    dt = np.diff(t)
    d_pred = -0.5 * (Om[1:] + Om[:-1]) * dt * RAD_PER_GYR
    res = np.mod(np.diff(ph) - d_pred + np.pi / 2, np.pi) - np.pi / 2
    return 0.5 * (t[1:] + t[:-1]), np.degrees(res), dt


def phase_check_summary(df):
    """Median per-step phase residual and the systematic Omega_b error it corresponds to."""
    _, res, dt = phase_check(df)
    med = float(np.median(res))
    sem = float(1.2533 * np.std(res, ddof=1) / np.sqrt(len(res)))     # standard error of a median
    scale = float(np.mean(dt) * RAD_PER_GYR * 180.0 / np.pi)          # deg per (km/s/kpc) of offset
    return dict(median_deg=med, sem_deg=sem, rms_deg=float(np.sqrt(np.mean(res**2))),
                dOmega=med / scale, dOmega_err=sem / scale)


def plot(df):
    set_style()
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7.8, 9.0), sharex=True,
                                     gridspec_kw=dict(height_ratios=[2.2, 1.25, 1.0], hspace=0.32))

    a1.fill_between(df.t_Gyr, df.Omega_p - df.Omega_p_err, df.Omega_p + df.Omega_p_err,
                    color=BLUE_RAMP[0], lw=0, label=r"$\pm 1\sigma$ (method)")
    a1.plot(df.t_Gyr, df.Omega_p, "-o", color=BLUE_RAMP[5], lw=1.8, ms=3.4, label=r"$\Omega_b$ (Fourier method)")
    tt = TRUTH_TIME * nbody.TIME_UNIT_GYR
    if df.t_Gyr.min() - 1e-6 <= tt <= df.t_Gyr.max() + 1e-6:
        row = df.iloc[int(np.argmin(np.abs(df.t_Gyr - tt)))]
        a1.plot(row.t_Gyr, row.Omega_p, "*", color=ACCENT, ms=16, mec="white", mew=0.8, zorder=6,
                label=rf"$t_0$ = {TRUTH_TIME:g} (mock tests): {row.Omega_p:.2f}")
    a1.set_ylabel(r"$\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    a1.legend(loc="upper right", fontsize=8.5)
    a1.set_title(f"Isolated N-body bar, {rel(nbody.SNAPSHOT_DIR)}: pattern speed and bar evolution",
                 color=INK_2, fontsize=10)

    a2.plot(df.t_Gyr, df.A2_max, "-", color=BLUE_RAMP[5], lw=1.8, label=r"$A_2^{\rm max}$")
    a2.set_ylabel(r"$A_2^{\rm max}$", color=BLUE_RAMP[5])
    b2 = a2.twinx()
    b2.plot(df.t_Gyr, df.R_bar_A2_kpc, "-", color=ACCENT, lw=1.5, label=r"$R_{\rm bar}$ ($A_2$ half max)")
    b2.plot(df.t_Gyr, df.R_CR_kpc, "--", color=INK_MUTED, lw=1.3, label=r"$R_{\rm CR}$")
    b2.set_ylabel("radius [kpc]")
    b2.set_ylim(0, None)
    a2.set_zorder(b2.get_zorder() + 1)
    a2.patch.set_visible(False)
    a2.set_ylim(0.30, None)        # headroom below the A2 curve for the legend
    b2.set_ylim(0, None)
    h = a2.get_legend_handles_labels()[0] + b2.get_legend_handles_labels()[0]
    a2.legend(h, [x.get_label() for x in h], loc="lower left", fontsize=8.5, ncol=3,
              columnspacing=1.4, handlelength=1.8)

    tm, res, _ = phase_check(df)
    s = phase_check_summary(df)
    a3.axhline(0, color=INK_MUTED, lw=0.7)
    a3.axhspan(s["median_deg"] - s["sem_deg"], s["median_deg"] + s["sem_deg"], color=ACCENT, alpha=0.3, lw=0)
    a3.axhline(s["median_deg"], color=ACCENT, lw=1.2)
    a3.plot(tm, res, "o", color=BLUE_RAMP[4], ms=3.2)
    a3.set_ylim(-90, 90)
    a3.set_ylabel("phase residual\nper step [deg]")
    a3.set_xlabel(r"$t$ [Gyr]")
    a3.set_title(rf"check: measured $\Delta\psi_2$ minus $-\Omega_b\,\Delta t$ per step (mod $180^\circ$); "
                 rf"median {s['median_deg']:+.1f} $\pm$ {s['sem_deg']:.1f}$^\circ$ "
                 rf"= $\Delta\Omega_b$ {s['dOmega']:+.2f} $\pm$ {s['dOmega_err']:.2f}",
                 color=INK_2, fontsize=8.5)
    for ax in (a1, a2, a3):
        ax.set_xlim(df.t_Gyr.min(), df.t_Gyr.max())
    return savefig(fig, "nbody_pattern_speed_evolution")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tmin", type=float, default=2.0, help="first time (default 2)")
    ap.add_argument("--tmax", type=float, default=8.0, help="last time (default 8)")
    ap.add_argument("--step", type=float, default=None, help="time step (default: every snapshot)")
    ap.add_argument("--model-time", action="store_true", help="--tmin/--tmax/--step are in model units, not Gyr")
    ap.add_argument("--plot-only", action="store_true", help="re-plot the saved CSV without reading snapshots")
    args = ap.parse_args()

    if args.plot_only:
        print(f"plotting {rel(OUT_CSV)}")
        plot(pd.read_csv(OUT_CSV).sort_values("t_model"))
        return

    unit = 1.0 if args.model_time else nbody.TIME_UNIT_GYR
    available = nbody.available_times()
    want = available[(available >= args.tmin / unit - 1e-9) & (available <= args.tmax / unit + 1e-9)]
    if args.step is not None:
        want = want[np.isclose((want - want[0])[:, None] % (args.step / unit), 0, atol=1e-6).any(1)]
    if not len(want):
        raise SystemExit(f"no snapshots in range; available model times {available.min()}-{available.max()}")
    print(f"{len(want)} snapshots, model t = {want[0]:g} to {want[-1]:g} "
          f"({want[0] * nbody.TIME_UNIT_GYR:.2f} to {want[-1] * nbody.TIME_UNIT_GYR:.2f} Gyr)", flush=True)

    rows, t0 = [], time.time()
    done = {}
    if OUT_CSV.exists():      # resumable: keep rows already measured
        done = {round(r.t_model, 6): r._asdict() for r in pd.read_csv(OUT_CSV).itertuples(index=False)}
    for i, t in enumerate(want):
        if round(float(t), 6) in done:
            rows.append(done[round(float(t), 6)])
            print(f"[{i + 1}/{len(want)}] t = {t:g}: cached", flush=True)
            continue
        rows.append(analyse(t))
        r = rows[-1]
        print(f"[{i + 1}/{len(want)}] t = {t:g} ({r['t_Gyr']:.2f} Gyr): "
              f"Omega_b = {r['Omega_p']:6.2f} +/- {r['Omega_p_err']:.2f} km/s/kpc, "
              f"A2_max = {r['A2_max']:.3f}, R_bar = {r['R_bar_A2_kpc']:.2f} kpc, "
              f"R_CR = {r['R_CR_kpc']:.2f} kpc  ({time.time() - t0:.0f} s elapsed)", flush=True)
        pd.DataFrame(rows).to_csv(OUT_CSV, index=False)      # flush after every snapshot

    df = pd.DataFrame(rows).sort_values("t_model").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {rel(OUT_CSV)}")
    plot(df)
    s = phase_check_summary(df)
    print(f"\nOmega_b over {df.t_Gyr.min():.2f}-{df.t_Gyr.max():.2f} Gyr: "
          f"{df.Omega_p.max():.2f} down to {df.Omega_p.min():.2f} km/s/kpc "
          f"(median quoted error {df.Omega_p_err.median():.2f})")
    print(f"phase check: median per-step residual {s['median_deg']:+.2f} +/- {s['sem_deg']:.2f} deg "
          f"(rms {s['rms_deg']:.1f}), i.e. a systematic Omega_b offset of "
          f"{s['dOmega']:+.2f} +/- {s['dOmega_err']:.2f} km/s/kpc")


if __name__ == "__main__":
    main()

"""Do the library's kinematic moments change with the potential, and is that change visible above the data errors?

chi2 against the data mixes two things: how much the PREDICTED moments move as the potential changes, and how
large the data errors are. This separates them using the saved library coefficients of every design point
(coef_lib in the sfratio_NNN.npz of the re-score).

For each design point i, the library coefficients are rotated from the trial bar frame back to the reference frame
of the data fit (the exact c' = T c of test1_nbody_qmc_focused_fourier_stars.rotation), and compared with the
library of the best point:
    D_i = chi2(library_i vs library_best | data covariance)   (scheme B terms, Hartlap, same as the likelihood)
D_i is how distinguishable potential i's prediction is from the best one, in units of the data errors; it does not
involve the data values. Its profile in |Omega_b| (smallest D_i per bin, tempered as D / 2 sqrt(2 N_data)) is the
best the data could ever do. Then the two ingredients separately, for <v_R> m = 2 in rings with 1 < R < 4 kpc,
|z| < 1: the spread of the library amplitude across the design (km/s, rotation invariant) and the bootstrap error
of the data coefficients (km/s).

Usage: python diagnose_library_sensitivity.py
Figure: figures/test2_library_sensitivity_osarg_vs_mock.png
"""
import json

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP
import test1_nbody_qmc_focused_fourier_stars as fs
from test2b_qmc_distance_error_run import load_design

USB = "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs"
LOCAL = "/Users/hanyuan/mira_dr4_runs"
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
STEM = "nbody_mock_selection_sigmaDM0.271_N50k"
CASES = [  # label, results run (design, evals), trajectory/score dir, variant, mock, mask, colour, marker
    ("OSARG data, Hunter+24", "osarg_zhang2024_mw_hunter24", f"{LOCAL}/osarg_zhang2024_mw_hunter24", "sfratio",
     "osarg_zhang2024_observed_R8.hdf5", None, INK_2, "D"),
    ("t = 7 mock, Hunter+24 (same design)", "test2_t7_sf_12p5pct_mwH24_osargdesign",
     f"{USB}/test2_t7_sf_12p5pct_mwH24_osargdesign", "sfratio", f"{STEM}.hdf5", None, ACCENT, "s"),
    ("t = 7 mock, own family", "test2_t7_sf_12p5pct", f"{USB}/test2_t7_sf_12p5pct", "sfratio", f"{STEM}.hdf5",
     None, BLUE_RAMP[5], "o"),
    ("t = 7 mock, own family, OSARG coverage + N", "test2_t7_sf_12p5pct_osarg_cover", f"{USB}/test2_t7_sf_12p5pct",
     "sfratio_osarg_cover", f"{STEM}.hdf5", f"results/rv_selection_masks/{STEM}__osarg_cover.npy", BLUE_RAMP[3], "^"),
]
RINGS = [r for r in range(fs.N_RING_DEFAULT) if 1.0 <= 0.5 * (fs.RING_EDGES[r] + fs.RING_EDGES[r + 1]) <= 4.0]


def data_fit(mock, mask):
    with h5py.File(MOCK_DIR / mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    if mask is not None:
        w = w[np.load(ROOT / mask)]
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)),
                                                      size=fs.N_BOOT).astype(float)
    return fs.data_fit(w, BAR_ANGLE_DEG, boot), len(w)


def main():
    set_style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.5, 5.6), gridspec_kw=dict(wspace=0.25))
    edges = np.arange(15.0, 45.01, 1.25)
    mid = 0.5 * (edges[1:] + edges[:-1])
    rows = []
    for lab, run, sdir, var, mock, mask, c, mk in CASES:
        meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
        design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
        design = design.set_index("idx")
        (cd0, cov0, n_stars, ok), n = data_fit(mock, mask)
        ev = pd.read_csv(RESULTS_DIR / run / "evals.csv").set_index("idx")
        ev = ev[ev.kind == "qmc"]
        lib = {}
        for i in ev.index:
            cl = np.load(f"{sdir}/{var}_{i:03d}.npz")["coef_lib"].reshape(fs.NB, -1)
            T = fs.rotation(BAR_ANGLE_DEG - float(design.loc[i, "phi_b_deg"]))      # trial frame -> data frame
            lib[i] = (cl @ T.T).reshape(fs.NB, 6, fs.NK)
        ibest = int(ev[f"logLs_{B}_{var}"].idxmax())
        D, nd = {}, None
        for i, cl in lib.items():
            valid = ok & np.all(np.isfinite(cl), axis=(1, 2)) & np.all(np.isfinite(lib[ibest]), axis=(1, 2))
            D[i], nd_i = fs.scheme_chi2(lib[ibest], cl, cov0, valid, fs.SCHEMES[B])
            nd = nd_i if i == ibest else nd
        Ds = pd.Series(D) / (2 * np.sqrt(2 * nd))
        om = design.loc[Ds.index, "Omega_abs"]
        prof = np.array([Ds[(om >= a) & (om < b)].min() if ((om >= a) & (om < b)).any() else np.nan
                         for a, b in zip(edges[:-1], edges[1:])])
        a1.plot(mid, -prof, "-", marker=mk, ms=4.5, color=c, lw=1.5, label=lab)
        # the two ingredients for <v_R> m = 2, 1 < R < 4, |z| < 1
        blk = [r * fs.NL for r in RINGS]
        amp = np.array([[np.hypot(*lib[i][b, 0, 1:3]) for b in blk] for i in lib])          # (points, rings)
        spread = np.nanpercentile(amp, 84, axis=0) - np.nanpercentile(amp, 16, axis=0)
        err = np.array([np.sqrt(0.5 * (cov0[b, 1, 1] + cov0[b, 2, 2])) for b in blk])
        dat = np.array([np.hypot(*cd0[b, 0, 1:3]) for b in blk])
        lib_med = np.nanmedian(amp, axis=0)
        rows.append(dict(case=lab, stars=n, N_data=nd, best_Omega=float(design.loc[ibest, "Omega_abs"]),
                         lib_spread_kms=float(np.median(spread)), data_err_kms=float(np.median(err)),
                         spread_over_err=float(np.median(spread / err)), data_amp_kms=float(np.median(dat)),
                         lib_amp_kms=float(np.median(lib_med)),
                         D_at_dOmega5=float(np.nanmin(Ds[np.abs(om - om[ibest]) > 5])) if (np.abs(om - om[ibest]) > 5).any() else np.nan))
        a2.plot([0.5 * (fs.RING_EDGES[r] + fs.RING_EDGES[r + 1]) for r in RINGS], spread / err, "-", marker=mk,
                ms=5, color=c, lw=1.5, label=lab)
        print(f"{lab}: {n:,} stars, N_data {nd}, best point idx {ibest} (Omega {design.loc[ibest, 'Omega_abs']:.1f})")
    t = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print(t.round(2).to_string(index=False))
    t.to_csv(RESULTS_DIR / "library_sensitivity_osarg_vs_mock.csv", index=False)
    for lev in (0.5, 2.0, 4.5):
        a1.axhline(-lev, color=INK_2, lw=0.6, ls=":")
    a1.set(xlim=(15, 45), ylim=(-40, 1), xlabel=r"$|\Omega_b|$  [km s$^{-1}$ kpc$^{-1}$]",
           ylabel=r"$-\min D$ in the bin  (tempered, $D/2\sqrt{2N}$)")
    a1.set_title("how distinguishable the predicted moments are from the best model's\n"
                 "(library vs library, in units of the data errors; no data values)", color=INK, fontsize=10.5)
    a1.legend(loc="lower left", frameon=False, fontsize=8.5)
    a2.axhline(1, color=INK_2, lw=0.8, ls=":")
    a2.set(xlabel="R [kpc]  (|z| < 1)", ylabel=r"16–84% spread of library $\langle v_R\rangle$ m=2 amplitude / data error")
    a2.set_title(r"$\langle v_R\rangle$ quadrupole: how much the prediction moves across the design," + "\nin units "
                 "of the data's bootstrap error", color=INK, fontsize=10.5)
    a2.legend(loc="upper right", frameon=False, fontsize=8.5)
    savefig(fig, "test2_library_sensitivity_osarg_vs_mock")


if __name__ == "__main__":
    main()

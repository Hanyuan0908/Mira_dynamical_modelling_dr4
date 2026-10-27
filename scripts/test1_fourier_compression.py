"""Test 1b: compress the polar moment maps into Fourier series in phi, per ring and |z| layer.

For each ring (R 0.5 kpc), |z| layer and field (<v_R>, <v_phi>, <v_z sign z>, sigma_R, sigma_phi, sigma_z), the 15-deg
sector values (folded phi in [0, pi), sectors with N >= N_MIN) are replaced by Fourier coefficients up to m_max (even,
full-plane m: 1, cos 2phi, sin 2phi, ..., cos m_max phi, sin m_max phi), fitted by generalised least squares with the
data's bootstrap errors D. The same linear operator is applied to the data and the library, so
    chi2(m_max) = dm^T D^-1 F (F^T D^-1 F)^-1 F^T D^-1 dm,   dm = m_data - m_lib,
the part of the current sector chi^2 that lies in modes m <= m_max; N_data per block = number of coefficients. With
12 sectors and m_max = 12, F is square and chi2 equals the current sector chi^2 exactly (cos 12phi vanishes at the
sector centres and is dropped). Blocks with fewer usable sectors than coefficients keep their sectors.

No new integrations: data sector moments and bootstrap errors (test1_nbody_moments.boot_errors, same seed), saved
library sector moments of the global search (000-064), the pilot offsets and accuracy-noise twin, and the diag ridge
points. Also the non-axisymmetric part alone, chi2F2to<m> = chi2F<m> - chi2F0 (modes 2..m with the error-weighted ring
average projected out; N_data = N(m) - N(0)). Outputs: results/test1_nbody/fourier_compression.json, columns
chi2F<m>_<set>, chi2F2to<m>_<set> and tempered logLs_F<m>_<set>, logLs_F2to<m>_<set> in qmc_evals.csv,
figures/test1_nbody_fourier_compression.png.
"""
import json

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, N_MIN, SEED
from test1_nbody_qmc import OUT_CSV, ARR_DIR, MOCK, TRUTH

NR, NP, NZ = len(ssl.CELL_R) - 1, len(ssl.CELL_PHI) - 1, len(ssl.CELL_Z) - 1
PHI_C = 0.5 * (ssl.CELL_PHI[1:] + ssl.CELL_PHI[:-1])
M_LIST = (0, 2, 4, 6, 8, 10, 12)
SETS = {"means": [0, 1, 2], "meansdisp": [0, 1, 2, 3, 4, 5], "Rphi_meansdisp": [0, 1, 3, 4]}
FIELDS = ["<v_R>", "<v_phi>", "<v_z>", "sigma_R", "sigma_phi", "sigma_z"]
# per-field mode schemes at m_max = 4: {field: keep the ring average (m = 0)?}; m = 2, 4 always kept (cos and sin)
SCHEMES = {"F4_inplane": {0: False, 1: True, 3: True, 4: True},          # <v_R> m 2,4; <v_phi>, sigma_R, sigma_phi m 0,2,4
           "F4_inplane_nos0": {0: False, 1: True, 3: False, 4: False}}   # same, dispersions without their ring average
MODELS = [("S 1.25", "pilot_09.npz"), ("S 0.8", "pilot_08.npz"), ("f 1.4", "pilot_07.npz"), ("f 0.7", "pilot_06.npz"),
          ("ridge S 0.82, f 1.50", "diag_SF_S0.82_f1.5.npz"), ("ridge S 1.25, f 0.70", "diag_SF_S1.25_f0.7.npz"),
          ("truth, accuracy 3e-6 (integrator noise only)", "pilot_01.npz")]


def basis(m_max, phi):
    cols = [np.ones_like(phi)]
    for m in range(2, m_max + 1, 2):
        cols += [np.cos(m * phi), np.sin(m * phi)]
    F = np.column_stack(cols)
    return F[:, np.linalg.norm(F, axis=0) > 1e-8]


def chi2_fourier(obs, lib_m1, lib_s, m_max, fields):
    """chi^2 and N_data of library vs data with each (ring, layer, field) block projected on modes <= m_max."""
    N, m1, s, em, es = obs
    good = (N >= N_MIN) & np.all(np.isfinite(lib_m1), 1) & np.all(em > 0, 1) & np.all(es > 0, 1)
    D = np.column_stack([m1 - lib_m1, s - lib_s]).reshape(NR, NP, NZ, 6)
    E = np.column_stack([em, es]).reshape(NR, NP, NZ, 6)
    good = good.reshape(NR, NP, NZ)
    chi, nd = 0.0, 0
    for r in range(NR):
        for z in range(NZ):
            av = good[r, :, z]
            n = int(av.sum())
            if n == 0:
                continue
            F = basis(m_max, PHI_C[av])
            if F.shape[1] > n or np.linalg.matrix_rank(F) < F.shape[1]:
                F = np.eye(n)
            for f in fields:
                w = 1.0 / E[r, av, z, f]**2
                b = F.T @ (w * D[r, av, z, f])
                chi += float(b @ np.linalg.solve(F.T @ (w[:, None] * F), b))
                nd += F.shape[1]
    return chi, nd


def scheme_chi2(per_field, scheme):
    """chi^2 (or N_data) of a scheme from per-field values {(field, m_max): x} at m_max = 0 and 4."""
    return sum(per_field[(f, 4)] - (0 if keep0 else per_field[(f, 0)]) for f, keep0 in scheme.items())


def lib(fn):
    a = np.load((OUT_CSV.parent if fn.startswith(("pilot", "diag")) else ARR_DIR) / fn)
    return a["lib_m1"], a["lib_s"]


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    df = pd.read_csv(OUT_CSV).sort_values("idx").reset_index(drop=True)
    obs_t = boot_errors(ssl.rotate_z(stars_mw, np.radians(TRUTH["phi_b_deg"])), np.random.default_rng(SEED))
    out = {"models": {}, "data_minus_truth": {}, "n_data_truth": {}}

    # ---- single offsets and ridge points (Omega_b, phi_b at the truth): realised and noise-free chi^2 vs m_max ----
    truth = lib("pilot_00.npz")
    zero_err = (obs_t[0], truth[0], truth[1], obs_t[3], obs_t[4])      # "data" = truth library -> noise-free
    for sname, fields in SETS.items():
        for m in M_LIST:
            c_t, nd = chi2_fourier(obs_t, *truth, m, fields)
            out["n_data_truth"].setdefault(sname, {})[m] = nd
            out["data_minus_truth"].setdefault(sname, {})[m] = c_t
            for lab, fn in MODELS:
                L = lib(fn)
                c, _ = chi2_fourier(obs_t, *L, m, fields)
                c0, _ = chi2_fourier(zero_err, *L, m, fields)
                d = out["models"].setdefault(lab, {}).setdefault(sname, {})
                d[m] = dict(delta_chi2=c - c_t, dlogL=-(c - c_t) / (2 * np.sqrt(2 * nd)), noise_free=c0)

    # ---- per field, m_max = 0 and 4: N_data at the truth, and Delta chi^2 of the models per field ----
    nd_pf = {(f, m): chi2_fourier(obs_t, *truth, m, [f])[1] for f in range(6) for m in (0, 4)}
    ct_pf = {(f, m): chi2_fourier(obs_t, *truth, m, [f])[0] for f in range(6) for m in (0, 4)}
    out["per_field"] = {"N_data": {f"{FIELDS[f]} m<={m}": nd_pf[(f, m)] for f in range(6) for m in (0, 4)},
                        "truth_chi2": {f"{FIELDS[f]} m<={m}": ct_pf[(f, m)] for f in range(6) for m in (0, 4)},
                        "models": {}}
    for lab, fn in MODELS[:6]:
        L = lib(fn)
        c_pf = {(f, m): chi2_fourier(obs_t, *L, m, [f])[0] for f in range(6) for m in (0, 4)}
        out["per_field"]["models"][lab] = {f"{FIELDS[f]} m<={m}": c_pf[(f, m)] - ct_pf[(f, m)]
                                           for f in range(6) for m in (0, 4)}
        for sc, spec in SCHEMES.items():
            nd, dchi = scheme_chi2(nd_pf, spec), scheme_chi2(c_pf, spec) - scheme_chi2(ct_pf, spec)
            out["models"][lab].setdefault(sc, {})[4] = dict(delta_chi2=dchi, dlogL=-dchi / (2 * np.sqrt(2 * nd)))

    # ---- global search: all 65 points, every m_max; same bootstrap errors per phi as test1_nbody_moments ----
    cache = {}
    cols = {f"chi2F{m}_{s}": [] for m in M_LIST for s in ("means", "meansdisp")}
    pf = {(f, m): [] for f in range(6) for m in (0, 4)}
    for _, row in df.iterrows():
        phi = round(float(row.phi_b_deg), 6)
        if phi not in cache:
            cache = {phi: boot_errors(ssl.rotate_z(stars_mw, np.radians(phi)), np.random.default_rng(SEED))}
        L = lib(f"{int(row.idx):03d}.npz")
        for m in M_LIST:
            for s in ("means", "meansdisp"):
                cols[f"chi2F{m}_{s}"].append(chi2_fourier(cache[phi], *L, m, SETS[s])[0])
        for (f, m) in pf:
            pf[(f, m)].append(chi2_fourier(cache[phi], *L, m, [f])[0])
    for k, v in cols.items():
        df[k] = v
    t = df.kind == "truth"
    check = float(df.loc[t, "chi2F12_meansdisp"].iloc[0] - df.loc[t, "chi2b_meansdisp"].iloc[0])
    out["check_m12_equals_sector_chi2_meansdisp"] = check
    out["qmc"] = {}
    variants = [(str(m), m, None) for m in M_LIST] + [(f"2to{m}", m, 0) for m in (4, 6)]
    for key, m, lo in variants:
        for s in ("means", "meansdisp"):
            nd = out["n_data_truth"][s][m] - (0 if lo is None else out["n_data_truth"][s][lo])
            if lo is not None:
                df[f"chi2F{key}_{s}"] = df[f"chi2F{m}_{s}"] - df[f"chi2F{lo}_{s}"]
            df[f"logLs_F{key}_{s}"] = -0.5 * df[f"chi2F{key}_{s}"] / np.sqrt(2 * nd)
            dl = (df[f"logLs_F{key}_{s}"] - df.loc[t, f"logLs_F{key}_{s}"].iloc[0])[~t]
            out["qmc"].setdefault(s, {})[key] = dict(N_data=nd, best=float(dl.max()), best_idx=int(df.idx[dl.idxmax()]),
                                                   n_gt_m2=int((dl > -2).sum()), n_gt_m4p5=int((dl > -4.5).sum()),
                                                   truth_rank=int((dl > 0).sum()) + 1)
    # hybrid: means in modes 0..4, dispersions in modes 2..4 only (their ring averages misfit even at the truth)
    nd = (out["n_data_truth"]["meansdisp"][4] - out["n_data_truth"]["meansdisp"][0]
          + out["n_data_truth"]["means"][0])
    df["chi2F4m_F2to4d"] = df["chi2F4_meansdisp"] - df["chi2F0_meansdisp"] + df["chi2F0_means"]
    df["logLs_F4m_F2to4d"] = -0.5 * df["chi2F4m_F2to4d"] / np.sqrt(2 * nd)
    dl = (df["logLs_F4m_F2to4d"] - df.loc[t, "logLs_F4m_F2to4d"].iloc[0])[~t]
    out["qmc"]["hybrid_F4m_F2to4d"] = dict(N_data=nd, best=float(dl.max()), best_idx=int(df.idx[dl.idxmax()]),
                                           n_gt_m2=int((dl > -2).sum()), n_gt_m4p5=int((dl > -4.5).sum()))
    for sc, spec in SCHEMES.items():
        nd = scheme_chi2(nd_pf, spec)
        df[f"chi2_{sc}"] = scheme_chi2({k: np.array(v) for k, v in pf.items()}, spec)
        df[f"logLs_{sc}"] = -0.5 * df[f"chi2_{sc}"] / np.sqrt(2 * nd)
        dl = (df[f"logLs_{sc}"] - df.loc[t, f"logLs_{sc}"].iloc[0])[~t]
        out["qmc"][sc] = dict(N_data=nd, best=float(dl.max()), best_idx=int(df.idx[dl.idxmax()]),
                              n_gt_m2=int((dl > -2).sum()), n_gt_m4p5=int((dl > -4.5).sum()))
    df.to_csv(OUT_CSV, index=False)
    with open(OUT_CSV.parent / "fourier_compression.json", "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- print ----
    print(f"check: chi2F12_meansdisp - chi2b_meansdisp at the truth = {check:.3e}")
    for sname in SETS:
        print(f"\n=== {sname} ===   m_max: " + "  ".join(f"{m:>16d}" for m in M_LIST))
        print(f"{'N_data':38s}" + "  ".join(f"{out['n_data_truth'][sname][m]:16d}" for m in M_LIST))
        print(f"{'chi2/N_data, truth vs data':38s}" + "  ".join(
            f"{out['data_minus_truth'][sname][m] / out['n_data_truth'][sname][m]:16.2f}" for m in M_LIST))
        for lab, _ in MODELS:
            d = out["models"][lab][sname]
            print(f"{lab:38s}" + "  ".join(f"{d[m]['delta_chi2']:+6.0f} ({d[m]['dlogL']:+5.2f}) {d[m]['noise_free']:4.0f}"
                                           for m in M_LIST))
        nd24 = out["n_data_truth"][sname][4] - out["n_data_truth"][sname][0]
        c24 = out["data_minus_truth"][sname][4] - out["data_minus_truth"][sname][0]
        print(f"  modes 2..4 only (N_data {nd24}, truth chi2/N {c24 / nd24:.2f}):  " + "  ".join(
            f"{lab.split(' (')[0]} {out['models'][lab][sname][4]['delta_chi2'] - out['models'][lab][sname][0]['delta_chi2']:+.0f}"
            f" ({-(out['models'][lab][sname][4]['delta_chi2'] - out['models'][lab][sname][0]['delta_chi2']) / (2 * np.sqrt(2 * nd24)):+.2f})"
            for lab, _ in MODELS[:6]))
    print("\nper field: truth chi2/N_data at m = 0 and in m = 2..4; Delta chi2 of each model at m = 0 / in m = 2..4")
    P = out["per_field"]
    for f in range(6):
        n0, n4 = P["N_data"][f"{FIELDS[f]} m<=0"], P["N_data"][f"{FIELDS[f]} m<=4"]
        c0, c4 = P["truth_chi2"][f"{FIELDS[f]} m<=0"], P["truth_chi2"][f"{FIELDS[f]} m<=4"]
        print(f"  {FIELDS[f]:10s} truth m0 {c0 / n0:.2f} (N {n0})  m2-4 {(c4 - c0) / (n4 - n0):.2f} (N {n4 - n0}) | " + "  ".join(
            f"{lab.split(' (')[0]} {d[f'{FIELDS[f]} m<=0']:+5.0f}/{d[f'{FIELDS[f]} m<=4'] - d[f'{FIELDS[f]} m<=0']:+5.0f}"
            for lab, d in P["models"].items()))
    for sc in SCHEMES:
        print(f"scheme {sc} (N_data {scheme_chi2(nd_pf, SCHEMES[sc])}): " + "  ".join(
            f"{lab} {out['models'][lab][sc][4]['delta_chi2']:+.0f} ({out['models'][lab][sc][4]['dlogL']:+.2f})"
            for lab, _ in MODELS[:6]) + f"  | global search: {out['qmc'][sc]}")
    for s in ("means", "meansdisp"):
        print(f"\nglobal search, {s}: " + "  ".join(
            f"\n   F{m}: best {q['best']:+.2f} (#{q['best_idx']}), n>-2 {q['n_gt_m2']}, n>-4.5 {q['n_gt_m4p5']}, rank {q['truth_rank']}"
            for m, q in out["qmc"][s].items()))

    # ---- figure ----
    colors = [INK, ACCENT, "#e87ba4", BLUE_RAMP[2], BLUE_RAMP[5], "#1baf7a", INK_MUTED]
    styles = ["-", "-", "-", "-", "--", "--", ":"]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), gridspec_kw=dict(wspace=0.3))
    s = "meansdisp"
    nd = np.array([out["n_data_truth"][s][m] for m in M_LIST])
    ax = axes[0]
    for (lab, _), c, ls in zip(MODELS, colors, styles):
        ax.plot(M_LIST, [out["models"][lab][s][m]["noise_free"] for m in M_LIST], ls, marker="o", ms=4, color=c,
                lw=1.6, label=lab)
    ax.set_title("Noise-free separation from the truth model\nΣ (Δ/err)² of library − truth library, modes ≤ m_max",
                 color=INK_2, fontsize=9.5)
    ax.set_ylabel("Σ (Δ/err)²")
    ax.legend(fontsize=7.5, loc="upper left")
    ax = axes[1]
    for (lab, _), c, ls in zip(MODELS[:-1], colors, styles):
        ax.plot(M_LIST, [out["models"][lab][s][m]["dlogL"] for m in M_LIST], ls, marker="o", ms=4, color=c, lw=1.6,
                label=lab)
    ax.axhline(0, color=INK_MUTED, lw=0.8)
    ax.axvline(12, color=INK_MUTED, lw=0.8, ls=":")
    ax.text(11.8, ax.get_ylim()[0] * 0.97, "= current\nsector χ²", ha="right", va="bottom", fontsize=7.5, color=INK_2)
    ax.set_title("Tempered ΔlogL vs the truth on the data\n−Δχ² / (2 √(2 N_data(m_max)))", color=INK_2, fontsize=9.5)
    ax.set_ylabel("ΔlogL (tempered)")
    ax = axes[2]
    inc = np.diff([0.0] + [out["data_minus_truth"][s][m] for m in M_LIST]) / np.diff(np.r_[0, nd])
    ax.plot(M_LIST, inc, "o-", color=INK, lw=1.6, ms=5, label="truth model vs data")
    ax.axhline(1, color=INK_MUTED, lw=0.8, ls="--")
    ax.set_title("χ² per added coefficient, truth model vs data\n(m_max = 0: the ring averages; 1 = pure noise)",
                 color=INK_2, fontsize=9.5)
    ax.set_ylabel("Δχ² / ΔN_data")
    for ax in axes:
        ax.set_xlabel("m_max (even m only; 12 = all 12 sectors)")
        ax.set_xticks(M_LIST)
    fig.suptitle("Fourier compression of the polar moment maps (means + dispersions of v_R, v_φ, v_z; all R, all |z|; "
                 "Ω_b, φ_b at the truth)", y=1.02, fontsize=11)
    savefig(fig, "test1_nbody_fourier_compression")


if __name__ == "__main__":
    main()

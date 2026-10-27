"""Test 1b: do we need both cos and sin Fourier terms, or only the ones allowed by the bar's symmetry?

A steady bar response that is symmetric under (y -> -y, t -> -t) in the bar frame has <v_R> odd about the bar axis
(sin m phi terms only) and <v_phi>, sigma_R, sigma_phi even (cos m phi terms only). Proposed scheme (per ring and
|z| layer, GLS fit of the 15-deg sectors with the data errors, as in test1_fourier_compression.py):
  <v_R>: m = 2, 4 (ring average fitted and projected out);  <v_phi>, sigma_R, sigma_phi: m = 0, 2, 4.
Compared here:
  both       cos and sin of every m (what the proposed scheme used so far)
  symmetric  only the symmetry-allowed terms: <v_R> sin 2phi, sin 4phi; the others 1, cos 2phi, cos 4phi
  forbidden  the rest (both minus symmetric; GLS-nested, so chi2 differences are exact projections)
for the S / bar_frac, bar-angle and pattern-speed offsets of the pilot (other parameters at the truth) and for the
global search (columns chi2_/logLs_fourier_symm_vR_sin24_vphi_sigR_sigphi_cos024 in qmc_evals.csv).
Also measured: the fraction of the m = 2 power of the ground-truth library in the forbidden term.
Numbers: results/test1_nbody/fourier_symmetry_terms.json
"""
import json

import numpy as np
import pandas as pd
import h5py

import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, N_MIN, SEED
from test1_nbody_qmc import OUT_CSV, MOCK, TRUTH
from test1_fourier_compression import lib, NR, NP, NZ, PHI_C

FIELDS = {0: "<v_R>", 1: "<v_phi>", 3: "sigma_R", 4: "sigma_phi"}
BOTH = [(0, "c"), (2, "c"), (2, "s"), (4, "c"), (4, "s")]
SCHEMES = {   # field: (terms, nuisance terms fitted and projected out)
    "both": {0: (BOTH[1:], [(0, "c")]), 1: (BOTH, []), 3: (BOTH, []), 4: (BOTH, [])},
    "symmetric": {0: ([(2, "s"), (4, "s")], [(0, "c")]), 1: ([(0, "c"), (2, "c"), (4, "c")], []),
                  3: ([(0, "c"), (2, "c"), (4, "c")], []), 4: ([(0, "c"), (2, "c"), (4, "c")], [])},
}
COLUMN = "fourier_symm_vR_sin24_vphi_sigR_sigphi_cos024"
PILOT = [("S 1.25", 9, 27.0), ("S 0.8", 8, 27.0), ("f 1.4", 7, 27.0), ("f 0.7", 6, 27.0),
         ("φ_b 17°", 4, 17.0), ("φ_b 37°", 5, 37.0), ("Ω_b 20", 2, 27.0), ("Ω_b 32", 3, 27.0)]


def design(phi, terms):
    return np.column_stack([np.ones_like(phi) if m == 0 else (np.cos(m * phi) if k == "c" else np.sin(m * phi))
                            for m, k in terms])


def gls(F, w, d):
    b = F.T @ (w * d)
    return float(b @ np.linalg.solve(F.T @ (w[:, None] * F), b))


def chi2_scheme(obs, lib_m1, lib_s, scheme):
    """chi^2, N_data of a per-field term scheme (blocks with too few sectors keep their sectors)."""
    N, m1, s, em, es = obs
    good = ((N >= N_MIN) & np.all(np.isfinite(lib_m1), 1) & np.all(em > 0, 1) & np.all(es > 0, 1)).reshape(NR, NP, NZ)
    D = np.column_stack([m1 - lib_m1, s - lib_s]).reshape(NR, NP, NZ, 6)
    E = np.column_stack([em, es]).reshape(NR, NP, NZ, 6)
    chi, nd = 0.0, 0
    for r in range(NR):
        for z in range(NZ):
            av = good[r, :, z]
            n = int(av.sum())
            if n == 0:
                continue
            for f, (terms, nuis) in scheme.items():
                w, d = 1.0 / E[r, av, z, f]**2, D[r, av, z, f]
                F = design(PHI_C[av], nuis + terms)
                full = F.shape[1] <= n and np.linalg.matrix_rank(F) == F.shape[1]
                F = F if full else np.eye(n)
                c0 = gls(design(PHI_C[av], nuis), w, d) if nuis else 0.0
                chi += gls(F, w, d) - c0
                nd += F.shape[1] - len(nuis)
    return chi, nd


def main():
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    obs = {phi: boot_errors(ssl.rotate_z(stars, np.radians(phi)), np.random.default_rng(SEED)) for phi in (27.0, 17.0, 37.0)}
    truth = lib("pilot_00.npz")
    out = {"pilot": {}, "n_data": {}}
    ref = {sc: chi2_scheme(obs[27.0], *truth, SCHEMES[sc]) for sc in SCHEMES}
    for sc in SCHEMES:
        out["n_data"][sc] = ref[sc][1]
    out["n_data"]["forbidden"] = ref["both"][1] - ref["symmetric"][1]
    out["truth_chi2_per_N"] = {sc: ref[sc][0] / ref[sc][1] for sc in SCHEMES}
    out["truth_chi2_per_N"]["forbidden"] = (ref["both"][0] - ref["symmetric"][0]) / out["n_data"]["forbidden"]
    for lab, i, phi in PILOT:
        L = lib(f"pilot_{i:02d}.npz")
        c = {sc: chi2_scheme(obs[phi], *L, SCHEMES[sc])[0] - ref[sc][0] for sc in SCHEMES}
        c["forbidden"] = c["both"] - c["symmetric"]
        out["pilot"][lab] = {sc: dict(delta_chi2=v, dlogL=-v / (2 * np.sqrt(2 * out["n_data"][sc]))) for sc, v in c.items()}

    # fraction of the ground-truth library's m = 2 power in the forbidden term, per field (|z| < 1, all rings)
    N, m1, s, em, es = obs[27.0]
    good = ((N >= N_MIN) & np.all(em > 0, 1)).reshape(NR, NP, NZ)
    val = np.column_stack([truth[0], truth[1]]).reshape(NR, NP, NZ, 6)
    err = np.column_stack([em, es]).reshape(NR, NP, NZ, 6)
    out["truth_library_m2_forbidden_power_fraction"] = {}
    for f, name in FIELDS.items():
        pa = pf = 0.0
        for r in range(NR):
            for z in range(2):
                av = good[r, :, z]
                if av.sum() < 5:
                    continue
                F = design(PHI_C[av], BOTH)
                w = 1 / err[r, av, z, f]**2
                c = np.linalg.solve(F.T @ (w[:, None] * F), F.T @ (w * val[r, av, z, f]))
                a2c, a2s = c[1], c[2]
                allowed, forbid = (a2s, a2c) if f == 0 else (a2c, a2s)
                pa += allowed**2
                pf += forbid**2
        out["truth_library_m2_forbidden_power_fraction"][name] = pf / (pa + pf)

    # global search column for the symmetric scheme
    df = pd.read_csv(OUT_CSV).sort_values("idx").reset_index(drop=True)
    cache, chis, chis_both = {}, [], []
    for _, row in df.iterrows():
        phi = round(float(row.phi_b_deg), 6)
        if phi not in cache:
            cache = {phi: boot_errors(ssl.rotate_z(stars, np.radians(phi)), np.random.default_rng(SEED))}
        L = lib(f"{int(row.idx):03d}.npz")
        chis.append(chi2_scheme(cache[phi], *L, SCHEMES["symmetric"])[0])
        chis_both.append(chi2_scheme(cache[phi], *L, SCHEMES["both"])[0])
    df[f"chi2_{COLUMN}"] = chis
    df[f"logLs_{COLUMN}"] = -0.5 * df[f"chi2_{COLUMN}"] / np.sqrt(2 * out["n_data"]["symmetric"])
    t = df.kind == "truth"
    out["check_both_equals_existing_column"] = float(
        chis_both[int(np.flatnonzero(t)[0])] - df.loc[t, "chi2_fourier_vR_m24_vphi_sigR_sigphi_m024"].iloc[0])
    df.to_csv(OUT_CSV, index=False)
    with open(OUT_CSV.parent / "fourier_symmetry_terms.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print("check (both == existing proposed-scheme column at the truth):", f"{out['check_both_equals_existing_column']:.2e}")
    print("N_data:", out["n_data"], " truth chi2/N:", {k: round(v, 2) for k, v in out["truth_chi2_per_N"].items()})
    print("fraction of the truth library's m = 2 power in the forbidden term (|z| < 1):",
          {k: round(v, 3) for k, v in out["truth_library_m2_forbidden_power_fraction"].items()})
    print(f"\n{'offset':10s} " + "  ".join(f"{sc:>22s}" for sc in ("both", "symmetric", "forbidden")))
    for lab, d in out["pilot"].items():
        print(f"{lab:10s} " + "  ".join(f"{d[sc]['delta_chi2']:+8.0f} ({d[sc]['dlogL']:+6.2f})    " for sc in ("both", "symmetric", "forbidden")))


if __name__ == "__main__":
    main()

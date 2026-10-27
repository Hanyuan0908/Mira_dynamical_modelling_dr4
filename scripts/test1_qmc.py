"""Test 1 quick emulator run: steady-state log-likelihood of the perfect mock (50k stars, R < 8 kpc)
at the truth and at a scrambled Sobol design of N_QMC points in (|Omega_b|, phi_b, ln S, ln bar_frac).
Rows are appended to results/test1_qmc/evals.csv as they finish, and the run resumes if restarted.
Per-star log p(v_i|x_i) of every point is saved to results/test1_qmc/logp/NNN.npy (mock star order),
so subsets (e.g. radial ranges) can be analysed later; points finished before this was added are
recomputed at the end and checked against their logged logL.
"""
import csv
import time
from pathlib import Path

import numpy as np
import h5py
from scipy.stats import qmc

from common import ROOT, MOCK_DIR
import steady_state_likelihood as ssl

OUT_DIR = ROOT / "results" / "test1_qmc"
OUT_CSV = OUT_DIR / "evals.csv"
LOGP_DIR = OUT_DIR / "logp"
MOCK = MOCK_DIR / "mock_perfect_R8_N50k.hdf5"
N_QMC = 64
SEED = 2026
TRUTH = dict(Omega_abs=36.27, phi_b_deg=27.0, S=1.0, bar_frac=1.0)
# prior box: |Omega_b| [km/s/kpc], phi_b [deg], ln S, ln bar_frac
BOX = dict(Omega_abs=(28.0, 52.0), phi_b_deg=(15.0, 45.0), lnS=(np.log(0.8), np.log(1.7)),
           lnf=(np.log(0.5), np.log(1.6)))
FIELDS = ["idx", "kind", "Omega_abs", "phi_b_deg", "S", "bar_frac", "logL", "logp_median", "tau_median",
          "t_integrate_s", "t_likelihood_s"]


def design():
    u = qmc.Sobol(d=4, scramble=True, seed=SEED).random(N_QMC)
    lo = np.array([v[0] for v in BOX.values()])
    hi = np.array([v[1] for v in BOX.values()])
    p = lo + u * (hi - lo)
    pts = [dict(kind="truth", **TRUTH)]
    pts += [dict(kind="qmc", Omega_abs=a, phi_b_deg=b, S=float(np.exp(c)), bar_frac=float(np.exp(d)))
            for a, b, c, d in p]
    return pts


def main():
    LOGP_DIR.mkdir(parents=True, exist_ok=True)
    done = set()
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            done = {int(r["idx"]) for r in csv.DictReader(f)}
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    pts = design()
    t_start = time.time()
    with open(OUT_CSV, "a", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=FIELDS)
        if not done:
            wr.writeheader()
        for i, p in enumerate(pts):
            if i in done:
                continue
            res, logp = ssl.evaluate(stars, p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"])
            np.save(LOGP_DIR / f"{i:03d}.npy", logp.astype(np.float32))
            row = dict(idx=i, **p, **{k: res[k] for k in FIELDS if k in res})
            wr.writerow(row)
            f.flush()
            el = time.time() - t_start
            print(f"[{i + 1}/{len(pts)}] {p['kind']:5s} Omega={p['Omega_abs']:.2f} phi={p['phi_b_deg']:.1f} "
                  f"S={p['S']:.3f} f={p['bar_frac']:.3f} logL={res['logL']:.1f} "
                  f"({res['t_integrate_s'] + res['t_likelihood_s']:.0f} s; elapsed {el / 60:.1f} min)", flush=True)

    # backfill per-star log p for points evaluated before it was saved
    with open(OUT_CSV) as f:
        logged = {int(r["idx"]): float(r["logL"]) for r in csv.DictReader(f)}
    for i, p in enumerate(pts):
        if (LOGP_DIR / f"{i:03d}.npy").exists():
            continue
        res, logp = ssl.evaluate(stars, p["Omega_abs"], p["phi_b_deg"], p["S"], p["bar_frac"])
        np.save(LOGP_DIR / f"{i:03d}.npy", logp.astype(np.float32))
        print(f"backfill [{i}] logL {res['logL']:.3f} vs logged {logged[i]:.3f} "
              f"(diff {res['logL'] - logged[i]:+.3f})", flush=True)


if __name__ == "__main__":
    main()

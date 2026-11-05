"""Is the phase drift caused by stars being in the wrong PLACE, or having the wrong VELOCITY?

A period-luminosity distance error does two things at once, and they are usually conflated:

  position   the star is placed at the wrong distance along its line of sight;
  velocity   its tangential velocity is v_t = 4.74 * mu * d, so a distance wrong by 12.5% makes
             the inferred v_t wrong by 12.5% as well. v_los is untouched, so the velocity VECTOR
             is also rotated, not just rescaled.

Both feed the m = 2 phase of <|v_R|/v_tot>. Separating them decides what a mechanism has to
explain. The two are isolated by rescaling the proper motions:

  full      position at d_obs, proper motions as observed          (what the pipeline does)
  pos only  position at d_obs, proper motions * d_true/d_obs       (true velocity, wrong place)
  vel only  position at d_true, proper motions * d_obs/d_true      (right place, wrong velocity)

The selection mask is taken from the full case and reused for all three, so the three variants
contain exactly the same stars and differ only in where they are put and how fast they move.

Usage: python diagnose_drift_position_vs_velocity.py [--workers 4] [--n-real 4]
Table: results/bar_angle_drift_pos_vs_vel.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import h5py
from scipy import stats

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy
from bar_angle_vecsum_verification import phase
from scan_bar_angle_windows import NBODY, AURIGA
import observe, selection

FRAC = 0.125
WIN = "kin0.4"


def variants(gname, seed):
    """m = 2 phase of <|v_R|/v_tot> for the three error decompositions, plus the no-error case."""
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{FRAC * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*W.T)
    d0 = o["d"]
    d1 = observe.scatter_distance_modulus(d0, FRAC * 5 / np.log(10),
                                          np.random.default_rng(seed))
    r = d0 / d1                                   # true / observed

    def gc(d, pml, pmb):
        return np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, pml, pmb, o["vlos"]))

    # one selection mask, from the full (realistic) case, reused everywhere
    w_full = gc(d1, o["pml"], o["pmb"])
    S = selection.evaluate(tab, o["l"], o["b"], d1, np.hypot(w_full[:, 0], w_full[:, 1]))
    keep = np.random.default_rng(seed + 500_000).random(len(S)) < np.minimum(S, 1.0)

    out = {}
    for name, w in (("full", w_full),
                    ("pos", gc(d1, o["pml"] * r, o["pmb"] * r)),
                    ("vel", gc(d0, o["pml"] / r, o["pmb"] / r)),
                    ("none", gc(d0, o["pml"], o["pmb"]))):
        out[name] = phase(w[keep], "absvR", WIN) + BAR_ANGLE_DEG
    return gname, out


def job(a):
    return variants(*a)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-real", type=int, default=4)
    ap.add_argument("--snaps", nargs="+", default=AURIGA)
    args = ap.parse_args()

    jobs = [(g, 11 + s) for g in args.snaps for s in range(args.n_real)]
    print(f"{len(jobs)} runs over {len(args.snaps)} galaxies", flush=True)
    acc = {}
    with ProcessPoolExecutor(args.workers) as ex:
        for k, (g, ph) in enumerate(ex.map(job, jobs), 1):
            acc.setdefault(g, []).append(ph)
            if k % 8 == 0:
                print(f"  {k}/{len(jobs)}", flush=True)

    rows = []
    for g, lst in acc.items():
        m = {k: float(np.mean([x[k] for x in lst])) for k in lst[0]}
        rows.append(dict(galaxy=g, phase_none=m["none"], phase_full=m["full"],
                         drift_full=m["full"] - m["none"],
                         drift_pos=m["pos"] - m["none"],
                         drift_vel=m["vel"] - m["none"]))
    D = pd.DataFrame(rows).sort_values("drift_full").reset_index(drop=True)
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    h = [int(g[2:]) for g in D.galaxy]
    D["R_CR"] = Q.loc[h, "R_CR"].values
    D["Omega_p"] = Q.loc[h, "Omega_p"].values
    D.to_csv(RESULTS_DIR / "bar_angle_drift_pos_vs_vel.csv", index=False)

    print("\nwhere the 12.5% drift comes from (deg, <|v_R|/v_tot>, window kin0.4)\n")
    print(f"  {'galaxy':8s} {'R_CR':>6s} | {'full':>7s} = {'position':>9s} + {'velocity':>9s} "
          f"{'(sum)':>8s}")
    print("  " + "-" * 60)
    for _, r in D.iterrows():
        print(f"  {r.galaxy:8s} {r.R_CR:6.2f} | {r.drift_full:+7.2f} = {r.drift_pos:+9.2f} + "
              f"{r.drift_vel:+9.2f} {r.drift_pos + r.drift_vel:+8.2f}")
    print("  " + "-" * 60)
    print(f"  {'mean |.|':8s} {'':6s} | {D.drift_full.abs().mean():7.2f}   "
          f"{D.drift_pos.abs().mean():9.2f}   {D.drift_vel.abs().mean():9.2f}")

    print()
    for c in ("drift_full", "drift_pos", "drift_vel"):
        rr = stats.pearsonr(D.R_CR, D[c])
        ss = stats.spearmanr(D.R_CR, D[c])
        print(f"  corr(R_CR, {c:10s}) = {rr.statistic:+.2f} (p = {rr.pvalue:.4f}), "
              f"rho = {ss.statistic:+.2f}")
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_angle_drift_pos_vs_vel.csv')}")


if __name__ == "__main__":
    main()

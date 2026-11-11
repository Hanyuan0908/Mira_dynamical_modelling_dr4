"""A run's own scrambled Sobol design, for searches that do not re-use the focused design (e.g. the real data).

Same conventions as the focused design: |Omega_b| and phi_b uniform, S and bar_frac log-uniform, which is also the
flat prior of the emulator (test1_gp_emulator works in (Omega, phi, ln S, ln f)). Every row is kind "qmc"; there
is no truth row.

Usage: python make_sobol_design.py --run osarg_zhang2024_mw_hunter24 --n 256 --seed 2041 \
           --omega 25 45 --phi 20 30 --S 0.5 1.5 --bar-frac 0.5 1.5
       python make_sobol_design.py --run osarg_zhang2024_mw_hunter24 --append --n 128 --seed 2042 \
           --omega 15 25 --phi 20 30 --S 0.5 1.5 --bar-frac 0.5 1.5
--append extends an existing design into an adjacent slab: the new points get the next idx and batch, the points
already evaluated are kept, and the recorded box becomes the union. Give the slab as many points per unit volume
as the existing design has, so the combined design stays evenly filled.
Output: results/<run>/design.csv and results/<run>/design_box.json (the box in the emulator's coordinates)
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import qmc

from common import ROOT, rel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--n", type=int, default=256, help="number of points (a power of 2 keeps Sobol balanced)")
    ap.add_argument("--seed", type=int, default=2041)
    ap.add_argument("--omega", type=float, nargs=2, required=True, metavar=("LO", "HI"))
    ap.add_argument("--phi", type=float, nargs=2, required=True, metavar=("LO", "HI"))
    ap.add_argument("--S", type=float, nargs=2, required=True, metavar=("LO", "HI"))
    ap.add_argument("--bar-frac", type=float, nargs=2, required=True, metavar=("LO", "HI"))
    ap.add_argument("--append", action="store_true", help="add the points to the run's existing design")
    args = ap.parse_args()
    box = dict(Omega_abs=tuple(args.omega), phi_b_deg=tuple(args.phi), lnS=tuple(np.log(args.S)),
               lnf=tuple(np.log(args.bar_frac)))
    lo, hi = np.array([v[0] for v in box.values()]), np.array([v[1] for v in box.values()])
    u = qmc.Sobol(4, scramble=True, seed=args.seed).random(args.n)
    p = lo + u * (hi - lo)
    out = ROOT / "results" / args.run
    out.mkdir(parents=True, exist_ok=True)
    old, meta = None, dict(n=args.n, seed=args.seed, slabs=[])
    if args.append:
        old = pd.read_csv(out / "design.csv")
        meta = json.load(open(out / "design_box.json"))
        meta.setdefault("slabs", [dict(box=meta["box"], n=meta["n"], seed=meta["seed"])])
        box = {k: (min(box[k][0], meta["box"][k][0]), max(box[k][1], meta["box"][k][1])) for k in box}
    first, batch = (0, 1) if old is None else (int(old.idx.max()) + 1, int(old.batch.max()) + 1)
    df = pd.DataFrame(dict(idx=first + np.arange(args.n), kind="qmc", batch=batch, Omega_abs=p[:, 0],
                           phi_b_deg=p[:, 1], S=np.exp(p[:, 2]), bar_frac=np.exp(p[:, 3])))
    meta["slabs"].append(dict(box={k: [float(a), float(b)] for k, (a, b) in zip(box, zip(lo, hi))},
                              n=args.n, seed=args.seed, first_idx=first))
    if old is not None:
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(out / "design.csv", index=False)
    meta.update(box={k: list(map(float, v)) for k, v in box.items()}, n=len(df),
                sampling="scrambled Sobol; Omega, phi uniform; S, bar_frac log-uniform")
    json.dump(meta, open(out / "design_box.json", "w"), indent=2)
    print(f"wrote {rel(out / 'design.csv')}: {len(df)} points ({args.n} new from idx {first}, seed {args.seed})")
    print(df.describe().loc[["min", "max"], ["Omega_abs", "phi_b_deg", "S", "bar_frac"]].round(3).to_string())


if __name__ == "__main__":
    main()

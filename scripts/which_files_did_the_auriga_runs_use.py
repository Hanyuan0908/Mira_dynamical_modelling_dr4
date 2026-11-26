"""Which trial-family files did the cluster Auriga runs integrate in? (2026-11-26)

The 13 cluster-built Auriga families (copied to the portable drive) carry the bar in their "axi" files
(check_auriga_family_axisymmetry.py). To test whether the Test 2 runs used them, restart saved orbits from their own
first snapshot (traj_NNN.npy holds the bar-frame phase space at t = 0, 25, 50, ... Myr) and integrate 100 Myr in
  (a) the cluster family files, exactly as the run did (steady_state_likelihood.bar_aligned_potential), and
  (b) a correct family: Au<N>_baryon_axi re-expanded here with symmetry="axisymmetric" (as potential_family.build_files
      does, laptop AGAMA), everything else identical.
Compared with the saved snapshots 1-4 (each run's own cadence: batch 1 saved every 25 Myr, batch 2 every 50 Myr):
the median position offset over 400 stars. Design points with bar_frac far
from 1 (where (a) and (b) differ most), and one near 1 as a control.
Usage: python which_files_did_the_auriga_runs_use.py [--halo 9]  ->  results/which_files_auriga_runs_used_au<N>.csv
"""
import argparse
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import agama

import galaxies
import steady_state_likelihood as ssl
import test1_nbody_qmc as old
from build_au18_potential import CYL_PARAMS
from common import RESULTS_DIR
from test2b_qmc_distance_error_run import scaled_design

USB = Path("/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4")
N_STARS, N_STEP = 400, 4


def run(pot, ic, Omega_abs, n_saved):
    dt = ssl.T_GYR / (n_saved - 1) / ssl.TIME_UNIT_GYR       # the saved file's own cadence (25 or 50 Myr)
    orbs = agama.orbit(potential=pot, ic=ic, time=N_STEP * dt, trajsize=N_STEP + 1, Omega=-abs(Omega_abs),
                       accuracy=old.ACCURACY, **ssl.ORBIT_KW)
    return np.stack(orbs[:, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--halo", type=int, default=9)
    ap.add_argument("--points", type=int, nargs="+", default=[65, 98, 120, 132])
    args = ap.parse_args()
    agama.setUnits(mass=1, length=1, velocity=1)
    g, pd_ = f"au{args.halo}", USB / "mock" / f"au{args.halo}_potential"
    run_dir = USB / "runs" / f"test2_{g}_sf_12p5pct"
    design, _ = scaled_design(g)
    design = design.set_index("idx")
    # (b) a correct family in a scratch folder: same files, axisymmetric baryons rebuilt here
    tmp = Path(tempfile.mkdtemp(prefix=f"{g}_correct_"))
    for f in pd_.iterdir():
        (tmp / f.name).symlink_to(f)
    (tmp / f"Au{args.halo}_baryon_axi.ini").unlink()
    grid = {k: v for k, v in CYL_PARAMS.items() if k not in ("type", "symmetry", "mmax")}
    agama.Potential(type="CylSpline", potential=agama.Potential(file=str(pd_ / f"Au{args.halo}_baryon_full.ini")),
                    symmetry="axisymmetric", **grid).export(str(tmp / f"Au{args.halo}_baryon_axi.ini"))
    galaxies.galaxy(g)
    rows = []
    for i in args.points + [int((design.bar_frac - 1).abs().idxmin())]:
        p = design.loc[i]
        saved = np.load(run_dir / f"traj_{i:03d}.npy", mmap_mode="r")
        sel = np.linspace(0, saved.shape[0] - 1, N_STARS).astype(int)
        ic = np.asarray(saved[sel, 0], dtype=np.float64)
        ref = np.asarray(saved[sel, 1:N_STEP + 1, :3], dtype=np.float64)
        r = dict(idx=i, Omega_abs=p.Omega_abs, phi_b=p.phi_b_deg, S=p.S, bar_frac=p.bar_frac)
        for lab, d in (("cluster_files", pd_), ("correct_family", tmp)):
            galaxies.GALAXIES[g]["pot_dir"] = d
            ssl._CACHE.pop(g, None)
            tr = run(ssl.bar_aligned_potential(p.bar_frac, p.S, g), ic, p.Omega_abs, saved.shape[1])
            off = np.linalg.norm(tr[:, 1:, :3] - ref, axis=2)
            for k in range(N_STEP):
                r[f"{lab}_median_offset_kpc_step{k + 1}"] = float(np.median(off[:, k]))
        r["saved_cadence_myr"] = 1000 * ssl.T_GYR / (saved.shape[1] - 1)
        rows.append(r)
        print(f"idx {i:3d} (Omega {p.Omega_abs:5.1f}, phi {p.phi_b_deg:4.1f}, S {p.S:.2f}, bar_frac {p.bar_frac:.2f}): "
              f"median |x - saved| after 1-4 saved steps of {r['saved_cadence_myr']:.0f} Myr: cluster files "
              + " ".join(f"{r[f'cluster_files_median_offset_kpc_step{k + 1}']:.1e}" for k in range(N_STEP))
              + " | correct family " + " ".join(f"{r[f'correct_family_median_offset_kpc_step{k + 1}']:.1e}" for k in range(N_STEP)),
              flush=True)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"which_files_auriga_runs_used_{g}.csv", index=False)


if __name__ == "__main__":
    main()

"""Direct chi2 (12 jittered runs each) at the typical parameters of chain windows of the OSARG Sormani+22 posterior:
the two |Omega_b| peaks (argument "peaks") or the valley between them ("valley"). See check_bimodality_sormani22.py."""
import sys, json, time
sys.path.insert(0, "/Users/hanyuan/Library/CloudStorage/Dropbox/python_script/mira_dynamics_dr4/scripts")
import numpy as np, pandas as pd, h5py
import test1_gp_emulator as ge
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
from rescore_empirical_selection import data_lbdR, library_fit_empirical
from chi2_jitter_small_perturbations import RANDOM_WIDTH
from check_bimodality_sormani22 import sample, RUN, MOCK, B
ge.use_run(RUN)
ev = pd.read_csv(RESULTS_DIR / RUN / "evals.csv"); tr = ev[ev.kind == "qmc"].reset_index(drop=True)
y = tr[f"logLs_{B}_sfratio"].values; y = y - y.max(); u = ge.to_unit(tr)
emu = ge.Emulator(u, y); ch = sample(emu, u, y); P = ge.from_unit(ch)
nd = int(tr[f"ndata_{B}_sfratio"].mode().iloc[0]); temp = 2 * np.sqrt(2 * nd)
pts = {}
WINDOWS = {"peaks": (("peak ~30", 28.0, 32.0), ("peak ~39", 38.0, 41.0)), "valley": (("valley ~34.5", 33.5, 35.5),)}
for tag, lo, hi in WINDOWS[sys.argv[1] if len(sys.argv) > 1 else "peaks"]:
    m = (P[:, 0] >= lo) & (P[:, 0] < hi)
    med = np.median(P[m], axis=0)
    pts[tag] = med
    print(f"{tag}: {m.mean():.2f} of the chain; median parameters Omega {med[0]:.2f}, phi {med[1]:.2f}, S {med[2]:.3f}, f {med[3]:.3f}; emulator logL {float(emu.predict(ge.to_unit(pd.DataFrame([med], columns=ge.KEYS))[0])[0]):.2f}", flush=True)
meta = json.load(open(RESULTS_DIR / RUN / "run_meta.json")); sigma_dm = float(meta["sigma_DM"])
mock = MOCK_DIR / MOCK; dat = data_lbdR(mock)
with h5py.File(mock, "r") as f:
    w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
cd0, Cd, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot); cd0 = cd0.reshape(fs.NB, -1)
nr, spec = fs.SCHEMES[B]; idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
hart = (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1)
W = {b: np.linalg.inv(np.linalg.cholesky(Cd[b][np.ix_(idx, idx)])) for b in (r * fs.NL + z for r in range(nr) for z in range(fs.NL)) if ok[b]}
ssl.ACCURACY = old.ACCURACY; rng = np.random.default_rng(2125); res = {}
for tag, p in pts.items():
    vals = []
    for r in range(12):
        q = dict(zip(ge.KEYS, p))
        if r:
            q = {k: v + RANDOM_WIDTH[k] * rng.standard_normal() for k, v in q.items()}
        pot = ssl.bar_aligned_potential(q["bar_frac"], q["S"], "sormani22")
        _, traj = ssl.integrate(pot, ssl.rotate_z(w, np.radians(q["phi_b_deg"])), -abs(q["Omega_abs"]))
        cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, ::2]), q["phi_b_deg"], sigma_dm, dat, ok)
        cl = cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - q["phi_b_deg"]).T
        blocks = [b for b in W if np.all(np.isfinite(cl[b][idx]))]
        vals.append(hart * sum(float(np.sum((W[b] @ (cd0[b][idx] - cl[b][idx])) ** 2)) for b in blocks))
    res[tag] = np.array(vals)
    print(f"{tag}: chi2 {res[tag].mean():.1f} +- {res[tag].std(ddof=1):.1f} (sem {res[tag].std(ddof=1)/np.sqrt(12):.1f}); runs " + " ".join(f"{v:.0f}" for v in vals), flush=True)
tag_run = sys.argv[1] if len(sys.argv) > 1 else "peaks"
if tag_run == "peaks":
    a, b = res["peak ~30"], res["peak ~39"]
    d, se = b.mean() - a.mean(), np.sqrt(a.var(ddof=1) / 12 + b.var(ddof=1) / 12)
    print(f"peak ~39 minus peak ~30: dchi2 {d:+.1f} +- {se:.1f} = {-d / temp:+.2f} +- {se / temp:.2f} tempered logL")
json.dump({k: v.tolist() for k, v in res.items()} | {"params": {k: v.tolist() for k, v in pts.items()}},
          open(RESULTS_DIR / RUN / f"bimodality_{tag_run}_direct.json", "w"), indent=2)

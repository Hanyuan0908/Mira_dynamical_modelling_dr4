"""Bar pattern speed of every barred Auriga z = 0 galaxy, by the Dehnen et al. (2023) method.

Dehnen, Semczuk & Schoenrich (2023, MNRAS 518, 2712) estimate the pattern speed from a single
snapshot, using the continuity equation applied to the m = 2 Fourier moments: no time differencing,
so it works on one z = 0 output. Their FourierMethod also returns the radial range it used, which is
a useful consistency check against the m = 2 amplitude profile.

The method assumes an anticlockwise-rotating disc. prepare_au18_mwframe puts Auriga discs in the
MW-like frame, which rotates clockwise, so y and v_y are mirrored before the call - exactly as
measure_nbody_pattern_speed_evolution does for the isolated N-body run.

Corotation is where Omega(R) from the enclosed mass crosses the pattern speed; R_CR / R_bar is the
usual fast/slow classification (fast if < 1.4).

Usage: python measure_auriga_pattern_speed.py [--halos 2 9 ...]
Table: results/auriga_pattern_speed.csv
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import auriga as au
from common import RESULTS_DIR, EXTERNAL_DIR, rel
from prepare_au18_mwframe import shrinking_sphere, measure_bar, RMAX, R_DISC
from survey_auriga_bars import ROOT

PATTERN_SPEED_DIR = os.environ.get("PATTERN_SPEED_DIR") or str(EXTERNAL_DIR / "patternSpeed-main")
if not hasattr(np, "asfarray"):            # removed in NumPy 2.0; patternSpeed.py still calls it
    np.asfarray = lambda a, dtype=np.float64: np.asarray(a, dtype=dtype)
sys.path.insert(0, PATTERN_SPEED_DIR)
import patternSpeed  # noqa: E402

G = 4.300917270e-6     # kpc (km/s)^2 / Msun


def circular_frequency(x, pot, rmax=25.0, n=160):
    """Omega(R) from the snapshot's OWN potential, so every mass component is included.

    The stellar mass alone is a large underestimate inside 25 kpc, so an Omega(R) built from it
    never rises to meet a fast bar's pattern speed and corotation comes out undefined. Arepo stores
    the total potential at each particle, which includes dark matter and gas; differentiating the
    azimuthally averaged profile gives the true circular frequency,
        v_c^2 = R dPhi/dR,   Omega = v_c / R.
    """
    R = np.hypot(x[:, 0], x[:, 1])
    thin = np.abs(x[:, 2]) < 1.0
    edges = np.linspace(0.5, rmax, n + 1)
    Rc = 0.5 * (edges[1:] + edges[:-1])
    prof = np.full(n, np.nan)
    for i in range(n):
        m = thin & (R >= edges[i]) & (R < edges[i + 1])
        if m.sum() > 50:
            prof[i] = np.mean(pot[m])
    ok = np.isfinite(prof)
    Rc, prof = Rc[ok], prof[ok]
    dphi = np.gradient(prof, Rc)
    vc2 = Rc * dphi
    good = vc2 > 0
    return Rc[good], np.sqrt(vc2[good]) / Rc[good]


def one(halo, root):
    hd = Path(root) / f"halo_{halo}"
    hdr = au.read_header(hd / "snapdir_127")
    h = au.main_halo(hd / "groups_127", hdr)
    st = au.read_particles(hd / "snapdir_127", 4,
                           ["Velocities", "Masses", "GFM_StellarFormationTime", "Potential"],
                           h["pos"], RMAX, hdr)
    keep = st["GFM_StellarFormationTime"] > 0
    st = {k: v[keep] for k, v in st.items()}
    x = st["Coordinates"].astype(np.float64)
    v = st["Velocities"].astype(np.float64)
    m = st["Masses"].astype(np.float64) * 1e10          # Msun
    x -= shrinking_sphere(x, m)
    s = np.linalg.norm(x, axis=1) < R_DISC
    v -= np.average(v[s], axis=0, weights=m[s])
    L = np.sum(m[s, None] * np.cross(x[s], v[s]), axis=0)
    R = au.rotation_to(L, np.array([0.0, 0.0, -1.0]))    # MW-like: clockwise
    x, v = x @ R.T, v @ R.T
    bar = measure_bar(x[:, 0], x[:, 1], x[:, 2], m)
    # mirror y so the method sees an anticlockwise disc
    pat, _ = patternSpeed.FourierMethod.patternSpeed(x[:, 0], -x[:, 1], v[:, 0], -v[:, 1], mu=m)
    Om = float(pat["Om"])
    Rc, Omega = circular_frequency(x, st["Potential"].astype(np.float64))
    f = Omega - Om
    i = np.flatnonzero(np.diff(np.sign(f)) != 0)
    Rcr = float(Rc[i[0]] - f[i[0]] * (Rc[i[0] + 1] - Rc[i[0]]) / (f[i[0] + 1] - f[i[0]])) \
        if len(i) else np.nan
    return dict(halo=halo, Omega_p=abs(Om), Omega_p_err=float(pat["Om_e"]),
                fit_R0=float(pat["R0"]), fit_R1=float(pat["R1"]),
                A2_max=float(bar["A2_max"]), R_peak=float(bar["R_peak"]),
                R_bar=float(bar["R_bar"]), R_CR=Rcr,
                R_CR_over_R_bar=Rcr / bar["R_bar"] if bar["R_bar"] > 0 else np.nan,
                M_star_disc=float(m[np.hypot(x[:, 0], x[:, 1]) < R_DISC].sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--halos", nargs="+", type=int, default=None)
    ap.add_argument("--root", default=ROOT)
    args = ap.parse_args()
    if args.halos:
        halos = args.halos
    else:
        s = pd.read_csv(RESULTS_DIR / "auriga_bar_survey.csv")
        halos = sorted(s[s.barred & (s.R_peak < 5) & s.R_bar.notna() & (s.R_bar < 8)].halo)
        print(f"barred halos with a well-defined bar: {halos}\n")
    rows = []
    for h in halos:
        try:
            r = one(h, args.root)
            rows.append(r)
            print(f"  Au{h:<3d} Omega_p {r['Omega_p']:6.2f} +/- {r['Omega_p_err']:4.2f}  "
                  f"(fit {r['fit_R0']:.1f}-{r['fit_R1']:.1f} kpc)  A2 {r['A2_max']:.3f}  "
                  f"R_bar {r['R_bar']:4.2f}  R_CR {r['R_CR']:5.2f}  "
                  f"R_CR/R_bar {r['R_CR_over_R_bar']:4.2f}", flush=True)
        except Exception as e:
            print(f"  Au{h:<3d} FAILED: {e.__class__.__name__}: {e}", flush=True)
    D = pd.DataFrame(rows).sort_values("Omega_p")
    D.to_csv(RESULTS_DIR / "auriga_pattern_speed.csv", index=False)
    print(f"\nwrote {rel(RESULTS_DIR / 'auriga_pattern_speed.csv')}")
    print(f"  Omega_p spans {D.Omega_p.min():.1f}-{D.Omega_p.max():.1f} km/s/kpc over {len(D)} galaxies")
    print(f"  fast bars (R_CR/R_bar < 1.4): {int((D.R_CR_over_R_bar < 1.4).sum())} of {len(D)}")


if __name__ == "__main__":
    main()

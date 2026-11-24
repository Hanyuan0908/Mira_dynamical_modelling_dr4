"""Why do the degeneracies differ between the Milky Way family and the mocks' families? Gradients of the predictions.

For a (star set, family) pair at a reference point (|Omega_b|, phi_b, S, bar_frac), the predicted scheme B coefficients
are differentiated by central differences, one parameter at a time:
  |Omega_b| +-10%, phi_b +-3 deg, S +-0.15, bar_frac +-0.25,
each library built from a 4 Gyr integration at 50 Myr (75 phases; noise floor ~0.65 tempered, library_noise_floor_test)
so the differences are not swamped by Monte Carlo noise. The differences are whitened with the star set's data
covariance (scheme B terms, the blocks valid at the reference), giving one vector g_i per parameter (per unit of the
parameter). Then:
  - cos(g_i, g_j): parameters whose predicted changes look alike are degenerate;
  - |g_i|: how strongly the kinematics react (tempered chi2 per unit squared, |g|^2 / (2 sqrt(2N)));
  - the Fisher matrix F = G^T G / (2 sqrt(2N)) (the tempered likelihood), its correlation matrix, and the widths of S and
    bar_frac marginal, with phi_b fixed, and with phi_b and |Omega_b| fixed;
  - the noise of a difference: a +0.8% |Omega_b| step (should be ~0 for a smooth model; its size is the noise).
Cells: t7 stars and OSARG stars in t7's own family; OSARG stars in Sormani+22; the steady-state mock in Hunter+24; OSARG
stars in the two hybrids (rotation curve of one family, bar of the other).

Usage: python fisher_gradients_by_family.py
Output: results/fisher_gradients_by_family.json and a printed table
"""
import json
import time

import numpy as np
import agama

from common import RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from potential_family import TIME_UNIT_GYR, ORBIT_KW
from rescore_empirical_selection import library_fit_empirical
from cross_response_stars_vs_potential import B
from library_noise_floor_test import load_stars
import cross_response_stars_vs_potential as cr

T_LONG, N_LONG = 4.0, 161
CELLS = [("t7cov", "nbody", 26.0), ("osarg", "nbody", 26.0), ("osarg", "sormani22", 32.5), ("h24mock", "mw", 32.5),
         ("osarg", "hyb_t7axi_mwbar", 26.0), ("osarg", "hyb_mwaxi_t7bar", 32.5)]
NAMES = ["Omega", "phi", "S", "f"]


def family_potential(family, f, S):
    if family in ("nbody", "mw", "sormani22"):
        return ssl.bar_aligned_potential(f, S, family)
    # hybrids at (bar_frac, S): rescale the non-axisymmetric part as the families do
    from potential_family import files as family_files
    from galaxies import galaxy as galaxy_info
    axi_gal, bar_gal = ("mw", "nbody") if family == "hyb_mwaxi_t7bar" else ("nbody", "mw")
    fixed = (agama.Potential(str(galaxy_info("mw")["pot_dir"] / "fixed.ini")) if axi_gal == "mw"
             else agama.Potential(str(galaxy_info("nbody")["pot_dir"] / "dm_spherical.ini")))
    Fa, Fb = family_files(axi_gal), family_files(bar_gal)
    rot = 0.0 if bar_gal == "mw" else np.radians(BAR_ANGLE_DEG)
    bary = agama.Potential(agama.Potential(file=str(Fa["baryon_axi"])),
                           agama.Potential(file=str(Fb["baryon_full"]), scale=[f * S, S], rotation=rot),
                           agama.Potential(file=str(Fb["baryon_axi"]), scale=[-f * S, S]))
    return agama.Potential(fixed, agama.Potential(type="CylSpline", potential=bary, mmax=8, symmetry="none", **ssl.NBODY_GRID))


def library(w, dat, sigma_dm, ok, family, om, phi, S, f):
    pot = family_potential(family, f, S)
    orbs = agama.orbit(potential=pot, ic=ssl.rotate_z(w, np.radians(phi)), time=T_LONG / TIME_UNIT_GYR, trajsize=N_LONG,
                       Omega=-abs(om), accuracy=old.ACCURACY, **ORBIT_KW)
    times = np.asarray(orbs[0, 0]).astype(np.float32) * TIME_UNIT_GYR
    traj = np.stack(orbs[:, 1]).astype(np.float32)
    cols = np.arange(0, N_LONG, 2)
    cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, cols]), phi, sigma_dm, dat, ok, times=times[cols])
    return (cl.reshape(fs.NB, -1) @ fs.rotation(BAR_ANGLE_DEG - phi).T).reshape(fs.NB, 6, fs.NK)   # data frame


def main():
    out = {}
    nr, spec = fs.SCHEMES[B]
    idx = [fi * fs.NK + k for fi, terms in spec.items() for k in terms]
    hart = np.sqrt((fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1))
    for stars, family, om0 in CELLS:
        t0 = time.time()
        w, dat, sigma_dm, phi0 = load_stars(stars)
        boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
        cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
        ssl.ACCURACY = old.ACCURACY
        ref = dict(om=om0, phi=phi0, S=1.0, f=1.0)
        steps = dict(om=0.10 * om0, phi=3.0, S=0.15, f=0.25)
        L = {"ref": library(w, dat, sigma_dm, ok, family, **ref)}
        for k, h in steps.items():
            for sgn in (+1, -1):
                p = dict(ref)
                p[k] += sgn * h
                L[f"{k}{'+' if sgn > 0 else '-'}"] = library(w, dat, sigma_dm, ok, family, **p)
        p = dict(ref)
        p["om"] *= 1.008
        L["noise"] = library(w, dat, sigma_dm, ok, family, **p)
        blocks = [b for b in range(nr * fs.NL) if ok[b] and all(np.all(np.isfinite(l[b].reshape(-1)[idx])) for l in L.values())]
        Wb = {b: hart * np.linalg.inv(np.linalg.cholesky(cov0[b][np.ix_(idx, idx)])) for b in blocks}
        vec = lambda c: np.concatenate([Wb[b] @ c[b].reshape(-1)[idx] for b in blocks])
        n = len(blocks) * len(idx)
        temp = 2 * np.sqrt(2 * n)
        G = np.column_stack([(vec(L[f"{k}+"]) - vec(L[f"{k}-"])) / (2 * h) for k, h in steps.items()])
        noise = vec(L["noise"]) - vec(L["ref"])
        F = G.T @ G / temp
        cos = (G.T @ G) / np.outer(np.linalg.norm(G, axis=0), np.linalg.norm(G, axis=0))
        Cov = np.linalg.inv(F)
        corr = Cov / np.sqrt(np.outer(np.diag(Cov), np.diag(Cov)))
        sd = np.sqrt(np.diag(Cov))
        sub = lambda keep: np.sqrt(np.diag(np.linalg.inv(F[np.ix_(keep, keep)])))   # others fixed
        sS_phi, sf_phi = sub([0, 2, 3])[1], sub([0, 2, 3])[2]
        sS_phiom, sf_phiom = sub([2, 3])
        # change of the prediction for each full step (tempered), and the noise of one difference (tempered)
        Dstep = {k: float(np.sum((vec(L[f"{k}+"]) - vec(L[f"{k}-"])) ** 2) / 4 / temp) for k in steps}
        key = f"{stars} in {family}"
        out[key] = dict(ref=ref, steps=steps, N=n, cos=cos.tolist(), corr=corr.tolist(), sd=sd.tolist(),
                        sd_S_phi_fixed=float(sS_phi), sd_f_phi_fixed=float(sf_phi), sd_S_phi_om_fixed=float(sS_phiom),
                        sd_f_phi_om_fixed=float(sf_phiom), D_half_step=Dstep, noise_D=float(noise @ noise / temp))
        print(f"== {key} (N {n}; {time.time() - t0:.0f} s); change of the prediction for one step, tempered: "
              + ", ".join(f"{k} {v:.2f}" for k, v in Dstep.items()) + f"; noise of a difference {noise @ noise / temp:.2f}")
        print("   cos between the gradients: " + ", ".join(f"{NAMES[i]}-{NAMES[j]} {cos[i, j]:+.2f}" for i in range(4) for j in range(i + 1, 4)))
        print("   Fisher correlations:       " + ", ".join(f"{NAMES[i]}-{NAMES[j]} {corr[i, j]:+.2f}" for i in range(4) for j in range(i + 1, 4)))
        print(f"   Fisher sd: Omega {sd[0]:.2f} ({100 * sd[0] / om0:.1f}%), phi {sd[1]:.2f} deg, S {sd[2]:.3f}, f {sd[3]:.3f}; "
              f"with phi fixed: S {sS_phi:.3f}, f {sf_phi:.3f}; with phi and Omega fixed: S {sS_phiom:.3f}, f {sf_phiom:.3f}", flush=True)
        json.dump(out, open(RESULTS_DIR / "fisher_gradients_by_family.json", "w"), indent=2)


if __name__ == "__main__":
    main()

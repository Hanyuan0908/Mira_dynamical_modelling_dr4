"""Stars against potential: which one makes the predicted kinematics react weakly to the pattern speed?

Three star sets, all with OSARG's star count and observed spatial coverage:
  t7cov    the t = 7 N-body mock (12.5% distance error), thinned to OSARG's coverage and N (mask osarg_cover)
  osarg    the real OSARG data (Zhang et al. 2024; 9.6%)
  h24mock  the OSARG stars phase-mixed in the Hunter+24 family at |Omega_b| 32.5, phi_b 25, S = bar_frac = 1 (9.6%;
           make_mw_hunter24_mock_from_osarg.py)
are each integrated in two families at their own bar (S = bar_frac = 1):
  nbody    the t = 7 snapshot's family, around its own pattern speed (Omega_ref 26.0)
  mw       the Hunter+24 family, around Omega_ref 32.5
Both references put corotation at about 1.6 bar lengths (compare_bar_strength_and_resonances.py: 1.65 and 1.59).
Two hybrid families separate the rotation curve (and so the resonances) from the bar:
  hyb_mwaxi_t7bar   Hunter+24's axisymmetric potential + t = 7's non-axisymmetric part (Omega_ref 32.5)
  hyb_t7axi_mwbar   t = 7's axisymmetric potential + Hunter+24's (Portail+17) non-axisymmetric part (Omega_ref 26.0)
The non-axisymmetric part of a family is baryon_full - baryon_axi, exactly what the family rescales with bar_frac.
  sormani22         the Sormani+22 Milky Way family (build_sormani22_potential.py), Omega_ref 32.5 like Hunter+24
|Omega_b| is varied in FRACTIONAL steps of Omega_ref with everything else fixed. phi_b is the star set's own
(27 deg for t7cov, 25 deg for the other two). Scoring is the pipeline's (1 Gyr, 50 Myr, sfratio with K = 10 against
the star set's own data, scheme B). Recorded per step, tempered by sqrt(2 N):
  D      chi2(library(Omega) vs library(Omega_ref) | the star set's data covariance): how much the PREDICTION changes
  data   chi2(data vs library(Omega)) - chi2(data vs library(Omega_ref)): how much the DATA notice
A step of +-0.8% (about 0.25 km/s/kpc) measures the library noise floor of each cell.

Usage: python cross_response_stars_vs_potential.py [--cells t7cov:nbody osarg:mw ...]
Output: results/cross_response_<stars>_<family>.npz and a printed table
"""
import argparse
import time

import numpy as np
import h5py

from common import ROOT, MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
from rescore_empirical_selection import data_lbdR, library_fit_empirical

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
STARS = {
    "t7cov": dict(mock="nbody_mock_selection_sigmaDM0.271_N50k.hdf5",
                  keep="results/rv_selection_masks/nbody_mock_selection_sigmaDM0.271_N50k__osarg_cover.npy", phi=27.0),
    "osarg": dict(mock="osarg_zhang2024_observed_R8.hdf5", keep=None, phi=25.0),
    "h24mock": dict(mock="mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5", keep=None, phi=25.0),
}
FAMILY = {"nbody": 26.0, "mw": 32.5, "sormani22": 32.5, "hyb_mwaxi_t7bar": 32.5, "hyb_t7axi_mwbar": 26.0}
# Rotation-curve swaps under the Sormani+22 bar, with the simulation's rotation curve RESCALED to the Milky Way's
# magnitude (user, 2026-11-25): axisymmetric part = the snapshot's spherical halo + its baryons re-expanded here with
# symmetry axisymmetric (the cluster-built *_axi files of t = 2, 3, 4, 5.5 are not axisymmetric), both multiplied by
# k^2 so that the mean v_c over 2-6 kpc equals Sormani+22's; bar = Sormani+22's baryon_full - baryon_axi.
# Omega_ref 32.5 as for the sormani22 cell.
SNAP_HYB = {"hyb_t2axi_s22bar": "nbody_t2", "hyb_t5.5axi_s22bar": "nbody_t5.5", "hyb_t7axi_s22bar": "nbody"}
FAMILY.update({k: 32.5 for k in SNAP_HYB})
# Sormani+22 with a compact central mass (Hernquist, a = 0.3 kpc) added to its axisymmetric part, mass chosen so that
# v_c(1 kpc) / <v_c>(2-6 kpc) reaches the target (S22: 0.91; the rescaled N-body snapshots: ~1.11), then everything
# axisymmetric rescaled by k^2 to keep <v_c>(2-6 kpc) at Sormani+22's; the bar is Sormani+22's, unchanged.
CORE_HYB = {"s22_core_r1.00": (0.3, 1.0, 1.00), "s22_core_r1.11": (0.3, 1.0, 1.11),   # (scale radius, R_match, target)
            "s22_core_a1_R2_r1.09": (1.0, 2.0, 1.09)}   # extended: v_c(2 kpc)/<v_c>(2-6) as the rescaled t = 2 (1.092)
FAMILY.update({k: 32.5 for k in CORE_HYB})


def _vc_mean(pot, R=np.linspace(2, 6, 41), nphi=36):
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    X = np.column_stack([(R[:, None] * np.cos(ph)).ravel(), (R[:, None] * np.sin(ph)).ravel(), np.zeros(R.size * nphi)])
    F = pot.force(X)
    FR = (F[:, 0] * np.cos(np.tile(ph, R.size)) + F[:, 1] * np.sin(np.tile(ph, R.size))).reshape(R.size, -1).mean(1)
    return float(np.mean(np.sqrt(-R * FR)))


def snapshot_hybrid(family):
    """Snapshot rotation curve (rescaled to Sormani+22's mean v_c over 2-6 kpc) + the Sormani+22 bar."""
    import agama
    import os
    import tempfile
    from potential_family import files as family_files
    from galaxies import galaxy as galaxy_info
    pd = galaxy_info(SNAP_HYB[family])["pot_dir"]
    Ff, Fs = family_files(SNAP_HYB[family]), family_files("sormani22")
    tmp = tempfile.mkdtemp(prefix="hyb_")
    axi_file = os.path.join(tmp, "baryon_axi_rebuilt.ini")
    agama.Potential(type="CylSpline", potential=agama.Potential(file=str(Ff["baryon_full"])), symmetry="axisymmetric",
                    mmax=0, **ssl.NBODY_GRID).export(axi_file)
    s22_axi = agama.Potential(agama.Potential(str(galaxy_info("sormani22")["pot_dir"] / "fixed.ini")),
                              agama.Potential(file=str(Fs["baryon_axi"])))
    raw = agama.Potential(agama.Potential(str(pd / "dm_spherical.ini")), agama.Potential(file=axi_file))
    k2 = (_vc_mean(s22_axi) / _vc_mean(raw)) ** 2
    axi = agama.Potential(agama.Potential(file=str(pd / "dm_spherical.ini"), scale=[k2, 1.0]),
                          agama.Potential(file=axi_file, scale=[k2, 1.0]))
    bar = agama.Potential(agama.Potential(file=str(Fs["baryon_full"])), agama.Potential(file=str(Fs["baryon_axi"]), scale=[-1.0, 1.0]))
    pot = agama.Potential(agama.Potential(file=str(pd / "dm_spherical.ini"), scale=[k2, 1.0]),
                          agama.Potential(type="CylSpline", mmax=8, symmetry="none", **ssl.NBODY_GRID,
                                          potential=agama.Potential(agama.Potential(file=axi_file, scale=[k2, 1.0]), bar)))
    ph = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    def m2(p, R):
        X = np.column_stack([R * np.cos(ph), R * np.sin(ph), np.zeros_like(ph)])
        F = p.force(X); FR = F[:, 0] * np.cos(ph) + F[:, 1] * np.sin(ph)
        return 2 * abs(np.mean(FR * np.exp(-2j * ph)))
    s22 = ssl.bar_aligned_potential(1.0, 1.0, "sormani22")
    print(f"   hybrid {family}: k^2 = {k2:.4f}; mean v_c 2-6 kpc: raw {_vc_mean(raw):.1f}, rescaled {_vc_mean(axi):.1f}, "
          f"hybrid {_vc_mean(pot):.1f}, Sormani+22 {_vc_mean(s22_axi):.1f}; bar m=2 force hybrid / Sormani+22 at R 2,3,4,5: "
          + " ".join(f"{m2(pot, R) / m2(s22, R):.3f}" for R in (2, 3, 4, 5)), flush=True)
    return pot


def core_hybrid(family):
    """Sormani+22 + a compact central mass, rescaled to Sormani+22's <v_c>(2-6 kpc); bar unchanged."""
    import agama
    from scipy.optimize import brentq
    from potential_family import files as family_files
    from galaxies import galaxy as galaxy_info
    Fs = family_files("sormani22")
    fixed, baxi = str(galaxy_info("sormani22")["pot_dir"] / "fixed.ini"), str(Fs["baryon_axi"])
    s22_axi = agama.Potential(agama.Potential(fixed), agama.Potential(file=baxi))
    a_core, r_match, target = CORE_HYB[family]
    def ratio(M, k2=1.0):
        p = agama.Potential(agama.Potential(file=fixed, scale=[k2, 1.0]), agama.Potential(file=baxi, scale=[k2, 1.0]),
                            agama.Potential(type="Dehnen", gamma=1.0, scaleRadius=a_core, mass=M * k2))
        return float(np.sqrt(-r_match * p.force([[r_match, 0, 0]])[0][0])) / _vc_mean(p), p
    M = brentq(lambda m: ratio(m)[0] - target, 0.0, 1e12, xtol=1e6)
    k2 = (_vc_mean(s22_axi) / _vc_mean(ratio(M)[1])) ** 2
    r, axi = ratio(M, k2)
    bar = agama.Potential(agama.Potential(file=str(Fs["baryon_full"])), agama.Potential(file=baxi, scale=[-1.0, 1.0]))
    pot = agama.Potential(agama.Potential(file=fixed, scale=[k2, 1.0]),
                          agama.Potential(type="Dehnen", gamma=1.0, scaleRadius=a_core, mass=M * k2),
                          agama.Potential(type="CylSpline", mmax=8, symmetry="none", **ssl.NBODY_GRID,
                                          potential=agama.Potential(agama.Potential(file=baxi, scale=[k2, 1.0]), bar)))
    vc = lambda p, R: float(np.sqrt(-R * p.force([[R, 0, 0]])[0][0]))
    print(f"   {family}: central mass {M:.3g} Msun, scale {a_core} kpc (x k^2 {k2:.4f} = {M * k2:.3g}); v_c({r_match:g})/<v_c>(2-6) {r:.3f}; "
          f"v_c at 1/2/5 kpc {vc(axi, 1):.0f} {vc(axi, 2):.0f} {vc(axi, 5):.0f} (S22 {vc(s22_axi, 1):.0f} {vc(s22_axi, 2):.0f} "
          f"{vc(s22_axi, 5):.0f}); <v_c>(2-6) {_vc_mean(pot):.1f} (S22 {_vc_mean(s22_axi):.1f})", flush=True)
    return pot


def family_potential(family):
    """The family at S = bar_frac = 1, bar along x; the hybrids swap the axisymmetric part and the bar."""
    if family in CORE_HYB:
        return core_hybrid(family)
    if family in ("nbody", "mw", "sormani22"):
        return ssl.bar_aligned_potential(1.0, 1.0, family)
    if family in SNAP_HYB:
        return snapshot_hybrid(family)
    import agama
    from potential_family import files as family_files
    from galaxies import galaxy as galaxy_info
    axi_gal, bar_gal = ("mw", "nbody") if family == "hyb_mwaxi_t7bar" else ("nbody", "mw")
    fixed = (agama.Potential(str(galaxy_info("mw")["pot_dir"] / "fixed.ini")) if axi_gal == "mw"
             else agama.Potential(str(galaxy_info("nbody")["pot_dir"] / "dm_spherical.ini")))
    Fa, Fb = family_files(axi_gal), family_files(bar_gal)
    rot = 0.0 if bar_gal == "mw" else np.radians(BAR_ANGLE_DEG)
    bary = agama.Potential(agama.Potential(file=str(Fa["baryon_axi"])),
                           agama.Potential(file=str(Fb["baryon_full"]), rotation=rot),
                           agama.Potential(file=str(Fb["baryon_axi"]), scale=[-1.0, 1.0]))
    return agama.Potential(fixed, agama.Potential(type="CylSpline", potential=bary, mmax=8, symmetry="none",
                                                  **ssl.NBODY_GRID))
FRACS = np.array([-0.3, -0.2, -0.1, -0.05, -0.008, 0.0, 0.008, 0.05, 0.1, 0.2, 0.3, 0.4])
FIELDS = {0: "<v_R>", 1: "<v_phi>", 3: "sig_R^2", 4: "sig_phi^2"}


def field_chi2(a, b, cov, f, valid):
    nr, spec = fs.SCHEMES[B]
    idx = [f * fs.NK + k for k in spec[f]]
    tot = 0.0
    for r in range(nr):
        for z in range(fs.NL):
            i = r * fs.NL + z
            if not valid[i]:
                continue
            d, C = (a[i] - b[i]).reshape(-1)[idx], cov[i][np.ix_(idx, idx)]
            tot += (fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1) * float(d @ np.linalg.solve(C, d))
    return tot


def run_cell(stars, family):
    st = STARS[stars]
    mock = MOCK_DIR / st["mock"]
    with h5py.File(mock, "r") as f:
        w = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
        sigma_dm = float(f.attrs["sigma_DM"])
    dat = data_lbdR(mock)
    if st["keep"]:
        keep = np.load(ROOT / st["keep"])
        w, dat = w[keep], tuple(np.asarray(a)[keep] for a in dat)
    phi = st["phi"]
    boot = np.random.default_rng(fs.SEED).multinomial(len(w), np.full(len(w), 1 / len(w)), size=fs.N_BOOT).astype(float)
    cd0, cov0, _, ok = fs.data_fit(w, BAR_ANGLE_DEG, boot)
    T = fs.rotation(phi - BAR_ANGLE_DEG)
    cd = (cd0.reshape(fs.NB, -1) @ T.T).reshape(fs.NB, 6, fs.NK)
    cov = np.einsum("ij,bjk,lk->bil", T, cov0, T)
    ssl.ACCURACY = old.ACCURACY
    pot = family_potential(family)
    ic = ssl.rotate_z(w, np.radians(phi))
    om_ref = FAMILY[family]
    omegas = om_ref * (1 + FRACS)
    libs, chi2 = [], []
    print(f"{stars} in {family}: {len(w):,} stars, sigma_DM {sigma_dm:.4f}, phi_b {phi}, Omega_ref {om_ref}", flush=True)
    for om in omegas:
        t0 = time.time()
        _, traj = ssl.integrate(pot, ic, -abs(om))
        cl, _ = library_fit_empirical(np.ascontiguousarray(traj[:, ::2]), phi, sigma_dm, dat, ok)
        cl = cl.reshape(fs.NB, 6, fs.NK)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        x2, nd = fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[B])
        libs.append(cl)
        chi2.append(x2)
        print(f"   Omega {om:6.2f}: chi2/N {x2 / nd:.3f}  ({time.time() - t0:.0f} s)", flush=True)
    libs, chi2 = np.array(libs), np.array(chi2)
    iref = int(np.where(FRACS == 0)[0][0])
    D, Df = [], []
    for cl in libs:
        valid = ok & np.all(np.isfinite(cl), axis=(1, 2)) & np.all(np.isfinite(libs[iref]), axis=(1, 2))
        D.append(fs.scheme_chi2(libs[iref], cl, cov, valid, fs.SCHEMES[B])[0])
        Df.append([field_chi2(libs[iref], cl, cov, f, valid) for f in FIELDS])
    np.savez(RESULTS_DIR / f"cross_response_{stars}_{family}.npz", fracs=FRACS, omegas=omegas, chi2=chi2, D=np.array(D),
             D_field=np.array(Df), fields=list(FIELDS.values()), ndata=nd, coef_lib=libs, coef_data=cd, cov_data=cov,
             phi=phi, sigma_dm=sigma_dm, n_stars=len(w))
    temp = 2 * np.sqrt(2 * nd)
    print(f"   {stars} in {family}: N_data {nd}; tempered by {temp:.1f}")
    print("   frac     " + " ".join(f"{x:+6.3f}" for x in FRACS))
    print("   model D  " + " ".join(f"{x / temp:6.2f}" for x in D))
    print("   data     " + " ".join(f"{(x - chi2[iref]) / temp:6.2f}" for x in chi2), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", nargs="+", default=[f"{s}:{f}" for s in STARS for f in FAMILY])
    args = ap.parse_args()
    for c in args.cells:
        run_cell(*c.split(":"))


if __name__ == "__main__":
    main()

"""Steady-state likelihood of a stellar sample in a trial barred potential.

log L(theta) = sum_i log p_theta(v_i | x_i), where p_theta(v|x) is estimated from the time-sampled
orbits of the sample itself, integrated in the trial potential (a discrete Schwarzschild model with
the stars' own orbits as the library and equal population weights):

  - frame: the potential's bar lies along the x-axis; a trial bar angle phi_b rotates the stars
    (MW-like frame, bar at -phi_b) by +phi_b so that their assumed bar lies along x;
  - orbits are integrated in the bar frame, agama.orbit(..., Omega=Omega_b), Omega_b < 0 (clockwise);
  - the library keeps samples with t >= T_MIN and is folded four-fold with the symmetries of motion
    in a rotating, even-m, z-symmetric potential: (x,y,z) -> (-x,-y,z), (x,y,-z) and both
    (velocities transformed alike);
  - selection: each orbit is weighted by 1 / tau_i, tau_i being the fraction of time it spends
    inside the selected region (here cylindrical R < R_SEL);
  - leave-one-orbit-out: star i's own orbit (all four images) is excluded from p(v_i | x_i);
  - p(v|x): conditional KDE over the K_NN nearest library points in position, Gaussian spatial
    weights with bandwidth r_K / 2, and a velocity kernel whose covariance is the local weighted
    velocity covariance times n_eff^(-2/7) (Scott's rule in 3D).
Also available (evaluate(..., extended=True)): a per-star local-Gaussian p(v|x) = N(v | mu_i, C_i) from
the same neighbours (mu_i, C_i saved for post-hoc variants), and a moment chi^2 comparing the
observed and library means / dispersions in bar-frame cells (R, phi mod pi, |z|) with errors
sigma/sqrt(N) and sigma/sqrt(2N). Trial potentials: for both galaxies the static barred baryons are
re-expanded into one mmax-8 CylSpline (au18: disc grid; nbody: NBODY_GRID, finer at the centre for
the bulge cusp: bar-force error 0.2% rms, 1.6x faster than the exact composite); DM (+BH) kept as is.
"""
import time

import numpy as np
import agama
from scipy.spatial import cKDTree

from common import BAR_ANGLE_DEG
from potential_family import make_static_fast, TIME_UNIT_GYR, POT_DIR, files as family_files
from galaxies import galaxy as galaxy_info

T_GYR = 1.0           # integration time
N_SNAP = 41           # snapshots over [0, T] (every 25 Myr)
T_MIN_GYR = 0.25      # library uses t >= T_MIN (after initial phase mixing)
R_SEL = 8.0           # Test 1 selection: cylindrical R < 8 kpc
K_NN = 1000           # neighbours per star for p(v|x)
ACCURACY = 1e-6
CHUNK = 4000

CELL_R = np.arange(0.0, 8.0 + 1e-9, 0.5)            # moment cells (bar frame, folded)
CELL_PHI = np.radians(np.arange(0, 181, 15))
CELL_Z = np.array([0.0, 0.5, 1.0, 2.0, np.inf])
BOX_XY = np.arange(-8.0, 8.0 + 1e-9, 0.5)            # library moment boxes (bar frame, Cartesian)
BOX_Z = np.arange(-4.0, 4.0 + 1e-9, 0.5)
PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))   # (v_R, v_phi, v_z) second-moment pairs
NBODY_GRID = dict(gridSizeR=36, gridSizez=36, Rmin=0.03, Rmax=50.0, zmin=0.02, zmax=50.0)
_CACHE = {}


def bar_aligned_potential(bar_frac, S, galaxy="au18"):
    """Trial potential (static, bar along x)."""
    if galaxy not in _CACHE:
        pd = galaxy_info(galaxy)["pot_dir"]
        _CACHE[galaxy] = (agama.Potential(agama.Potential(str(POT_DIR / "dm.ini")),
                                          agama.Potential(str(POT_DIR / "bh.ini"))) if galaxy == "au18"
                          else agama.Potential(str(pd / "dm_spherical.ini")))
    if galaxy == "au18":     # merged mmax-8 CylSpline for the baryons
        return make_static_fast(bar_frac, S, phi0=np.radians(BAR_ANGLE_DEG), mmax=8, symmetry="none",
                                pot_dm_bh=_CACHE[galaxy])
    F = family_files(galaxy)
    bary = agama.Potential(agama.Potential(file=str(F["baryon_axi"])),
                           agama.Potential(file=str(F["baryon_full"]), scale=[bar_frac * S, S],
                                           rotation=np.radians(BAR_ANGLE_DEG)),
                           agama.Potential(file=str(F["baryon_axi"]), scale=[-bar_frac * S, S]))
    return agama.Potential(_CACHE[galaxy], agama.Potential(type="CylSpline", potential=bary, mmax=8,
                                                           symmetry="none", **NBODY_GRID))


def cell_moments(w, weight=None):
    """Weighted N, means and dispersions of (v_R, v_phi, v_z sign z) in folded bar-frame cells
    (R, phi mod pi, |z|); arrays of shape (n_cells,) and (n_cells, 3). Dispersions use the unbiased
    (reliability-weights) variance, W^2 / (W^2 - sum w^2) times the 1/W one, i.e. N / (N - 1) for unit weights."""
    weight = np.ones(len(w)) if weight is None else weight
    R = np.hypot(w[:, 0], w[:, 1])
    phi = np.mod(np.arctan2(w[:, 1], w[:, 0]), np.pi)
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5] * np.sign(w[:, 2])])
    i, j, k = (np.digitize(R, CELL_R) - 1, np.digitize(phi, CELL_PHI) - 1, np.digitize(np.abs(w[:, 2]), CELL_Z) - 1)
    ok = (i >= 0) & (i < len(CELL_R) - 1) & (j >= 0) & (j < len(CELL_PHI) - 1)
    shape = (len(CELL_R) - 1, len(CELL_PHI) - 1, len(CELL_Z) - 1)
    flat = np.ravel_multi_index((i[ok], j[ok], k[ok]), shape)
    n = int(np.prod(shape))
    W = np.bincount(flat, weights=weight[ok], minlength=n)
    W2 = np.bincount(flat, weights=weight[ok]**2, minlength=n)
    N = np.bincount(flat, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        m1 = np.column_stack([np.bincount(flat, weights=weight[ok] * v[ok, a], minlength=n) for a in range(3)]) / W[:, None]
        m2 = np.column_stack([np.bincount(flat, weights=weight[ok] * v[ok, a]**2, minlength=n) for a in range(3)]) / W[:, None]
        unbiased = W**2 / (W**2 - W2)
    return N, m1, np.sqrt(np.maximum((m2 - m1**2) * unbiased[:, None], 0))


def box_sums(w, weight=None):
    """Moment sums in 0.5 kpc Cartesian boxes (BOX_XY x BOX_XY x BOX_Z): raw count N, weighted mass S0,
    S1 = sum w v and S2 = sum w v_a v_b for v = (v_R, v_phi, v_z) (v_phi > 0 along rotation; pairs in
    PAIRS order RR, pp, zz, Rp, Rz, pz). Sums, so they can be rebinned exactly; means = S1/S0,
    covariances = S2/S0 - mean_a mean_b."""
    weight = np.ones(len(w)) if weight is None else weight
    R = np.hypot(w[:, 0], w[:, 1])
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                         -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, w[:, 5]])
    ix, iy, iz = (np.digitize(w[:, 0], BOX_XY) - 1, np.digitize(w[:, 1], BOX_XY) - 1, np.digitize(w[:, 2], BOX_Z) - 1)
    nx, nz = len(BOX_XY) - 1, len(BOX_Z) - 1
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < nx) & (iz >= 0) & (iz < nz)
    flat = np.ravel_multi_index((ix[ok], iy[ok], iz[ok]), (nx, nx, nz))
    n, wt, vv = nx * nx * nz, weight[ok], v[ok]
    out = dict(N=np.bincount(flat, minlength=n), S0=np.bincount(flat, weights=wt, minlength=n),
               S1=np.column_stack([np.bincount(flat, weights=wt * vv[:, a], minlength=n) for a in range(3)]),
               S2=np.column_stack([np.bincount(flat, weights=wt * vv[:, a] * vv[:, b], minlength=n) for a, b in PAIRS]))
    return {k: a.reshape((nx, nx, nz) + a.shape[1:]) for k, a in out.items()}


def moment_chi2(obs, lib, rmin=0.0, disp=False, n_min=20):
    """chi^2 of library vs observed cell moments (observational errors only)."""
    N, mo, so = obs
    _, ml, sl = lib
    Rc = np.repeat(0.5 * (CELL_R[1:] + CELL_R[:-1]), (len(CELL_PHI) - 1) * (len(CELL_Z) - 1))
    good = (N >= n_min) & (Rc >= rmin) & np.all(np.isfinite(ml), 1) & np.all(so > 0, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.sum((mo - ml)[good]**2 / (so[good]**2 / N[good, None]))
        if disp:
            c += np.sum((so - sl)[good]**2 / (so[good]**2 / (2 * N[good, None])))
    return float(c), int(good.sum())


def rotate_z(w, angle):
    c, s = np.cos(angle), np.sin(angle)
    out = np.array(w, dtype=np.float64, copy=True)
    for i in (0, 3):
        out[:, i], out[:, i + 1] = c * w[:, i] - s * w[:, i + 1], s * w[:, i] + c * w[:, i + 1]
    return out


def integrate(pot, ic, Omega_b):
    T = T_GYR / TIME_UNIT_GYR
    orbs = agama.orbit(potential=pot, ic=ic, time=T, trajsize=N_SNAP, Omega=Omega_b, accuracy=ACCURACY,
                       verbose=False)
    times = orbs[0, 0] * TIME_UNIT_GYR
    traj = np.stack(orbs[:, 1]).astype(np.float32)          # (N, N_SNAP, 6), bar frame
    return times, traj


def build_library(times, traj):
    """Folded library (positions, velocities, orbit id, weight) and tau per orbit."""
    n = traj.shape[0]
    R = np.hypot(traj[:, 1:, 0], traj[:, 1:, 1])            # t > 0
    tau = np.clip(np.mean(R < R_SEL, axis=1), 1.0 / (traj.shape[1] - 1), 1.0)
    keep = times >= T_MIN_GYR - 1e-9
    w = traj[:, keep].reshape(-1, 6).astype(np.float64)
    oid = np.repeat(np.arange(n), keep.sum())
    images = [w]
    for sx, sz in ((-1, 1), (1, -1), (-1, -1)):              # rotation by pi, z-mirror, both
        m = w.copy()
        m[:, [0, 1, 3, 4]] *= sx
        m[:, [2, 5]] *= sz
        images.append(m)
    lib = np.concatenate(images)
    lib_oid = np.tile(oid, 4)
    lib_w = (1.0 / tau)[lib_oid]
    return lib, lib_oid, lib_w, tau


def log_p_v_given_x(stars, lib, lib_oid, lib_w, k=K_NN, return_local=False):
    """Leave-one-orbit-out conditional KDE log p(v_i | x_i) for every star (star i = orbit i).
    With return_local: also the local-Gaussian log N(v_i | mu_i, C_i) and mu_i, C_i."""
    tree = cKDTree(lib[:, :3])
    out = np.empty(len(stars))
    gauss = np.empty(len(stars))
    mus = np.empty((len(stars), 3))
    covs = np.empty((len(stars), 3, 3))
    for s in range(0, len(stars), CHUNK):
        q = stars[s:s + CHUNK]
        d, idx = tree.query(q[:, :3], k=k + 150, workers=-1)
        own = lib_oid[idx] == np.arange(s, s + len(q))[:, None]
        d = np.where(own, np.inf, d)
        order = np.argsort(d, axis=1)[:, :k]                     # k nearest non-own neighbours
        d = np.take_along_axis(d, order, 1)
        idx = np.take_along_axis(idx, order, 1)
        hx = d[:, -1:] / 2
        W = np.exp(-0.5 * (d / hx)**2) * lib_w[idx]
        W /= W.sum(1, keepdims=True)
        V = lib[idx, 3:]
        mu = np.einsum("nk,nkj->nj", W, V)
        dV = V - mu[:, None]
        C = np.einsum("nk,nki,nkj->nij", W, dV, dV)
        neff = 1.0 / np.sum(W**2, 1)
        H = C * (neff ** (-2.0 / 7.0))[:, None, None]
        Hinv = np.linalg.inv(H)
        _, logdet = np.linalg.slogdet(H)
        u = V - q[:, None, 3:]
        maha = np.einsum("nki,nij,nkj->nk", u, Hinv, u)
        logk = -0.5 * maha - 0.5 * (3 * np.log(2 * np.pi) + logdet)[:, None]
        mx = logk.max(1, keepdims=True)
        out[s:s + len(q)] = mx[:, 0] + np.log(np.sum(W * np.exp(logk - mx), 1))
        if return_local:
            r = q[:, 3:] - mu
            _, logdetC = np.linalg.slogdet(C)
            gauss[s:s + len(q)] = -0.5 * np.einsum("ni,nij,nj->n", r, np.linalg.inv(C), r) \
                - 0.5 * (3 * np.log(2 * np.pi) + logdetC)
            mus[s:s + len(q)], covs[s:s + len(q)] = mu, C
    return (out, gauss, mus, covs) if return_local else out


def evaluate(stars_mw, Omega_abs, phi_b_deg, S, bar_frac, galaxy="au18", extended=False, keep_traj=False):
    """log L for one trial parameter set; stars_mw: (N, 6) in the MW-like frame.
    Returns (summary, per-star KDE log p), or with extended=True (summary, arrays) where the summary
    also has the local-Gaussian logL and moment chi^2 values, and arrays holds per-star KDE and
    Gaussian log p, mu, C, and the observed and library cell moments."""
    t0 = time.time()
    pot = bar_aligned_potential(bar_frac, S, galaxy)
    stars = rotate_z(stars_mw, np.radians(phi_b_deg))         # assumed bar -> x-axis
    times, traj = integrate(pot, stars, -abs(Omega_abs))
    t1 = time.time()
    lib, lib_oid, lib_w, tau = build_library(times, traj)
    res = log_p_v_given_x(stars, lib, lib_oid, lib_w, return_local=extended)
    logp = res[0] if extended else res
    summary = dict(logL=float(np.sum(logp)), logp_median=float(np.median(logp)),
                   tau_median=float(np.median(tau)), tau_min=float(tau.min()))
    if not extended:
        t2 = time.time()
        summary.update(t_integrate_s=round(t1 - t0, 1), t_likelihood_s=round(t2 - t1, 1))
        return summary, logp
    _, gauss, mu, C = res
    obs_m, lib_m = cell_moments(stars), cell_moments(lib, lib_w)
    for name, rmin, disp in (("chi2_means", 0.0, False), ("chi2_means_R1", 1.0, False),
                             ("chi2_meansdisp_R1", 1.0, True)):
        summary[name], summary[name + "_ncells"] = moment_chi2(obs_m, lib_m, rmin, disp)
    summary["logL_gauss"] = float(np.sum(gauss))
    t2 = time.time()
    summary.update(t_integrate_s=round(t1 - t0, 1), t_likelihood_s=round(t2 - t1, 1))
    arrays = dict(logp_kde=logp.astype(np.float32), logp_gauss=gauss.astype(np.float32),
                  mu=mu.astype(np.float32), C=C.astype(np.float32),
                  obs_N=obs_m[0], obs_m1=obs_m[1], obs_s=obs_m[2],
                  lib_W=lib_m[0], lib_m1=lib_m[1], lib_s=lib_m[2],
                  tau=tau.astype(np.float32),
                  **{f"box_{k}": v for k, v in box_sums(lib, lib_w).items()})
    if keep_traj:
        arrays.update(times=times, traj=traj)
    return summary, arrays

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
"""
import time

import numpy as np
import agama
from scipy.spatial import cKDTree

from common import BAR_ANGLE_DEG
from potential_family import make_static_fast, TIME_UNIT_GYR, POT_DIR

T_GYR = 1.0           # integration time
N_SNAP = 41           # snapshots over [0, T] (every 25 Myr)
T_MIN_GYR = 0.25      # library uses t >= T_MIN (after initial phase mixing)
R_SEL = 8.0           # Test 1 selection: cylindrical R < 8 kpc
K_NN = 1000           # neighbours per star for p(v|x)
ACCURACY = 1e-6
CHUNK = 4000

_DM_BH = None


def bar_aligned_potential(bar_frac, S):
    """Trial potential (static, bar along x) with merged mmax-8 CylSpline for the baryons."""
    global _DM_BH
    if _DM_BH is None:
        _DM_BH = agama.Potential(agama.Potential(str(POT_DIR / "dm.ini")),
                                 agama.Potential(str(POT_DIR / "bh.ini")))
    return make_static_fast(bar_frac, S, phi0=np.radians(BAR_ANGLE_DEG), mmax=8, symmetry="none",
                            pot_dm_bh=_DM_BH)


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


def log_p_v_given_x(stars, lib, lib_oid, lib_w, k=K_NN):
    """Leave-one-orbit-out conditional KDE log p(v_i | x_i) for every star (star i = orbit i)."""
    tree = cKDTree(lib[:, :3])
    out = np.empty(len(stars))
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
    return out


def evaluate(stars_mw, Omega_abs, phi_b_deg, S, bar_frac):
    """log L for one trial parameter set; stars_mw: (N, 6) in the MW-like frame."""
    t0 = time.time()
    pot = bar_aligned_potential(bar_frac, S)
    stars = rotate_z(stars_mw, np.radians(phi_b_deg))         # assumed bar -> x-axis
    times, traj = integrate(pot, stars, -abs(Omega_abs))
    t1 = time.time()
    lib, lib_oid, lib_w, tau = build_library(times, traj)
    logp = log_p_v_given_x(stars, lib, lib_oid, lib_w)
    t2 = time.time()
    return dict(logL=float(np.sum(logp)), logp_median=float(np.median(logp)),
                tau_median=float(np.median(tau)), tau_min=float(tau.min()),
                t_integrate_s=round(t1 - t0, 1), t_likelihood_s=round(t2 - t1, 1)), logp

"""Selection function of a mock, and the library weighting it implies.

`make_mock_selection.py` stores the fitted selection function in the mock file as three tables, so that
S can be evaluated anywhere:

    S(l, b, d) = S_sky[HEALPix(l, b)] * S_dist[d bin] * S_R[R bin],   R = R(l, b, d),

with S = 0 outside the observed R < R_max of the mock, and each particle kept with probability min(S, 1).
S is defined on **observed** quantities, so to evaluate it on an orbit library the library has to be observed
the same way the data were: rotated back to the MW frame, converted to (l, b, d), smeared in distance modulus
with the mock's sigma_DM, and converted back.

Library weighting (the plan in PROGRESS, "Next steps"), generalising the R < R_SEL indicator that
`steady_state_likelihood.build_library` and `test1_nbody_qmc_focused_fourier_stars.library_fit` use:

    weight(sample j of orbit i) = S(x_j) / tau_i,     tau_i = <S(x(t))> along orbit i.

S(x_j) makes the library follow the same selection as the data; 1/tau_i undoes the bias from having drawn the
orbits themselves from a selected sample. With S replaced by the indicator R < R_SEL this is exactly the
current weighting, so the two are the same scheme with different S.
"""
import numpy as np
import h5py
import healpy as hp

import observe
import steady_state_likelihood as ssl

TAU_FLOOR = 1e-3          # orbits that essentially never fall in the selection get this, not 1/0


def load_tables(mock_path):
    """The mock's selection tables, or None if it has no selection function."""
    with h5py.File(mock_path, "r") as f:
        if "S_sky" not in f:
            return None
        return dict(S_sky=f["S_sky"][:], dist_edges=f["dist_edges"][:], S_dist=f["S_dist"][:],
                    R_edges=f["R_edges"][:], S_R=f["S_R"][:],
                    nside=int(f.attrs["healpix_nside"]), nest=bool(f.attrs.get("healpix_nest", False)),
                    R_max=float(f.attrs["R_max_obs"]), sigma_dm=float(f.attrs["sigma_DM"]))


def evaluate(tab, l, b, d, R):
    """Selection probability min(S, 1) at observed (l, b, d) with observed Galactocentric R; 0 for R >= R_max."""
    pix = hp.ang2pix(tab["nside"], np.asarray(l), np.asarray(b), nest=tab["nest"], lonlat=True)
    di = np.clip(np.digitize(d, tab["dist_edges"]) - 1, 0, len(tab["S_dist"]) - 1)
    ri = np.clip(np.digitize(R, tab["R_edges"]) - 1, 0, len(tab["S_R"]) - 1)
    S = tab["S_sky"][pix] * tab["S_dist"][di] * tab["S_R"][ri]
    return np.where(R < tab["R_max"], np.minimum(S, 1.0), 0.0)


def observed_library(traj, times, phi_deg, sigma_dm, seed=11):
    """Library samples (t >= T_MIN) observed with distance errors, returned in the bar frame.

    Returns the observed bar-frame samples, their observed (l, b, d_obs, R_obs), and the number of
    snapshots kept per orbit (so that per-orbit quantities can be formed by reshaping)."""
    keep = times >= ssl.T_MIN_GYR - 1e-9
    n_keep = int(keep.sum())
    w = traj[:, keep].reshape(-1, 6).astype(np.float64)
    mw = ssl.rotate_z(w, -np.radians(phi_deg))                 # bar frame -> MW frame
    g = observe.to_galactic(*mw.T)
    d_obs = observe.scatter_distance_modulus(g["d"], sigma_dm, np.random.default_rng(seed))
    obs = np.column_stack(observe.to_galactocentric(g["l"], g["b"], d_obs, g["pml"], g["pmb"], g["vlos"]))
    R_obs = np.hypot(obs[:, 0], obs[:, 1])
    return ssl.rotate_z(obs, np.radians(phi_deg)), (g["l"], g["b"], d_obs, R_obs), n_keep


def sf_weights(tab, lbdR, n_orbits, n_keep):
    """Sample weights S(x_j) / tau_i and the per-orbit tau_i = <S> along the orbit."""
    S = evaluate(tab, *lbdR)
    tau = np.clip(S.reshape(n_orbits, n_keep).mean(axis=1), TAU_FLOOR, None)
    return S / np.repeat(tau, n_keep), tau


def library_fit_sf(traj, times, ok, tab, phi_deg, seed=11, fs=None):
    """Per-star Fourier fit of the library weighted by the selection function.

    Drop-in replacement for `test1_nbody_qmc_focused_fourier_stars.library_fit` for a mock that has one."""
    if fs is None:
        import test1_nbody_qmc_focused_fourier_stars as fs
    obs_bar, lbdR, n_keep = observed_library(traj, times, phi_deg, tab["sigma_dm"], seed)
    wt, tau = sf_weights(tab, lbdR, traj.shape[0], n_keep)
    b, X, v, m = fs.prepare(obs_bar)
    return fs.fit(b, X, v, wt[m], ok, leverage=False), tau


# ---------------------------------------------------------------------------------------------------
# Empirical selection function: S from the data / library density ratio, with no external knowledge
# ---------------------------------------------------------------------------------------------------
# The weighting above needs S, which a real survey does not hand you. But if the trial potential is right,
# the only reason the observed sample and the (re-observed) orbit library differ in (l, b, d) is the
# selection, so S can be read off their ratio. Fitted in the same separable form the mock uses,
# S = s_sky[HEALPix] * s_dist[d] * s_R[R], by iterative proportional fitting against the *data's* three
# margins, with the library as the base population instead of the simulation's full particle set.
#
# This costs no discriminating power: the likelihood fits <v>(phi) by least squares, so it is already
# independent of the density, and reweighting by a function of position cannot change the library's mean
# velocity at a given position. The risk is the opposite one - S is refitted for every trial potential, so
# it can absorb genuine model error. Whether it does is an empirical question (see PROGRESS).

N_IPF = 10


def _bin_indices(lbdR, nside, d_edges, R_edges, nest=False):
    l, b, d, R = lbdR
    return (hp.ang2pix(nside, np.asarray(l), np.asarray(b), nest=nest, lonlat=True),
            np.clip(np.digitize(d, d_edges) - 1, 0, len(d_edges) - 2),
            np.clip(np.digitize(R, R_edges) - 1, 0, len(R_edges) - 2))


def _ipf(index, weights, targets, n_iter=N_IPF):
    """Weighted iterative proportional fitting of S = prod_k f_k[index_k] so that the weighted base
    population reproduces each target margin. Generalises make_mock_selection.fit_selection to weights."""
    factors = [np.ones(len(t)) for t in targets]
    for _ in range(n_iter):
        for k, (ix, t) in enumerate(zip(index, targets)):
            others = weights * np.prod([f[i] for j, (f, i) in enumerate(zip(factors, index)) if j != k], axis=0)
            have = np.bincount(ix, weights=others, minlength=len(t))
            factors[k] = np.divide(t, have, out=np.zeros(len(t)), where=have > 0)
    return factors


def empirical_tables(data_lbdR, lib_lbdR, lib_w, nside=32, d_bin=0.5, d_max=20.0, R_bin=0.5, R_max=8.0,
                     nest=False, n_iter=N_IPF):
    """Selection tables fitted to the data's (sky, distance, R) margins using the library as the base.

    Bins default to coarser than the mock's generating grid (0.5 kpc rather than 0.25) because the data
    side is only ~50k stars. Returns the same dict `evaluate` takes."""
    d_edges = np.arange(0.0, d_max + 1e-9, d_bin)
    R_edges = np.arange(0.0, R_max + 1e-9, R_bin)
    npix = hp.nside2npix(nside)
    di = _bin_indices(data_lbdR, nside, d_edges, R_edges, nest)
    li = _bin_indices(lib_lbdR, nside, d_edges, R_edges, nest)
    inside = np.asarray(data_lbdR[3]) < R_max
    t_sky = np.bincount(di[0][inside], minlength=npix).astype(float)
    t_dist = np.bincount(di[1][inside], minlength=len(d_edges) - 1).astype(float)
    t_R = np.bincount(di[2][inside], minlength=len(R_edges) - 1).astype(float)
    ok = np.asarray(lib_lbdR[3]) < R_max
    w = np.where(ok, lib_w, 0.0)
    n = inside.sum()
    s_sky, s_dist, s_R = _ipf(li, w, [t * n / max(t.sum(), 1) for t in (t_sky, t_dist, t_R)], n_iter)
    return dict(S_sky=s_sky, dist_edges=d_edges, S_dist=s_dist, R_edges=R_edges, S_R=s_R,
                nside=nside, nest=nest, R_max=R_max, sigma_dm=np.nan)


def evaluate_raw(tab, l, b, d, R):
    """S without the min(S, 1) cap - correct for an empirical S, whose normalisation is arbitrary."""
    pix = hp.ang2pix(tab["nside"], np.asarray(l), np.asarray(b), nest=tab["nest"], lonlat=True)
    di = np.clip(np.digitize(d, tab["dist_edges"]) - 1, 0, len(tab["S_dist"]) - 1)
    ri = np.clip(np.digitize(R, tab["R_edges"]) - 1, 0, len(tab["S_R"]) - 1)
    return np.where(R < tab["R_max"], tab["S_sky"][pix] * tab["S_dist"][di] * tab["S_R"][ri], 0.0)


# ---------------------------------------------------------------------------------------------------
# Direct density ratio: S = rho_obs(l, b, d) / rho_lib(l, b, d), no functional form assumed
# ---------------------------------------------------------------------------------------------------
# The simplest reading of "make the library look like the data on the sky": bin both in (l, b, d) and take
# the ratio cell by cell. No separable form, no fitting - S is a deterministic function of the data and the
# library. The method is density-free by construction (fs.fit regresses <v>(phi)), so forcing the spatial
# distributions to agree costs no discriminating power; what it buys is the *within-block* distribution,
# which matters because the |z| layers are 1 kpc wide and the selection varies inside them.
#
# Resolution is a trade: a 50k-star sample gives a median of ~6 stars per cell at nside 8 with 1 kpc distance
# bins, so the ratio is Poisson-noisy. Vertical resolution is what matters most - a HEALPix pixel subtends
# 7.3 deg at nside 8 (1.0 kpc at 8 kpc, comparable to a whole layer) and 3.7 deg at nside 16 (0.5 kpc).

def ratio_tables(data_lbdR, lib_lbdR, lib_w, nside=16, d_bin=1.0, d_max=20.0, R_max=8.0, nest=False):
    """S = rho_obs / rho_lib per (HEALPix pixel, distance bin) cell, normalised to the data count.

    Cells the library never reaches get S = 0 (they cannot contribute anyway); cells the data never reach
    also get S = 0, which is the correct statement that the survey does not sample them."""
    d_edges = np.arange(0.0, d_max + 1e-9, d_bin)
    nd = len(d_edges) - 1
    npix = hp.nside2npix(nside)

    def cells(lbdR):
        l, b, d, R = lbdR
        pix = hp.ang2pix(nside, np.asarray(l), np.asarray(b), nest=nest, lonlat=True)
        di = np.clip(np.digitize(d, d_edges) - 1, 0, nd - 1)
        return pix * nd + di, np.asarray(R) < R_max

    dc, din = cells(data_lbdR)
    lc, lin = cells(lib_lbdR)
    n = npix * nd
    n_obs = np.bincount(dc[din], minlength=n).astype(float)
    w_lib = np.bincount(lc[lin], weights=np.asarray(lib_w)[lin], minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(w_lib > 0, n_obs / np.maximum(w_lib, 1e-30), 0.0)
    return dict(kind="ratio", S=S, nside=nside, nest=nest, d_edges=d_edges, nd=nd, R_max=R_max,
                n_cells_used=int(np.sum((n_obs > 0) & (w_lib > 0))),
                frac_data_in_empty_lib_cells=float(n_obs[w_lib <= 0].sum() / max(n_obs.sum(), 1)))


def evaluate_ratio(tab, l, b, d, R):
    pix = hp.ang2pix(tab["nside"], np.asarray(l), np.asarray(b), nest=tab["nest"], lonlat=True)
    di = np.clip(np.digitize(d, tab["d_edges"]) - 1, 0, tab["nd"] - 1)
    return np.where(np.asarray(R) < tab["R_max"], tab["S"][pix * tab["nd"] + di], 0.0)

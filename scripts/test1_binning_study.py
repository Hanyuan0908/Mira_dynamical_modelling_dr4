"""Test 1b: does re-binning the moment maps tighten the bar length / strength constraints?

No new integrations. The data and every saved orbit library (global-search points 000-064 and the four
S / bar_frac variants of run_SF_box_models.py) are compared in bins built as unions of the saved 0.5 kpc
library boxes (box -> bin by its centre; the same mapping for data and library), folded into y >= 0, z >= 0
(the library is four-fold symmetric; the stars are folded the same way, v_z -> v_z sign z). Only boxes
wholly inside R < 8 kpc and |z| < 4 kpc are used.

Binnings (|z| layers 0-0.5-1-2-4 kpc in all):
  polar15 / polar30 / polar45   R 0.5 kpc x phi 15 / 30 / 45 deg (polar15 ~ the current likelihood cells)
  polarR1_15                    R 1 kpc x phi 15 deg
  sq0.5 / sq1                   0.5 / 1 kpc squares
  adaptive100                   R 0.5 kpc rings, 12/6/4/3/2/1 sectors so that a bin holds >= ~100 stars
  fourier024 / fourier02        polar15 sectors of each ring compressed to their m = 0, 2, 4 (or 0, 2)
                                Fourier coefficients (rings with all 12 sectors >= N_MIN stars; other rings
                                keep their sectors)
Statistics per bin: <v_R>, <v_phi>, <v_z sign z>, sigma_R, sigma_phi, sigma_z, and the star count (library:
1/tau-weighted mass normalised to the data total). Errors: N_BOOT bootstrap resamplings of the stars
(the same resamplings for every point). chi^2 over bins with >= N_MIN stars; tempered
logL = -chi^2 / (2 sqrt(2 N_data)), N_data = (statistics) x (bins used at the truth).
Output: results/test1_nbody/binning_study.json
"""
import json

import numpy as np
import h5py
import pandas as pd

import steady_state_likelihood as ssl
from test1_nbody_qmc import ARR_DIR, MOCK, TRUTH, OUT_CSV

N_MIN = 20
N_BOOT = 100
SEED = 0
E_XY, E_Z = ssl.BOX_XY, ssl.BOX_Z
XC, ZC = 0.5 * (E_XY[1:] + E_XY[:-1]), 0.5 * (E_Z[1:] + E_Z[:-1])
X, Y, Z = np.meshgrid(XC, XC, ZC, indexing="ij")
RB, PB = np.hypot(X, Y), np.degrees(np.arctan2(Y, X))
QUAD = (Y > 0) & (Z > 0) & (RB + 0.25 * np.sqrt(2) < ssl.R_SEL)
IZ = np.digitize(Z, [0.0, 0.5, 1.0, 2.0, 4.0]) - 1
VARIANTS = [(1.25, 1.0, "S 1.25"), (0.8, 1.0, "S 0.8"), (1.0, 1.4, "f 1.4"), (1.0, 0.7, "f 0.7")]
STATS = ["<v_R>", "<v_phi>", "<v_z>", "sigma_R", "sigma_phi", "sigma_z", "counts"]
SETS = {"means": [0, 1, 2], "meansdisp": [0, 1, 2, 3, 4, 5], "counts": [6], "meansdisp+counts": list(range(7))}


def compress(*keys):
    lab = np.full(X.shape, -1)
    _, inv = np.unique(np.stack([np.asarray(k)[QUAD] for k in keys], 1), axis=0, return_inverse=True)
    lab[QUAD] = inv.ravel()
    return lab.ravel()


def polar(dR, dphi):
    return compress(np.floor(RB / dR), np.floor(PB / dphi), IZ)


def square(d):
    return compress(np.floor(X / d), np.floor(Y / d), IZ)


def polar_adaptive(nbox, n_target):
    ring = np.floor(RB / 0.5)
    nsec = np.ones(X.shape)
    for r in np.unique(ring[QUAD]):
        for k in range(4):
            sel = QUAD & (ring == r) & (IZ == k)
            n = nbox.reshape(X.shape)[sel].sum()
            nsec[sel] = max([m for m in (12, 6, 4, 3, 2, 1) if n / m >= n_target] + [1])
    return compress(ring, IZ, nsec, np.floor(PB / (180.0 / nsec)))


def fourier_projection(lab, modes):
    """Matrix P (rows x bins of lab = polar15) and, per row, the bins it uses: rings (R 0.5, |z| layer) whose
    12 sectors all exist are replaced by their Fourier coefficients for the given modes; others keep sectors."""
    nb = lab.max() + 1
    ring, layer, sec = (np.zeros(nb, int) for _ in range(3))
    for b in range(nb):
        i = np.flatnonzero(lab == b)[0]
        ring[b], layer[b], sec[b] = int(RB.ravel()[i] // 0.5), IZ.ravel()[i], int(PB.ravel()[i] // 15)
    rows = []
    for key in sorted(set(zip(ring, layer))):
        bins = np.flatnonzero((ring == key[0]) & (layer == key[1]))
        if len(bins) == 12:
            phi = np.radians(15 * sec[bins] + 7.5)
            for m in modes:
                for f in ((np.ones_like(phi),) if m == 0 else (np.cos(m * phi), np.sin(m * phi))):
                    row = np.zeros(nb)
                    row[bins] = f / 12 * (1 if m == 0 else 2)
                    rows.append(row)
        else:
            for b in bins:
                row = np.zeros(nb)
                row[b] = 1
                rows.append(row)
    return np.array(rows)


def stats_from_sums(N, S1, S2, unbiased=False):
    """unbiased=True for star counts (data): variance times N / (N - 1). Library box sums keep 1/W (their
    effective number of orbits per cell is ~10^3, so the bias is < 10^-3 relative)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        m = S1 / N[:, None]
        s = np.sqrt(np.maximum((S2 / N[:, None] - m**2) * ((N / (N - 1))[:, None] if unbiased else 1), 0))
    return np.column_stack([m, s, N])


def data_sums(sl, v, nb, wt):
    ok = sl >= 0
    s, w, vv = sl[ok], wt[ok], v[ok]
    N = np.bincount(s, weights=w, minlength=nb)
    S1 = np.column_stack([np.bincount(s, weights=w * vv[:, a], minlength=nb) for a in range(3)])
    S2 = np.column_stack([np.bincount(s, weights=w * vv[:, a]**2, minlength=nb) for a in range(3)])
    return N, S1, S2


def library_stats(arr, lab, nb):
    ok = lab >= 0
    l = lab[ok]
    S0 = np.bincount(l, weights=arr["box_S0"].ravel()[ok], minlength=nb)
    S1 = np.column_stack([np.bincount(l, weights=arr["box_S1"].reshape(-1, 3)[ok, a], minlength=nb) for a in range(3)])
    S2 = np.column_stack([np.bincount(l, weights=arr["box_S2"].reshape(-1, 6)[ok, a], minlength=nb) for a in range(3)])
    return stats_from_sums(S0, S1, S2)


def fold_stars(stars_mw, phi_deg):
    w = ssl.rotate_z(stars_mw, np.radians(phi_deg))
    f = w[:, 1] < 0
    w[np.ix_(f, [0, 1, 3, 4])] *= -1
    f = w[:, 2] < 0
    w[np.ix_(f, [2, 5])] *= -1
    R = np.hypot(w[:, 0], w[:, 1])
    v = np.column_stack([(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R, -(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R,
                         w[:, 5]])
    ix, iy, iz = (np.digitize(w[:, 0], E_XY) - 1, np.digitize(w[:, 1], E_XY) - 1, np.digitize(w[:, 2], E_Z) - 1)
    ok = (ix >= 0) & (ix < len(XC)) & (iy >= 0) & (iy < len(XC)) & (iz >= 0) & (iz < len(ZC))
    box = np.full(len(w), -1)
    box[ok] = np.ravel_multi_index((ix[ok], iy[ok], iz[ok]), X.shape)
    return box, v


class Binning:
    def __init__(self, name, lab, P=None):
        self.name, self.lab, self.nb = name, lab, lab.max() + 1
        self.P = np.eye(self.nb) if P is None else P
        self.uses = (np.abs(self.P) > 0).astype(float)

    def project(self, A):
        """P @ A, NaN for every row that uses a bin with a non-finite value."""
        out = self.P @ np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
        out[self.uses @ (~np.isfinite(A)).astype(float) > 0] = np.nan
        return out

    def observed(self, box, v, boot_w):
        sl = np.where(box >= 0, self.lab[np.maximum(box, 0)], -1)
        main = stats_from_sums(*data_sums(sl, v, self.nb, np.ones(len(v))), unbiased=True)
        boots = np.array([self.project(stats_from_sums(*data_sums(sl, v, self.nb, w), unbiased=True)) for w in boot_w])
        N = main[:, 6]
        rowN = np.array([N[r != 0].min() for r in self.P])          # fewest stars among the bins a row uses
        return dict(d=self.project(main), err=np.nanstd(boots, 0), rowN=rowN, binN=N, Ntot=N[N >= N_MIN].sum())

    def chi2(self, obs, lib):
        """chi^2 per statistic (7,) and the number of rows used."""
        scale = obs["Ntot"] / lib[obs["binN"] >= N_MIN, 6].sum()     # library mass -> data counts
        L = lib.copy()
        L[:, 6] *= scale
        mdl = self.project(L)
        used = (obs["rowN"] >= N_MIN) & np.all(np.isfinite(mdl), 1) & np.all(obs["err"] > 0, 1)
        return np.sum(((obs["d"] - mdl)[used] / obs["err"][used])**2, axis=0), int(used.sum())


def main():
    with h5py.File(MOCK, "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    rng = np.random.default_rng(SEED)
    boot_w = rng.multinomial(len(stars_mw), np.full(len(stars_mw), 1 / len(stars_mw)), size=N_BOOT).astype(float)
    box_t, v_t = fold_stars(stars_mw, TRUTH["phi_b_deg"])
    nbox = np.bincount(box_t[box_t >= 0], minlength=X.size)
    p15 = polar(0.5, 15)
    binnings = [Binning("polar15", p15), Binning("polar30", polar(0.5, 30)), Binning("polar45", polar(0.5, 45)),
                Binning("polarR1_15", polar(1.0, 15)), Binning("sq0.5", square(0.5)), Binning("sq1", square(1.0)),
                Binning("adaptive100", polar_adaptive(nbox, 100)),
                Binning("fourier024", p15, fourier_projection(p15, (0, 2, 4))),
                Binning("fourier02", p15, fourier_projection(p15, (0, 2)))]

    df = pd.read_csv(OUT_CSV).sort_values("idx").reset_index(drop=True)
    points = [(f"{int(r.idx):03d}.npz", r.phi_b_deg, r.kind) for r in df.itertuples()]
    points += [(f"maps_S{S:g}_f{fr:g}.npz", TRUTH["phi_b_deg"], lab) for S, fr, lab in VARIANTS]
    res = {b.name: {} for b in binnings}
    obs_cache = {}
    for fn, phi, kind in points:
        key = round(float(phi), 6)
        if key not in obs_cache:
            box, v = fold_stars(stars_mw, phi)
            obs_cache = {key: [b.observed(box, v, boot_w) for b in binnings]}
        with np.load(ARR_DIR / fn) as a:
            arr = {k: a[k] for k in ("box_S0", "box_S1", "box_S2")}
        for b, obs in zip(binnings, obs_cache[key]):
            c, n = b.chi2(obs, library_stats(arr, b.lab, b.nb))
            res[b.name][fn] = dict(chi2=c.tolist(), n_rows=n, kind=kind)

    summary = {}
    n_stars = len(stars_mw)
    for b, obs in zip(binnings, [b.observed(box_t, v_t, boot_w[:10]) for b in binnings]):
        r = res[b.name]
        t = r["000.npz"]
        out = dict(n_rows_truth=t["n_rows"], n_bins=int(b.nb),
                   frac_stars_used=float(obs["Ntot"] / n_stars), sets={})
        for sname, idx in SETS.items():
            nd = len(idx) * t["n_rows"]
            temper = np.sqrt(2 * nd)
            c_t = sum(t["chi2"][i] for i in idx)
            dl = {fn: -(sum(x["chi2"][i] for i in idx) - c_t) / (2 * temper) for fn, x in r.items()}
            qmc = np.array([dl[fn] for fn in r if r[fn]["kind"] == "qmc"])
            out["sets"][sname] = dict(
                N_data=nd, sqrt_2N=float(temper), chi2_per_N_truth=float(c_t / nd),
                delta_chi2={lab: float(sum(r[f"maps_S{S:g}_f{fr:g}.npz"]["chi2"][i] for i in idx) - c_t)
                            for S, fr, lab in VARIANTS},
                dlogL={lab: float(dl[f"maps_S{S:g}_f{fr:g}.npz"]) for S, fr, lab in VARIANTS},
                qmc_best=float(qmc.max()), qmc_n_gt_m2=int((qmc > -2).sum()), qmc_n_gt_m4p5=int((qmc > -4.5).sum()))
        summary[b.name] = out
    with open(OUT_CSV.parent / "binning_study.json", "w") as fh:
        json.dump(dict(summary=summary, per_point=res, n_boot=N_BOOT, n_min=N_MIN), fh, indent=1)

    for sname in SETS:
        print(f"\n=== {sname} ===")
        print(f"{'binning':12s} {'rows':>5s} {'N_data':>6s} {'stars':>6s} {'chi2/N':>6s} {'sqrt2N':>6s} | "
              + " ".join(f"{l:>14s}" for _, _, l in VARIANTS) + " | best  n>-2 n>-4.5")
        for name, o in summary.items():
            s = o["sets"][sname]
            print(f"{name:12s} {o['n_rows_truth']:5d} {s['N_data']:6d} {o['frac_stars_used']:6.3f} "
                  f"{s['chi2_per_N_truth']:6.2f} {s['sqrt_2N']:6.1f} | "
                  + " ".join(f"{s['delta_chi2'][l]:+6.0f} ({s['dlogL'][l]:+5.2f})" for _, _, l in VARIANTS)
                  + f" | {s['qmc_best']:+5.2f} {s['qmc_n_gt_m2']:3d} {s['qmc_n_gt_m4p5']:4d}")


if __name__ == "__main__":
    main()

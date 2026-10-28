"""Test 2b, first look: what Mira-like distance errors do to the per-star Fourier coefficients the likelihood compares.

Data: the perfect 50k N-body mock and the same stars with sigma_DM = 0.265 mag (make_mock_distance_error.py), both
rotated by the true bar angle (27 deg) and fitted exactly as in test1_nbody_qmc_focused_fourier_stars.py (300 bootstraps).
Model: the truth library (global search point 000), fitted on the blocks usable in each data set, (a) as it is (the
"do nothing" case) and (b) forward-modelled: every library sample with true R < 8 kpc (the mock's selection) is rotated
back to the MW-like frame, observed from the Sun, given the same distance-modulus error, converted back and rotated to the
bar frame; samples that land at observed R > 8 kpc drop out of the grid, as the smeared data do. (b) tests the
comparison in observed space when the orbits are right; the library here starts from the true positions, which real
data do not provide.
Numbers: chi2 / N of each data set against the truth model (and of the smeared data against the forward-modelled
truth model) for the likelihood schemes (joint covariance, Hartlap),
per field, and the shift of the smeared coefficients from the perfect ones in units of the errors.
Figures: figures/test1_nbody_distance_error_fourier_coefficients_perfect_vs_smeared_absz_{0to1,1to2}.png
Numbers: results/test1_nbody_qmc_focused/distance_error_effect_on_fourier_coefficients.json
"""
import json

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
import test1_nbody_qmc_focused_fourier_stars as fs
import test1_nbody_qmc as old
import observe
from test1_nbody_qmc_focused_run import library_samples
from plot_fourier_coefficients_data_vs_models import as_kms, FIELDS, TERMS, LAYERS, RC

SMEARED = MOCK_DIR / "nbody_mock_distance_error_sigmaDM0.265_R8_N50k.hdf5"
NK = fs.NK


def stars(path):
    with h5py.File(path, "r") as f:
        return np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])


def forward_smeared_library_fit(traj, times, ok, sigma_dm, phi_deg, seed=1):
    """Per-star Fourier fit of the library after observing every sample (true R < 8 kpc) with distance errors."""
    w, wt, _ = library_samples(times, traj)
    sel = np.hypot(w[:, 0], w[:, 1]) < ssl.R_SEL
    w, wt = w[sel], wt[sel]
    mw = ssl.rotate_z(w, -np.radians(phi_deg))
    g = observe.to_galactic(*mw.T)
    d_obs = observe.scatter_distance_modulus(g["d"], sigma_dm, np.random.default_rng(seed))
    obs = np.column_stack(observe.to_galactocentric(g["l"], g["b"], d_obs, g["pml"], g["pmb"], g["vlos"]))
    b, X, v, m = fs.prepare(ssl.rotate_z(obs, np.radians(phi_deg)))
    return fs.fit(b, X, v, wt[m], ok, leverage=False)


def main():
    set_style()
    rng = np.random.default_rng(fs.SEED)
    boot_w = rng.multinomial(50000, np.full(50000, 1 / 50000), size=fs.N_BOOT).astype(float)
    traj = np.load(old.ARR_DIR / "traj_000.npy")
    times = np.load(old.ARR_DIR / "000.npz")["times"]
    sets = {}
    for tag, path in (("perfect", old.MOCK), ("smeared", SMEARED)):
        cd, cov, n, ok = fs.data_fit(stars(path), old.TRUTH["phi_b_deg"], boot_w)
        cl = fs.library_fit(traj, times, ok)
        valid = ok & np.all(np.isfinite(cd), axis=(1, 2)) & np.all(np.isfinite(cl), axis=(1, 2))
        sets[tag] = dict(cd=cd, cov=cov, n=n, cl=cl, valid=valid, ok=ok)
    sigma_dm = float(h5py.File(SMEARED, "r").attrs["sigma_DM"])
    cl_fwd = forward_smeared_library_fit(traj, times, sets["smeared"]["ok"], sigma_dm, old.TRUTH["phi_b_deg"])
    sets["smeared_vs_forward_model"] = dict(sets["smeared"], cl=cl_fwd,
                                            valid=sets["smeared"]["valid"] & np.all(np.isfinite(cl_fwd), axis=(1, 2)))
    res = dict(sigma_DM=sigma_dm, chi2_per_N_vs_truth_model={})
    schemes = dict(fs.SCHEMES)
    for f, name in enumerate(fs.FIELD_NAMES):
        schemes[f"field_{name}_m024"] = (fs.N_RING_DEFAULT, {f: fs.ALL})
        schemes[f"field_{name}_m0"] = (fs.N_RING_DEFAULT, {f: [0]})
        schemes[f"field_{name}_m24"] = (fs.N_RING_DEFAULT, {f: fs.NO_M0})
    for name, sc in schemes.items():
        row = {}
        for tag, s in sets.items():
            chi, nd = fs.scheme_chi2(s["cd"], s["cl"], s["cov"], s["valid"], sc)
            row[tag] = [round(chi / nd, 3), nd]
        res["chi2_per_N_vs_truth_model"][name] = row
    p, s = sets["perfect"], sets["smeared"]
    res["stars_per_layer"] = {k: [int(x) for x in sets[k]["n"].reshape(fs.NR, fs.NL).sum(0)] for k in ("perfect", "smeared")}
    res["error_ratio_smeared_over_perfect_median"] = {
        name: round(float(np.nanmedian(np.sqrt(np.diagonal(s["cov"], axis1=1, axis2=2)[:, f * NK:(f + 1) * NK]
                                               / np.diagonal(p["cov"], axis1=1, axis2=2)[:, f * NK:(f + 1) * NK]))), 3)
        for f, name in enumerate(fs.FIELD_NAMES)}

    shape = (fs.NR, fs.NL, 6, NK)
    for layer, tag, label in LAYERS:
        fig, axes = plt.subplots(len(FIELDS), len(TERMS), figsize=(17.5, 12.5), gridspec_kw=dict(hspace=0.42, wspace=0.32))
        for r, (name, f, kind) in enumerate(FIELDS):
            mod, _ = as_kms(s["cl"].reshape(shape)[:, layer, f], kind=kind)
            fwd, _ = as_kms(cl_fwd.reshape(shape)[:, layer, f], kind=kind)
            series = {}
            for key, st in (("perfect", p), ("smeared", s)):
                e = np.sqrt(np.diagonal(st["cov"], axis1=1, axis2=2)).reshape(shape)[:, layer, f]
                series[key] = as_kms(st["cd"].reshape(shape)[:, layer, f], e, kind) + (st["valid"].reshape(fs.NR, fs.NL)[:, layer],)
            for k, term in enumerate(TERMS):
                ax = axes[r, k]
                ax.axvspan(6.0, 8.0, color=INK_MUTED, alpha=0.12, lw=0)
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                use = series["smeared"][2]
                ax.plot(RC[use], mod[use, k], color=ACCENT, lw=1.8, label="truth model (#0)" if (r, k) == (0, 0) else None)
                ax.plot(RC[use], fwd[use, k], color=BLUE_RAMP[5], lw=1.8, ls="--",
                        label="truth model, observed with the same distance errors" if (r, k) == (0, 0) else None)
                notes = []
                for key, mfc, off in (("perfect", "white", -0.05), ("smeared", INK, 0.05)):
                    val, err, ok = series[key]
                    ax.errorbar(RC[ok] + off, val[ok, k], yerr=err[ok, k], fmt="o", ms=3.4, mfc=mfc, mec=INK, color=INK,
                                ecolor=INK_2, elinewidth=0.8, capsize=0, zorder=4,
                                label=(f"data, {key}" + (" (σ_DM = 0.265 mag)" if key == "smeared" else ""))
                                if (r, k) == (0, 0) else None)
                    mod_k, _ = as_kms(sets[key]["cl"].reshape(shape)[:, layer, f], kind=kind)
                    r6 = ok & (np.arange(fs.NR) < fs.N_RING_DEFAULT)
                    notes.append(f"{np.mean(((val[r6, k] - mod_k[r6, k]) / err[r6, k])**2):.2f}")
                val, err, ok = series["smeared"]
                r6 = ok & (np.arange(fs.NR) < fs.N_RING_DEFAULT)
                notes.append(f"{np.mean(((val[r6, k] - fwd[r6, k]) / err[r6, k])**2):.2f}")
                ax.text(0.02, 0.97, f"χ²/N  perfect: {notes[0]}\nsmeared: {notes[1]}\nsmeared, observed model: {notes[2]}",
                        transform=ax.transAxes,
                        fontsize=6.0, va="top", color=INK_2, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
                ax.set_xlim(0, 8)
                ax.tick_params(labelsize=7.5)
                if r == 0:
                    ax.set_title(term, fontsize=10)
                if k == 0:
                    ax.set_ylabel(f"{name}\n" + ("km/s" if kind == "mean" else "σ [km/s]"), fontsize=9.5)
                if r == len(FIELDS) - 1:
                    ax.set_xlabel("R [kpc] (observed)", fontsize=9)
        fig.legend(loc="upper center", ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.965))
        fig.suptitle(f"Distance errors (σ_DM = 0.265 mag, 12.2%): Fourier coefficients per 0.3 kpc ring, {label}, bar frame "
                     f"at the true angle; σ terms: m = 0 is σ₀, others c_m / 2σ₀\nχ²/N per panel against the truth model "
                     f"(diagonal errors, R < 6 kpc) for the perfect data, the smeared data, and the smeared data against the "
                     f"model observed with the same errors", fontsize=10.5, y=0.997)
        savefig(fig, f"test1_nbody_distance_error_fourier_coefficients_perfect_vs_smeared_{tag}")
    out = ROOT / "results" / "test1_nbody_qmc_focused" / "distance_error_effect_on_fourier_coefficients.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

"""Why are |Omega_b| and phi_b well constrained but the bar length S and strength bar_frac not? (focused global search)

From the GP-emulator posterior of one logLs column (same emulator, seeds and sampler settings as test1_gp_emulator.py):
  1. correlation matrix of (|Omega_b|, phi_b, ln S, ln bar_frac);
  2. marginal 1-sigma widths (half the 16-84% range) against conditional ones (Delta logL = -0.5 along each axis through
     the emulator maximum, the other three fixed). A parameter that trades off against another has conditional << marginal;
  3. principal axes of the (ln S, ln bar_frac) posterior: the best- and the worst-constrained combination;
  4. what the data pin down. The bar part of a trial potential is bar_frac * Phi_bar(r / S): AGAMA's scale = [f S, S]
     multiplies the potential by f S / S = f and stretches it by S (checked numerically here). Its m = 2 amplitude in the
     plane is therefore A2(R) = bar_frac * A2_0(R / S), and the bar's tangential-to-radial force ratio
     Q_T(R) = 2 A2(R) / (R |dPhi_0/dR|) (the axisymmetric part does not depend on S or bar_frac). The fractional 1-sigma
     width of A2(R) over the posterior samples, as a function of R, shows at which radius the bar strength is measured
     and how well.
Usage: python analyse_bar_length_strength_degeneracy.py [--column logLs_...] [--csv path/to/evals.csv]
Figure: figures/test1_nbody_qmc_focused_bar_length_strength_degeneracy_<column>.png
Numbers: results/test1_nbody_qmc_focused/bar_length_strength_degeneracy_<column>.json
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt
import emcee
import agama
from scipy.optimize import minimize, brentq

from common import ROOT, BAR_ANGLE_DEG, set_style, savefig, INK_2, ACCENT, BLUE_RAMP, CMAP_SEQ
import steady_state_likelihood  # noqa: F401  (AGAMA units, as for the trial potentials)
from potential_family import files as family_files
from galaxies import galaxy as galaxy_info
import test1_gp_emulator as g

KEYS = ["Omega_abs", "phi_b_deg", "S", "bar_frac"]
LABELS = ["|Ω_b|", "φ_b", "ln S", "ln bar_frac"]
R_GRID = np.round(np.arange(0.2, 8.0 + 1e-9, 0.1), 6)
N_PHI = 64


def posterior(column, csv):
    """Emulator, its maximum and the posterior chain in the unit cube (as in test1_gp_emulator.main)."""
    g.use_run("nbody_qmc_focused")
    g.RUN["csv"] = csv
    df, _ = g.load_table(column)
    tr = df[df.kind == "qmc"]
    u, y = g.to_unit(tr), tr["logL"].values
    emu = g.Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r

    def logpost(q):
        return -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0]
    nw = 32
    p0 = np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((nw, 4)), 0, 1)
    sampler = emcee.EnsembleSampler(nw, 4, logpost)
    sampler.random_state = np.random.RandomState(2).get_state()
    sampler.run_mcmc(p0, 3000, progress=False)
    return emu, best, sampler.get_chain(discard=1000, thin=5, flat=True), len(tr)


def to_log_params(p):
    """(N, 4) physical (Omega, phi, S, f) -> (Omega, phi, ln S, ln f)."""
    return np.column_stack([p[:, 0], p[:, 1], np.log(p[:, 2]), np.log(p[:, 3])])


def conditional_widths(emu, best):
    """Half-width of the Delta logL = -0.5 interval along each axis through the maximum, in (Omega, phi, ln S, ln f)."""
    top = -best.fun
    out = []
    for j in range(4):
        def f(t):
            q = best.x.copy()
            q[j] = t
            return emu.predict(q)[0] - (top - 0.5)
        ends, hit = [], []
        for edge in (0.0, 1.0):
            if f(edge) < 0:
                ends.append(brentq(f, min(edge, best.x[j]), max(edge, best.x[j])))
                hit.append(True)
            else:
                ends.append(edge)
                hit.append(False)
        lo, hi = (g.LO[j] + e * (g.HI[j] - g.LO[j]) for e in ends)     # box is in (Omega, phi, ln S, ln f)
        out.append(dict(half_width=float((hi - lo) / 2), closed=bool(all(hit))))
    return out


def bar_profiles():
    """A2_0(R): m = 2 amplitude of the unscaled bar potential in the plane; |dPhi_0/dR| of the axisymmetric potential."""
    F = family_files("nbody")
    bar = agama.Potential(agama.Potential(file=str(F["baryon_full"]), rotation=np.radians(BAR_ANGLE_DEG)),
                          agama.Potential(file=str(F["baryon_axi"]), scale=[-1.0, 1.0]))
    axi = agama.Potential(agama.Potential(file=str(F["baryon_axi"])),
                          agama.Potential(str(galaxy_info("nbody")["pot_dir"] / "dm_spherical.ini")))
    phi = np.linspace(0, 2 * np.pi, N_PHI, endpoint=False)

    def a2(pot, R):
        x = np.column_stack([np.repeat(R, N_PHI) * np.tile(np.cos(phi), len(R)),
                             np.repeat(R, N_PHI) * np.tile(np.sin(phi), len(R)), np.zeros(len(R) * N_PHI)])
        P = pot.potential(x).reshape(len(R), N_PHI)
        return np.abs(2.0 / N_PHI * np.sum(P * np.exp(-2j * phi), axis=1))
    R_fine = np.linspace(0.01, 12.0, 1200)
    A2_0 = a2(bar, R_fine)
    FR = np.abs(axi.force(np.column_stack([R_GRID, np.zeros_like(R_GRID), np.zeros_like(R_GRID)]))[:, 0])
    # check of the scaling convention used by ssl.bar_aligned_potential
    S, f = 1.3, 0.8
    scaled = agama.Potential(agama.Potential(file=str(F["baryon_full"]), scale=[f * S, S], rotation=np.radians(BAR_ANGLE_DEG)),
                             agama.Potential(file=str(F["baryon_axi"]), scale=[-f * S, S]))
    Rc = np.array([1.0, 3.0, 5.0])
    check = float(np.max(np.abs(a2(scaled, Rc) / (f * np.interp(Rc / S, R_fine, A2_0)) - 1)))
    return R_fine, A2_0, FR, check


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--column", default="logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024")
    ap.add_argument("--csv", default=str(ROOT / "results" / "test1_nbody_qmc_focused" / "evals.csv"))
    args = ap.parse_args()
    set_style()
    emu, best, chain, n_train = posterior(args.column, args.csv)
    post = g.from_unit(chain)
    lp = to_log_params(post)
    q = np.percentile(post, [16, 50, 84], axis=0)
    res = dict(column=args.column, csv=args.csv, n_train=n_train,
               posterior={k: [round(float(a), 4) for a in q[:, i]] for i, k in enumerate(KEYS)})

    # 1-2. correlations, marginal vs conditional widths
    res["correlation"] = np.round(np.corrcoef(lp.T), 3).tolist()
    marg = (np.percentile(lp, 84, axis=0) - np.percentile(lp, 16, axis=0)) / 2
    cond = conditional_widths(emu, best)
    res["widths"] = {lab: dict(marginal=round(float(m), 4), conditional=round(c["half_width"], 4),
                               conditional_closed_in_box=c["closed"], ratio=round(float(m / c["half_width"]), 2))
                     for lab, m, c in zip(LABELS, marg, cond)}
    res["widths_note"] = ("|Omega_b| in km/s/kpc, phi_b in deg, ln S and ln bar_frac dimensionless (= fractional); "
                          "conditional = Delta logL = -0.5 along one axis through the emulator maximum")

    # 3. principal axes in (ln S, ln f)
    w, V = np.linalg.eigh(np.cov(lp[:, 2:].T))
    narrow, wide = V[:, 0], V[:, 1]
    res["lnS_lnf_principal_axes"] = dict(
        narrow_sigma=round(float(np.sqrt(w[0])), 4), wide_sigma=round(float(np.sqrt(w[1])), 4),
        narrow_combination=f"ln bar_frac + {narrow[0] / narrow[1]:.2f} ln S",
        wide_direction_slope_dlnf_dlnS=round(float(wide[1] / wide[0]), 2))

    # 4. what radius the bar strength is measured at
    R_fine, A2_0, FR, check = bar_profiles()
    S, f = post[:, 2], post[:, 3]
    lnA2 = np.log(f[:, None] * np.interp(R_GRID[None, :] / S[:, None], R_fine, A2_0))
    widthA2 = (np.percentile(lnA2, 84, axis=0) - np.percentile(lnA2, 16, axis=0)) / 2
    QT = 2 * np.exp(lnA2) / (R_GRID * FR)
    QT_true = 2 * np.interp(R_GRID, R_fine, A2_0) / (R_GRID * FR)
    i_star = int(np.argmin(widthA2))
    R_star = float(R_GRID[i_star])
    res["scale_convention_check_max_rel_err"] = check
    res["bar_strength_vs_R"] = dict(R_best=R_star, frac_width_at_R_best=round(float(widthA2[i_star]), 4),
                                    frac_width_at={f"{r:g}": round(float(np.interp(r, R_GRID, widthA2)), 4)
                                                   for r in (1, 2, 3, 4, 5, 6)},
                                    QT_true_at_R_best=round(float(QT_true[i_star]), 4),
                                    QT_posterior_at_R_best=[round(float(x), 4) for x in
                                                            np.percentile(QT[:, i_star], [16, 50, 84])],
                                    QT_true_max=[round(float(QT_true.max()), 4), float(R_GRID[np.argmax(QT_true)])])

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), gridspec_kw=dict(wspace=0.3))
    ax = axes[0]
    H, xe, ye = np.histogram2d(S, f, bins=[np.geomspace(0.6, 1.65, 45), np.geomspace(0.45, 2.0, 45)])
    ax.pcolormesh(xe, ye, H.T, cmap=CMAP_SEQ, vmin=0.5, rasterized=True)
    Ss = np.geomspace(0.6, 1.65, 200)
    for R0, ls in ((2.0, ":"), (R_star, "-"), (5.5, "--")):
        target = np.interp(R0, R_fine, A2_0)
        ax.plot(Ss, target / np.interp(R0 / Ss, R_fine, A2_0), color=ACCENT if R0 == R_star else INK_2, lw=1.4, ls=ls,
                label=f"same bar m = 2 potential as the truth at R = {R0:.1f} kpc")
    ax.plot(1, 1, marker="*", color=ACCENT, ms=13, mec="white", ls="none", label="truth")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.6, 1.65)
    ax.set_ylim(0.45, 2.0)
    ax.set_xticks([0.6, 0.8, 1.0, 1.25, 1.6])
    ax.set_xticklabels(["0.6", "0.8", "1.0", "1.25", "1.6"])
    ax.set_yticks([0.5, 0.7, 1.0, 1.4, 2.0])
    ax.set_yticklabels(["0.5", "0.7", "1.0", "1.4", "2.0"])
    ax.set_xlabel("bar length S")
    ax.set_ylabel("bar strength bar_frac")
    ax.set_title(f"posterior samples; corr(ln S, ln f) = {res['correlation'][2][3]:+.2f}", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right", frameon=False)

    ax = axes[1]
    ax.plot(R_GRID, widthA2, color=BLUE_RAMP[5], lw=2, label="bar m = 2 potential at R (= bar force ratio Q_T(R))")
    for val, lab, ls in ((marg[2], "bar length S", "--"), (marg[3], "bar strength bar_frac", ":"),
                         (marg[0] / q[1, 0], "|Ω_b|", "-.")):
        ax.axhline(val, color=INK_2, lw=1.1, ls=ls, label=f"{lab}: {100 * val:.0f}%")
    ax.axvline(R_star, color=ACCENT, lw=1)
    ax.text(R_star + 0.1, widthA2[i_star] * 0.6, f"R = {R_star:.1f} kpc: {100 * widthA2[i_star]:.1f}%", color=ACCENT,
            fontsize=8.5)
    ax.set_ylim(0, max(0.45, float(widthA2.max()) * 1.05))
    ax.set_xlim(0, 8)
    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("fractional 1σ posterior width")
    ax.set_title("how well each quantity is measured", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")

    ax = axes[2]
    lo, mid, hi = np.percentile(QT, [16, 50, 84], axis=0)
    ax.fill_between(R_GRID, lo, hi, color=BLUE_RAMP[1], lw=0, label="posterior 16–84%")
    ax.plot(R_GRID, mid, color=BLUE_RAMP[5], lw=1.5, label="posterior median")
    ax.plot(R_GRID, QT_true, color=ACCENT, lw=1.5, ls="--", label="truth")
    ax.axvline(R_star, color=ACCENT, lw=1)
    ax.set_xlim(0, 8)
    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("bar force ratio Q_T(R) = 2 A₂(R) / (R |dΦ₀/dR|)")
    ax.set_title("bar strength profile", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle(f"Bar length vs strength: what the data constrain ({args.column}, {n_train} Sobol points, GP-emulator "
                 f"posterior)", fontsize=10.5, y=1.0)
    savefig(fig, f"test1_nbody_qmc_focused_bar_length_strength_degeneracy_{args.column}")
    out = ROOT / "results" / "test1_nbody_qmc_focused" / f"bar_length_strength_degeneracy_{args.column}.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

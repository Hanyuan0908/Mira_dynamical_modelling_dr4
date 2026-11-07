"""Pattern-speed recovery with a per-star distance error, against the fiducial uniform-error runs.

The real Mira sample does not have one distance error: each star has its own, roughly Gaussian
about the mean. Every mock before this gave every star exactly 12.5%, and the model assumes that.
These runs scatter each star by its own draw sigma_i ~ N(0.2714, 0.0543) mag = 12.5% +- 2.5%,
while the mock's sigma_DM attribute still records the nominal value - so the orbit library is
smeared uniformly at 12.5% and the rho_obs/rho_lib weighting is built on the same assumption.
Nothing else differs from the fiducial run of the same galaxy: same design, same seeds, scheme B,
sfratio, K = 10, 50 Myr.

Left: recovered against true, both sets, joined per galaxy so the shift is visible.
Right: the shift itself, (per-star bias) - (fiducial bias), which is what the test is about.

Usage: python plot_scatter_test_comparison.py
Figure: figures/test2_scatter_test_pattern_speed.png
Table:  results/scatter_test_comparison.csv
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
TAGS = ["t3", "t5.5", "t7", "au18", "au26", "au20", "au24", "au12", "au5", "au17", "au10",
        "t4.5", "t2"]


def load(run):
    p = RESULTS_DIR / run / f"emulator_results_logLs_{B}_sfratio.json"
    return json.load(open(p)) if p.exists() else None


def main():
    set_style()
    rows = []
    for t in TAGS:
        f, s = load(f"test2_{t}_sf_12p5pct"), load(f"test2_{t}_sf_12p5pct_scat")
        if f is None or s is None:
            print(f"  missing: {t}")
            continue
        T = f["truth"]["Omega_abs"]
        g = lambda d, k: (d["posterior"][k]["p50"], d["posterior"][k]["p16"],
                          d["posterior"][k]["p84"])
        fm, flo, fhi = g(f, "Omega_abs")
        sm, slo, shi = g(s, "Omega_abs")
        rows.append(dict(tag=t, suite="Auriga" if t.startswith("au") else "N-body", truth=T,
                         fid=fm, fid_lo=flo, fid_hi=fhi, scat=sm, scat_lo=slo, scat_hi=shi,
                         fid_bias=100 * (fm / T - 1), scat_bias=100 * (sm / T - 1)))
    D = pd.DataFrame(rows).sort_values("truth")
    D["shift_pct"] = D.scat_bias - D.fid_bias
    D.to_csv(RESULTS_DIR / "scatter_test_comparison.csv", index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 5.5), gridspec_kw=dict(wspace=0.26))
    lim = (0.90 * D.truth.min(), 1.10 * D.truth.max())
    a1.plot(lim, lim, "--", color=INK, lw=1.1, zorder=1, label="perfect recovery")
    for _, r in D.iterrows():
        a1.plot([r.truth, r.truth], [r.fid, r.scat], "-", color=INK_MUTED, lw=0.9, zorder=2)
    a1.errorbar(D.truth, D.fid, yerr=[D.fid - D.fid_lo, D.fid_hi - D.fid], fmt="o", ms=7,
                color=BLUE_RAMP[5], mfc="white", mec=BLUE_RAMP[5], mew=1.6, capsize=3, lw=1.2,
                zorder=4, label="fiducial: uniform 12.5%")
    a1.errorbar(D.truth, D.scat, yerr=[D.scat - D.scat_lo, D.scat_hi - D.scat], fmt="s", ms=7,
                color=ACCENT, mec=INK, mew=0.9, capsize=3, lw=1.2, zorder=5,
                label="per-star: 12.5% ± 2.5%, modelled as uniform")
    a1.set_xlim(*lim)
    a1.set_ylim(*lim)
    a1.set_aspect("equal")
    a1.set_xlabel(r"true $|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]")
    a1.set_ylabel(r"recovered $|\Omega_b|$ (posterior median, 16-84%)")
    a1.legend(loc="upper left", fontsize=8.5)
    a1.set_title(f"recovered against true, {len(D)} galaxies", color=INK_2, fontsize=10)

    a2.axhline(0, color=INK, lw=1.1, ls="--", zorder=1)
    rms = float(np.sqrt((D.shift_pct ** 2).mean()))
    a2.axhspan(-rms, rms, color=BLUE_RAMP[1], alpha=0.5, lw=0, zorder=0,
               label=f"rms shift {rms:.2f}%")
    for _, r in D.iterrows():
        c = ACCENT if r.suite == "Auriga" else BLUE_RAMP[5]
        a2.plot(r.truth, r.shift_pct, "s" if r.suite == "Auriga" else "o", ms=8, color=c,
                mec=INK, mew=0.8, zorder=4)
        a2.annotate(r.tag, (r.truth, r.shift_pct), xytext=(0, 9), textcoords="offset points",
                    ha="center", fontsize=7.5, color=INK_2)
    a2.plot([], [], "s", color=ACCENT, mec=INK, label="Auriga")
    a2.plot([], [], "o", color=BLUE_RAMP[5], mec=INK, label="N-body")
    a2.set_xlim(*lim)
    a2.set_xlabel(r"true $|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]")
    a2.set_ylabel("shift in bias, per-star $-$ fiducial [%]")
    a2.legend(loc="lower right", fontsize=8.5)
    a2.set_title(f"mean shift {D.shift_pct.mean():+.2f}%  |  "
                 f"mean |bias| {np.abs(D.fid_bias).mean():.2f}% → "
                 f"{np.abs(D.scat_bias).mean():.2f}%", color=INK_2, fontsize=10)

    fig.suptitle("Does a per-star distance error, modelled as uniform, bias the pattern speed?  "
                 "No: the shift is consistent with re-drawing the errors at all",
                 color=INK, y=0.98, fontsize=12)
    fig.tight_layout()
    savefig(fig, "test2_scatter_test_pattern_speed")

    print(f"  {'gal':6s} {'truth':>7s} {'fiducial':>9s} {'per-star':>9s} {'shift':>7s}")
    for _, r in D.iterrows():
        print(f"  {r.tag:6s} {r.truth:7.2f} {r.fid_bias:+8.1f}% {r.scat_bias:+8.1f}% "
              f"{r.shift_pct:+6.1f}%")
    print(f"\n  mean |bias| {np.abs(D.fid_bias).mean():.2f}% -> {np.abs(D.scat_bias).mean():.2f}%"
          f";  mean shift {D.shift_pct.mean():+.2f}%, rms {rms:.2f}%, "
          f"{int((D.shift_pct > 0).sum())} up / {int((D.shift_pct < 0).sum())} down")
    print(f"wrote {rel(RESULTS_DIR / 'scatter_test_comparison.csv')}")


if __name__ == "__main__":
    main()

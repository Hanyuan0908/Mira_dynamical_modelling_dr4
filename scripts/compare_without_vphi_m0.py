"""Posteriors with and without the <v_phi> m = 0 term (rescore_without_vphi_m0.py + test1_gp_emulator.py).

For every run with both emulator results: |Omega_b|, phi_b, S and bar_frac (median, 16-84%), tempered as the pipeline
(sqrt(2 N) with each scheme's own N). Truth: |Omega_b| from the simulation table's finite-difference column for the
N-body snapshots and the run's measured value for Auriga (as the presentation figures); 32.5 for the steady-state
Hunter+24 mock; S = bar_frac = 1 and phi_b = 27 deg (25 for the steady-state mock); none for OSARG.
Summary over the 18 mock galaxies: mean |bias|, median 16-84% width, and how often the truth lies inside 16-84%.

Usage: python compare_without_vphi_m0.py
Output: results/without_vphi_m0_summary.csv, figures/test2_pattern_speed_and_bar_recovery_without_vphi_m0.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from plot_presentation_figures import TRUTH_TABLE

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
B0 = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m24"
GAL = ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25", "t7",
       "au5", "au12", "au14", "au20", "au22", "au28"]
EXTRA = [("t7, OSARG coverage + N", "test2_t7_sf_12p5pct_osarg_cover", "sfratio_osarg_cover", "t7"),
         ("t7 in Hunter+24", "test2_t7_sf_12p5pct_mwH24_osargdesign", "sfratio", "t7"),
         ("steady-state mock, Hunter+24", "mw_hunter24_mock_from_osarg_ics", "sfratio", "h24mock"),
         ("OSARG, Hunter+24", "osarg_zhang2024_mw_hunter24", "sfratio", None)]
KEYS = ["Omega_abs", "phi_b_deg", "S", "bar_frac"]


def truth(g, tab):
    if g is None:
        return None
    if g == "h24mock":
        return dict(Omega_abs=32.5, phi_b_deg=25.0, S=1.0, bar_frac=1.0)
    if g.startswith("au"):
        om = float(json.load(open(RESULTS_DIR / f"test2_{g}_sf_12p5pct" / "run_meta.json"))["omega_truth"])
    else:
        t = float(g[1:])
        om = float(tab.loc[np.isclose(tab.snap_time, t), "pattern_speed_fd"].iloc[0])
    return dict(Omega_abs=om, phi_b_deg=27.0, S=1.0, bar_frac=1.0)


def post(run, scheme, var):
    p = RESULTS_DIR / run / f"emulator_results_logLs_{scheme}_{var}.json"
    return json.load(open(p))["posterior"] if p.exists() else None


def main():
    set_style()
    tab = pd.read_csv(TRUTH_TABLE)
    cases = [(g, f"test2_{g}_sf_12p5pct", "sfratio", g) for g in GAL] + EXTRA
    rows = []
    for lab, run, var, g in cases:
        tr = truth(g, tab)
        for scheme, tag in ((B, "scheme B"), (B0, "without vphi m=0")):
            p = post(run, scheme, var)
            if p is None:
                print(f"missing: {run} {tag}")
                continue
            r = dict(case=lab, run=run, scheme=tag, mock=lab in GAL)
            for k in KEYS:
                q = p[k]
                r[f"{k}_p50"], r[f"{k}_p16"], r[f"{k}_p84"] = q["p50"], q["p16"], q["p84"]
                if tr is not None:
                    t = tr[k]
                    r[f"{k}_truth"] = t
                    r[f"{k}_bias"] = q["p50"] - t
                    r[f"{k}_inside"] = q["p16"] <= t <= q["p84"]
            r["Omega_width_pct"] = 100 * (r["Omega_abs_p84"] - r["Omega_abs_p16"]) / r["Omega_abs_p50"]
            if tr is not None:
                r["Omega_bias_pct"] = 100 * r["Omega_abs_bias"] / tr["Omega_abs"]
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "without_vphi_m0_summary.csv", index=False)
    m = t[t.mock]
    print("\n18 mock galaxies (12 N-body + 6 Auriga), median width / mean |bias| / truth inside 16-84%:")
    for tag, s in m.groupby("scheme", sort=False):
        line = [f"Omega {s.Omega_width_pct.median():.1f}% / {s.Omega_bias_pct.abs().mean():.2f}% / {s.Omega_abs_inside.mean():.2f}"]
        for k in ("S", "bar_frac", "phi_b_deg"):
            w = (s[f"{k}_p84"] - s[f"{k}_p16"]).median()
            line.append(f"{k} {w:.3f} / {s[f'{k}_bias'].abs().mean():.3f} (mean {s[f'{k}_bias'].mean():+.3f}) / {s[f'{k}_inside'].mean():.2f}")
        print(f"  [{tag:17s}] " + "; ".join(line))
    print("\nper run: |Omega_b| median (16-84%) [bias %] | S | bar_frac | phi_b, scheme B -> without vphi m=0")
    for lab in t.case.unique():
        a = t[(t.case == lab) & (t.scheme == "scheme B")]
        b = t[(t.case == lab) & (t.scheme == "without vphi m=0")]
        if a.empty or b.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        f = lambda r, k: f"{r[k + '_p50']:.2f} ({r[k + '_p16']:.2f}-{r[k + '_p84']:.2f})"
        bias = lambda r: f" [{r['Omega_bias_pct']:+.1f}%]" if "Omega_bias_pct" in r and pd.notna(r.get("Omega_bias_pct")) else ""
        print(f"  {lab:30s} Om {f(a, 'Omega_abs')}{bias(a)} -> {f(b, 'Omega_abs')}{bias(b)} | S {a.S_p50:.2f} -> {b.S_p50:.2f} "
              f"| f {a.bar_frac_p50:.2f} -> {b.bar_frac_p50:.2f} | phi {a.phi_b_deg_p50:.1f} -> {b.phi_b_deg_p50:.1f}")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw=dict(wspace=0.28))
    for tag, col, off, mk in (("scheme B", INK_MUTED, -0.25, "o"), ("without vphi m=0", ACCENT, 0.25, "s")):
        s = m[m.scheme == tag]
        x = s.Omega_abs_truth
        axes[0].errorbar(x + off, s.Omega_bias_pct, yerr=[100 * (s.Omega_abs_p50 - s.Omega_abs_p16) / x,
                                                         100 * (s.Omega_abs_p84 - s.Omega_abs_p50) / x],
                         fmt=mk, color=col, ms=5, capsize=0, lw=1.2, label=tag)
        for ax, k in ((axes[1], "S"), (axes[2], "bar_frac")):
            ax.errorbar(x + off, s[f"{k}_p50"], yerr=[s[f"{k}_p50"] - s[f"{k}_p16"], s[f"{k}_p84"] - s[f"{k}_p50"]],
                        fmt=mk, color=col, ms=5, capsize=0, lw=1.2, label=tag)
    axes[0].axhline(0, color=INK, lw=0.8, ls="--")
    for ax in axes[1:]:
        ax.axhline(1, color=INK, lw=0.8, ls="--")
    axes[0].set(xlabel=r"true $|\Omega_b|$", ylabel=r"$|\Omega_b|$ bias [%] (median, 16–84%)")
    axes[1].set(xlabel=r"true $|\Omega_b|$", ylabel="bar length S (truth 1)")
    axes[2].set(xlabel=r"true $|\Omega_b|$", ylabel="bar strength bar_frac (truth 1)")
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle("18 mock galaxies at 12.5% distance error: scheme B against scheme B without the ⟨v_φ⟩ m = 0 term "
                 "(points offset ±0.25 for legibility)", color=INK, fontsize=11)
    savefig(fig, "test2_pattern_speed_and_bar_recovery_without_vphi_m0")


if __name__ == "__main__":
    main()

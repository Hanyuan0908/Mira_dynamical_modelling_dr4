"""Recovered against true pattern speed for the Test 2 searches.

One point per (run, library variant): the posterior median of |Omega_b| from the GP emulator, with the 16-84%
range as the error bar, against the galaxy's measured pattern speed. The three variants differ only in how the
orbit library is weighted and whether it is observed with the data's distance error:

  asis      library as integrated, weighted 1/tau (tau = time fraction at R < R_SEL); ignores the smearing
  observed  every library sample re-observed with the same sigma_DM, same 1/tau weights
  sf        re-observed and weighted S(x_j)/tau_i with tau_i = <S> along the orbit (selection.py)

Left: recovered against true, with the 1:1 line. Right: the fractional bias, which is what the left panel is
hard to read at this dynamic range.

Usage: python plot_test2_pattern_speed_recovery.py [--runs a b c] [--tag NAME]
Figure: figures/test2_pattern_speed_recovery<tag>.png
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt

from common import ROOT, RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from nbody import SNAPSHOT_DIR

SCHEME_A = "fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"
SCHEME_B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
RUNS = ["test2_t2_sf_10pct", "test2_t3_sf_10pct", "test2_t4_sf_10pct", "test2_t5.5_sf_10pct",
        "test2_t7_sf_10pct", "test2_t7_nosf_10pct"]
VARIANTS = ["asis", "observed", "sfratio", "sfemp", "sf"]
STYLE = {"asis": (BLUE_RAMP[5], "o", "no selection weighting"),
         "observed": (ACCENT, "s", "smeared, no selection"),
         "sfratio": ("#b07aa1", "^", r"S estimated: $\rho_{\rm obs}/\rho_{\rm lib}(l,b,d)$"),
         "sfemp": ("#d4a017", "v", "S estimated: separable fit"),
         "sf": ("#1baf7a", "D", "S known exactly (mock header)")}


TABLE = SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv"   # ships with the snapshots, in mock/Nbody


def tabulated(t_model, col="pattern_speed_dehnen"):
    """The simulation directory's own pattern speed at this snapshot time."""
    import pandas as _pd
    df = _pd.read_csv(TABLE)
    hit = df[np.isclose(df.snap_time, float(t_model))]
    return float(hit[col].iloc[0]) if len(hit) else np.nan


def collect(runs, scheme, inflate=False, truth="measured"):
    """(run, variant) -> truth, median, p16, p84, and whether the mock has a selection function."""
    out = []
    for r in runs:
        meta_p = RESULTS_DIR / r / "run_meta.json"
        if not meta_p.exists():
            print(f"  skipping {r}: no run_meta.json")
            continue
        meta = json.load(open(meta_p))
        for f in sorted((RESULTS_DIR / r).glob("emulator_results_*.json")):
            d = json.load(open(f))
            col = d["column"][len("logLs_"):]
            if not col.startswith(scheme + "_"):
                continue
            v = col[len(scheme) + 1:]
            if v not in VARIANTS:
                continue
            p = d["posterior"]["Omega_abs"]
            f = 1.0
            if inflate:
                import pandas as _pd
                ev = _pd.read_csv(RESULTS_DIR / r / "evals.csv")
                col = f"chi2N_{scheme}_{v}"
                if col in ev:
                    f = float(np.sqrt(ev[col].min()))
            gal = meta.get("galaxy", "nbody")
            auriga = str(gal).startswith("au")
            om = float(meta["omega_truth"])
            if truth != "measured" and not auriga:
                col = "pattern_speed_fd" if truth == "tabulated_fd" else "pattern_speed_dehnen"
                om = tabulated(7.0 if gal == "nbody" else float(gal.split("_t")[1]), col)
            # The Auriga haloes are separate z = 0 galaxies, not snapshots of the N-body run, so
            # there is no simulation table to read a finite-difference pattern speed from: their
            # truth is our own Dehnen-method measurement (results/auN_ground_truth.json, which is
            # what omega_truth holds). That is a DIFFERENT truth convention from the N-body points
            # on the same axes, and the convention is not cosmetic -- PROGRESS 2026-11-02/03 shows
            # the sf/sfratio ordering reversing between them. Flagged in the legend and the title.
            out.append(dict(run=r, variant=v, truth=om, med=p["p50"],
                            lo=p["p50"] - f * (p["p50"] - p["p16"]), hi=p["p50"] + f * (p["p84"] - p["p50"]),
                            sf="sf" in meta["variants"], loo=d["loo_rms"], inflation=f,
                            auriga=auriga, galaxy=gal))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=RUNS)
    ap.add_argument("--tag", default="")
    ap.add_argument("--truth", default="tabulated_fd", choices=["measured", "tabulated", "tabulated_fd"],
                    help="which Omega_b to call the truth: our own measurement (run_meta.json), or the "
                         "simulation directory's own table (pattern_speed_dehnen / _fd columns)")
    ap.add_argument("--inflate", action="store_true",
                    help="scale the error bars by sqrt(chi2_min/N), the usual correction for a model that "
                         "does not fit within its errors; the tempering assumes chi2/N = 1")
    ap.add_argument("--scheme", default=SCHEME_B, help="chi2 scheme to plot (default: scheme B, no sigma m=0)")
    ap.add_argument("--variants", nargs="+", default=None,
                    help="library weightings to show (default: all). e.g. --variants sfratio")
    args = ap.parse_args()
    set_style()
    if args.variants:
        VARIANTS[:] = [v for v in VARIANTS if v in args.variants]
    pts = collect(args.runs, args.scheme, args.inflate, args.truth)
    if not pts:
        raise SystemExit("no emulator results found")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.2, 5.4), gridspec_kw=dict(wspace=0.26))
    lim = (0.92 * min(p["truth"] for p in pts), 1.14 * max(p["truth"] for p in pts))
    a1.plot(lim, lim, color=INK, lw=1.1, ls="--", zorder=1, label="perfect recovery")
    a2.axhline(0, color=INK, lw=1.1, ls="--", zorder=1)

    # small horizontal offsets so the two t = 7 runs and the three variants stay readable
    off = {"asis": -0.70, "observed": -0.35, "sfratio": 0.0, "sfemp": 0.35, "sf": 0.70}
    for v in VARIANTS:
        sel = [p for p in pts if p["variant"] == v]
        if not sel:
            continue
        c, mk, lab = STYLE[v]
        for ax, yfun, efun in ((a1, lambda p: p["med"], lambda p: (p["med"] - p["lo"], p["hi"] - p["med"])),
                               (a2, lambda p: 100 * (p["med"] / p["truth"] - 1),
                                lambda p: (100 * (p["med"] - p["lo"]) / p["truth"],
                                           100 * (p["hi"] - p["med"]) / p["truth"]))):
            for p in sel:
                x = p["truth"] + off[v] + (0.9 if not p["sf"] else 0.0)
                lo, hi = efun(p)
                # Auriga: same colour and variant marker, but a dark edge and a larger symbol, so
                # the two suites are separable without spending the marker shape on it.
                ax.errorbar(x, yfun(p), yerr=[[lo], [hi]], fmt=mk, color=c,
                            ms=9 if p.get("auriga") else 7,
                            mfc=c if p["sf"] else "white",
                            mec=INK if p.get("auriga") else c,
                            mew=1.9 if p.get("auriga") else 1.4,
                            lw=1.3, capsize=3, zorder=4 if p.get("auriga") else 3)
        a1.plot([], [], mk, color=c, label=lab)

    for p in pts:
        if not p["sf"] and p["variant"] == "asis":
            a1.annotate("t = 7, no selection\n(open symbols)", (p["truth"] + 0.9, p["med"]),
                        xytext=(14, -26), textcoords="offset points", fontsize=8, color=INK_2,
                        arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.7))
    a1.set_xlim(*lim)
    a1.set_ylim(*lim)
    a1.set_aspect("equal")
    a1.set_xlabel(r"true $|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]")
    a1.set_ylabel(r"recovered $|\Omega_b|$ (posterior median, 16-84%)")
    if any(p.get("auriga") for p in pts):
        a1.plot([], [], "o", color="none", mec=INK, mew=1.9, ms=9,
                label="Auriga z = 0 (dark edge)")
        a1.plot([], [], "o", color="none", mec=INK_MUTED, mew=1.4, ms=7,
                label="N-body snapshots")
    a1.legend(loc="upper left", fontsize=8.5)
    a1.set_title("recovered against true pattern speed", color=INK_2, fontsize=10)

    a2.set_xlim(*lim)
    a2.set_xlabel(r"true $|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]")
    a2.set_ylabel(r"bias in $|\Omega_b|$ [%]")
    a2.set_title("fractional bias (filled = selection mock, open = no selection)", color=INK_2, fontsize=10)
    for v in VARIANTS:
        s = [p for p in pts if p["variant"] == v]
        if not s:
            continue
        bias = lambda q: np.mean([abs(100 * (x["med"] / x["truth"] - 1)) for x in q])
        nb = [x for x in s if not x.get("auriga")]
        au = [x for x in s if x.get("auriga")]
        lab = f"{v}: mean |bias| {bias(s):.1f}%"
        if nb and au:        # the two suites do not share a truth convention; quote them apart
            lab = f"{v}: N-body {bias(nb):.1f}%, Auriga {bias(au):.1f}%"
        a2.plot([], [], STYLE[v][1], color=STYLE[v][0], label=lab)
    a2.legend(loc="lower right", fontsize=8.5)

    sch = "scheme B: sigma ring averages dropped" if args.scheme == SCHEME_B else \
          ("scheme A: all m = 0, 2, 4" if args.scheme == SCHEME_A else args.scheme)
    extra = r" | error bars $\times\sqrt{\chi^2_{\min}/N}$" if args.inflate else ""
    extra += {"measured": "", "tabulated": " | truth = simulation table (Dehnen col.)",
              "tabulated_fd": " | truth = simulation table (finite-difference col.)"}[args.truth]
    if any(p.get("auriga") for p in pts) and args.truth != "measured":
        extra += "; Auriga truth = own Dehnen measurement"
    sig = sorted({round(float(json.load(open(RESULTS_DIR / r / "run_meta.json"))["sigma_DM"])
                        * np.log(10) / 5 * 100, 1) for r in args.runs
                  if (RESULTS_DIR / r / "run_meta.json").exists()})
    lvl = f"{sig[0]:g}%" if len(sig) == 1 else "/".join(f"{v:g}%" for v in sig)
    fig.suptitle(f"Test 2: pattern-speed recovery at {lvl} distance error - {sch}{extra}",
                 fontsize=11.5, y=0.98)
    savefig(fig, f"test2_pattern_speed_recovery{args.tag}")

    print(f"\n{'run':22s} {'variant':9s} {'truth':>7s} {'median':>8s} {'bias':>7s} {'16-84%':>15s}")
    for p in sorted(pts, key=lambda q: (q["truth"], q["variant"])):
        print(f"{p['run']:22s} {p['variant']:9s} {p['truth']:7.2f} {p['med']:8.2f} "
              f"{100 * (p['med'] / p['truth'] - 1):+6.1f}% {p['lo']:7.2f}-{p['hi']:.2f}")


if __name__ == "__main__":
    main()

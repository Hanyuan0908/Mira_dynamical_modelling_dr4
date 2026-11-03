"""Sort figures/ into a readable layout: headline figures on top, everything else filed by purpose.

figures/ had grown to ~390 loose PNGs, dominated by per-design-point emulator diagnostics (91 slices,
91 leave-one-out, 69 corners), which buried the dozen figures anyone actually looks at. This moves the
rest into `figures/archive/<theme>/` and leaves the headline set at the top level.

Nothing is deleted and nothing is renamed, so every figure is still found by its own name; only the
directory changes. Re-running a plotting script writes to the top level again (common.savefig always
writes to FIG_DIR), so run this afterwards to re-file it - it is idempotent.

  python organise_figures.py --dry-run      # list what would move
  python organise_figures.py                # move
  python organise_figures.py --undo         # put everything back at the top level
"""
import argparse
import re
import shutil

from common import FIG_DIR, rel

KEEP_DIRS = {"key_results", "key_results_12p5pct", "test2_schemeB", "test2_schemeB_selection",
             "mock_validation", "setup_checks", "archive"}

# Headline figures, kept at the top level. Exact names.
HEADLINE = [
    # campaign summaries
    "test2_error_level_comparison.png",
    "test2_error_level_comparison_sfratio.png",
    "test2_pattern_speed_recovery_schemeB_fdtruth.png",
    "test2_pattern_speed_recovery_12p5pct_schemeB.png",
    "test2_snapshot_coverage.png",
    # which snapshot, and what the mocks look like
    "nbody_pattern_speed_evolution.png",
    "test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.217.png",
    "test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.271.png",
    # the bar-angle result
    "test2_bar_angle_vecsum_drift.png",
    "test2_bar_angle_vecsum_transfer.png",
    # why bar length and strength are weak
    "test1_nbody_qmc_focused_bar_length_strength_degeneracy_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.png",
]

# (archive subfolder, regex). First match wins, so order matters.
RULES = [
    ("emulator_diagnostics", r"emulator_(slices|loo)_"),
    ("emulator_diagnostics", r"logl_corner"),
    ("posteriors",           r"emulator_corner"),
    ("../mock_validation",   r"mock_selection_vs_dr3_miras|spatial_mock_vs_dr3|mock_spatial|mock_vR|"
                             r"_overview|dr3_mira_|mock_grid"),
    ("bar_angle",            r"bar_angle|estimator_(scan|transfer)|vecsum"),
    ("au18",                 r"^au18_"),
    ("../setup_checks",      r"_check\.png$|potential_density|rotation_curve|force_check|bar_fourier|"
                             r"potential_family|stellar_density|pattern_speed\.png$|density_stationarity"),
    ("test1",                r"^test1|^test2b_|^mock_perfect"),
    ("pattern_speed_variants", r"^test2_pattern_speed_recovery"),
    ("../setup_checks",      r"^nbody_|^au18_"),
    ("test2_other",          r"^test2_"),
]


def classify(name):
    for sub, pat in RULES:
        if re.search(pat, name):
            return sub
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--undo", action="store_true", help="move every archived figure back to the top level")
    args = ap.parse_args()
    arch = FIG_DIR / "archive"

    if args.undo:
        n = 0
        for f in sorted(arch.rglob("*.png")):
            dest = FIG_DIR / f.name
            print(f"  {rel(f)} -> {rel(dest)}")
            if not args.dry_run:
                shutil.move(str(f), str(dest)); n += 1
        print(f"restored {n} figures")
        return

    moves, kept = [], []
    for f in sorted(FIG_DIR.glob("*.png")):
        if f.name in HEADLINE:
            kept.append(f.name)
            continue
        moves.append((f, (arch / classify(f.name) / f.name).resolve()))

    by_sub = {}
    for _, dest in moves:
        by_sub[dest.parent.name] = by_sub.get(dest.parent.name, 0) + 1
    print(f"top level: {len(kept)} headline figures kept, {len(moves)} to file\n")
    for sub, n in sorted(by_sub.items(), key=lambda kv: -kv[1]):
        print(f"  {sub:32s} {n:4d}")
    missing = [h for h in HEADLINE if not (FIG_DIR / h).exists()]
    if missing:
        print("\n  headline figures not present (not an error, just not made yet):")
        for m in missing:
            print(f"    {m}")
    if args.dry_run:
        print("\n(dry run; nothing moved)")
        return
    for src, dest in moves:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dest))
    print(f"\nmoved {len(moves)} figures into {rel(arch)}")


if __name__ == "__main__":
    main()

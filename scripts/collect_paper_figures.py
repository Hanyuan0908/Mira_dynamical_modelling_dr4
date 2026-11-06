"""Maintain figures/paper_figures/ - the results worth publishing, copied from across figures/.

Everything in paper_figures/ is a copy, so the originals stay where their scripts write them and
nothing breaks if the folder is reorganised. Re-running rebuilds the standing set; --add puts a new
figure in without touching the rest, which is how figures get added as they are judged worth keeping.

  python collect_paper_figures.py                                  # rebuild the standing set
  python collect_paper_figures.py --add figures/foo.png            # add one, keeping its name
  python collect_paper_figures.py --add figures/foo.png --as 06_bar_angle --section test2_12p5pct
  python collect_paper_figures.py --list                           # what is in there now
"""
import argparse
import shutil

from common import FIG_DIR, rel

DEST = FIG_DIR / "paper_figures"
A = FIG_DIR / "archive"

# (section, destination name, source) - source relative to figures/
STANDING = [
    ("test1_perfect_mock", "01_posterior_corner_schemeA",
     "archive/posteriors/test1_nbody_qmc_focused_emulator_corner_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.png"),
    ("test1_perfect_mock", "02_posterior_corner_schemeB",
     "archive/posteriors/test1_nbody_qmc_focused_emulator_corner_logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024.png"),
    ("test1_perfect_mock", "03_fourier_coefficients_data_vs_models_absz_0to1",
     "archive/test1/test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_absz_0to1.png"),
    ("test1_perfect_mock", "03_fourier_coefficients_data_vs_models_absz_1to2",
     "archive/test1/test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_absz_1to2.png"),
    ("test2_12p5pct", "01_pattern_speed_recovery_12_snapshots",
     "test2_pattern_speed_recovery_12p5pct_12pt.png"),
    ("test2_12p5pct", "02_pattern_speed_10pct_vs_12p5pct",
     "test2_error_level_comparison_sfratio.png"),
    ("test2_12p5pct", "03_mock_vs_dr3_miras",
     "test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.271.png"),
    ("test2_12p5pct", "04_all_mocks_vs_dr3_miras",
     "mock_validation/test2_mock_grid_12.5pct.png"),
    ("test2_12p5pct", "05_snapshot_coverage",
     "test2_snapshot_coverage.png"),
    # --- Auriga: the 7 galaxies of batch 1, scored in the adopted configuration (2026-11-05) ---
    ("test2_12p5pct", "06_pattern_speed_recovery_12pt_plus_13_auriga",
     "test2_pattern_speed_recovery_12p5pct_12pt_plus_auriga.png"),
    # --- bar angle. The calibration is a BAND because the compression depends on a bar property
    # that is still unidentified (PROGRESS 2026-11-04/05); the N-body set is held out as a blind
    # check of it. The amplitude-correction and 13-galaxy band figures are deliberately NOT here:
    # both were retracted when the bar-strength dependence failed on the 25-galaxy set.
    ("bar_angle", "01_auriga_calibration_band_nbody_blind",
     "test2_bar_angle_auriga_calibration.png"),
    ("bar_angle", "02_kinematic_window_comparison",
     "test2_bar_angle_window_comparison.png"),
    ("bar_angle", "03_transfer_25_galaxies",
     "test2_bar_angle_transfer_12p5.png"),
    ("bar_angle", "04_outlier_velocity_maps_au5_au20",
     "test2_smeared_velocity_maps_auriga_outliers.png"),
    ("bar_angle", "05_bias_vs_distance_error_au18_self03",
     "test2_bar_angle_bias_vs_distance_error.png"),
    ("bar_angle", "06_auriga_band_self03",
     "test2_bar_angle_auriga_band.png"),
    ("bar_angle", "07_offset_vs_bar_property",
     "test2_bar_angle_offset_vs_bar_property.png"),
    ("test2_12p5pct", "07_bar_strength_profile_recovery",
     "test2_bar_strength_profile_recovery.png"),
    ("emulator_validation", "01_leave_one_out", "emulator_validation/loo_sfratio.png"),
    ("emulator_validation", "02_pooled_draws_K10",
     "emulator_validation/pooled_draws_test2_t7_sf_12p5pct.png"),
]
# Known-S (`sf`) and the separable fit (`sfemp`) are deliberately NOT here. The selection function a
# real survey can build is the density ratio rho_obs/rho_lib, so that is the only one these figures
# show; the others live in figures/key_results_12p5pct/ and figures/test2_schemeB_selection/.

# (section, destination subfolder, source folder)
# corner plots are gathered by a glob over figures/, since the emulator writes them at the top level
STANDING_GLOBS = [
    ("test2_12p5pct", "corner_sfratio",
     "test2_{t}_sf_12p5pct_emulator_corner_logLs_"
     "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sfratio.png",
     ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25", "t7"]),
    ("test2_12p5pct", "corner_sfratio_auriga",
     "test2_{t}_sf_12p5pct_emulator_corner_logLs_"
     "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sfratio.png",
     ["au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24",
      "au26", "au27", "au28"]),
    ("emulator_validation", "loo_per_run_auriga",
     "test2_{t}_sf_12p5pct_emulator_loo_logLs_"
     "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sfratio.png",
     ["au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24",
      "au26", "au27", "au28"]),
    # The kinematic prior on phi_b: the same likelihood with phi_b pinned by the <|v_R|/v_tot>
    # calibration. Blue = flat prior, orange = with the prior, grey = truth. For an Auriga galaxy
    # the prior is leave-one-out, since the calibration is built from those 13.
    ("test2_12p5pct", "corner_kinematic_prior",
     "test2_{t}_sf_12p5pct_corner_kinprior_logLs_"
     "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sfratio.png",
     ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25",
      "t7", "au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24",
      "au26", "au27", "au28"]),
    ("emulator_validation", "loo_per_run",
     "test2_{t}_sf_12p5pct_emulator_loo_logLs_"
     "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sfratio.png",
     ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25", "t7"]),
]
STANDING_DIRS = []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--add", metavar="PATH", help="a figure to add (path, or a name to search for)")
    ap.add_argument("--as", dest="name", help="destination name (default: the source's own name)")
    ap.add_argument("--section", default="test2_12p5pct", help="subfolder of paper_figures/")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for f in sorted(DEST.rglob("*.png")):
            print(f"  {f.relative_to(DEST)}")
        return

    if args.add:
        src = FIG_DIR.parent / args.add
        if not src.exists():
            src = next(iter(sorted(FIG_DIR.rglob(args.add if args.add.endswith(".png")
                                                 else f"*{args.add}*.png"))), None)
        if src is None or not src.exists():
            raise SystemExit(f"not found: {args.add}")
        d = DEST / args.section
        d.mkdir(parents=True, exist_ok=True)
        out = d / ((args.name + ".png") if args.name else src.name)
        shutil.copy2(src, out)
        print(f"added {rel(src)} -> {rel(out)}")
        return

    n = 0
    for sec, name, src in STANDING:
        s = FIG_DIR / src
        if not s.exists():
            print(f"  missing: {src}")
            continue
        (DEST / sec).mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, DEST / sec / f"{name}.png")
        n += 1
    for sec, sub, tmpl, tags in STANDING_GLOBS:
        (DEST / sec / sub).mkdir(parents=True, exist_ok=True)
        for t in tags:
            hits = sorted(FIG_DIR.rglob(tmpl.format(t=t)))
            if not hits:
                print(f"  missing: {tmpl.format(t=t)}")
                continue
            shutil.copy2(hits[0], DEST / sec / sub / f"{t}.png")
            n += 1
    print(f"paper_figures/: {n} figures in place")


if __name__ == "__main__":
    main()

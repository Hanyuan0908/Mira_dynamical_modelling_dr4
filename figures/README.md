# figures/

Headline figures sit at the top level; everything else is filed under `archive/` by purpose.
Nothing has been renamed or deleted - only moved - so any figure named in `PROGRESS.md` is still
found by its own name, one level down.

`scripts/organise_figures.py` maintains this. The scripts that produce one figure per galaxy now
write directly into `mock_validation/` or `setup_checks/` (`savefig` takes a "/" in the name), so
they no longer pile up at the top level. Anything else that lands at the top level is re-filed by
running the organiser; it is idempotent, and `--undo` flattens everything back.

## Top level - the figures worth looking at first

| file | what it shows |
|---|---|
| `test2_error_level_comparison.png` | **pattern speed recovery at 10% against 12.5% distance error**, S known exactly |
| `test2_error_level_comparison_sfratio.png` | the same with S *estimated from the data* - the realistic case |
| `test2_pattern_speed_recovery_schemeB_fdtruth.png` | the 10% campaign's recovery figure |
| `test2_pattern_speed_recovery_12p5pct_schemeB.png` | the same for 12.5% |
| `test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.217.png` | the 10% mock against the DR3 Miras (sky, distance, R, face-on) |
| `test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.271.png` | the same for the 12.5% mock |
| `nbody_pattern_speed_evolution.png` | Omega_b against snapshot time - how the five test galaxies were chosen |
| `test2_bar_angle_vecsum_drift.png` | bar-angle estimator: distance-error drift, amplitude-weighted rings |
| `test2_bar_angle_vecsum_transfer.png` | whether that calibration transfers between snapshots |
| `test2_snapshot_coverage.png` | where the test galaxies sit on the run's pattern-speed track |
| `test1_nbody_qmc_focused_bar_length_strength_degeneracy_*.png` | why S and bar_frac are weakly constrained |

## Collection folders - scripts write straight into these, never to the top level

| folder | what |
|---|---|
| `mock_validation/` | every mock against the DR3 Miras: `make_mock_selection.py`, `plot_mock_vs_dr3_mira_spatial.py`, `plot_mock_overview.py`. One figure per galaxy per error level, so this grows with every new snapshot - which is why it has its own folder. |
| `setup_checks/` | per-galaxy potential, density, rotation-curve, force and bar-Fourier checks: `prepare_*_mwframe.py`, `check_au18_potential.py`, `test_potential_family.py` |

## Animations

| file | what it shows |
|---|---|
| `animations/vR_quadrupole_correct_vs_wrong_pattern_speed.{mp4,gif}` | the central idea: bar-frame ⟨v_R⟩ quadrupole of the 50k-star mock integrated for 1 Gyr in the correct potential (left) and with a 27% too fast pattern speed (right) |
| `animations/vR_quadrupole_correct_vs_wrong_bar_strength.{mp4,gif}` | the same with the bar at half strength (right) |
| `animations/*_t{0000,0500,1000}Myr.png` | stills at 0, 0.5 and 1 Gyr |

Made by `scripts/make_animation_vR_quadrupole_correct_vs_wrong_potential.py --wrong pattern_speed|bar_strength`.

## Campaign deliverables

| folder | what |
|---|---|
| `key_results/` | the 10% campaign: pattern speed, bar angle, estimator comparison, with a README |
| `key_results_12p5pct/` | the 12.5% campaign: summary figures and corner plots per snapshot, with a README |
| `test2_schemeB/` | scheme A against scheme B, the comparison that settled on B |
| `test2_schemeB_selection/` | how S is obtained (none / ratio / separable / true), 20 corner plots |

## archive/

| folder | n | what |
|---|---|---|
| `emulator_diagnostics/` | 213 | GP emulator 1D slices, leave-one-out checks, and the evaluated-point logL corners - one set per run and column. Progress checks, not results. |
| `posteriors/` | 69 | every per-run emulator corner plot. The ones that matter are copied into the `key_results*` folders. |
| `test1/` | 27 | Test 1b: square and polar maps, Fourier coefficients and compression, the per-star fit check, Test 2b distance-error studies |
| `bar_angle/` | 10 | bar-angle estimators: inversion curves, calibration universality, transfer tests |
| `au18/` | 8 | the Auriga 18 galaxy, superseded by the isolated N-body run on 2026-10-27 |
| `pattern_speed_variants/` | 7 | earlier recovery figures: scheme A, inflated errors, and the other truth conventions |

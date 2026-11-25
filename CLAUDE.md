# Working rules for this project

These rules were set by the user on 2026-11-25. They apply to every session and every agent working in this
repository.

## 1. Evidence: never claim anything without showing it

- **Never lie.** Every statement of fact or result comes with the numbers or the figure that support it, shown to the
  user in the same reply. Give the numbers themselves, not just a pointer to a file, and give figures as
  repo-relative paths.
- **Guesses and hypotheses are allowed, but label them** ("Hypothesis, untested: …"), and never let them stand in for
  results.
- **Do not make up mechanisms or arguments to explain a result.** If the cause is not established, say so and propose
  the test that would establish it.
- **When a claim turns out wrong, withdraw it explicitly.** Do this in the reply and in `PROGRESS.md`; mark the entry
  as WITHDRAWN rather than deleting it.

## 2. Bookkeeping: PROGRESS.md is the record

- **Log every test, check or decision in `PROGRESS.md`** as a dated entry. Record:
  - what was asked;
  - the scripts that reproduce it, run from `scripts/` with the `astro312` environment;
  - the figures, with full repo-relative paths;
  - the key numbers.
- Keep the figure index and "Next steps → Current status and next steps" up to date.
- Then commit and push.
- Record user decisions as "Decision (user)". Examples: the Sormani+22 Milky Way family, and the temperature T1.

## 3. Catch-up after every context compaction (and when resuming the project)

Before doing anything else:

1. **Re-read every user message,** including corrections and decisions:
   `python scripts/extract_user_messages.py <session transcript .jsonl>`. Transcripts are in
   `~/.claude/projects/<project>/`.
2. **Read `PROGRESS.md`:**
   - "Conventions";
   - the 2026-11-05 entry "THE PIPELINE, written down";
   - every entry since the previous compaction;
   - "Next steps → Current status and next steps".
3. **Read this file and the memory index** with the feedback files it points to.
4. **Read the scripts:** the docstring summaries, newest first, and the full source of any script about to be run or
   edited.
5. **Tell the user briefly** what was re-read and where the project stands, with numbers.

A local SessionStart hook (`.claude/settings.local.json`, matcher `compact`) injects this reminder automatically. It
is not in git.

## 4. Deliver what was asked

- **Give the user the exact figure they asked for.** A "corner plot" is a full corner plot, one per case if several
  are asked for, not a set of 1D marginals. Use square-pixel maps, not polar ones.
- **Give figures descriptive file names.** No codes or abbreviations.
- **Do not launch heavy runs nobody asked for.** Propose them, with an estimated time, and let the user decide.

## 5. Standing technical decisions (see PROGRESS.md for the evidence)

- **Environment:** `astro312` conda. Watch memory to avoid running out of it. Wrap long runs in `caffeinate -i -s`.
  In zsh, never pass CLI options through one unquoted variable.
- **Trial potentials:** the DM halo is spherical in both the truth and the trial potentials. The Milky Way family is
  **Sormani+22** (galaxy `sormani22`; decision of 2026-11-24, restored on 2026-11-25 after a brief switch to
  Hunter+24, which is kept for comparison only). The OSARG prior box is being widened to bar_frac 0.1–2.0 and
  S 0.4–1.7 (run `osarg_zhang2024_sormani22_widebar`).
- **Bulky outputs** of new runs go on the portable drive:
  `/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs/` (`MIRA_DR4_RUNS`). The drive
  is exFAT, so no symlinks can be made on it.
- **Likelihood:**
  - scheme B (`fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024`);
  - the `sfratio` selection;
  - K = 10 distance draws;
  - libraries sampled every 50 Myr over 0.25–1 Gyr.
- **Temperature:** T1 = σ₁, the std of χ² − χ²(best run) over tiny random perturbations at the best fit. Measure it
  per Fourier grid (decision of 2026-11-25). Do not use σ₂. Do not use √(2N) as a likelihood. Do not calibrate T on
  the mock suite.
- **Project goals:** measure the bar's pattern speed, and also its strength and length. Priors on the bar shape are
  a fallback, not the answer.

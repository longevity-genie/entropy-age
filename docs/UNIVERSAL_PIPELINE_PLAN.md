# Plan: one universal sample-processing path with best practices built in (2026-11-24)

Anton asked for the code to be more or less universal in how samples are processed, with best practices enforced.
Today the array age-slope series share one runner (`src/entropy_age/experiments/series.py`, one description file per
dataset in `config/series/*.toml`, 22 files); cultures, sequencing and the older array analyses do not.

## Steps

1. **Best-practice guards in the shared runner** (enforced, not documented only). The runner refuses, with a message
   naming the series and the rule: arrays of more than one platform in a group (exists); a group whose arrays fail
   the 90 % detection cut (exists, through `read_qc`); more than one array per donor in a group (exists: first kept,
   now also counted and reported); a CpG filled or imputed on the drift path (exists: complete CpGs only). It adds:
   a `guards` table per series (arrays read, passing detection, kept per donor, platforms, CpGs complete, chips, the
   share of age variance explained by chip) written beside `gaps`; and a flag `age_on_chip` when chip explains more
   than 90 % of age variance, in which case the chip-held fit is reported as unidentifiable rather than given.
2. **Clock step on a light preprocessing path.** A second SeSAMe output per series (`prep` without pOOBAH masking:
   quality mask, channel inference, dye bias, noob), written beside the strict betas. The runner scores the pinned
   biolearn clocks (Horvath 2013, skin & blood, Hannum, PhenoAge, DunedinPACE, Zhang) per array on that path, reports
   per clock and array the share of the clock's absolute weight that is measured, refuses a clock below a set share
   (default 0.9), uses shared-probe clocks on EPIC v2 (only CpGs present on the array; no imputation), and reports age
   acceleration as the residual on age within the series (with the series' held covariates). No absolute clock ages
   are compared across platforms. Mouse and mammalian arrays get no human clocks.
3. **Culture runner.** Same description-file pattern for passaged or arrested cultures: time axis = population
   doublings, passage or days (named in the file); line or donor intercepts; culture resampling within line (the
   existing `culture_bootstrap`); technical replicates collapsed to one array per culture and time point. Migrate
   `bj_telomerase`, `astrocyte_passage` (+ hTERT), `arrest_timecourse` behind it.
4. **Sequencing runner.** For RRBS/WGBS/single-cell count data: coverage filter (default at least 10 reads in every
   library of a group) in place of detection; strands kept separate unless the file says merged; complete CpGs only;
   promoter CpGs from the species' annotation; animal or donor resampling. Migrate `mouse_atlas_rrbs`, `jena_rrbs`,
   `wgbs_nuclei`, `meer_tissues`, `mouse_rrbs_tissues` behind it.
5. **Equivalence tests before retiring old code.** For each migrated analysis, a test (or a recorded run) shows the
   runner reproduces the old module's housekeeping gap within 0.01 SD on the same data and choices, or the difference
   is explained and recorded (e.g. a deliberate fix). The first such check is running: BXD liver through the runner's
   mouse path (R156) against the `bxd_liver` module's -0.30.
6. **Older array analyses** (headline blood/cortex `tissue_contrast`, Huang `mouse_tissues`, `bxd_liver`, GTEx,
   SATSA) move behind the runner only after 1-5; their bespoke parts (Lunnon same-donor contrasts, within-person
   pairs) stay as add-ons that read the runner's matrices.
7. **Separate package later.** Once the runners are proven equivalent, extract them (with the rebuild script and the
   readers) into a small library other projects can reuse; not before, so the code is not frozen while it changes.

## Order and cost

1 and 2 first (they change what every series reports); then 5 on the BXD control; then 3 and 4 with their
migrations and equivalence checks; 6 and 7 later. Each step lands with tests and a CHANGELOG line; reruns go through
Slurm.

## Open questions for the reviewer

- Is "chip explains > 90 % of age variance" the right identifiability rule, or should it be rank-based (age not
  estimable once chip is in the design)?
- Is a 0.9 measured-weight threshold for clocks defensible, and is "shared-probe clock without refitting" the right
  reading of the 450K/EPIC v2 compatibility paper (bioRxiv 10.1101/2024.10.17.618672)?
- Light path: is noob + dye bias without pOOBAH the right minimal SeSAMe recipe for clocks?
- Anything in the plan that would silently change existing results, or a best practice it misses (e.g. sex checks,
  sample swaps, cross-reactive probes, SNP probes)?

## Revised after Codex's review (docs/reviews/2026-11-24-codex-pipeline-plan.md): implement with these changes

1. **First: freeze baselines and fix identity, QC and identifiability.** Keep one physical array key through betas,
   QC and chip (done 2026-11-24: rebuild and QC reader keep the same first array). Fail on zero probe overlap (done).
   Identifiability by design rank: in the full fit and in every resample, age must raise the rank of the nuisance
   design; residual df from rank; > 90 % age variance explained by chip is a warning, reported with residual age
   variation. A per-series `guards` table: arrays read, passing detection, kept per donor, platforms, CpGs complete,
   chips, rank check, attrition by group. Sex inferred from X/Y against metadata and SNP-probe fingerprints for swaps
   and duplicates go into the rebuild script (SeSAMe `inferSex`, SNP betas) and are backfilled per series.
2. **Migrations keep each design.** Cultures (trajectories, lines, hTERT arms, arrest's month-zero reference) and
   sequencing (count-preserving; the RRBS modules' presence filter and fill are reproduced, not replaced by
   complete-only, unless a recorded decision changes them) move behind runners only with a full equivalence check of
   every reported effect, interval and table, not one gap.
3. **Clocks are validated before they are enabled.** Published scoring with intercepts and each clock's
   transformation and missing-input policy; DunedinPACE through its author-supported implementation (EPIC v2
   included); masked vs unmasked detection as a sensitivity path, not a light path that drops pOOBAH; coverage of
   platform-absent, masked and failed inputs reported separately; scores validated on matched DNA across platforms
   where available. The colleagues' "do not drop CpGs" is honoured by scoring with each clock's own policy, not by
   leaving failed probes unmasked.
4. Package extraction last.

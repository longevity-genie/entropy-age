# Plan: what it takes to make entropy-age an Aging Cell submission (2026-11-24)

Anton (Zulip, 2026-11-24): "make a plan on it, review it with codex and execute". Starting point: v0.3 is robust
(raw arrays, donor intervals, positive control, two reviews) but reads as "active promoters stay unmethylated",
with small effects, no consequence for aging, an unresolved mechanism and unexplained heterogeneity across tissues.

## Workstreams, in order

### A. Consequence: does the shield matter for clocks and outcomes?

A1. **Validated clock scoring first** (from the pipeline plan, step 3): published scoring with intercepts and each
clock's transformation and missing-input policy; DunedinPACE through its author-supported code; per array the share
of each clock's weight measured, with platform-absent, masked and failed inputs counted separately. Test against
published example outputs where available.

A2. **Where clocks sit**: rerun `clock_cpgs` and `clock_weights` on current code for blood and cortex (and the new
raw-array tissues): the share of each clock's CpGs and weight in shielded promoters against matched expectation
(same island x level strata), with gene-resampling and donor intervals. Question: do clocks avoid the shield, or
weight shielded CpGs up?

A3. **Incremental value for outcomes (pre-specified)**. One primary outcome per cohort, fixed before running:
EPIC-Italy incident cancer more than 2 years after the draw (rebuilt from IDATs); SJLIFE death by the data freeze
(GSE314261; raw files checked first); Vitality 90+ death (GSE68194, Cox). Models: clock acceleration alone (each
clock), then clock plus the shield share F (and look-alike drift B). Report the change in likelihood (LR test),
AUC/C-index with bootstrap intervals over people, and the F and B coefficients. Covariates as in each module (age,
sex, cell fractions, chip). A fixed-effects meta-analysis of the F and B log-odds/log-hazards across cohorts, with
the heterogeneity statistic. Only the pre-specified outcomes count as primary; the rest are labelled exploratory.

### B. Mechanism: one discriminating test of the division link

B0. **Audit the single-cell brain atlas (GSE292803) design before analysis**: cell types, ages, and whether there
are several animals per cell type and age. If each cell type and age is one pooled sample, donor intervals are
impossible and the atlas can only give cell-type point estimates; then it is reported as such, not as a test.

B1. If replicates exist: per cell type (dividing OPCs and microglia against non-dividing neurons), the housekeeping
gap of the age slope with coverage-matched CpGs and animal resampling, and the dividing-minus-postmitotic contrast
resampling animals together. Prediction of the division hypothesis: larger gap in OPCs/microglia than in neurons
of the same mice. The sorted-nuclei sequencing result (neurons vs oligodendrocyte lineage, no difference) is stated
beside it, whichever way B1 goes.

### C. Explain the heterogeneity across series (meta-regression)

C1. Collect every series' housekeeping gap (relative entropy, chip held where available, all CpGs) with its donor
interval into one table (about 30 human and mouse series). Random-effects meta-analysis of the gap, and
meta-regression on pre-specified moderators: age span covered (years, and as a share of lifespan), median age,
tissue turnover class, platform, and the tissue's mean promoter expression where known. Question: do age span and
turnover explain why breast, myometrium and young mouse liver show none and blood, colon, heart and lung do?

### D. Reframe and clean the text

D1. Headline: drift is constrained by promoter activity across tissues and species; essentiality adds little beyond
expression; the size of the constraint follows age span and turnover (from C). Title proposed to Anton, his call.
D2. Remove the remaining references to the withdrawn study-wide FDR (e.g. the EPIC-Italy passage), rerun the
claims table against the new numbers, regenerate S5, and rebuild figures from the final tables.
D3. A short "limitations" paragraph: mostly blood and brain arrays; GTEx deposited only; mechanism unresolved.

## Order, compute and stops

A1 -> A2 and C1 in parallel (cheap, tables already exist for C1) -> A3 (Slurm; SJLIFE is large) -> B0 -> B1 if the
design allows -> D. Each step lands with tests, a run-log entry and a CHANGELOG line, and an independent review of
new code before its numbers enter the paper. Stop and report to Anton if A3 shows no incremental value (then the
paper is Epigenetics & Chromatin, and the plan's D1 changes) or B0 shows no replicates.

## Questions for the reviewer

- Is the incremental-value design (A3) sound, and what would a reviewer of Aging Cell require beyond it?
- Is a meta-regression across heterogeneous series (C1) defensible, and how should within-study dependence (several
  tissues of one study) be handled?
- Is anything here likely to be read as fishing, and how to pre-register against it?
- What is missing for Aging Cell that this plan does not address?

## B0 result (2026-11-24)

GSE292803 holds 36 cell types x 3 ages (8 weeks, 9 and 18 months), one pooled sample per cell type and age (108
allc files, all fetched). No animal replicates: B1 cannot resample animals. The testable form is across cell types:
the age-slope gap per cell type from three pooled ages, then glial and dividing types (OPC, microglia, two astrocyte
types, oligodendrocytes; 5-6 types) against 30 neuronal types, with class labels permuted across cell types and a
cell-type bootstrap. That is weaker than an animal-level test and is reported as such.

## Revised after Codex's review (docs/reviews/2026-11-24-codex-aging-cell-plan.md): execute with changes

1. **Freeze claims and a protocol.** A dated protocol (docs/PROTOCOL_2026-11-24.md) listing every result already
   inspected and fixing, for each remaining analysis, eligibility, endpoint, lag, clock, score, preprocessing,
   covariates, contrasts, multiplicity and the decision threshold; every locked analysis is reported. Reconcile the
   text: remaining study-wide FDR passages, preprocessing descriptions, interval labels (gene vs donor), S4 and S5.
2. **Donor-level inference for the central results**: the same-sample expression matching (monocytes, fibroblasts,
   epidermis) and the division contrasts (six-donor division vs damage, colon vs hippocampus off the floor, WGBS
   oligodendrocyte vs neuron) get donor/animal-resampling intervals; muscle stem-cell intervals labelled as
   conditional on the six mice. Search for one donor-replicated lineage study (promoter methylation, expression and
   lineage abundance over matched doublings) as the discriminating test of maintenance vs selection.
3. **Harmonised forest plot** of every series' housekeeping gap with provenance; paired contrasts where donors or
   mice are shared (joint resampling); no post-hoc moderator claim. Atlas trajectories descriptive only.
4. **Validate one clock** (published scoring, author implementation where it exists) before any outcome
   reanalysis; outcome models compare clock, then + B, then + F given B on identical observations; EPIC-Italy
   breast >2 y as the explicit follow-up; SJLIFE descriptive and off the critical path; Vitality cohort-specific.
5. **Rebuild figures and text** on the frozen tables; submit on the biological result. Dropped: the claim that age
   span and turnover explain the heterogeneity; the cross-cohort fixed-effect average of outcomes.

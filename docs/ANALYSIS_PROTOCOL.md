# Analysis protocol (living; started 2026-11-24)

The protocol says how each analysis is done. It is not frozen: when a method changes, every result it affects is
recomputed with the new method and the text follows the new numbers (Anton, 2026-11-24: "Do not freeze the protocol,
just recompute old results if we decide to change it and they can be affected"). Two rules keep this honest: every
analysis listed here is reported whatever its result, and a change of method is recorded below with its date, its
reason and the results it moved, so no method is chosen after seeing which one gives a better number.

## Results already inspected (not confirmatory)

Every result in docs/RUN_LOG.md up to R158, including: the blood and cortex gaps; the expression-matched gaps in
monocytes (MESA), cultured fibroblasts (Sturm) and epidermis; the division contrasts (R114 colon vs hippocampus, R136
division vs damage in six donors, R144 WGBS oligodendrocyte vs neuron, R140 muscle stem cells); all array and
sequencing series of R137-R158; the species analysis; EPIC-Italy R96 (deposited and rebuilt), SJLIFE and Vitality 90+
outcome models. None of these is re-labelled as a confirmation.

## Locked analyses

### P1. Donor-level intervals for the expression matching

- Data: MESA monocytes (deposited), Sturm fibroblasts (rebuilt), epidermis (as in the paper).
- Statistic: housekeeping gap (relative entropy slope) matched on island x young level, then also on the gene's
  own expression decile in the same samples; the share of the gap removed.
- Inference: 200 resamples of donors (fibroblasts: of cultures within line), strata and SD fixed from the full
  sample; an interval only when at least half the resamples are identifiable (series rule). Reported beside the
  existing gene-resampling intervals.
- Claim kept if the donor interval of the share removed excludes zero; otherwise the text says "not resolved over
  donors".

### P2. Division contrasts with donor or animal resampling

- R136 (division vs damage, six donors): donors resampled with their arms together.
- R114 / R141c (colon vs hippocampus off the floor, 10-90 % band, all-CpG scale): mice resampled within tissue,
  contrast from independent resamples.
- R144 (WGBS oligodendrocyte lineage vs neurons, same donors): donors resampled with both fractions.
- R140 (muscle stem cells): intervals labelled as conditional on the six mice; a mouse-level leave-one-out reported.
- The division claim is stated as supported only where its donor/animal interval excludes zero.

### P3. Forest plot of every series

- Rows: every series' housekeeping gap (relative entropy slope, all CpGs, chip intercepts where estimable) with
  donor/animal interval, platform, processing (rebuilt/deposited), n, age span; one row per independent cohort and
  tissue. Paired contrasts only where donors or mice are shared, resampled jointly. No meta-regression and no
  moderator claim.

### P4. One validated clock, then outcome models

- Clock: PhenoAge (published scoring with intercept; validated against the authors' example output), and
  DunedinPACE through its author-supported implementation if it runs on our arrays. No other clock is used for the
  primary models.
- Outcome (primary, one): EPIC-Italy breast cancer diagnosed more than 2 years after the draw, women; rebuilt
  arrays; the matched case-control design of R96 (age within 3 years, same chip first), matched sets resampled.
- Models on identical observations: clock acceleration (with age, cell fractions, chip in the matching), then + B
  (housekeeping look-alike drift), then + F given B. Likelihood-ratio tests; odds ratios per SD with set-resampling
  intervals. Decision: B or F is reported as adding to the clock if its LR p < 0.05 and its interval excludes 1.
- SJLIFE and Vitality 90+: descriptive, cohort-specific, no pooled estimate.

## Multiplicity

P1-P3 are estimation with intervals, no tests. P4 has one primary test (B added to the clock, EPIC-Italy breast > 2
y) and one secondary (F given B); everything else is exploratory and labelled so.

## P2 outcome (2026-11-24)

Existing outputs already carry the required resampling: colon minus hippocampus resamples mice (all CpGs −0.11,
−0.22 to −0.08); off the floor each tissue resamples mice and the intervals separate (colon −0.30, −0.61 to −0.23;
hippocampus +0.00, −0.05 to +0.05); R136 division vs damage holds leaving out each of the six donors (F2 −0.36 to
−0.46); R144 WGBS resamples donors with both fractions (oligodendrocyte minus neuron +0.03, −0.08 to +0.12,
unresolved); muscle stem cells: −0.056 conditional on six mice, −0.03 to −0.09 leaving out one mouse. Text updated to
the current tables; no new runs needed.

## Change 1 (2026-11-24, before any data were downloaded): P5, maintenance vs selection in single-cell clones

Data: Flunkert et al. 2018 (GSE112812, GSE112873): two human fetal fibroblast strains, single-cell clones grown about
20 doublings, and their parental cultures; 450k IDATs, rebuilt with the shared recipe. Only non-irradiated clones in
the primary analysis; irradiated clones exploratory.
Statistic: per clone, each CpG's relative entropy from the parental culture of its strain (the parent's mean level
as the reference); the housekeeping gap of that per-clone distance against matched CpGs (island x parental level),
averaged over clones of a strain; intervals by resampling clones within strain.
Prediction: replication-coupled maintenance keeps the gap inside a lineage, so clones show a negative gap of the
size seen per doubling in cultures; selection among lineages cannot act inside one clone, so a clone's drift shows
no gap. Decision: "maintenance supported" if the clone interval lies below zero; "selection not excluded" if it
covers zero. The parental bulk cultures' own drift (where several passages exist) is reported beside it.

## Change 2 (2026-11-24, Anton; before any clock was scored): a clock panel and clock disagreement

Every clock in biolearn's gallery at the pinned commit is scored on each rebuilt series, with biolearn's own
implementation and its documented missing-CpG policy. A clock enters the panel of a series only without quality
loss: at least 90 % of its absolute weight measured (not imputed) in at least 90 % of the arrays, and it runs on
that platform without error. Per sample: each clock's age acceleration (residual on age, sex, cell fractions and chip
in each series); the panel's disagreement = the SD across the panel's chronological and phenotypic clocks of their
standardised accelerations.
Primary analyses stay as in P4 (PhenoAge, EPIC-Italy breast > 2 y). Added, exploratory and labelled so: (i) each panel
clock in the same nested models; (ii) disagreement as a predictor beside the mean acceleration; (iii) whether the
shield share F and look-alike drift B correlate with disagreement. All panel clocks are reported, none selected.

## Change 3 (2026-11-24): P5's decision rule replaced after review (docs/reviews/2026-11-24-codex-fibroblast-clones.md)

The clone endpoint mixes the founder cell's inherited state (which prior selection may have shaped), which founders
grew out, and change during expansion; the two parental bulk arrays cannot give the parent's cell-to-cell
distribution. A gap inside clones therefore cannot discriminate maintenance from selection, and the rule "maintenance
supported if the clone interval lies below zero" is withdrawn. The result is reported as: housekeeping promoters of
single-cell clones stay closer to their parent than matched CpGs; maintenance and selection before or during
cloning both remain possible. Discriminating data (sister branches of the same founders over doublings, with lineage
abundance in pooled cultures) do not exist in open repositories.

## Change 4 (2026-11-24): arrays with failed signal leave QC; clocks get their own unmasked rebuild

*Intensity.* Detection alone let 127 of EPIC-Italy's 839 detection-passing arrays through at about a tenth of the
series' median signal; they sit on 15-16 chips and carry 31 % of breast cases against 16 % of controls, and under exact
chip matching they produced a drift (B) odds ratio of 4.5 per SD that vanished on chips with normal signal
(docs/reviews/2026-11-24-exact-chip-B.md). From now on an array also fails QC when its mean intensity is below a
quarter of the median of the arrays passing detection in its series (`reprocessed.MIN_INTENSITY_SHARE`; flagged
`low_intensity`). Across every rebuilt series this removes the 127 EPIC-Italy arrays and no other. Every result on
the rebuilt EPIC-Italy arrays is recomputed (incident and matched cancer models, the clock models).

*Clock track.* SeSAMe's design-flag mask (Q) removes 20 % of GrimAge's weight on every array. Clocks are scored on a
second rebuild without Q (SeSAMe CDPB: detection masking kept), `data/interim/sesame_clock`; drift analyses keep QCDPB.
Masked and unmasked scores of the same clocks are compared and reported.

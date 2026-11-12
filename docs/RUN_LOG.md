# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

## R30 — clocks of our own (RM11), 2026-11-12

- Commit 275984f; job 6432498 failed after 8 s: one stratum (open sea, lowest decile) had 512 shielded CpGs and
  373 background CpGs. Fixed: both pools take the scarcer count per stratum. Rerun as R31.

## R29 — shield against clock acceleration, coupling-free measures (RM11), 2026-11-12

- Commit 03333c2; `figures shield_clocks --scale paper --seed 42`; 8 cores, 64 GB; under a minute.
- Result, as measured (blood, 5 components, housekeeping; partial r with acceleration, age and technical
  covariates out): B +0.21 to +0.58 for StocH, StocZ, DunedinPACE, epiTOC1; P at equal B −0.15 to −0.29 for the
  chronological and stochastic clocks, −0.79 epiTOC1, −0.63 cortical; P at equal B ≈ 0 for PhenoAge (−0.07),
  DunedinPACE (0.00), Zhang mortality (−0.04). Cortex mostly n.s. (epiTOC1 −0.51).

## R28 — shield against clock acceleration per person (RM11), 2026-11-12

- Commit 3f3008b; job 6432348, 27 s, 6.6 GB. P and F correlated negatively with almost every clock's
  acceleration: built in (P subtracts look-alike drift, clocks are made of drifting CpGs). Superseded by R29.

## R27 — hippocampus replication of the accessibility result (RM19), 2026-11-12

- Commit 60e4191; `figures atac_hippocampus --scale paper --seed 42`; 8 cores, 32 GB.
- Job 6432066 completed in 6 min, 4.4 GB at peak. Subclasses with 50+ nuclei in 20+ donors: astrocytes 34,
  microglia 40, OPC 38, oligodendrocytes 40, subiculum 20; 13,561 genes.
- Result, as measured: signs depend on the adjustment (oligodendrocytes core-essential raw +0.24, promoter
  fraction −0.23, both −0.01; OPC raw −0.18, both +0.14). With both adjustments: housekeeping +0.05 (astrocytes),
  +0.15 (OPC), +0.10 (subiculum), p ≤ 0.007; microglia and oligodendrocytes ≈ 0.

## R26 — hippocampus replication of the accessibility result (RM19), 2026-11-12

- Commit e084027; `figures atac_hippocampus --scale paper --seed 42`; 8 cores, 32 GB.
- Job 6432044 failed after 3 min: no barcode matched. The metadata's cell ids carry a library tag for some donors
  (`hc77_deep_...`) that the 10x barcodes lack. Fixed (match on the barcode alone; one tag per donor, no
  repeats); rerun as R27.

## R25 — clonal haematopoiesis and the shield (RM16), 2026-11-12

- Commit de1ddda; `figures chip --scale paper --seed 42`; 8 cores, 32 GB.
- Job 6431885 completed in 1 min 46 s, 0.9 GB at peak. 64 patients, 101,667 CpGs; any CH 29, DNMT3A 21, TET2 12
  against 35 without CH.
- Result, as measured: no effect on B, F or S in any arm, raw or with 5 components (all |F| ≤ 0.022, p ≥ 0.15).

## R24 — cortical cell types with the promoter read fraction as a covariate (RM19), 2026-11-12

- Commit acb8254; `figures atac_cortex --scale paper --seed 42` on R23's cached pseudobulk; 8 cores, 32 GB.
- Job 6431737 completed in 9 min, 1.3 GB at peak.
- Result, as measured (promoter fraction and 5 components): core-essential +0.07 to +0.23 SD in 6 of 7 cell types
  (p ≤ 0.055; inhibitory neurons +0.11, p 0.009); housekeeping +0.04 to +0.16 in excitatory neurons, microglia,
  OPCs and vascular cells (p ≤ 0.013), −0.02 to +0.03 in astrocytes, inhibitory neurons, oligodendrocytes (n.s.).

## R23 — accessibility spread per cortical cell type (RM19), 2026-11-12

- Commit 8c0b788; `figures atac_cortex --scale paper --seed 42`; 8 cores, 48 GB.
- Job 6431650 completed in 11 min, 7.6 GB at peak; pseudobulk cached. 121-348 donors per cell type, 4,173
  housekeeping and 656 core-essential promoter peaks.
- Result, as measured: essential promoters spread apart with age more than matched promoters in every cell type
  (raw housekeeping +0.004 to +0.17 SD, core-essential +0.02 to +0.16; with 5 components core-essential +0.08 to
  +0.23, all p < 0.03). Superseded for reading by R24, which adds the promoter read fraction as a covariate.

## R22 — look-alike drift and relative gap against lifespan (RM20), 2026-11-12

- Commit 58efc63; `figures species --scale paper --seed 42`; 8 cores, 40 GB.
- Job 6430954 waited on the CPU cap, then completed in 61 min, 25.8 GB at peak. The `gap_sd` fits equal R12's.
- Result, as measured (per log maximum lifespan, weight and tissue fixed): log look-alike drift rate −1.42
  (entropy, housekeeping), −1.15 (dispersion), p < 10⁻⁶; relative gap (gap / look-alike rate) −2.42 (entropy,
  housekeeping, p 0.04), −0.83 (core-essential, p 0.07), +0.46 and +0.49 (dispersion, p 0.62, 0.34). Groups with a
  non-positive look-alike rate are left out of these two outcomes (122-125 of 139 for entropy, 75-82 for
  dispersion).

## R20 — headline gaps against probe design and unreliable probes (RM18), 2026-11-12

- Commit 2f4654e (fixture fix c729b3c); `figures probe_robustness --scale paper --seed 42`; 8 cores, 64 GB.
- Job 6430935 waited on the CPU cap, then completed in 41 min, 5.6 GB at peak.
- Result, as measured: with probe type in the strata, blood lifespan gaps shrink by about 15 % (housekeeping
  entropy −0.348 → −0.295); cortex and same-donor gaps move by under 0.01; the Zhou mask moves no gap by more
  than 0.01; no gap changes sign. Set CpGs are type I in 42-45 %, all promoter CpGs in 35 %.

## R19 — unbiased bootstrap intervals, tissues and mouse (RM9), 2026-11-12

- Commit 6995d10; `figures tissues tissues_cells mouse --scale paper --seed 42`; 8 cores, 96 GB.
- Job 6430533 waited on the CPU cap, then completed in 2 h 52 min, 33 GB at peak. Compute nodes have no git,
  so the log carries no commit line; the commit is in the job comment.
- Result, as measured: intervals widen (median width x1.16 tissues, x1.18 tissues_cells, x1.24 mouse). Every
  human interval that excluded zero still does (20 of 21, 73 of 84). Mouse: 11 of 42 contrasts exclude zero
  (15 before); lost are core-essential and non-essential contrasts of renewing or slow against post-mitotic
  tissue; housekeeping entropy, renewing − post-mitotic −0.16 (−0.27 to −0.05), slow − post-mitotic −0.18
  (−0.28 to −0.07), holds. The lead's reading is in `docs/extension.md`.

## R18 — chromatin-writer syndromes in blood, with age (RM14), 2026-11-11

- Commit b7cf41c; `figures syndromes --scale paper --seed 42`; 4 cores (lowered while queued), 8 GB.
- Job 6430424 waited on the CPU cap, then completed in 16 min, 1.3 GB at peak. GSE97362: 3 without an age and 33
  outside the controls' 2-19 years left out; GSE74432 has no ages (no age covariate).
- Result, as measured (housekeeping, raw / 5 components): drift B up in Sotos +0.018/+0.021, NSD1 variants
  +0.013/+0.012, Weaver +0.011/+0.010, Kabuki +0.009/+0.011, CHARGE +0.007/+0.005; KMT2D and CHD7 variants ≈ 0.
  Shield share F: Sotos −0.051/−0.043, NSD1 variants −0.067/−0.039, Weaver −0.063/−0.037; Kabuki, CHARGE, KMT2D,
  CHD7 variants within ±0.016. The lead's reading is in `docs/extension.md`.

## R17 — chromatin-writer syndromes in blood (RM14), 2026-11-11

- Commit e4ec847; `fetch syndromes`, `figures syndromes --scale paper --seed 42`; 8 cores, 32 GB.
- Job 6430339 waited on the CPU cap, then completed in 26 min, 1.3 GB at peak.
- Superseded by R18: one GSE97362 sample has no age, and the age covariate was dropped for the whole series
  instead of the sample. Tables not committed.

## R16 — which way essential CpGs move under interventions (RM12), 2026-11-11

- Commit d2d7dcb; `figures mmc_interventions --scale paper --seed 42`; 8 cores, 16 GB.
- Job 6430333 completed in 76 s, 3.5 GB at peak. B and F rows equal R15's.
- Result, as measured (signed shift S, same-chip): Tet3 heterozygous core-essential +0.012 (cortex, p 3·10⁻⁴),
  +0.014 (striatum, p 2·10⁻⁴); GHRKO kidney core-essential +0.028 (p 8·10⁻⁶); CR liver housekeeping +0.005
  (p 0.001). The lead's reading is in `docs/extension.md`.

## R15 — interventions against same-chip controls (RM12), 2026-11-11

- Commit 0cb79ae; `figures mmc_interventions --scale paper --seed 42`; 8 cores, 16 GB.
- Job 6430322 completed in 73 s, 3.5 GB at peak. Same-chip controls exist for 13 arms; partial reprogramming
  has too few.
- Result, as measured (housekeeping F, same-chip): Tet3 heterozygous cortex +0.17 (p 0.006), striatum +0.25
  (p 9·10⁻⁶), 8 against 8; GHRKO liver +0.01 (p 0.67). The lead's reading is in `docs/extension.md`.

## R12 — shielding across mammals (RM10), 2026-11-11

- Commit fa6a796; `fetch mammals` (4.5 GB, 4.7 min) then `figures species --scale paper --seed 42`; 8 cores, 64 GB.
- Job 6429555 completed in 59 min, 25.6 GB at peak (the next run fits in 32 GB). 139 species-tissue groups of 65
  species, 11,529 samples.
- Result, as measured: housekeeping gap per log maximum lifespan +0.029 (entropy, p 0.14), +0.036 (dispersion,
  p 0.11), weight and tissue held fixed. The lead's reading is in `docs/extension.md`.

## R13 — mouse interventions across tissues (RM12), 2026-11-11

- Commit b1c9ae0; `figures mmc_interventions --scale paper --seed 42` on R12's download, started when R12 ended.
- Job 6429673 completed in 74 s, 3.5 GB at peak. 13 study-tissue-arm comparisons with at least 5 animals per arm.
- Result, as measured (housekeeping F): CR in liver −0.008 (p 0.77); Tet3 heterozygous cortex +0.16, striatum
  +0.13 (p < 10⁻⁶; controls pooled across studies); GHRKO liver −0.10 (p 0.01). The lead's reading is in
  `docs/extension.md`.

## R14 — directed demethylation outside promoters (RM15), 2026-11-11

- Commit f84d2a7; `fetch annotation` (solo-WCGW tags, 10 s, 1.9 GB) then `figures intergenic --scale paper --seed 42`;
  8 cores, 32 GB asked; first in priority, by antonkulaga's word.
- Job 6430056 completed in 41 s, 2.9 GB at peak. 130,991 open-sea CpGs outside promoters in Hannum (13,103
  solo-WCGW), 131,254 in Jaffe.
- Result, as measured: blood solo-WCGW −0.21 percentage points per decade (p 4·10⁻⁵) against −0.08 for other
  open-sea CpGs; cortex none. The lead's reading is in `docs/extension.md`.

## R11 — the mouse figures with corrected intervals (RM2, RM3), 2026-11-11

- Commit f52c6e5; `figures interventions bxd rrbs_tissues --scale paper --seed 42`; 8 cores, 24 GB asked.
- Job 6429326 completed in 54.4 min, 10.6 GB at peak. Every aging-gap point estimate now lies inside its interval;
  the intervals widened toward zero and none of the essential-gene intervals crossed it.

## R10 — six mouse tissues with numeric ages (RM2), 2026-11-11

- Commit d3e052d; `fetch rrbs-tissues` then `figures rrbs_tissues --scale paper --seed 42`; 8 cores, 64 GB asked.
- Job 6429231 failed after 54 min, in the figure: matplotlib refused a negative error bar, because the aging-gap
  bootstrap was biased toward zero and left point estimates outside their intervals (fixed in the code, see
  the changelog). The fetch and parse completed (2.2 GB peak); the tables (`tissues`, `aging_gaps`,
  `age_effects`) were written before the figure. 131,756-384,453 promoter CpGs per tissue, 54-72 mice each.
- Rerun of the figure in R11.

## R9 — the shield in BXD mouse liver (RM3), 2026-11-11

- Commit e8b851d; `figures bxd --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set.
- Job 6429222 completed in 40 s, 1.0 GB at peak (a rerun fits in 2 GB). 4,221 promoter probes, 321 females, 49 in
  the young chow reference.
- Result, as measured (housekeeping F): per 10 months −0.28 (p < 10⁻¹⁰), −0.055 with five composition components
  (p 10⁻⁴); HFD −0.105 (p 3·10⁻⁵), −0.017 with components (p 0.26). The lead's reading is in `docs/extension.md`.

## R8 — interventions in mouse blood RRBS, composition check (RM3), 2026-11-11

- Commit d5c7433; `figures interventions --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set. Waited
  on the account's CPU cap before it started.
- Job 6429130 completed in 37.5 min, 10.6 GB at peak (a rerun fits in 16 GB).
- Result, as measured (housekeeping, shield share F): CR in B6D2F1, one sample series, +0.025 (p 0.012) and +0.025
  with five composition components (p 0.002); CR in C57BL/6, two series, +0.027 (p 0.002) and −0.016 with
  components (p 0.06). The lead's reading is in `docs/extension.md`.

## R7 — interventions in mouse blood RRBS (RM3), 2026-11-11

- Commit 30ed4b8; `fetch rrbs` (243 blood tables, about 3 GB, and the mm10 tracks) then
  `figures interventions --scale paper --seed 42`; 8 cores, 64 GB asked; job comment set.
- Job 6429040 completed in 53 min: the fetch and parse 15.9 min at 19.7 GB peak, the figure 37.4 min at 9.5 GB. A
  rerun of the figure fits in 16 GB. 542,008 promoter CpGs covered at 10x or more in 90 % of the 243 samples.
- Result, as measured (housekeeping): CR moves the shield share F from −0.060 to −0.032 in C57BL/6 (p 0.002) and by
  +0.025 in B6D2F1 (p 0.012), with drift B changed by −7 % and −3 % (n.s.). Outputs superseded by R8, which adds a
  composition check.

## R6 — clocks from protected and unprotected CpGs (RM8), 2026-11-11

- Commit 7d40536; `figures chromatin_clocks --scale paper --seed 42`; 8 cores, 32 GB asked; job comment set. Waited
  on the account's shared CPU cap before it started.
- Job 6428565 completed in 39.4 min, 2.8 GB at peak (the next run fits in 6 GB). Pools of 33,786 (active in
  both), 19,330 (Polycomb in both) and 149,740 (any) shared CpGs, 4,000 per clock, 5 draws.
- Result, as measured: no pool transfers well between tissues (out-of-tissue r 0.06-0.70); the Polycomb pool's
  brain-minus-blood offset in Lunnon varies least between draws (SD 4-6 years against 20-42). The lead's reading
  is in `docs/extension.md`.

## R5 — published clocks on the protection map (RM8), 2026-11-11

- Commit 4787a9b; `figures clocks --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set. Waited about an
  hour on the account's shared CPU cap (`AssocGrpCpuLimit`) before it started.
- Job 6428487 completed in 11.8 min, 3.2 GB at peak (the next run fits in 6 GB). 16 clocks, 64,900 CpGs per cohort.
- Result, as measured: clock CpGs switch TSS activity between blood and cortex at 0.15-0.33 against 0.12 for
  background; the tissue offset does not follow that share (Spearman 0.11, p 0.68, 16 clocks). The lead's reading
  is in `docs/extension.md`.

## R4 — the chromatin-matched gap with sorted brain chromatin (RM7), 2026-11-11

- Commit b0184a3; `figures chromatin --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set.
- Job 6428474 completed in 27.8 min, 5.2 GB at peak (the next run fits in 10 GB). It rewrote R3's outputs with the
  sorted-ATAC matchings added; R3's own numbers are unchanged in them.
- Result, as measured (housekeeping, brain): sorted neuron and glia accessibility removes 30-37 % of the gap,
  with the bulk state 32-49 %. The lead's reading is in `docs/extension.md`.

## R3 — the chromatin-matched gap (RM7), 2026-11-11

- Commit 105cc56; `fetch annotation` (re-prepared the probe table with positions, 5 s, 1.8 GB), then
  `figures chromatin --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set.
- Job 6428419 completed in 20.5 min, 5.2 GB at peak: the next run fits in 10 GB. About 109,000 promoter CpGs per
  tissue in Hannum, Jaffe and Lunnon's five tissues.
- Result, as measured (housekeeping, dispersion slope in SDs): blood −0.289 matched on island and level, −0.034
  matched on its own chromatin (E062) as well; cortex −0.075 and −0.058. The lead's reading is in
  `docs/extension.md`. Tables in `data/output/chromatin_match/`, figure `docs/figures/chromatin_match.png`.

## R2 — the per-person score (RM6), 2026-11-11

- Commit 2fe148a (module `entropy_age.protection`, figure `per_person`); `--scale paper --seed 42`; 8 cores, 48 GB
  asked; job comment set by the account's rule.
- Job 6428234 completed in 13 min (88 CPU-minutes, 6.8 cores busy on average), 5.9 GB at peak: the next run fits in
  12 GB. Scored 632 Hannum and 159 Jaffe samples (references: 24 and the Jaffe controls aged 20-35), and 71-75
  Lunnon donors per region.
- Result, as measured (housekeeping, no components): P per decade −0.0018 (SE 0.0003) in blood, −0.00036
  (SE 0.00013) in cortex; reliability 0.997 and 0.955; blood's F below the same donor's region in 86-91 % of
  donors. The lead's reading is in `docs/extension.md`. Tables in `data/output/per_person/`, figure
  `docs/figures/per_person.png`.

## R1 — GSE179847 (RM5), 2026-11-10

- Commit b7036e1 (cohort `sturm`, figure `passage_sturm`); `--scale paper --seed 42`; 8 cores, 64 GB asked as the
  first measure of the parse's memory.
- Jobs 6427994 (parse) and 6427995 (figure, after it): the parse failed in 1 s, before any work: the job script
  sourced the shell profile under `set -u`, which the system profile does not survive. Fixed (`set -u` after the
  profile); the figure job cancelled.
- Resubmitted as 6427996 (parse) and 6427997 (figure, after it): both completed. The parse took 28 s and
  4.4 GB at peak; the figure 43 s and 1.9 GB (the next run of either fits in 8 GB). The job log printed
  "commit unknown", because the CPU nodes have no git; the clone was at ffbc2c5 (code as b7036e1).
- Result, as measured (housekeeping set, the five baseline lines): coupling change per 10 PD −0.055, interval
  −0.21 to +0.10; mean coupling −0.27. The lead's reading is in `docs/extension.md`. Tables in
  `data/output/passaged_cells_sturm/`, figure `docs/figures/passaged_cells_sturm.png`.

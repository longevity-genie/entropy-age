# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

## R136 — senescence by division and by damage in six donors (RM40, GSE227160), 2026-11-22

- Commit 287ec10; `figures senescence_donors --scale paper --seed 42`; job 6592188 after the SeSAMe rebuild of 48
  EPIC arrays (47 pass). Two fibroblast donors and four stromal-cell donors, each with young, middle-passage, two
  replicative-senescent, two irradiated and two doxorubicin cultures; every donor has one array on each of 8 chips.
  One model per CpG: methylation on donor, chip and treatment. 89,960 promoter CpGs.
- **Squared drift from young, noise removed (matched CpGs of the housekeeping set).** Replicative senescence 0.0141
  (0.0139 to 0.0165 leaving one donor out); irradiation 0.0006; doxorubicin −0.0015; middle passage −0.0003.
- **Share of squared drift in replicative senescence.** Housekeeping −0.39 (leave-one-donor-out −0.36 to −0.46);
  core-essential −0.35; non-essential +0.51. In the other arms the share divides two near-zero values.
- **Not read.** The gap of |coefficient| is negative in all four arms (housekeeping −0.064 to −0.117 SD), drift or
  not: lower array noise at these promoters, not a shield.
- Reading: repeats R111 (share −0.58, one cell type, chip as unit) in six donors of two cell types with donor and
  chip held: drift and the shield come with division to senescence; arrest by damage shows no drift.

## R135 — astrocytes immortalised with telomerase (RM40, GSE226079), 2026-11-22

- Commit d0abc68 plus the module commit after it; `figures astrocyte_htert --scale paper --seed 42`, then with
  `ENTROPY_AGE_CHIP=1`; job 6592161 after the SeSAMe rebuild of 70 EPIC arrays. hTERT cultures with doublings: 43 in
  three replicate lines (15, 14, 14), 29.5 to 72.2 doublings, 12 chips; 0.81 of the variance of doublings is left
  within line and chip. 87,104 promoter CpGs measured in all. Slopes on doublings with a line intercept.
- **Gap per doubling.** Housekeeping −0.083 SD (−0.092 to −0.074); chip held −0.084 (−0.092 to −0.074).
  Core-essential −0.041 (−0.060 to −0.023); chip held −0.042. Non-essential +0.141; chip held +0.142.
- Reading: a second immortalised cell type carries the gap per doubling, with chip held, a little below primary
  astrocytes (−0.102, R133) and fibroblasts. Three replicate lines of one source; no arrested control.

## R114c and R141c — mouse tissues off the 0 % floor (review item; Anton asked for the follow-up analyses outside blood and brain), 2026-11-24

- Commit 3060c3f; `figures jena_rrbs` and `figures mouse_atlas_rrbs`, tables `gaps_off_floor` (CpGs with young level
  10-90 %, own strata, mouse bootstrap). Gaps below are rescaled to the all-CpG SD of each tissue, as the review
  required for comparing bands.
- **Colon and hippocampus (R114c).** Housekeeping: colon −0.247 (−0.514 to −0.189; 35,349 CpGs), hippocampus +0.006
  (−0.098 to +0.096; 40,712). Core-essential −0.267 (−0.623 to −0.186) and −0.076 (−0.302 to +0.122).
- **Four tissues of GSE225166 (R141c).** Housekeeping: lung −0.078 (−0.191 to +0.003; 10,316), heart +0.053 (−0.035
  to +0.123; 11,262), cortex +0.017 (−0.235 to +0.163; 5,062), liver 290 CpGs (no estimate).
- Reading: in colon the gap is not a floor effect; in hippocampus all of it sits at CpGs near 0 %. The colon-over-
  hippocampus contrast is therefore stronger off the floor. In the four-tissue series no tissue has a resolved gap
  off the floor, so it adds little either way.

## R141b — mouse tissues with complete CpGs only (review item), 2026-11-23

- Commit 3ec7b03; job 6593358. R141 had kept CpGs present in 90 % of libraries and filled the rest; this run uses
  CpGs covered at least 10 times in every library. CpGs left: lung 28,727, heart 30,312, cortex 12,335, liver 518.
- Housekeeping, relative entropy: lung −0.086 (−0.139 to −0.015); heart −0.046 (−0.126 to +0.008); cortex −0.011
  (−0.083 to +0.132); liver +0.099 (−0.125 to +0.316) on 518 CpGs, not an estimate.
- Reading: R141's heart and liver gaps do not survive the nothing-filled rule; only lung is resolved. Complete
  coverage in every library is a harsh filter for RRBS; slopes that tolerate missing libraries without filling would
  keep more CpGs and are the proper next step. R141's numbers are superseded.

## Reruns after the interval fix (991f93d), 2026-11-23

- A code review found that the donor and culture resampling intervals of R138, R140 and R142 divided each resample
  by its own SD, which pulls them toward zero. They are now scaled by the full sample's SD; point estimates are
  unchanged. The intervals given in the entries below for these three runs are superseded by:
- **R142, neurons (82 donors), housekeeping.** Without chip −0.051, donors −0.091 to −0.008; with chip −0.030, donors
  −0.102 to +0.037. Controls (28): −0.030, −0.085 to +0.001; with chip not estimable.
- **R138, culture bootstrap, housekeeping.** Telomerase fibroblasts −0.160 to −0.072; mortal cells on shared chips
  −0.036 to −0.011; earlier batch −0.002 to +0.158; astrocytes −0.128 to −0.052; astrocytes with telomerase −0.089
  to −0.080.
- **R140, single muscle stem cells, divergence.** Housekeeping −0.056 (−0.094 to −0.027); core-essential −0.074
  (−0.141 to −0.010); with expression matched −0.031 (−0.062 to +0.005) and −0.058 (−0.113 to +0.015).

## R142 — the shared array age-series runner, regression on two series, 2026-11-23

- Commit 6a3bc03; `figures series_neurons`, `figures series_brain_regions` (description files in `config/series`);
  job 6593085. Gaps equal those of `neuron_series` to three decimals (housekeeping −0.051, −0.030 with chip;
  controls −0.030, +0.024).
- **New: donor-resampling intervals.** Neurons, every donor: −0.060 to −0.007 without chip; −0.047 to +0.022 with
  chip. The gene-resampling interval (−0.040 to −0.020 with chip) had made the chip-held neuronal gap look resolved;
  it is not. Paper corrected.

## R141 — tissues of the same mice at three ages by RRBS (GSE225166), 2026-11-23

- Commit after 83d0444; `figures mouse_atlas_rrbs --scale paper --seed 42`; job 6593084. Liver 32 mice (3x12, 15x9,
  24x11; 32,524 promoter CpGs), heart 30 (86,787), cortex 31 (73,648), lung 30 (101,665); 8 lanes, 5 to 8 of them
  holding two or more ages; sex and lane held; 200 resamples of mice. The ten other tissues have 5 to 10 mice and
  are not analysed.
- **Housekeeping, relative entropy.** Liver −0.104 (−0.218 to −0.048); lung −0.104 (−0.154 to −0.009); heart −0.079
  (−0.182 to −0.009); cortex −0.051 (−0.139 to +0.025). Core-essential −0.105, −0.105, −0.083, −0.077 (cortex
  interval to +0.006). Non-essential +0.10 in heart and lung.
- Reading: with lane held, heart carries a gap of the size seen in dividing tissues' lower range; cortex is the
  smallest and unresolved. The order liver = lung > heart > cortex is the one of the consortium array (R118).

## R140 — between-cell drift in single muscle stem cells (GSE121437), 2026-11-23

- Commit after 304388a (plain between-cell variance, promoters matched on island status, young level and coverage
  fifth); `figures musc_single_cell --scale paper --seed 42`; job 6592905. scBS-seq of 275 cells (136 young from
  mice Y2, Y7, Y8; 139 old from O1, O5, O8); 26,144 promoters with at least 4 reads in at least 20 cells per age;
  16,200 with expression from the count table of the same cells. 200 resamples of cells within mouse.
- **Divergence (old minus young between-cell variance), set minus matched, SD.** Housekeeping −0.056 (−0.072 to
  −0.020), each mouse left out −0.024 to −0.085; core-essential −0.074 (−0.109 to −0.008), left out −0.033 to
  −0.076; non-essential +0.055 (−0.016 to +0.105). With the expression decile matched: housekeeping −0.031 (−0.047
  to +0.004), core-essential −0.058 (−0.085 to +0.011).
- **Between-cell variance in young cells.** Higher at essential promoters than matched (housekeeping +0.113,
  core-essential +0.088), +0.047 and +0.037 with expression matched; not interpreted: read sampling within one cell
  differs between CpG-dense and CpG-poor promoters.
- **Shift of the mean.** No difference (housekeeping −0.013, −0.015 to +0.017).
- Caveats: three mice per age; old cells carry about 19 more reads per promoter per cell than young (mean 60),
  which lowers their measured variance; promoters are matched on mean coverage, not on that difference. Not in the
  paper until the coverage difference is controlled (downsample old cells to the young depth).
- Reading: suggestive that age adds less cell-to-cell variation at essential promoters within one sorted cell
  type, about half of it explained by expression.

## R139 — cortex accessibility with a measure free of counting noise (RM19), 2026-11-23

- Commit 16492c1; `figures atac_cortex --scale paper --seed 42` on the cached pseudobulk; job 6592847. Seven cell
  types, 121 to 348 donors each; table `gaps` now carries `measure`.
- **Spread of log counts** (as R23, R24): essential promoters spread more with age (core-essential +0.07 to +0.23 SD
  with both adjustments; housekeeping −0.02 to +0.16).
- **Excess over counting noise** (age slope of the squared coefficient of variation beyond Poisson): core-essential
  −0.011 to +0.012, housekeeping −0.008 to +0.012 across cell types and adjustments; three of 14 with p < 0.05
  (excitatory neurons −0.004, inhibitory neurons −0.011, vascular cells +0.012), no common sign.
- Reading: the earlier "essential promoters spread more" was a property of the log-count measure, whose noise
  depends on read depth. With counting noise removed, accessibility shows no gap in either direction.

## R137 — three brain regions of about 300 donors (RM40, GSE349373), 2026-11-23

- Commit 90f59f5; `figures brain_regions --scale paper --seed 42`; job 6592740 after the SeSAMe rebuild of 904 EPIC
  arrays (job 6592305). Prefrontal cortex, amygdala, hippocampus; controls, MDD and PTSD about 100 each per region;
  one array per donor and region; Jaffe's adult control donors left out; sex and diagnosis held.
- **Every donor, housekeeping, without / with chip.** PFC −0.086 / −0.081 (−0.089 to −0.069; 294); amygdala −0.065 /
  −0.070 (−0.081 to −0.058; 297); hippocampus −0.079 / −0.069 (−0.081 to −0.056; 300). Core-essential with chip
  −0.059, −0.043, −0.053. 10-90 % with chip: −0.061, −0.055, −0.045.
- **Controls, with chip.** −0.043 (100), −0.036 (99), −0.078 (98); off the floor the amygdala controls are not
  resolved (−0.012).
- 0.73 to 0.79 of the variance of age left within diagnosis and chip (0.36 to 0.49 in controls).
- Reading: the largest brain replication; the gap has the size of the discovery cohort (−0.06) in three regions,
  holds with chip intercepts and off the floor.

## R134b — sorted neurons, one array per donor (GSE112179), 2026-11-22

- Commit after 54ca722 (neuron_series keeps one array per donor); job 6592267. The series holds 100 arrays of 83
  donors (9 donors with replicates); R134 below counted arrays as donors and is superseded by this run.
- **82 donors (28 controls), sex and diagnosis held.** Housekeeping −0.051 (−0.060 to −0.042); with chip −0.030
  (−0.040 to −0.020). Core-essential −0.043 (−0.060 to −0.027); with chip −0.009 (−0.028 to +0.014). 10-90 %:
  housekeeping −0.078, with chip −0.043. 0.62 of the variance of age left within diagnosis and chip.
- **Controls alone (28).** Housekeeping −0.030 (−0.041 to −0.018); with chip +0.024 with 0.18 of the variance of age
  left: not estimable with chip held.
- Reading unchanged: neurons carry a small housekeeping gap that survives chip intercepts.

## R134 — sorted cortical neurons, third series (RM40, GSE112179), 2026-11-22

- Commit d0abc68; `figures neuron_series --scale paper --seed 42`; job 6592157 after the SeSAMe rebuild of 100 EPIC
  arrays (99 pass detection). NeuN+ nuclei, frontal cortex, ages 23 to 77, 13 chips; 33 controls on 11 chips.
- **Every donor, sex and diagnosis held.** Housekeeping: relative entropy −0.054 (−0.063 to −0.046); with chip
  −0.034 (−0.043 to −0.024). Core-essential −0.046 (−0.062 to −0.030); with chip −0.017 (−0.036 to +0.004).
  10-90 %: housekeeping −0.090, with chip −0.051. Dispersion slope: −0.036, with chip +0.032. 0.61 of the variance
  of age is left within diagnosis and chip.
- **Controls alone (33).** Housekeeping −0.028 (−0.040 to −0.016); with chip +0.023, but only 0.18 of the variance
  of age is left within chips: not estimable with chip held.
- Reading: neurons carry a small housekeeping gap that survives chip intercepts (about half the size in unsorted
  nuclei, −0.05, R109). The first sorted series (R109) had −0.02 and −0.008 within chips. "Absent where cells do not
  divide" therefore holds for arrested cultures, not for neurons; the dispersion measure does not show the gap with
  chip held.

## R133 — a second cell type in culture: passaged fetal astrocytes (RM40, GSE202554), 2026-11-20

- Commit 82e3311; `figures astrocyte_passage --scale paper --seed 42`, then with `ENTROPY_AGE_CHIP=1`; job 6558106
  (test partition) after the SeSAMe rebuild of 63 EPIC arrays. Primary human fetal astrocytes, unsorted cultures
  with doublings recorded: 45 (donor 1: 37, 0 to 11.7 doublings; donor 2: 8, 4.2 to 15.4), 8 chips each spanning
  at least 6 doublings; 0.68 of the variance of doublings is left within donor and chip. 89,883 promoter CpGs
  measured in all. Slopes on doublings with a donor intercept.
- **Gap per doubling.** Housekeeping: relative entropy −0.102 SD (−0.111 to −0.094), entropy slope −0.173,
  dispersion slope −0.078; with chip intercepts −0.101 (−0.110 to −0.092), −0.167, −0.074. Core-essential −0.072
  (−0.083 to −0.061); chip held −0.070. Non-essential +0.342 (+0.272 to +0.411); chip held +0.333.
- Reading: a second cell type in culture gives the fibroblast number (Sturm −0.102, R124; hTERT BJ −0.136, R130),
  over a span of only 15 doublings and with chip held. Two donors, most arrays from one; no arrested control; no
  expression from these cultures.

## R132 — a third cortex series: bulk orbitofrontal cortex (RM40, GSE254293), 2026-11-20

- Commit b79d37c; `figures orbitofrontal --scale paper --seed 42`, job 6556959 (test partition) after the SeSAMe
  rebuild of 96 EPIC arrays. Bulk orbitofrontal cortex, 90 people after duplicates and the detection cut, 24 to 84
  years, 12 chips (chip explains 0.02 of the variance of age over everyone; 0.25 in the 34 with status 0, 0.09 in
  the 56 with status 1). The record codes diagnosis 0 and 1 without naming them.
- **Housekeeping gap, everyone with status held** (84,572 CpGs): relative entropy −0.017 SD (−0.030 to −0.008);
  with chip intercepts −0.010 (−0.022 to −0.001). Dispersion slope −0.037 (−0.050 to −0.027); chip held −0.022
  (−0.034 to −0.012). Core-essential: −0.024 (−0.043 to −0.003) and −0.016 (−0.037 to +0.010); dispersion slope
  −0.046 and −0.037.
- **By status.** Status 0 (34): relative entropy −0.024 (−0.035 to −0.014), chip held −0.005 (−0.015 to +0.006);
  dispersion slope −0.054 and −0.020. Status 1 (56): −0.004 (−0.020 to +0.008), chip held −0.018 (−0.032 to
  −0.006); dispersion slope −0.015 and −0.052. The two groups swap places when chip is held, so neither is read
  on its own.
- Reading: a third cortex series gives the same size as the others (Jaffe −0.06 deposited and rebuilt, White
  −0.005 to −0.02, R109): a small gap, about a tenth of the one in monocytes, epidermis or cultured fibroblasts,
  with an interval that just excludes zero. Bulk tissue, cell mix not held.

## R131 — four tissues of the same adults on shared chips (RM40, GSE280465), 2026-11-20

- Commit 9464bd5; `figures cross_tissue --scale paper --seed 42` (31 min), then with `ENTROPY_AGE_CHIP=1`; job
  6555050 (test partition) after the SeSAMe rebuild of the 163 EPIC v2 arrays. GSE280465 (Apsley et al. 2025):
  buccal swab, saliva, dried blood spot and blood mononuclear cells of adults aged 19 to 70. 17 of 27 chips carry
  all four tissues; in 88 to 100 % of donors the two tissues of a pair sit on one chip.
- **Saliva is out.** 16 of 37 saliva arrays fall below 90 % detection (median 0.93 against 0.98 elsewhere), which
  leaves 16 to 20 donors per saliva pair, under the minimum of 25. The cut was not loosened.
- **Between-person dispersion, housekeeping gap** (age and sex held; first, second, difference with donor
  bootstrap, 400 resamples): buccal −0.046, dried blood spot −0.071, +0.026 (−0.004 to +0.049), p 0.12, 33 donors;
  buccal −0.042, mononuclear cells −0.107, +0.064 (+0.012 to +0.095), p 0.015, 43 donors; dried blood spot −0.071,
  mononuclear cells −0.079, +0.008 (−0.032 to +0.039), p 0.80, 32 donors. Core-essential: no pair differs
  (p 0.17 to 0.86).
- **Age slopes, housekeeping gap.** Relative entropy: buccal −0.110 and −0.126, dried blood spot −0.070 and −0.065,
  mononuclear cells −0.085 and −0.091; no pair differs (p 0.34, 0.64, 0.53). Entropy slope: buccal −0.26, dried
  blood spot −0.16, mononuclear cells −0.27; dried blood spot against mononuclear cells +0.097 (+0.007 to +0.182),
  p 0.025, the others p 0.22 and 0.88.
- **Chip held** (chip intercepts in the dispersion model; the age slopes carry none, about 40 donors on 27 chips
  leave too little): buccal −0.049, dried blood spot −0.061, +0.012 (−0.029 to +0.027), p 0.96; buccal −0.048,
  mononuclear cells −0.089, +0.042 (−0.027 to +0.062), p 0.55; dried blood spot −0.062, mononuclear cells −0.101,
  +0.039 (−0.023 to +0.060), p 0.66.
- Reading: with tissue no longer tied to slide, every tissue shows the gap (−0.04 to −0.11 SD on dispersion,
  −0.07 to −0.13 on relative entropy). No difference between tissues holds once chip is in the model: the one
  seen without it (mononuclear cells more shielded than buccal cells, +0.064) keeps its direction and loses its
  interval. All tissues here renew, so the series does not test turnover, and with 32 to 43 donors a difference
  of the size seen in Lunnon (0.08 to 0.10 SD) is at the edge of what it can detect. Buccal swabs are a mixture
  of epithelial cells and leukocytes, not held.

## R130 — one fibroblast line with and without telomerase: shield per doubling (RM40, GSE268145), 2026-11-20

- Commits c82b42a and the column fix after it; `figures bj_telomerase --scale paper --seed 42`, then with
  `ENTROPY_AGE_CHIP=1`; jobs 6554733 and its rerun 6554813 (test partition), 7 min with the smoke run. GSE268145: BJ
  foreskin fibroblasts (one male line), trajectories A and B, mortal and with hTERT, EPIC, rebuilt from IDATs (job
  6554698). Slopes on population doublings with a trajectory intercept. Three groups, each alone: hTERT (20 arrays
  after the detection cut, 8 chips, doublings 2 to 38); mortal on the same chips (12 arrays, 6 chips, 4 to 42; four
  arrays fall just below 90 % detection); mortal of an earlier batch (23 arrays, 9 other chips, 2.3 to 52). The two
  deeply senescent arrays are left out.
- **hTERT cells.** Housekeeping gap per doubling: relative entropy −0.136 SD (−0.143 to −0.129), entropy slope
  −0.199, dispersion slope −0.122; core-essential −0.064 (−0.082 to −0.037); non-essential +0.137. With chip
  intercepts: −0.138, −0.198, −0.061; core-essential −0.067. Share of the variance of doublings left within
  trajectory and chip: 0.71.
- **Mortal cells on the same chips.** Housekeeping −0.026 (−0.036 to −0.018), entropy slope −0.049, dispersion
  slope −0.096; core-essential −0.025 (−0.044 to 0.000); chip held −0.025, −0.046, −0.045. Variance left: 0.45.
- **Mortal cells, earlier batch.** The measures disagree: relative entropy +0.087 (+0.080 to +0.095), dispersion
  slope +0.045, entropy slope −0.164; core-essential +0.057 and −0.123. Mean drift per doubling on relative
  entropy is ten times that of the other groups (0.010 against 0.001 and 0.000), and 69,658 CpGs are complete
  against 88,000. Variance left: 0.59.
- Reading: cells that divide without approaching senescence carry a shield per doubling as large as any measured
  in culture (Sturm −0.102, R124), so the shield does not need the approach to senescence. The mortal arms do not
  give a second estimate: one is small (12 arrays), and the other's sign depends on the measure. One cell line, no
  arrested control. No comparison between hTERT and mortal is claimed.

## R129 — epidermis of 117 women: age drift against expression of the same women (RM40, fetch order 4), 2026-11-20

- Commit 5f04674; `figures epidermis_expression --scale paper --seed 42`, then the same with `ENTROPY_AGE_CHIP=1`;
  job 6553578 (test partition), 9 min for the first paper pass, after the SeSAMe rebuild of both series and a smoke
  run. Holzscheck et al. 2020: suction-blister epidermis of healthy women. E-MTAB-8993 (EPIC, 65 women; 4 fall
  below 90 % detection, 61 kept, 23 to 79 years, 84,425 promoter CpGs measured in all) and E-MTAB-8992 (450k, 56
  other women, 21 to 76, 96,007 CpGs), never pooled. Expression: E-MTAB-8995, the depositors' Salmon TPM table
  (34,673 genes, ComBat-corrected on the linear scale; 6 % of values below zero set to zero), mean log2(TPM + 1)
  per array group; 27,616 genes with one approved symbol. Age slopes hold the depositors' batch.
- **Chip.** The IDAT files were renamed, and the dataset search recorded the chip as not recoverable. It is in each
  file's header (`scripts/idat_barcodes.py`): 9 EPIC chips and 5 450k chips, each spanning 27 to 54 years of age;
  chip explains 0.24 (EPIC) and 0.18 (450k) of the variance of age.
- **Gap, relative entropy (EPIC; 450k).** Housekeeping −0.106 (−0.112 to −0.097); −0.105 (−0.113 to −0.095). With
  the gene's expression decile: −0.031 (−0.037 to −0.022), 71 % removed; −0.031 (−0.039 to −0.022), 71 %. With
  keratinocyte chromatin state (Roadmap E057): −0.052 (51 %); −0.055 (47 %). With both: −0.019 (82 %); −0.021 (80 %).
  Core-essential −0.058 → −0.003 (−0.016 to +0.012); −0.054 → −0.001. Non-essential +0.218 → +0.088; +0.203 → +0.081.
  Entropy slope, housekeeping: −0.225 → −0.100; −0.237 → −0.099. Dispersion slope: −0.189 → −0.077; −0.191 → −0.068.
- **Outside every gene set** (9,761; 10,014 genes): slope relative to matched CpGs from +0.17 (lowest expression
  decile) to −0.14 (highest) in both groups; rank correlation −0.071 (p 2·10⁻¹²) and −0.086 (p 7·10⁻¹⁸); entropy
  slope −0.19 and −0.21.
- **Chip held.** Housekeeping −0.102 → −0.028 (72 %) on EPIC, −0.110 → −0.034 (70 %) on 450k; core-essential to
  zero in both; gradient −0.067 and −0.096.
- **Join check.** Across women, per-gene correlation of promoter methylation with expression: mean r −0.003
  (p 0.87) on EPIC and −0.028 (p 0.045) on 450k; with chip held −0.022 (p 0.15) and −0.024 (p 0.025). Weak, as in
  monocytes (−0.007); the join is carried by the gradient across genes.
- Reading: a third tissue and the first with both layers from the same people on rebuilt arrays in two
  independent groups of women. The two groups agree to the second decimal. Expression of the gene in the same
  tissue accounts for about 70 % of the housekeeping gap and all of the core-essential one; adding chromatin state
  brings it to about 80 %; a remainder of −0.02 SD stays, its interval excluding zero. Limits: one laboratory,
  women only, cell composition not held, chromatin reference from cultured keratinocytes.

## R127 — Lunnon's tissues normalised together from the raw intensities (RM43 item 7), 2026-11-20

- Commit 52ce0b3 and the reader fix after it; `ENTROPY_AGE_LUNNON=joint figures tissues` (job 6541264, test
  partition, 2 h 11 min with the build and smoke run) and `figures chromatin` (job 6541910, itmz, 30 min),
  `--scale paper --seed 42`. GEO's `GSE59685_signal_intensities.csv.gz` (methylated, unmethylated, detection p per
  array): methylated and unmethylated intensities quantile-normalised across all 531 arrays within probe design,
  beta = M / (M + U + 100), detection p > 0.01 masked. The deposited betas were dasen-normalised per tissue. Against
  them, per sample: r 0.999, mean absolute difference 0.01.
- **Same donors, blood minus brain region** (housekeeping dispersion gap; joint, then deposited): frontal cortex
  −0.105 (−0.127 to −0.081) against −0.092; superior temporal gyrus −0.104 (−0.127 to −0.082) against −0.096;
  entorhinal cortex −0.081 (−0.104 to −0.054) against −0.072; cerebellum −0.095 (−0.119 to −0.073) against −0.087.
- **Gap within each tissue** (island x level): blood −0.179 against −0.146; frontal cortex −0.069 against −0.050;
  superior temporal gyrus −0.068 against −0.044; entorhinal cortex −0.094 against −0.070.
- **Chromatin matching.** The region's own state removes 12 % (frontal, E073), 18 % (temporal, E072) and 18 %
  (entorhinal, E071) of the gap, against 5, 8 and 10 % on deposited betas; blood's state removes 35 to 40 % in the
  brain regions (25 to 27 %); in blood its own state removes 71 % (67 %).
- Reading: normalising each tissue separately did not make the blood-brain difference; under one normalisation it
  is slightly larger. This is not a rebuild (no IDATs exist) and slide still goes with tissue, so the same-donor
  contrast stays batch-confounded and is reported as such.

## R128 — the human groups of the species analysis with the human annotation (RM43 item 5), 2026-11-20

- Commit 0ccea8d; `ENTROPY_AGE_SPECIES_ANNOTATION=human figures species --scale paper --seed 42`, job 6540997 (test
  partition), 2.5 min. The species analysis reads promoters, islands and gene sets from the mouse genome for every
  species. Here the 13 human tissue groups (1,075 samples, 20 to 282 per tissue) are redone with the human
  annotation, on the 2,843 autosomal promoter probes the mammalian array shares with the 450k array (housekeeping
  401 CpGs, core-essential 44, non-essential 114; 4,640 probes under the mouse annotation). No chip terms, as R68c.
- Housekeeping, relative entropy: negative in 9 of 13 tissues under either annotation; median −0.049 SD (human)
  against −0.024 (mouse); correlation across tissues 0.81. Entropy slope: 11 of 13 under both; −0.097 against
  −0.113; correlation 0.92. Largest in skin (−0.205) and muscle (−0.147); positive in lung, spleen, adipose.
- Core-essential, relative entropy: 10 of 13 negative, median −0.041 (human, 44 CpGs) against 11 of 13 and −0.139
  (mouse); correlation 0.56. Too few CpGs under the human annotation to size it.
- Reading: the human rows of the species result do not depend on borrowing the mouse annotation; for
  housekeeping promoters the two agree tissue by tissue. Deposited betas, as the whole analysis.

## R126 — clonal haematopoiesis and the writer syndromes with the bisulfite-conversion index held (RM43 item 4), 2026-11-20

- Commit 6230082; `ENTROPY_AGE_BETAS=sesame ENTROPY_AGE_GCT=1 figures chip` and `figures syndromes`, `--scale paper
  --seed 42`, job 6540787 (itmz), 42 min with the smoke runs. Each array's SeSAMe conversion index (GCT) enters
  the group contrast as a covariate. An exclusion cut is not possible: every array of the three series is above
  1.2 (medians 1.59, 1.47, 1.39; about 1.07 in Jaffe and EPIC-Italy).
- **Clonal haematopoiesis (GSE315366).** The index does not differ by group (any CH 1.550 against 1.550, p 1.0;
  DNMT3A p 0.93; TET2 p 0.77), and every effect is the same to the third decimal. The nulls of R99 stand.
- **Writer syndromes: the index goes with case status.** GSE74432: Sotos 1.534 against 1.420 in controls (+0.93 SD,
  p 2·10⁻⁵), NSD1 variants +0.67 SD, Weaver +0.76 SD. GSE97362: Kabuki 1.258 against 1.413 (−2.30 SD, p 2·10⁻²⁰),
  KMT2D variants −0.94 SD, CHARGE −0.47 SD (p 0.04), CHD7 variants +0.08. Cases and controls were not converted
  alike.
- **Effects with the index held** (housekeeping F, 5 components; without → with): Sotos −0.066 → −0.068
  (p 2·10⁻¹³), NSD1 variants −0.053 → −0.058, Weaver −0.054 → −0.054 (p 0.003), CHARGE −0.019 → −0.019 (p 0.007),
  Kabuki −0.010 → −0.007 (p 0.47), KMT2D variants −0.001 → +0.000. Core-essential: Kabuki −0.022 → −0.027
  (p 0.009), Weaver −0.035 → −0.036.
- Reading: the syndrome effects do not move when the index is held, so incomplete conversion as such does not
  make them. But an index that differs by 0.5 to 2.3 SD between cases and controls says the groups were processed
  in different batches, and a linear covariate cannot remove everything a batch does. The syndrome results keep
  their place with that caveat stated; the comparison stays within each series.

## R116 — CR blood on the mouse recipe (RM43 item 4), 2026-11-20

- `ENTROPY_AGE_BETAS=sesame figures cr_blood_replication --scale paper --seed 42`, job 6537808 (itmz, queued 2 h,
  ran 3 h 19 min). GSE282661 rebuilt with SeSAMe's mouse steps (strain mask, detection over the species mask): all
  45 arrays pass, 66,406 promoter probes (67,065 under the first recipe, R98).
- CR against controls (8 against 8, two chips with both), within-chip refitting randomisation, 720 assignments,
  exact: core-essential F −0.077 (p 0.40), housekeeping F −0.036 (p 0.47); with 3 composition components −0.073
  (p 0.16) and −0.036 (p 0.22). TLS+CR: −0.040 (p 0.45) and −0.021 (p 0.57). Long-term TLS: +0.002 and −0.008
  (p ≥ 0.83). The 12-week contrast of the controls has no chip overlap and no p.
- The recipe changes nothing against R98 (−0.078, p 0.36; −0.037, p 0.46). R56's CR finding stays withdrawn.

## R124 — fibroblasts in culture: drift per doubling against expression of the same cultures (RM40 item 4), 2026-11-19

- Commit c6f938d; `ENTROPY_AGE_BETAS=sesame figures sturm_expression --scale paper --seed 42`, job 6539446 (test
  partition), 3 min with the smoke run. Sturm GSE179847 rebuilt from IDATs: 89 untreated cultures of 6 healthy
  lines at ambient oxygen (HC1 31, HC2 29, HC3 10, HC4 7, HC5 7, HC6 5; 0 to 79 doublings), 83,740 promoter CpGs
  measured in all. Per-CpG slopes on doublings with a line intercept. Expression: mean log2(TPM + 1) over the 77
  untreated healthy-line RNA-seq samples of GSE179848 at 21 % oxygen (28,633 genes).
- **Gap per doubling, relative entropy.** Housekeeping −0.102 SD (−0.109 to −0.094); with the gene's expression
  decile in the matching −0.028 (−0.036 to −0.021), 72 % removed. Core-essential −0.055 → −0.006 (−0.024 to
  +0.011). Non-essential +0.113 → +0.023 (−0.043 to +0.075). Entropy slope: housekeeping −0.208 → −0.053.
- **Outside every gene set** (10,662 genes): slope relative to matched CpGs by expression decile from +0.106
  (lowest) to −0.100 (highest); rank correlation −0.059 (p 1·10⁻⁹); entropy slope −0.255, dispersion slope −0.164.
- Reading: the monocyte pattern (R122) appears again where time is cell divisions in a dish: promoters of
  expressed genes drift less per doubling, and expression in the same cultures accounts for about three quarters
  of the housekeeping gap and all of the core-essential one.
- Against R94: there the per-sample coupling of the shield to matched drift was near zero on the rebuilt arrays.
  That statistic asks whether a culture with more drift has proportionally more shield; this one compares per-CpG
  slopes over all cultures. They are not the same question, and this run does not restore R94's claim.
- **With chip held (R125**, commit ce2fe91, `ENTROPY_AGE_CHIP=1`, job 6539692, 3 min). The 89 cultures sit on 50
  chips, 27 of them with two or more cultures; a third of the variance of doublings (0.325) is left within line and
  chip, so the slopes are estimable but on less. Relative entropy: housekeeping −0.102 (−0.109 to −0.094) → −0.024
  (−0.031 to −0.016) with the expression decile, 77 % removed; core-essential −0.053 → −0.002; non-essential +0.164
  → +0.056. Gradient outside the sets +0.141 to −0.106, rank correlation −0.061 (p 2·10⁻¹⁰). The result does not
  come from chip.
- Limits: HC1 and HC2 give two thirds of the cultures; expression is a culture-wide mean, not paired sample by
  sample; the chip-held slopes rest on a third of the spread in doublings.

## R123 — the same in CD4+ T cells (RM40 item 4), 2026-11-19

- Commit 31b8d82; `figures mesa_expression --scale paper --seed 42`, job 6539223 (test partition), 27 min with the
  smoke run. CD4+ T cells of 214 people (45-79 years), GSE56581 (450k, deposited) with GSE56580 (HT-12), joined on
  the participant id. The monocyte rows reproduce R122 to the last digit.
- Join check: mean r −0.0085 over 15,885 genes (51 % negative; permutation p 0.005).
- **Relative entropy:** there is almost no shield to explain. Housekeeping −0.015 SD (−0.024 to −0.006), −0.013
  with the expression decile (9 % removed), −0.010 with the chromatin state; core-essential −0.006 (−0.029 to
  +0.018). Outside the gene sets no gradient with expression (rank correlation +0.02; deciles −0.00 to +0.07).
- **Entropy slope:** housekeeping −0.219 → −0.069 with the expression decile (69 % removed, 66 to 72), −0.031 with
  the state, −0.019 with both; core-essential −0.091 → +0.007. Outside the gene sets the entropy slope falls with
  expression (+0.13 in the lowest decile to −0.18 in the highest; rank correlation −0.197, p 4·10⁻¹⁰²); the
  dispersion slope barely (−0.054).
- Reading: where T cells show a shield at all (entropy slope, which carries the floor effect of R68), expression
  in the same cells accounts for it as in monocytes. On the corrected measure T cells have next to none, in a
  series six times smaller and with a narrower age range; the monocyte gradient is not reproduced here and this
  cannot say whether that is the cell type or the sample size.

## R122 — the shield against expression measured in the same cells: MESA monocytes (RM40 item 4), 2026-11-19

- Commit 4ae98cb; `figures mesa_expression --scale paper --seed 42`, job 6539095 (test partition), 17 min with the
  smoke run. CD14+ monocytes of 1,202 people (44-83 years), methylation GSE56046 (450k, deposited, no IDATs) and
  expression GSE56045 (HT-12 array, deposited normalised values), joined one to one on the participant id; 109,137
  promoter CpGs, 20,716 genes with an expression probe. First use in this project of a second layer measured in
  the same samples as the methylation.
- **Join check.** Promoter methylation against expression across people, per gene, age and covariates held: mean
  r −0.0068 over 15,890 genes (52.5 % negative; permutation null 0.0000 ± 0.0007, p 0.005). Negative, as it must
  be, and very small: in one sorted cell type few promoters vary between people.
- **Matching on the gene's own expression decile.** Housekeeping gap on relative entropy −0.099 SD (−0.107 to
  −0.089) with island x level; −0.017 (−0.026 to −0.008) with expression added: 82 % removed (75 to 92). With the
  mononuclear chromatin state instead +0.001; with both +0.017 (+0.010 to +0.026), the sign reversed.
  Core-essential: −0.035 → +0.020 (−0.003 to +0.044) with expression. Entropy slope, housekeeping: −0.137 → −0.040
  (71 % removed). Non-essential promoters, which drift more than matched CpGs: +0.218 → +0.124 with expression
  (43 % removed), +0.073 with both.
- **Outside every gene set** (11,678 genes): the slope relative to the island x level stratum falls with
  expression decile, from +0.072 SD (lowest) through −0.10 (deciles 5 to 9) to −0.199 (highest); rank correlation
  across genes −0.117 (p 5·10⁻³⁷); dispersion slope +0.106 to −0.247 (−0.167).
- **Where the sets sit.** 75 % of housekeeping and 82 % of core-essential CpGs are in expression deciles 6 to 10,
  against 27 % of CpGs outside the sets; 86 % of non-essential CpGs are in deciles 1 to 3.
- Reading: in one cell type, with expression from the same samples, the gene's expression level accounts for four
  fifths of the housekeeping shield and all of the core-essential one, and genes in no set are shielded in step
  with their expression. This is the direct form of the activity reading (R108) and removes its weakest point,
  reference maps from other people's cells. Limits: deposited values on both layers; one array of each; the cell
  type is the one where turnover selection would act on precursors, not on the cells measured.

## R121 — blood with the neutrophil chromatin map (RM43 item 3), 2026-11-19

- Commits 5a835c6 and the figure fix after it. `ENTROPY_AGE_BLOOD_STATE=E030 figures equilibrium_drift`,
  `drift_hotspots`, `essentiality_gradient` (jobs 6538747, 6538749, 6538751; 35, 3 and 7 min) and `figures chromatin`
  with a new matching on both blood maps at once (job 6538876); `--scale paper --seed 42`, test partition, deposited
  Hannum betas (no IDATs). Blood rows only change; each is given as neutrophils (E030) against mononuclear cells
  (E062, the earlier run).
- **Gain and loss rates (R74).** Active-TSS island promoters, set minus matched, loss per decade: core-essential
  +0.0053 (+0.0026 to +0.0081) against +0.0063 (+0.0031 to +0.0097); housekeeping +0.0041 (+0.0022 to +0.0061)
  against +0.0015 (+0.0001 to +0.0034), on 18,137 CpGs instead of 7,256. Per-CpG resting level with the state in the
  strata: housekeeping −0.076 against −0.060 SD, core −0.066 against −0.062 (from −0.335 and −0.243 unmatched).
- **Drift hotspots (R72).** Same pattern: Polycomb-repressed 2.02 against 2.11 times expected (dispersion slope),
  bivalent TSS 1.78 against 1.99, active TSS 0.32 against 0.59, flanking active TSS 0.56 against 0.40.
- **Essentiality gradient (R58).** Per unit of essentiality with the state held: entropy slope +0.054 (0.041 to
  0.071) against +0.046; relative entropy +0.007 (p 0.28) against +0.004 (p 0.54); with expression and breadth
  added +0.023 against +0.022 (entropy), −0.009 against −0.007 (relative entropy).
- **Chromatin matching (R104).** Housekeeping gap on relative entropy: −0.147 unmatched; −0.013 with E062 (91 %
  removed), −0.025 with E030 (83 %), −0.006 with both maps held (96 %; interval reaches zero). Core-essential:
  −0.085; −0.001, −0.006, −0.001. Lunnon whole-blood dispersion, housekeeping: 67 %, 66 %, 78 % removed.
- Reading: no blood conclusion depends on which of the two blood maps is used. The neutrophil map leaves slightly
  more of the gap than the mononuclear one; holding both leaves none on relative entropy.

## R120 — Tet3 against Tet1 on shared CpGs and shared genes (RM43 item 8), 2026-11-19

- `ENTROPY_AGE_BETAS=sesame figures tet_shared --scale paper --seed 42`, commit 7537767, jobs 6538512 and 6538532
  (test partition), 38 s, 1.2 GB. Tet1: rebuilt 285k mouse arrays (6 mutants against 7 wild types per arm). Tet3:
  deposited 40k mammalian arrays, 8 heterozygous against 8 same-chip wild types per tissue.
- **The arrays share almost no CpGs.** By mm10 position 150 CpGs sit on both arrays, 42 at promoters; 37 promoter
  CpGs are measured in every animal of both series, 5 of them in the housekeeping set. On those 37 no score can be
  estimated (housekeeping F: Tet1 knockout +0.21, SE 0.26; Tet3 cortex +1.9, SE 0.68). The shared-CpG comparison the
  review asked for cannot be made with these arrays.
- **On shared genes** (each array's own CpGs, the gene sets cut to the 314 genes with promoter CpGs on both arrays):
  housekeeping F, Tet1 knockout −0.190 (SE 0.015, p 6·10⁻⁸), catalytic-dead −0.193 (p 3·10⁻⁵), on 1,034 set CpGs;
  Tet3 heterozygous cortex +0.169 (p 0.009), striatum +0.246 (p 5·10⁻⁶), on 553. Core-essential: Tet1 −0.20 and
  −0.22, Tet3 +0.25 and +0.53. On each array's full set: Tet1 −0.127 and −0.122, Tet3 +0.168 and +0.245.
- Reading: the opposite signs are not made by the two arrays covering different genes. What still differs: the
  CpGs within those genes, array chemistry and processing (rebuilt against deposited), laboratory, heterozygous
  against homozygous loss, and age. Both mutants raise the matched CpGs' deviation B (Tet1 +0.02 to +0.03, Tet3
  +0.005 to +0.009).

## R118 — the species shield with chip intercepts (RM43 item 5), 2026-11-19

- `ENTROPY_AGE_CHIP=1 figures species --scale paper --seed 42`, job 6538040 (test partition), 1 h 07 min. Deposited
  consortium betas (no IDATs). Chip intercepts in every age slope; 17 of 139 groups left out for too few animals per
  chip, leaving 122 groups of 57 species.
- Housekeeping relative-entropy gap negative in 86 % of groups (88.5 % of the same groups without chip terms),
  median −0.108 (−0.126), interval below zero in 79 groups and above in 3; correlation of the two fits across groups
  0.93. Core-essential: negative in 80 %, median −0.129. Dispersion slope: 81 %, median −0.111.
- Lifespan: the gap is still not tuned to maximum lifespan (core-essential entropy gap, p 0.93); the matched CpGs'
  own rate still falls with lifespan (p < 0.001).
- **Mouse tissues within chips** (numeric ages, 19 to 903 animals per tissue): lung −0.49 (−0.54 to −0.44), spleen
  −0.45, liver −0.35, blood −0.32, skin −0.22 (−0.32 to −0.11); whole brain −0.13, fibroblasts −0.12; muscle −0.08,
  cortex −0.06, heart −0.04, kidney −0.04 (intervals across zero); hippocampus +0.03, cerebellum +0.04, striatum
  +0.10 (+0.03 to +0.17).
- Reading: the species result does not come from chip going with age. And with chip held, mouse tissues on the
  consortium array order as colon and hippocampus do by RRBS (R114): large in renewing tissues and liver, absent in
  brain regions and heart. This is the within-chip turnover check asked for under RM40 item 1, on a second platform
  and a third dataset. (Kidney sits with the post-mitotic tissues; liver with the renewing ones.)

## R117, R119 — two checks from the methods review: complete CpGs only; SATSA on one platform (RM43 items 9, 6), 2026-11-19

- **R117** (`ENTROPY_AGE_BETAS=sesame ENTROPY_AGE_COMPLETE=1 figures tissues`, job 6538038, test partition, 57 min):
  the main tissue comparison on CpGs measured in every sample, nothing median-filled. Only the rebuilt cortex had
  missing values (94,513 → 89,936 CpGs); Hannum and Lunnon are complete as deposited. Every gap agrees with the filled
  run to the third decimal (cortex housekeeping relative entropy −0.082 and −0.082, dispersion −0.118 and −0.119;
  core-essential −0.073 and −0.099 unchanged). The median fill does not move the headline gaps.
- **R119** (`figures within_person_450k`, job 6538161, 5 min): deposited SATSA on its 450k draws only, 1,094 draws of
  447 people, 647 consecutive pairs (pooled run: 873 pairs), no chip term. Per-person shield share F, housekeeping,
  per decade within a person: −0.126 (−0.226 to −0.030), against −0.142 (−0.228 to −0.054) pooled; core-essential
  −0.053 (−0.161 to +0.048) against −0.073. Signed shift S, housekeeping −0.135 (−0.229 to −0.034). The within-person
  growth of the shield share does not depend on pooling the two arrays. The small CpG-level divergence gap does:
  housekeeping −0.017 (p 0.001) pooled, −0.003 (p 0.41) on 450k alone; the measurement-floor gap holds (−0.036).

## R114 — colon against hippocampus by RRBS in one mouse cohort (RM40), 2026-11-19

- `figures jena_rrbs --scale paper --seed 42`, job 6537454 (test partition), 46 min after a smoke pass in the same
  job; two earlier attempts failed in the readers. GSE233734 (colon, 82 mice) and GSE304754 (hippocampus, 83 mice):
  C57BL/6J/Ukj males, 16-18 per age at 3, 9, 15, 24 and 28 months; processed per-CpG tables; promoter CpGs covered
  10 times in 90 % of a tissue's mice: 439,323 and 462,430. 200 mouse resamples per tissue, tissues resampled
  independently.
- **Batch:** our audit of both series (`data/input/geo/audit/`): the hippocampus libraries ran on four flowcells,
  the same four that carry the colon libraries (a fifth holds colon only), and every flowcell holds every age bin
  in both tissues (R2 of age by flowcell 0.02 and 0.01). Age and tissue are both free of flowcell.
- **Gaps (relative entropy):** housekeeping, colon −0.189 (−0.294 to −0.170), hippocampus −0.078 (−0.094 to −0.065);
  core-essential −0.166 and −0.071; non-essential +0.537 and +0.155.
- **Colon minus hippocampus:** housekeeping −0.111 (−0.222 to −0.085), core-essential −0.095 (−0.196 to −0.068),
  non-essential +0.382 (+0.312 to +0.660); p 0.005 each (the floor of 200 resamples). Entropy slope and dispersion
  slope agree (housekeeping −0.144 and −0.179).
- **Reading:** the shield is more than twice as large in the renewing tissue as in the post-mitotic one, in a
  sequencing series where neither age nor tissue is tied to batch. This is the turnover ordering that was withdrawn
  for the mouse array (R92) because of chip confounding, now seen on clean data. Limits: two tissues, one strain,
  males; hippocampus is bulk tissue and contains dividing glia; a one-base uncertainty in the hippocampus
  coordinates could not be resolved by overlap (222,930 shared CpGs unshifted, 223,304 shifted) and does not matter
  for promoter membership; sporadic missing values (up to 10 % per CpG) are median-filled as in the other RRBS runs.

## R115 — glia against neurons in a second sorted series, GSE66351 (RM40), 2026-11-19

- `figures nuclei_paired_gse66351 --scale paper --seed 42`, job 6537455 (test partition), 15 min; two earlier
  attempts failed on a column name. Rebuilt from IDATs (189 of 190 arrays pass; 450k). 31 donors with both NeuN+ and
  NeuN− nuclei (16 controls, 15 Alzheimer's), 18-94 y, only 2 under 60; sex and diagnosis held; 92,035 complete
  CpGs; 400 paired donor resamples.
- Glia minus neurons, relative entropy: housekeeping −0.043 against +0.001, difference −0.044 (−0.108 to +0.015),
  p 0.11; core-essential −0.019 against +0.007, difference −0.026 (−0.088 to +0.036), p 0.38. Absolute slope:
  −0.035 and −0.032 (p 0.34, 0.48).
- Reading: the same direction and nearly the same size as R110 (−0.055) in an independent series on another array,
  not significant alone, as expected with an age slope resting on two young donors. Two series now agree that glia
  carry more of the housekeeping shield than neurons; neither is strong alone.

## R113 — the breast-cancer lead in a second series, GSE283951 (C12, RM40), 2026-11-19

- `figures breast_replication --scale paper --seed 42`, job 6536367 (test partition), 23 s. Rebuilt from IDATs (85
  of 87 arrays pass); 32 cases, 53 controls, mean age 54.0 and 52.8; 16 chips, 9 holding both groups; 82,067 complete
  promoter CpGs; scores against controls, centred within chip.
- Odds ratio per SD, age held: B housekeeping 1.11 (0.71 to 1.73), p 0.66; B core-essential 1.09 (0.70 to 1.70);
  F 0.88 (0.56 to 1.38); S housekeeping 0.73 (0.44 to 1.23). With mean methylation held: unchanged (B 1.11, 0.68 to
  1.79).
- Reading: no association. The interval holds both 1.0 and the 1.7-1.8 first seen in EPIC-Italy, so 32 cases neither
  replicate nor refute; the point estimate equals EPIC-Italy's on rebuilt arrays with adjustment (1.11). Blood here
  was not shown to be drawn before diagnosis, and cell fractions are not held. C12 stays weak, with no support from
  a second series.

## R112 — drift and the shield in non-dividing fibroblasts over six months, E-MTAB-15859 (RM40), 2026-11-19

- `figures arrest_timecourse --scale paper --seed 42`, commit after ccb2d35 (arrest_timecourse), job 6536205 (test
  partition), 3.5 min. Rebuilt from IDATs (79 of 79 arrays pass; EPIC v1 and v2 in one study). Used: untreated
  confluent fibroblasts on EPIC v1, three donors, six monthly arrays each, each donor within one chip; 76-80 thousand
  complete promoter CpGs. Donor 2 (EPIC v2, months split over three chips) left out.
- **Drift** (mean squared slope per month against all 720 reorderings of the months): donor 1 0.88 of the null
  median (p 0.67), donor 4 0.86 (p 0.64): none. Donor 3 1.96 (p 0.007), but in donor 3 month runs down the chip in
  order (correlation with array row 0.99; −0.15 and −0.14 in donors 1 and 4), so its trend cannot be told from
  array position.
- **Shield:** matched gaps of the absolute slope and of the relative-entropy slope are within ±0.02 SD in every
  donor and set, intervals across or touching zero (housekeeping: −0.010, +0.013, +0.004).
- **Reading:** six months without division leave promoter methylation where it was, within what six arrays per
  donor can see; the noise level is about 0.004 per month in beta, so a change of the size seen in replicative
  senescence (R111: about 0.1 in beta at mid-level CpGs) would have been plain. With R111 this says that promoter
  drift in fibroblast culture goes with division, not with time, and that where there is no drift there is no
  shield to measure. Limits: three donors, six arrays each, no replicates; two series from different laboratories
  are being read together; replicative and arrested cultures are not from the same experiment.

## R111 — replicative against X-ray and RAS senescence within chips, GSE197723 (RM40), 2026-11-19

- `figures senescence_modes --scale paper --seed 42`, commit d67772f, job 6536147 (test partition), 18 s; first pass
  6534668 without the noise floor (ccb2d35). Rebuilt from IDATs (313 of 316 arrays pass); fibroblasts on 12 chips,
  each with 2 untreated, 2 replicative-senescent, 2 X-ray and 2 RAS cultures; 86,856 complete promoter CpGs; reference
  per chip = its two untreated cultures; one value per chip and arm; exact sign-flip tests over 12 chips.
- **Uncorrected scores mislead.** Matched-CpG deviation B: replication 0.126, X-ray 0.032, RAS 0.039; shield share F:
  −0.33, −0.03, −0.08 (housekeeping). The low F of the arrested arms is array noise in the denominator.
- **With the noise floor subtracted** (squared deviations; floor from each chip's untreated pair): real drift B2 is
  0.090 (0.079 to 0.101) in replicative senescence and 0.000 to 0.001 in the X-ray and RAS arms, that is, none
  measurable. Shield in replicative senescence: P2 −0.052 (housekeeping) and −0.033 (core-essential), **F2 −0.58
  (−0.60 to −0.55) and −0.56 (−0.60 to −0.52), negative on 12 of 12 chips.** In the arrested arms P2 is 0.000 (±0.001)
  and F2 is undefined (intervals of several units).
- **Reading:** (1) cultures that reached senescence by dividing show a large shield on rebuilt, chip-balanced arrays:
  essential and housekeeping promoters carry less than half the squared drift of matched CpGs. This is the first
  culture result that is clean of chip confounding. (2) Cultures arrested by X-rays or RAS show no drift from their
  untreated chip-mates, so there is nothing for a shield to act on and the division-against-arrest comparison of
  shield share cannot be made in this series. (3) How long the arrested cultures were kept, their donors and the
  doublings of the replicative arm are not in the record, so "drift needs division" is not shown either: the arrested
  cultures may simply be younger in culture time. E-MTAB-15859 (confluent cultures held up to 6 months) is the series
  that can separate time from division.

## R110 — glia against neurons in the same donors (RM40), 2026-11-19

- `figures nuclei_paired --scale paper --seed 42`, commit ccb2d35, job 6534603 (test partition), 27 min after a
  smoke pass in the same job. GSE279509 rebuilt; donors with both fractions passing QC; CpGs complete in both; sex
  held; 400 donor resamples that take both of a donor's arrays.
- **Other glia (NeuN−/SOX10−) minus neurons, 65 donors, 18-92 y:** housekeeping relative entropy −0.084 against −0.029,
  difference **−0.055 (−0.121 to −0.007), p 0.025**; core-essential −0.054 against −0.038, difference −0.016 (−0.057
  to +0.017), p 0.27. Absolute slope of the mean: housekeeping −0.049 (−0.097 to −0.008), p 0.02; core-essential
  −0.038 (−0.083 to +0.004), p 0.09.
- **Oligodendrocyte lineage (SOX10+) minus neurons, 101 donors:** housekeeping −0.021 against −0.018, difference
  −0.003 (−0.044 to +0.037), p 0.75; core-essential −0.007, p 0.72. Absolute slope: −0.017 (p 0.045) and −0.015
  (p 0.05).
- Reading: one of the two glial fractions is more shielded than the same donors' neurons, for housekeeping genes, at
  p 0.025 among four primary tests; the oligodendrocyte lineage is not distinguishable from neurons on relative
  entropy. R109's ordering (other glia > oligodendrocyte lineage > neurons) is kept in sign, but only its extreme
  pair is separated. Suggestive, not established; a second sorted series (GSE107729 WGBS, GSE66351) is the check.

## R102, R105, R109b — cortex on rebuilt arrays in the main comparison; blood minus cortex matched; GSE252501 with depression held, 2026-11-19

- **R102** (`ENTROPY_AGE_BETAS=sesame figures tissues`, job 6520488, test partition, 1 h 02 min; Jaffe rebuilt, 208
  donors; Hannum and Lunnon deposited, flagged). Cortex housekeeping gap: relative entropy −0.058 → −0.082,
  dispersion −0.075 → −0.118; core-essential −0.072 → −0.073 and −0.081 → −0.099; non-essential 0.00 → +0.03 (p
  0.11). Blood minus cortex with the rebuilt cortex (mixed processing): housekeeping −0.066 (−0.145 to +0.005),
  core-essential −0.017 (−0.084 to +0.051): not resolved.
- **R105** (`figures tissues_matched`, job 6520487, 1 h 07 min; both deposited; ages 20.8-85.2, so 616 of 656 blood
  donors and all 213 cortex donors; 11 of 12 EpiDISH fractions held in blood). Blood minus cortex, housekeeping:
  relative entropy −0.067 (−0.134 to −0.014; unmatched −0.093), dispersion −0.103 (−0.169 to −0.038; unmatched
  −0.209); core-essential +0.001 and −0.023 (intervals across zero). Matching ages and holding cell fractions removes
  a third to a half of the difference; the rest is still a contrast of two studies.
- **R109b** (job 6521034, 14 min; commit b8f4e69): GSE252501 DLPFC, 52 donors without alcohol use disorder,
  depression diagnosis held: housekeeping −0.005 (−0.017 to +0.006), −0.013 (−0.024 to −0.003) with chip;
  core-essential −0.015 and −0.020 (intervals exclude zero). The 23 with neither diagnosis: housekeeping −0.018
  (−0.028 to −0.008), core-essential −0.015 (−0.028 to −0.001). Holding the diagnosis changes nothing. The sorted
  nuclei groups are as in R109.
- Reading: "blood more shielded than cortex" is at most −0.07 SD for housekeeping genes and absent for core-essential
  genes, in a cross-study contrast whose interval touches zero once cortex is rebuilt. Not a finding.

## R109 — the shield in two new brain series, sorted nuclei included (RM40, RM44), 2026-11-19

- `figures brain_replication --scale paper --seed 42`, commit b9b866f, job 6520486 (test partition), 13 min. Both
  series rebuilt from IDATs (QCDPB), complete CpGs only (81-90 thousand), gene-resampling intervals. Audits of both
  with the external agent's tool are in `data/input/geo/audit/`. GSE279509: 524 of 749 arrays pass 90 % detection
  (2 unreadable); in NeuN+ detection falls with age (r −0.23). GSE252501: 232 of 232 pass; Jaffe's donors excluded.
- **GSE279509 (Walker 2025), housekeeping relative-entropy gap, sex held / sex and chip held:**
  unsorted nuclei (148 donors, 18-95 y) −0.053 (−0.062 to −0.043) / −0.050; other glia, NeuN−/SOX10− (89) −0.077
  (−0.086 to −0.068) / −0.067; oligodendrocyte lineage, SOX10+ (127) −0.034 (−0.045 to −0.020) / −0.036; **neurons,
  NeuN+ (130) −0.021 (−0.032 to −0.009) / −0.008 (−0.019 to +0.003).** Core-essential: unsorted −0.054 / −0.044,
  other glia −0.055 / −0.050, oligodendrocyte lineage −0.023 / −0.026 (intervals reach zero), neurons −0.032 (to
  +0.001) / −0.003. Off the floor (10-90 %) the glial and unsorted gaps hold (housekeeping −0.04 to −0.09), the
  neuronal ones do not survive chip intercepts.
- **GSE252501 (White 2024), 52 DLPFC donors without alcohol use disorder, 25-78 y:** housekeeping −0.005 (−0.016 to
  +0.006) / −0.014 (−0.025 to −0.005); core-essential −0.015 / −0.019 (intervals exclude zero). Nucleus accumbens
  (51): −0.009 / −0.024. Same sign as Jaffe, a quarter of its size or less. This run did not hold the
  major-depression diagnosis that half of these donors carry; R109b (job 6521034) does.
- Dispersion slopes with chip intercepts change sign between groups (unsorted +0.064 off the floor) and are not
  interpreted: with 72-101 chips for 89-148 arrays the chip terms absorb most of the between-person spread.
- **Reading:** the cortex shield replicates in an independent series of unsorted nuclei at the size first seen in
  Jaffe. By cell type it is largest in the glial fractions, which can divide, and smallest in neurons, which cannot;
  within chips the neuronal gap is not distinguishable from zero. That is the ordering a division-dependent cause
  (selection among cells, or replication-coupled maintenance) predicts, and not what "a shield in neurons would
  exclude selection" hoped for. Limits: one series; intervals resample genes, not donors; no direct test of the
  neuron-glia difference yet (the same donors and chips carry both, so a paired donor bootstrap is possible and is
  the next run).

## R108 — the joint test in cortex (RM44), 2026-11-19

- `ENTROPY_AGE_BETAS=sesame figures joint_cortex --scale paper --seed 42`, commit b9b866f, job 6520485 (test
  partition), 9 min; first attempt 6519456 failed in the cortex loader (an array that failed QC). Jaffe rebuilt, 208
  donors; 94,513 promoter CpGs, 26,871 at 10-90 % young level, 15,385 at 20-80 %. Under the fullest matching 90-99 %
  of each set's CpGs still have a control.
- **Relative entropy, housekeeping:** all CpGs −0.082 (−0.094 to −0.070) with island and level; −0.045 with
  expression; −0.046 with chromatin state; −0.047 with sorted accessibility; **−0.019 (−0.030 to −0.005) with all
  three.** At 10-90 %: −0.091 → −0.032 (−0.061 to +0.005). At 20-80 %: −0.031 → +0.015 (−0.032 to +0.075).
- **Core-essential:** all CpGs −0.073 → −0.035 (−0.056 to −0.018); 10-90 % −0.106 → −0.030 (−0.092 to +0.026);
  20-80 % −0.066 → +0.007.
- **Dispersion slope:** all CpGs housekeeping −0.118 → −0.041 (−0.053 to −0.029), core-essential −0.099 → −0.059; off
  the floor the gap is absent already under the base matching (10-90 %: −0.027 and +0.014).
- **Reading:** held together, expression, chromatin state and accessibility remove about three quarters of the
  housekeeping gap and half of the core-essential gap; what is left (−0.02 to −0.04 SD) sits at CpGs near 0 %, and
  off the floor nothing distinguishable from zero remains. In cortex too, the supported statement is "constitutively
  active promoters drift less"; an essentiality-specific shield beyond activity is not shown. The three reference
  tracks come from other people than the methylation.

## R107 — Huang mouse tissues on SeSAMe's mouse recipe (RM43 item 4), 2026-11-19

- `ENTROPY_AGE_BETAS=sesame figures mouse`, commit 4fd2503 (with the paired-chip dispersion fix), job 6518430 (itmz),
  2 h 38 min; rebuild TQCDPB, 581 of 583 arrays. Housekeeping relative entropy, renewing − post-mitotic: −0.191
  (−0.507 to −0.066) unadjusted, +0.023 (−0.054 to +0.182) within chip; slow − post-mitotic −0.099 and +0.031. As R97
  and R92: the turnover ordering is not separable from chip in this series.

## R98 — CR blood (GSE282661) on R SeSAMe betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures cr_blood_replication --scale paper --seed 42`, commit d4997f3, job 6515170
  (itmz), 3 h 04 min. First recipe (QCDPB); all 45 arrays pass; 67,065 promoter probes.
- CR against controls, within-chip refitting randomisation (720 assignments, exact): core-essential F −0.078
  (p 0.36), housekeeping F −0.037 (p 0.46); deposited −0.09 (p 0.16) and p 0.23. TLS+CR and long-term TLS: all
  p ≥ 0.40. The 12-week aging contrast still has no chip overlap and no p.
- R56's CR finding stays withdrawn; the rebuild weakens it further. Paper and Table S1 name the rebuild.

## R103b — QC, batch and sample distributions over the rebuilt series (RM41, RM42), 2026-11-18

- `figures preprocessing_qc`, job 6518679 (test partition), 6 min; R103 (6517574) failed on mixed column types.
  15 series (GSE282661 still on the first recipe).
- **Batch share of the 10 leading components explained by chip:** Endicott GSE197512 0.79, Browder 0.61, EPIC-Italy
  0.55, Jaffe 0.51 (adult controls 0.31), GSE74432 0.41, CH 0.25, BXD 0.20, Huang 0.20, Sturm 0.15, others ≤ 0.16.
  Row and position ≤ 0.11 everywhere.
- **Design confounded with chip (share of variance):** Endicott population doublings 0.77, days in culture 0.81, cell
  line 0.81, batch date 1.00; Jaffe group 0.77 (full series); EPIC-Italy sex 0.62. Endicott's passage contrast is
  largely a chip contrast.
- **Sample distributions** (distance of each array's beta quantile curve from its series median; flagged beyond
  median + 5 MAD): Jaffe full 48 of 660, Browder 10 of 155, Tet1 7 of 26, EPIC-Italy 7 of 839, BXD 5 of 288,
  GSE74432 4 of 119; none in Endicott, Sturm, Huang, CR, CH, Dnmt3a, GSE97362. Flags are recorded, not yet used to
  exclude.
- **Rebuilt against deposited:** Huang identical (r 1.000); Endicott per-probe r 0.999; Sturm 0.976; Jaffe per-sample
  r 0.990 but per-probe r 0.815 (mean |difference| 0.05), so the authors' processing changed between-sample
  variation at single CpGs.

## R106 — BXD, Browder, Tet1, Dnmt3a on SeSAMe's per-platform mouse recipe (RM43 item 4), 2026-11-18

- Rebuilds with TQCDPB (MM285) and SHCDPB (Mammal40), detection over the species-mask probes, 90 % cut (4fd2503);
  `ENTROPY_AGE_BETAS=sesame figures bxd|reprogramming|tet1|dnmt3a`, job 6518427 (test), 66 s.
- **BXD:** 288 of 338 arrays pass, 230 females with 160 in the aging fit. Aging gap, relative entropy: housekeeping −0.279
  (−0.452 to −0.143), core-essential −0.282. Housekeeping F per 10 months −0.248, −0.109 with 5 components; HFD −0.086
  (p 0.001), **−0.035 with components (p 0.072)**: R100's "part survives composition" (p 0.011) does not hold under
  the mouse masks; unresolved. extension.md and Supplementary Note 3 corrected.
- **Browder:** spleen core-essential F +0.171 (p 0.007), housekeeping +0.109 (p 0.068); muscle housekeeping −0.107
  (p 0.066; deposited −0.115, p 0.019). Directions hold; muscle no longer below 0.05.
- **Tet1:** 19 of 24 BS arrays, F −0.12 to −0.17 (p < 0.001), as R99. **Dnmt3a:** F −0.54 to −0.65, as R99.

## R104 — chromatin matching with region-matched brain maps (RM43 item 3), 2026-11-18

- `figures chromatin --scale paper --seed 42`, commit 1316d9d, job 6517856 (test partition), 21 min. Lunnon superior
  temporal gyrus now matched with Roadmap E072 (inferior temporal lobe), entorhinal cortex with E071 (hippocampus
  middle), frontal cortex with E073; cerebellum left out (no Roadmap map); BA9 sorted ATAC used for frontal cortex
  and Jaffe only. Deposited betas, as before.
- Housekeeping dispersion gap removed by the region's own state: frontal cortex 5 % (E073), superior temporal gyrus
  8 % (E072; was 14 % with E073), entorhinal cortex 10 % (E071; was 8 %); by blood's E062 25-27 %. Blood (E062 91 %,
  E030 83 % on relative entropy) and Jaffe (E073 25 %, sorted ATAC 45 %, both 56 %) unchanged.
- Reading: a region-matched brain map explains no more of a brain region's gap than the prefrontal map did; blood's
  map still removes more there. In brain the reference state acts as a generic label of promoter activity (C2).
  Paper Results updated.

## R101 — muscle fibres on R SeSAMe betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures fibres --scale paper --seed 42`, commit f37945f (after the replicate-name fix),
  job 6517530 (test partition), 7.5 min; first attempt 6517106 failed on a replicate's name. 92 of 92 arrays mapped;
  12 people before training, 11 at 12 weeks; 88,247-88,549 CpGs.
- Whole muscle minus pure fibres, housekeeping: pre WM − T1 +0.008 (0.000 to +0.016), WM − T2 −0.006 (−0.013 to 0.000);
  12 weeks +0.007 and +0.002, intervals across zero (deposited: +0.004, −0.006, +0.001, +0.007). Still no difference
  between non-dividing fibres and whole muscle. Not in the paper's text; Table S1 names the rebuild.

## R100 — BXD liver and Browder cyclic OSKM on R SeSAMe betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures bxd|reprogramming --scale paper --seed 42`, commit fa2ce10, job 6516403 (test
  partition), 40 s.
- **BXD liver (GSE199979):** all 339 arrays pass; 3,927 promoter probes (deposited 4,221), 321 females. Aging gap,
  relative entropy: housekeeping −0.256 (−0.378 to −0.153), core-essential −0.277 (was −0.320, −0.315). Housekeeping
  F per 10 months −0.244 (was −0.283), −0.116 with 5 components (was −0.055). **HFD −0.097 (was −0.105), and −0.040
  (p 0.011) with components, where the deposited betas gave −0.017 (p 0.26):** part of the diet effect is not
  composition. extension.md and Supplementary Note 3 changed.
- **Browder OSKM (GSE190665):** 167 of 168 arrays in the rebuild; the same mice per tissue. Muscle F housekeeping
  −0.109 (p 0.043; was −0.115, p 0.019); spleen core-essential +0.158 (p 0.006; was +0.211, p 0.001), housekeeping
  +0.094 (p 0.075; was +0.104, p 0.046). Young − old skin −0.122 (p 0.004). Directions hold; spleen weakens.
- Manuscript: the BXD gap names the rebuild, and the RRBS-against-array turnover ranking is removed (RM41 item 9b:
  different platforms are reported side by side, not ranked).

## R99 — CH, Tet1 and Dnmt3a gain of function on R SeSAMe betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures chip|tet1|dnmt3a --scale paper --seed 42`, commit b25c777, job 6515945 (test
  partition), 2 min for all three.
- **CH (GSE315366):** all 64 patients pass; 87,669 CpGs (deposited 101,667). Still no effect: |F| ≤ 0.031, p ≥ 0.19
  for any CH, DNMT3A and TET2, raw or with 5 components (deposited |F| ≤ 0.022, p ≥ 0.15).
- **Tet1 (GSE271375):** 19 of 24 bisulfite arrays pass detection, 6 KO, 6 Hxd and 7 wild types (was 8 each).
  KO F core-essential −0.169/−0.128 (raw / 3 components), housekeeping −0.128/−0.094; Hxd −0.136/−0.134 and
  −0.123/−0.116; all F p ≤ 0.002. S −0.015 to −0.020 (p < 10⁻⁵). The shield strengthening holds (was F −0.11 to −0.17).
- **Dnmt3a W326R (GSE324160):** all 72 arrays pass. F core-essential −0.68, −0.54, −0.56 (bone marrow, intestine,
  liver; was −0.60, −0.53, −0.61), housekeeping −0.65, −0.54, −0.57; B +0.065 to +0.093.
- No claim changes. Manuscript Results and the supplement's data table name the rebuild.

## R97 — Huang mouse tissues on R SeSAMe betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures mouse --scale paper --seed 42`, job 6514013 (test partition), 2 h 05 min, on the
  chip-aware code of R92 (before 4711d89, so its within-chip dispersion intervals are degenerate and not reported).
- **The deposited Huang betas are R SeSAMe's own output.** On the 580 arrays both hold, all 279,644 shared probes agree
  to float32 rounding (max |difference| 3e-8) with identical missing cells. The rebuild adds 16,426 probes and drops 3
  arrays that fail the 90 % detection cut (young controls: two white-blood-cell, one hippocampus).
- Relative-entropy contrasts move by at most 0.02 (entropy slopes up to 0.04): housekeeping relative entropy, renewing − post-mitotic −0.191 (−0.508 to −0.066)
  unadjusted and +0.023 (−0.057 to +0.181) within chip; slow − post-mitotic −0.095 and +0.031. Per tissue, white
  blood cells −0.072 unadjusted, +0.028 within chip; kidney −0.162 (was −0.176). R92's reading holds on raw data.
- Claims unchanged beyond R92; E32 and Supplementary Note 3 name R97. Huang counts as rebuilt for RM41.

## R92 — mouse turnover contrast within chip (selection audit), 2026-11-18

- `figures mouse --scale paper --seed 42` on deposited betas, after the chip-aware code (3518b5f) and the bootstrap
  key fix; job 6512998 (test partition), 2 h 02 min. Adds `aging, within chip` (controls on chips shared by both
  ages, chip intercepts) and `aging, same tissues` (those tissues, no chip terms) beside the unadjusted `aging`.
- Unadjusted contrasts reproduce the published ones after the relative-entropy fix: housekeeping relative entropy,
  renewing − post-mitotic −0.186 (−0.451 to −0.066, p 0.005), slow − post-mitotic −0.099 (−0.195 to −0.019).
- **Within chip the ordering is gone:** renewing − post-mitotic +0.023 (−0.058 to +0.187, p 0.54), slow −
  post-mitotic +0.031 (−0.038 to +0.177, p 0.38); core-essential +0.007 and +0.010. Per tissue, white blood cells
  −0.075 → +0.039, liver −0.049 → +0.008, post-mitotic tissues −0.14 to +0.09. Colon, spleen, kidney and adipose
  drop out (fewer than 2 aged and 2 young on a shared chip).
- Same tissues without chip: renewing −0.130 (−0.935 to +0.075), slow −0.107 (−0.313 to +0.028). Losing tissues
  widens the interval; holding chip moves the estimate to zero. The slow − post-mitotic entropy-slope contrast
  (not relative entropy) holds on the same tissues (−0.190, p 0.005) but not within chip (−0.092, p 0.08).
- Defect: within-chip dispersion-slope intervals are degenerate (±1e12), cause not traced; not reported (RM41).
- Claims C2 and C3 and E32 changed: the mouse turnover ordering is not separable from chip. Manuscript Results,
  Discussion (i) and Supplementary Note 3 updated. R97 reruns the same on R SeSAMe betas.

## R87 — directed against noise pace of tissues across mammals (RM39 item 3), 2026-11-18

- Commit ab708cc; `figures species_pace --scale paper --seed 42`; job 6511070 (itmz), 41 min. 139 species-tissue
  groups of 65 species (`species_shield`'s), look-alike CpGs; directed pace = median |age slope of the mean|,
  noise pace = median age slope of dispersion, per year and per lifespan.
- Result, as measured: within the 9 species with >= 3 tissues, Spearman of directed against noise pace, median
  +0.18 (mouse +0.18 over 17 tissues, rat +0.48 over 13, human +0.07 over 13, naked mole rat −0.71 over 8), positive in
  7 of 9, Wilcoxon p 0.20; per lifespan the same. Tissue pairs seen in >= 5 species: 8; blood against skin ordered
  alike in 50 % (directed) and 42 % (noise) of 12 species; no pair ordered consistently.
- Reading: tissues that change fastest in a directed way are not clearly the noisiest, and no tissue order is
  shared across species, as R40 found in rodent RNA for noise. The data hold few species with several tissues.
  Table `data/output/species_pace/`.

## R84b, R85, R86, R88, R90 — the matched cancer design and four of the eight open ideas (RM38 item 4, RM39), 2026-11-18

- Jobs on the test partition (account cluster) and itmz, twins cancelling each other: R84b 6511660 (12 min, commit
  17f2439), R85 6511066 (28 min, b16f126), R86 6511661 (16 min, 4b2e7a4), R88 6511662 (10 min, 69d6d25), R90
  6511663 (11 min, 1d818af); all `--scale paper --seed 42`.
- **R84b, breast cancer as a matched case-control (EPIC-Italy, women):** 233 cases matched to 340 cancer-free
  controls (age within 3 years, same chip first; controls dealt in rounds), conditional logistic regression, 1,000
  set resamples. Drift at housekeeping look-alikes B: all cases OR 1.46 (1.19-1.90), with mean beta and 11 cell
  fractions 1.48 (1.16-2.14); diagnosed over 2 years after the draw (167 sets) 1.71 (1.35-2.48), adjusted 1.76
  (1.37-2.94); core-essential look-alikes 1.80 and 1.82. Shield share F 1.02-1.07, intervals across 1. Clocks:
  DunedinPACE 1.02-1.26, PhenoAge 1.11-1.27 (adjusted, over 2 years, 1.27, 1.03-1.69). The cohort-wide result (R46)
  holds in the design the cohort was sampled for; one cohort.
- **R85, sex across mammals (RM39 item 7):** 73 species-tissue groups of 38 species with both sexes. Female/male
  ratio of the look-alike drift rate, median log ratio +0.03 to +0.06 (relative entropy), share higher in females
  0.52-0.54, Wilcoxon p 0.77-0.90; shield gap female minus male −0.01 to −0.02, p 0.29-0.38 (core dispersion −0.035,
  p 0.07). No sex difference found across mammals, as in human tissues (R78).
- **R86, Braak order within Alzheimer's brains (RM39 item 2):** 109 Lunnon donors with three or four regions. The
  within-donor gradient of drift along entorhinal > temporal > frontal > cerebellum: mean −0.004 to −0.007 (SE 0.04),
  slope on Braak stage −0.008 to −0.017 per stage (p 0.43-0.68); with 5 components removed the same. Disorder does
  not follow the pathology's regional order.
- **R88, organs of one mouse at one age (RM39 item 1):** Huang controls, young and aged groups apart, 11 tissues, 51
  tissue pairs with >= 10 shared mice. Look-alike drift B is correlated across tissues of the same mouse, mean r
  0.15-0.19 by set (null mean 0.00, 95th centile 0.06, p 0.001); the shield share F, core-essential 0.19 (p
  0.001), housekeeping 0.12; signed shift S, little. Each mouse has a body-wide drift level at a fixed age, the
  clock-free counterpart of R41's human result. Chip is not held (R92 pending), and tissues of one mouse often share
  a chip, so part of this may be array batch.
- **R90, writer syndromes at noisy against directed aging CpGs (RM39 item 6):** noisy-aging and directed-aging CpGs
  from Hannum blood (top 5 %). Carriers gain extra variance at noisy-aging CpGs, beyond level- and island-matched
  CpGs, in NSD1 variants (+0.38 log variance ratio, permutation p 0.012), Weaver (EZH2, +0.58, p 0.016) and CHARGE
  (CHD7, +0.25, p 0.002); not in Sotos (+0.05, p 0.48), Kabuki or KMT2D variants (about 0). Directed-aging CpGs:
  no arm shifts them the way age does (p 0.40-0.89). Seven arms; CHARGE alone survives a Bonferroni cut, and the
  two NSD1 arms disagree. Where these variants act, they look like noise added at the CpGs aging makes noisy, not
  a faster directed clock.

## R96 — EPIC-Italy breast cancer on rebuilt arrays (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures cancer cancer_matched --scale paper --seed 42` (commits b36c116, 8c00678 fix:
  reference CpGs masked in every array left out of the cell-mix projection); job 6514018 (test partition), 24 min.
  GSE51032 rebuilt with R SeSAMe (6 of 845 arrays fail detection); mean beta and cell fractions from the rebuild too.
- Result, as measured (matched case-control among women, odds ratio per SD; deposited R84b → rebuilt). Drift B,
  housekeeping look-alikes, cases over 2 years after the draw: age only 1.71 → 1.84 (1.37-3.48, p 1e-4); with mean
  beta and 11 cell fractions 1.76 → 1.11 (0.69-2.18, p 0.66). All cases: 1.46 → 1.52 (age only), 1.48 → 1.15
  (adjusted). Shield share F adjusted: 1.02 → 0.69 (0.49-0.95, p 0.015). DunedinPACE 1.26 → 1.19, PhenoAge 1.27 →
  1.18, n.s. Sets 167 → 163.
- Reading: on raw-processed arrays the drift signal is there but cannot be separated from global methylation and
  cell mix; the processed-data result (adjusted 1.76) depended on the authors' normalisation. The F signal that
  appears on rebuilt data was not seen before and is not claimed. The lead is downgraded: unresolved, needs a
  second cohort and a decision on which adjustment fits the question.

## R91 — the CR blood series with a within-chip, refitting randomisation (selection audit), 2026-11-18

- Commit 7357864; `figures cr_blood_replication --scale paper --seed 42`; job 6511076 (itmz), 3 h 35 min (smoke
  6510231). Treatment reassigned within each array chip, each chip's treated count kept; the controls' reference,
  composition components and leave-one-out scores refitted per assignment; the observed assignment reproduces the
  reported effect (asserted). CR against controls 720 assignments (exact), TLS+CR 1,000 of 4,620 drawn, long-term
  TLS 84 (exact, one shared chip); controls at 12 weeks against Time 0 share no chip: no p.
- Result, as measured (3 components; effect, old fixed-score p → corrected p). CR: core-essential F −0.090, 0.0065
  → 0.161; housekeeping F −0.027, 0.026 → 0.235; B +0.002, 0.17-0.16; S +0.002-0.003, 0.30-0.33. TLS+CR: F −0.043
  and −0.015, p 0.53 and 0.67. Long-term TLS: all p > 0.12. Chip-intercept effects (descriptive): CR core F −0.075,
  housekeeping −0.028.
- Reading: no effect in this series survives the corrected test. R56's "CR's changes spare essential promoters"
  is withdrawn as a finding and kept only as a point estimate; the 12-week aging contrast stays unidentifiable.
  Claims, manuscript, supplement and extension.md corrected.

## R95 — writer syndromes and their concordance on rebuilt betas (RM41), 2026-11-18

- `ENTROPY_AGE_BETAS=sesame figures syndromes syndrome_concordance --scale paper --seed 42` (commit b36c116); job
  6513665 (test partition), 29 min. GSE74432 and GSE97362 from R SeSAMe (HM450; 2 and 1 arrays fail detection).
  Outputs `data/output/syndromes_rsesame/`, `syndrome_concordance_rsesame/` beside the deposited ones.
- Result, as measured (housekeeping, deposited → rebuilt). Look-alike drift B: Sotos +0.018 → +0.019, NSD1 variants
  +0.013 → +0.014, Weaver +0.012 → +0.014, Kabuki +0.009 → +0.009, CHARGE +0.007 → +0.007 (all p < 0.001 both);
  KMT2D and CHD7 variants about 0 both. Shield share F: Sotos −0.051 → −0.069, NSD1 variants −0.067 → −0.079,
  Weaver −0.063 → −0.078 (all p ≤ 0.02 both); the others within ±0.02, n.s. Noise concordance at noisy-aging
  CpGs: NSD1 variants +0.38 (p 0.012) → +0.42 (p 0.004), Weaver +0.58 (0.016) → +0.82 (0.012), CHARGE +0.25 (0.002)
  → +0.14 (0.064), Sotos, Kabuki, KMT2D and CHD7 variants about 0 both; directed concordance none either way.
- Reading: the syndrome results hold on rebuilt arrays, the NSD1 and EZH2 shield loss a little stronger. In R90 the
  arm that survived Bonferroni changes: on rebuilt betas it is NSD1 variants (p 0.004), and CHARGE no longer
  reaches 0.05. The concordance claim rests on NSD1 and EZH2, not CHARGE.

## R94 — the culture test on rebuilt betas: deposited against R SeSAMe on the same samples (RM41), 2026-11-18

- Commit 53b4db5 (after a filter that rebuilt a set per feature ran for hours, 6512015 cancelled); `figures
  passage_pipelines --scale paper --seed 42`; job 6513615 (test partition). Endicott GSE197512 and Sturm GSE179847
  each scored twice on exactly the samples R SeSAMe passes (Endicott 333 of 372: 33 dropped in the TERT arm, 4
  baseline and 2 reduced-oxygen samples of one foreskin line; Sturm 462 of 479), on the CpGs complete in both
  sources (no fill).
- Result, as measured (coupling = gap per unit of noise-free matched drift; its change per 10 doublings, mean over
  lines with t interval). Endicott, 7 lines, housekeeping: coupling −0.202 deposited, −0.003 rebuilt; change −0.113
  (−0.257, +0.032) deposited, −0.176 (−0.392, +0.040, p 0.09) rebuilt. Core-essential: coupling −0.229 → +0.133;
  change −0.086 → −0.191 (−0.395, +0.013, p 0.06). Sturm, 4 lines (HC6 falls below 5 samples): housekeeping coupling
  −0.492 → −0.368, change −0.009 (−0.085, +0.066) → −0.025 (−0.082, +0.033); core −0.526 → −0.381, change +0.003 →
  −0.005. Original estimates on all samples for reference: Endicott −0.139 (−0.334, +0.056), Sturm −0.055.
- Reading: processing, not sample selection, moves Endicott. On rebuilt betas the shield in Endicott's cultures is
  about zero at the reference (housekeeping) or reversed (core), and its change with doublings is more negative,
  the direction selection predicts, with intervals still including zero. Sturm keeps a shield on both versions and
  no change with doublings. The culture result no longer argues against selection; it is unresolved and
  processing-sensitive. Tables `data/output/passaged_cells_{deposited_subset,sesame}`,
  `passaged_cells_sturm_{deposited_subset,sesame}`.

## R93 — Jaffe cortex, paired pipeline pilot: deposited, methylprep and R SeSAMe on identical donors (RM41), 2026-11-18

- Raw step: `scripts/sesame_preprocess.R` (R 4.4.3, sesame 1.24.0, sesameData 1.24.0, ExperimentHub 2.14.0,
  preprocessCore 1.68.0 rebuilt with --disable-threading because the conda build's threads failed in dyeBiasNL;
  `prepSesame("QCDPB")`, auto platform HM450) on the IDATs of the 213 adult controls of `figures tissues`; job
  6511393 (test partition, 4 cores), 11 min. 5 donors fail detection < 90 % (0.82-0.88); bisulfite index GCT
  1.04-1.14. Outputs `data/interim/sesame/GSE74193_adult/` (betas, QC, session).
- Pilot: commit 80f95ad, `figures reprocessed_jaffe --scale paper --seed 42`; job 6511655, 34 min. 208 donors passing
  in all three sources, 89,446 CpGs measured without a mask in every donor under every source (no fill); the same
  100 donor resamples for each source; neuron fraction, plate, sex and race held as before.
- Result, as measured (matched gap in SD with 95 % interval; deposited / methylprep / R SeSAMe): housekeeping relative
  entropy −0.064 (−0.091, −0.034) / −0.086 (−0.136, −0.050) / −0.081 (−0.125, −0.046); core-essential −0.057 / −0.077 /
  −0.074; housekeeping spread −0.072 (−0.111, −0.028) / −0.124 / −0.119 (−0.153, −0.063); core-essential spread
  −0.053 / −0.100 / −0.099; non-essential +0.006 / +0.033 / +0.030 (relative entropy), intervals across 0. Paired
  differences: R SeSAMe minus methylprep at most 0.004 SD in every set and metric (intervals within ±0.01);
  R SeSAMe minus deposited, relative entropy −0.017 (−0.067, +0.012) housekeeping and core, spread −0.047
  (−0.081, −0.003) and −0.046 (−0.087, −0.001).
- Reading: the two raw routes agree, so methylprep's earlier rebuild (R61) stands as a SeSAMe-equivalent here.
  Rebuilding from IDATs does not weaken cortex's protection; it strengthens it on spread (by about two thirds) and
  leaves relative entropy within its interval. The deposited betas understate the cortex gap, not inflate it. The
  clock-shield association is not rerun here: it waits for the clock-scorer audit (RM41 handoff parts 3-4).

## R83 — per-CpG gain and loss with donor intervals (RM38 item 2), 2026-11-18

- Commit a276d20 (`figures per_cpg_interval --scale paper --seed 42`), ENTROPY_AGE_BOOT=100; 4 cores, 64 GB; job
  6506416, 5 h 32 min, 2.9 GB peak (smoke 6506415 passed). Hannum 656 donors, Jaffe 213; R67's per-CpG bounded
  relaxation fits; each resample's gap on the full sample's SD.
- Result, as measured (set minus look-alikes, SD, 95 % donor interval). Loss rate `k_loss`, island x level: blood
  housekeeping +0.114 (+0.051, +0.147), core-essential +0.073 (+0.020, +0.099), non-essential −0.079 (−0.090,
  −0.036); cortex housekeeping −0.020 (−0.034, +0.038), core −0.017 (−0.049, +0.045). With the CpG's own chromatin
  state also matched: blood housekeeping +0.020 (−0.003, +0.041), core +0.012 (−0.025, +0.037); cortex −0.026,
  −0.020, intervals across 0. Gain rate `k_gain`: within ±0.02 everywhere, intervals across 0. Resting level
  (equilibrium): lower at essential promoters in both tissues, island x level blood −0.34/−0.24, cortex −0.15/−0.08,
  with chromatin state −0.06 to −0.07 and −0.03 to −0.07, all excluding 0.
- Reading: per CpG, faster loss at essential promoters is a blood result, and it goes with the promoters' own
  chromatin state; in cortex no faster loss is seen per CpG. What holds in both tissues, beyond chromatin state, is
  a lower resting level. R67/R74's pooled rates (faster loss also in cortex, +0.006-0.007 per decade within active
  TSS) rest on the pooled fit and do not survive per CpG. Some point estimates sit near an interval's edge (blood
  equilibrium): the resampled fits identify a different set of CpGs each time, so these intervals are approximate.
  Table `data/output/per_cpg_interval/per_cpg_gaps.csv`.

## Audit — selection hypothesis, outside review of 2026-11-18 (R91, R92 pending)

- The review (another model, read-only, on 2261801) compared deposited GEO metadata with the saved tables. Its four
  scripts were rerun on this checkout (1d818af and later): every number reproduces to the last digits. Nothing in
  the 21 commits since changed the culture, CR, transplant, Huang or Tet3 tables.
- **CR blood, GSE282661 (R56).** All 45 labels agree. All 8 Time 0 animals sat on one chip, the 8 endpoint
  controls on two others: the 12-week aging contrast is confounded with chip and is now reported as unresolved
  (was "no drift to shield"). The permutation p shuffled labels on scores fitted against the observed controls;
  `cr_blood_replication` now reassigns treatment within chips (each chip's treated count kept), refits the
  reference, the components and the leave-one-out scores per assignment, and asserts that the observed
  assignment reproduces the reported effect (7357864). CR against controls has 720 assignments, enumerated; a
  comparison whose arms share no chip gets no p (`within_chip_overlap` 0). Smoke R91 (job 6510231, test partition);
  the paper-scale rerun is pending, so R56's p values (core-essential F p 0.006, housekeeping 0.026) stand
  uncorrected until it lands. Chip-adjusted housekeeping F (3 components), from the saved scores: −0.0282 against
  −0.0274 unadjusted (descriptive). The same fixed-score shuffle in `reprogramming_rrbs` is valid for its OSKM
  contrast and approximate for its age contrast; flagged in its docstring, not refitted.
- **Mouse turnover, Huang GSE279437 (R68 mouse).** Aged and young controls share no chip in spleen and adipose, one
  in colon and kidney (age explained by chip: R2 1.0, 1.0, 0.94, 0.87). The published class contrast (renewing −
  post-mitotic, housekeeping relative entropy −0.190, −0.471 to −0.065) is unadjusted. `mouse_tissues` now writes
  `chip_overlap` and adds `aging, within chip` (controls on shared chips, chip intercepts; identifiable only where
  >= 2 aged and 2 young share chips: white blood cells among the renewing tissues; colon, kidney, spleen and
  adipose drop out) and `aging, same tissues` (unadjusted, on those tissues), so chip adjustment and the change of
  target can be told apart (3518b5f). The bootstrap already resampled animals once for all tissues. The mouse
  outputs predate the relative-entropy fix (77aa712); R92 (job 6510291, requeue) reruns both. Claims corrected
  meanwhile: "five renewing tissues" was three.
- **Culture (Endicott, Sturm).** Inclusions match GEO exactly (174, 46). Housekeeping coupling change per 10
  doublings −0.139 (−0.334 to +0.056, 7 lines), −0.142 with batch intercepts; Sturm −0.055 (−0.213 to +0.102, 5
  lines); on the 4 lines a plate model can fit, −0.0025 (−0.089 to +0.084) unadjusted and +0.0026 (−0.065 to
  +0.070) adjusted, so plate moves nothing and the shift from −0.055 is HC6's removal; all 5 Sturm lines used a
  pooled noise floor. "Constant share" and "rejects selection" replaced by the slopes and intervals: zero is
  inside, but so are selection-sized changes.
- **Transplant (R57) and Tet3 (R15, R16).** Core-essential minus non-essential bottleneck effect −0.0013 (pair
  bootstrap −0.035 to +0.033); on common pairs and overlapping times +0.013: little evidence of a
  core-essential-specific effect, as stated. Tet3: 8 heterozygotes and 8 wild types per tissue, all 0.5 years,
  genotypes balanced on every chip; bulk cortex and striatum; GEO records no litters. "Littermates" became
  "same-chip, same-age wild types" and "in neurons" became "in bulk brain tissue".
- Text changed: `docs/CLAIMS.md` (C2, C3 renamed "not supported as the main cause", C4), `docs/evidence_matrix.csv`
  (E5 and E11 weakened, E6 scoped to GTEx, E32 added), `manuscript/paper.md`, `manuscript/supplement.md`,
  `docs/extension.md`, `docs/PAPER_DRAFT.md`, `docs/NOVELTY_REVIEW.md`, `docs/ROADMAP.md`. The Drive round v0.1
  Docs are not touched; these go into v0.2 (RM37).

## R89 — drift per division or per day, across every culture arm (RM39 item 5), 2026-11-18

- Commit (this one); `figures culture_rates`, login node, seconds: the saved per-sample scores of the RM5 runs
  (Endicott GSE197512, 354 samples, 13 experiments: passaging, mitomycin arrest, serum withdrawal, 3 % oxygen,
  TERT; Sturm GSE179847, 440 samples, 19 experiments: dexamethasone, galactose, 2-deoxyglucose,
  beta-hydroxybutyrate, SURF1). Look-alike drift ``B`` regressed jointly on doublings and days, one intercept per
  experiment. Reused, not rerun (Anton's rule): the scores are R1's and the RM5 runs'.
- Result, as measured (per 10 doublings / per 10 days, SE): Sturm (doublings and days r 0.73) housekeeping-matched
  B +0.0154 (0.0007) / +0.0012 (0.0002); core-essential-matched +0.0112 / +0.0009; non-essential-matched +0.0319 /
  +0.0013. Endicott (r 0.97; the arrest arms carry the split) +0.0148 (0.0022) / +0.0052 (0.0011); core +0.0103 /
  +0.0035; non-essential-matched 0.000 (0.004) / +0.0127 (0.0022).
- Reading: in cells, drift at unmethylated promoter strata rises with divisions at the same rate in both series
  (0.0148 and 0.0154 per 10 doublings). Time adds its own part: at about 3 days per doubling, divisions give
  about four times what time gives in Sturm, about as much in Endicott, whose arrested arms carry the split. Methylated
  strata (matched to non-essential promoters) drift per day in Endicott and per division in Sturm; not settled.
  The shield's relative gap is left out: its fall with doublings is the array-noise artefact found in RM5.
  Table `data/output/culture_rates/rates.csv`.

## R82 — placebo gene sets matched on expression (RM38 item 3), 2026-11-18

- Commit 6c22e54; `figures placebo_expression --scale paper --seed 42`; 4 cores, 64 GB; job 6506412, 32 min. R69's
  design (Hannum and Jaffe, every autosomal promoter CpG, the shield's matching), but each placebo set drawn decile
  by decile of the tissue's GTEx v8 median expression (Whole Blood; Brain - Frontal Cortex BA9) to the real set's
  counts; 200 draws per set and tissue. Balance: no decile ran short; mean log2(TPM + 1) real against placebo,
  blood housekeeping 3.59 / 3.66, core 3.85 / 3.79; cortex 4.54 / 4.59, 4.64 / 4.63.
- Result, as measured (real gap vs expression-matched placebo mean [95 %], SD): blood housekeeping relative
  entropy −0.138 vs −0.111 [−0.120, −0.103], dispersion −0.273 vs −0.212; blood core-essential relative entropy
  −0.104 vs −0.093 [−0.118, −0.069] (18 % of placebos below it), dispersion −0.214 vs −0.177. Cortex housekeeping
  relative entropy −0.053 vs −0.014 [−0.022, −0.006], dispersion −0.067 vs −0.015; core-essential −0.074 vs −0.011
  [−0.035, +0.011], dispersion −0.086 vs −0.011. Non-essential sets sit below their low-expression placebos
  (blood +0.118 vs +0.187).
- Reading: in blood, expression alone gives about four fifths of the housekeeping gap and all of the core-essential
  gap; what essentiality adds beyond expression is about −0.03 SD (housekeeping) and nothing (core). In cortex the
  gap is four to seven times what expression gives (−0.04 to −0.06 SD beyond it). With R80 (blood's gap is mostly
  the 0 % floor, cortex keeps its size) the robust shield is cortex's; blood's is mostly expressed, unmethylated
  promoters. R69's random placebos (gap 0) did not separate these. Table `data/output/placebo_expression/`.

## R80 — the shield away from the 0 % floor (RM38 item 1), 2026-11-18

- Commit 963a9ff (audit fixes in); `figures floor_check --scale paper --seed 42`, ENTROPY_AGE_BOOT=100; 4 cores,
  64 GB; job 6503584, 38 min (first submission 6502924 cancelled before the audit fixes). Hannum blood and Jaffe
  cortex, `tissue_contrast` features, covariates and matching (island x young-level decile), CpGs kept when their
  young level is in the band; 95 % donor bootstrap on the full sample's SD.
- Result, as measured (matched gap, SD; all levels → 10-90 % → 20-80 %):
  blood housekeeping, dispersion −0.29 (−0.39, −0.18) → −0.19 (−0.25, −0.12) → −0.08 (−0.12, −0.02); relative
  entropy −0.14 → −0.08 (−0.15, −0.03) → 0.00 (−0.05, 0.04). Blood core-essential, dispersion −0.17 → −0.07 (−0.15,
  0.05) → +0.03; relative entropy −0.08 → −0.03 → +0.01. Blood non-essential stays positive: dispersion +0.27 →
  +0.17 → +0.21. Cortex housekeeping, dispersion −0.07 → −0.09 (−0.17, 0.01) → −0.08 (−0.19, 0.03); core-essential
  −0.08 → −0.10 (−0.19, −0.01) → −0.06 (−0.19, 0.08); relative entropy alike. CpGs left at 20-80 %: 10 % of
  housekeeping, 8 % of core-essential (blood); 17 % and 13 % (cortex).
- Reading: in blood most of the shield is the floor. Away from it the housekeeping gap on spread shrinks to a third
  and the relative-entropy gap and the core-essential gap vanish; what survives in blood is the contrast with
  non-essential promoters, which drift more at every level. In cortex the gap keeps its size off the floor but its
  interval reaches 0 at 20-80 %. The species result (R59) held off the floor; human blood's does only in part.
  The paper's blood shield must be stated with this bound. Table `data/output/floor_check/gaps_by_band.csv`.

## R81 — outside audit: old against fixed relative entropy on the headline gaps, 2026-11-18

- An outside audit (2026-11-18, another model, read-only) found that `relative_entropy` clipped only the reference
  to [0.01, 0.99], so near 0 % a sample moving away could score closer; fixed in 77aa712 (one pseudo-count on both
  levels), with the bootstrap scale of floor_check, allelic_controls and sex_drift, the global FDR columns, NaN
  strata, rescale's duplicate check and the hotspot null (tests/test_audit_fixes.py).
- Commit 77aa712; `/data/uq2347/logs/kl_impact.py` (the old function inline against the fixed one); 4 cores, 64 GB;
  job 6503280, 1 min 25 s. Hannum and Jaffe, `tissue_contrast` features, covariates and matching, 200 permutations.
- Result, as measured (matched gap on relative entropy, SD, old → fixed): blood housekeeping −0.146 → −0.141,
  core-essential −0.084 → −0.082, non-essential +0.135 → +0.133; cortex housekeeping −0.056 → −0.056, core-essential
  −0.070 → −0.070, non-essential −0.000 → −0.001; every p unchanged. 15.6 % of blood housekeeping CpGs start below
  1 % (7.3 % of the rest), none in cortex.
- Reading: the fix moves the headline gaps by at most 0.005 SD; no claim changes. The global FDR rerun (963a9ff):
  8,532 tests, 4,207 survive, largest surviving p 0.024 as before. Reruns of R77-R80 on the fixed code follow;
  earlier relative-entropy tables stay as they are, with this measured shift as their correction bound.

## R79 — hotspots of rising disorder across mammals (RM35 item 3), 2026-11-17

- Commit eb834c1; `figures species_hotspots --scale paper --seed 42`; 4 cores, 96 GB; job 6495752, 1 h 36 min.
  139 species-tissue groups of GSE223748 (the `species_shield` groups), 27,889 autosomal mapped CpGs; top 5 % per
  group on relative-entropy slope and on dispersion slope; null: each group's labels shuffled within its young-level
  deciles, 200 draws.
- Result, as measured: hotspots are barely shared. CpGs hot in ≥ 10 groups: 6,926 against 6,585 expected (1.05×,
  relative entropy), 6,388 against 5,190 (1.23×, dispersion); hot in half the groups: 3 CpGs (relative entropy), 0
  (dispersion). Mean pairwise rank correlation of the slopes between groups 0.06 (relative entropy) and 0.03
  (dispersion); higher in the same tissue across orders (0.12) than in the same order across tissues (0.07). The
  CpGs hot in most groups are CpG islands of developmental Polycomb targets: Tlx3 (two CpGs, 66-72 groups), Hoxa11
  (71), Lhfpl4 (74, 65), Celf6 (67), Runx1 (66), Prdm13 (65), all low-methylated when young.
- Reading: across mammals the place of age drift is mostly specific to a species and tissue; tissue shares more
  than phylogeny. The small shared tail is the Polycomb developmental islands that universal clocks are built on
  (Lu 2023), the programme R72 found in human blood. No universal noise hotspot outside it. Tables
  `data/output/species_hotspots/`.
- Rerun on the audit fixes (963a9ff; null shuffled among measured CpGs, fixed relative entropy): job 6503583,
  1 h 41 min. Hot in ≥ 10 groups: 6,965 against 6,716 (1.04×, relative entropy), 6,388 against 5,190 (1.23×,
  dispersion); hot in half the groups: 2 CpGs (relative entropy), 0 (dispersion). Slope correlation 0.058 and 0.028;
  same tissue 0.12, same order 0.07. The reading stands.

## R77 — loci held at one half by allele as a noise floor (RM35 item 1), 2026-11-16

- Commit dc4a72a; `figures allelic_controls --scale paper --seed 42`, ENTROPY_AGE_BOOT=50; 8 cores, 64 GB; job
  6494447, 46 min (smoke 6494294 failed on a membership column mismatch, fixed e368864; 6494332 passed). Imprinted
  genes: geneimprint "Imprinted", 131 human genes (read 2026-11-16), island or shore CpGs at young level 0.25-0.75.
- Result, as measured (matched gap on relative entropy, SD, 95 % donor bootstrap; island relation x 20 level bins
  against 60,000 random autosomal CpGs): imprinted half-methylated CpGs in Hannum blood −0.51 (−0.56, −0.31),
  Lunnon blood men −0.20 (−0.26, −0.01), women −0.08 (−0.19, 0.11); cortex 0 (Jaffe −0.03, Lunnon −0.02/−0.02, all
  intervals across 0). X-island CpGs at one half in women: Lunnon blood −0.29 (−0.35, −0.09); Jaffe cortex −0.03
  (−0.16, 0.07), Lunnon frontal cortex +0.10 (−0.26, 0.22). Housekeeping in the same runs: Hannum −0.16, Jaffe −0.05.
- Rerun on the audit fixes (77aa712; resamples on the full sample's SD, fixed relative entropy): job 6503581,
  50 min. Imprinted, Hannum blood −0.51 (−0.71, −0.39); Jaffe cortex −0.03 (−0.17, 0.09); X at one half in women,
  Lunnon blood −0.29 (−0.55, −0.13), Lunnon frontal cortex +0.10 (−0.43, 0.32); housekeeping Hannum −0.15 (−0.22,
  −0.11), Jaffe −0.05 (−0.09, −0.03). Intervals now contain their estimates; the reading stands.
- Reading: loci a maintenance machinery holds at one half are protected from drift in blood, two to three times
  the essential-promoter shield, and not in cortex. The blood-brain difference of the shield is a property of the
  tissue's maintenance of held states, not of essential genes alone. Tables `data/output/allelic_controls/`.

## R78 — drift and the shield by sex (RM35 item 2), 2026-11-16

- Commit 2d9f86a; `figures sex_drift --scale paper --seed 42`, 100 resamples; 4 cores, 64 GB; job 6494333, 62 min.
  Hannum, Jaffe, Lunnon's five tissues, GTEx colon and lung (other GTEx tissues have < 20 donors of a sex).
- Result, as measured (women − men, 95 % donor bootstrap): the background drift rate (median slope of spread and of
  relative entropy) does not differ in any of nine tissues (p_boot 0.14-1.0). The housekeeping shield on relative
  entropy does not differ either (−0.09 to +0.07 SD, every interval across 0); on dispersion it is weaker in women in
  Jaffe cortex (+0.08, 0.01-0.10) and GTEx colon (+0.14, 0.02-0.24), nowhere on relative entropy.
- Rerun on the audit fixes (77aa712): job 6503582, 47 min. Relative entropy, women − men, housekeeping: −0.09 to
  +0.06 SD, every interval across 0; dispersion weaker in women in Jaffe cortex (+0.08, 0.01-0.15) and GTEx colon
  (+0.14, 0.01-0.30). The reading stands.
- Reading: no sex difference in drift or the shield; sex need not be split in the shield's analyses. Tables
  `data/output/sex_drift/`.

## R73 — metric audit leftovers: mouse interventions, BXD, RRBS, Meer (RM33), 2026-11-14

- Commit 40fd45d; `figures interventions bxd rrbs_tissues meer_tissues --scale paper --seed 42`; 4 cores, 64 GB;
  job 6480746, 2 h 23 min. Old measures reproduce.
- Result, as measured (housekeeping aging gap, SD, entropy → relative entropy): GSE80672 −0.350 → −0.304; BXD
  −0.524 → −0.320; RRBS blood −0.160 → −0.077, lung −0.268 → −0.210, kidney −0.327 → −0.184, liver −0.257 → −0.131,
  adipose −0.195 → −0.028, muscle −0.160 → −0.080; Meer cortex −0.092 → −0.064, heart −0.100 → −0.071, liver −0.356
  → −0.220, lung −0.290 → −0.138.

## R76 — input checksums (RM34 item 5), 2026-11-16

- `sha256sum` over 753 archive files under `data/input/` (192 GB); 1 core; job 6493412, 10 min. Table
  `docs/INPUT_CHECKSUMS.tsv`.

## R75 — RNA turnover test in Tabula Muris Senis bulk (RM4), 2026-11-16

- Commit aee749c; `figures rna_tms --scale paper --seed 42`; 8 cores, 64 GB; job 6493369, 11 min. GSE132040, 831
  samples of 17 organs aged 3-27 months; subsample 31 animals.
- Result, as measured (housekeeping, spread slope, covariates): organ intervals include 0 in 15 of 17 (small
  intestine −0.52 to −0.16, pancreas +0.08 to +0.46); class contrasts renewing − brain +0.20 (p 0.21), slow −
  brain +0.27 (p 0.02), others p ≥ 0.25. Tables `data/output/rna_tms/`.

## R66 — SATSA within-person, rebuilt from IDATs (RM28), 2026-11-14

- R66a: `scripts/reprocess_idats.py` on E-MTAB-7309, 1,072 450k arrays in batches of 100; 4 cores, 64 GB; job
  6479360, 12 h 22 min, killed out of memory while combining (methylprep's `beta_values.pkl` was complete). R66c
  (`--from-pickle`, commit after 43e9f2c; 128 GB, job 6481946, 4 min): 485,512 CpGs, 1,072 samples; then `figures
  satsa_idats --scale paper --seed 42`: 444 people, 58,340 CpGs, 628 pairs (873 with the deposited betas).
- Result, as measured (housekeeping, set − look-alikes, SD; deposited R52 → IDATs): within-person divergence
  −0.017 → −0.014 (p 0.016); directed change −0.002 (p 0.67) → −0.140 (p 0.001); floor −0.029 → −0.073;
  cross-sectional dispersion −0.12 → −0.36, |mean| +0.00 → −0.11, relative entropy −0.23 (IDATs only).

## R74 — gain and loss rates with chromatin state held (RM32), 2026-11-14

- Commit 599c2fa; `figures equilibrium_drift --scale paper --seed 42`; 4 cores, 32 GB; job 6480933, 27 min.
- Result, as measured (island promoters, set − look-alikes, loss per decade, 95 %): active TSS only, blood core
  +0.0063 (+0.0031, +0.0097), housekeeping +0.0015 (+0.0001, +0.0034); cortex core +0.0070 (+0.0041, +0.0097),
  housekeeping +0.0057 (+0.0034, +0.0083). Per-CpG with own state in strata: resting level blood −0.335 → −0.060
  (housekeeping), −0.243 → −0.062 (core); cortex −0.145 → −0.071, −0.077 → −0.030; per-CpG loss blood +0.114 →
  +0.020. Tables `data/output/equilibrium_drift/`.

## R72 — drift hotspots by region class (RM34 item 4), 2026-11-14

- Commit a3805b7; `figures drift_hotspots --scale paper --seed 42`; 4 cores, 96 GB (21 GB peak); job 6480722, 1 min
  16 s. All autosomal CpGs: Hannum 473,034, Jaffe 473,864; top 5 % by dispersion slope and by relative entropy.
- Result, as measured (observed / expected within young-level deciles, dispersion): blood 13_ReprPC 2.11 (z 72.8),
  10_TssBiv 1.99, 12_EnhBiv 1.50, 6_EnhG 0.20, 4_Tx 0.36, housekeeping 0.42, core-essential 0.46, non-essential
  1.62; cortex range 0.71 (core-essential) to 1.32 (ZNF/repeats), open sea 1.18, solo-WCGW 1.25. Table
  `data/output/drift_hotspots/enrichment`.

## R71 — one FDR over every reported p (RM34 item 3), 2026-11-14

- `figures global_fdr` run on the login node (reads tables only, seconds), after R70. 7,626 tests in 53 output
  folders; 3,932 at p < 0.05, 3,662 at global q < 0.05; largest surviving p 0.024. Lost in claim-bearing folders:
  CR housekeeping F/P in the second blood series, CR C57BL/6 core F, several MMC knockout contrasts, any-cancer and
  colorectal drift, GTEx smoking and death type, species all-CpG lifespan slopes. Tables `data/output/global_fdr/`.

## R69 — placebo gene sets (RM34 item 1), 2026-11-14

- Commit 9321cb6; `figures placebo_sets --scale paper --seed 42`; 4 cores, 64 GB; job 6480607, 34 min. Hannum 656
  donors, 192,719 promoter CpGs; Jaffe 213, 192,998; 200 placebo draws per set and tissue.
- Result, as measured (real gap vs placebo mean [95 %], SD): blood housekeeping relative entropy −0.143 vs +0.019
  [−0.001, +0.036], core −0.108 vs +0.014 [−0.026, +0.064]; cortex housekeeping −0.053 vs +0.013 [−0.001, +0.028],
  core −0.074 vs +0.014 [−0.022, +0.050]; no placebo below any essential gap on any metric. Tables
  `data/output/placebo_sets/`.

## R70 — Hannum shield with 12 cell fractions held (RM34 item 2), 2026-11-14

- Commit 176a250; `figures composition_gaps --scale paper --seed 42`; 4 cores, 32 GB; job 6480608, 8 min. 656
  donors, 108,992 CpGs; fractions by Houseman projection on EpiDISH cent12CT450k.
- Result, as measured (housekeeping gap, SD, without → with fractions): relative entropy −0.152 → −0.152; |mean|
  −0.159 → −0.177; entropy −0.348 → −0.359; dispersion −0.289 → −0.199; core-essential relative entropy −0.087 →
  −0.091; all p 0.001. Tables `data/output/composition_gaps/`.

## R68 — metric audit: relative entropy beside binary entropy (RM33), 2026-11-14

- Commit e3ad29c (`kl_slope`, `K`); `--scale paper --seed 42`. R68a `tissues`, `mouse` (8 cores, 96 GB; job
  6479900, 2 h 25 min); R68b `chromatin`, `essentiality_gradient`, `sorted_cells`, `reprocessed_jaffe`,
  `probe_robustness`, `gill_reprogramming` (4 cores; job 6479901, 2 h 08 min); R68c `species`, `species_floor`
  (job 6479902, 2 h 14 min); R68d `gtex` (job 6479903, 2 h 54 min). Old measures reproduce.
- Result, as measured (housekeeping gap, SD, entropy → relative entropy): blood −0.348 → −0.152, cortex −0.207 →
  −0.058; chromatin share cortex 0.43-0.69 → 0.25-0.56, blood 0.77-0.84 → 0.83-0.91; R58 blood remainder +0.022 →
  −0.007 (p 0.33), cortex +0.018 → +0.022; CD4+ T −0.219 → −0.015; mouse post-mitotic +0.03 to +0.09 (renewing
  −0.05 to −0.19); species 96 % → 89 % negative, R22 relative-gap lifespan slope −2.42 (p 0.04) → +1.43 (p 0.23).
  Verdicts in `docs/METRIC_AUDIT.md`.

## R67 — gain and loss rates per gene set (RM32), 2026-11-13

- Commit 9ffc595; `figures equilibrium_drift --scale paper --seed 42`; 4 cores, 32 GB; job 6479803, 16 min, 3.2 GB
  peak. Hannum 656 donors, 108,992 CpGs; Jaffe adult controls 213, 109,137 CpGs; 200 donor resamples, 1,000
  permutations. Figure redrawn locally from the saved table (x ticks thinned).
- Result, as measured (island promoters, set − look-alikes, per decade, 95 % interval): blood gain core −0.0012
  (−0.0017, −0.0008), housekeeping −0.0012; loss core +0.0025 (+0.0012, +0.0037), housekeeping +0.0018 (+0.0009,
  +0.0026); cortex gain −0.0002, loss core +0.0054 (+0.0035, +0.0072), housekeeping +0.0042 (+0.0027, +0.0058).
  Per CpG (matched SD): resting level blood −0.24 / −0.33, cortex −0.08 / −0.15; distance to it blood −0.28 /
  −0.42, cortex −0.03 / −0.11; loss rate blood +0.07 / +0.11, cortex −0.02 / −0.02; p ≤ 0.02 except cortex core
  loss (0.24). Fits identified: 76-84 % of CpGs. Tables `data/output/equilibrium_drift/`.

## R65 — Gill 2022 transient reprogramming from IDATs (RM17), 2026-11-13

- R65a: `scripts/reprocess_idats.py` on GSE165179_RAW.tar, 96 EPIC arrays in batches of 50; 4 cores; job 6479358,
  2 h 05 min; 865,859 CpGs.
- R65b: `figures gill_reprogramming --scale paper --seed 42` (commit bd5e0c9); 4 cores, 32 GB; job 6479825, 36 s.
  86,366 promoter CpGs; 13 reprogrammed, 21 failed, 33 controls over three donors; 5,000 within-donor permutations.
- Result, as measured (within-donor effect, permutation p): reprogrammed − failed, B housekeeping +0.045, core
  +0.036, non-essential +0.125 (0.0002 each); P housekeeping −0.006, core −0.004 (0.0002); F housekeeping −0.016
  (0.25), core −0.015 (0.35), non-essential −0.060 (0.014); E +0.005 / +0.003 (0.0002); Horvath skin & blood raw
  sum −0.36 (0.085). Failed − control: B housekeeping +0.010 (0.005), F −0.018 (0.14), non-essential F −0.11
  (0.0004). Tables `data/output/gill_reprogramming/`.

## R61 — Jaffe cortex reprocessed from IDATs (RM18 item 3), 2026-11-13

- R61a: `scripts/reprocess_idats.py` (methylprep 1.7.1, sesame=True; Python 3.11, pandas<2, numpy<2) on
  GSE74193_RAW.tar, 675 arrays in batches of 100; 4 cores; job 6469509, 7 h 17 min, 50 GB peak; 485,512 CpGs.
  A first run failed at once (pandas 2 removed `DataFrame.append`), fixed by the pins.
- R61b: `figures reprocessed_jaffe --scale paper --seed 42`; 4 cores, 64 GB; job 6479292, 7 min (figure legend
  redrawn in job 6479321). 213 adult controls, 93,509 CpGs in both versions.
- Result, as measured (matched gap, SD, deposited -> reprocessed): entropy housekeeping −0.218 -> −0.154,
  core-essential −0.133 -> −0.101, non-essential +0.103 -> +0.047; dispersion −0.078 -> −0.124, −0.071 -> −0.100,
  +0.040 -> +0.043; mean change −0.057 -> −0.118, −0.053 -> −0.096, −0.076 -> +0.004. Per-CpG correlation of the
  slopes: entropy 0.75, dispersion 0.56, mean change 0.82.

## R64 — published clocks' weights on protected CpGs (RM8), 2026-11-13

- Asked by antonkulaga. Commit bc1345f; `figures clock_weights --scale paper --seed 42`; 4 cores, 32 GB; job
  6472732, 21 s (figure redrawn in jobs 6472944, 6473274, 6474053). Hannum blood 656 donors, Jaffe adult cortex controls 213;
  12 clocks with 5+ protected CpGs; 1,000 resamples of each clock's CpGs.
- Result, as measured (blood; weight beyond slope and young level, 95 % CI): Horvath 22 % protected, 26 % of the
  signal, slope 0.64x, spread 0.33x, weight 1.30x, beyond slope 1.07 (0.63-1.73); PhenoAge 17 %, 9 %, 0.55x,
  0.29x, 1.10x, 0.88 (0.63-1.22); DamAge 1.66 (1.14-2.36, p 0.003); AdaptAge 0.51 (0.22-1.25); the other eight
  0.91-1.58, every interval including 1. Cortex alike (DamAge 1.69, 1.19-2.41).

## R63 — four mouse tissues 6-30 months, RRBS (RM2), 2026-11-13

- Commit 7cad5d2; `figures meer_tissues --scale paper --seed 42`; 4 cores, 48 GB; job 6469574, 44 min, 3.6 GB.
  GSE121141 (Meer 2018): male C57BL/6, 6/12/20/30 months, cortex 20, heart 21, liver 20, lung 20 mice; promoter
  CpGs covered 10+ times in 90 % of a tissue's mice (530,000-603,000); 200 resamples of mice.
- Result, as measured (housekeeping, SD, 95 % CI): entropy cortex −0.092 (−0.207 to −0.060), heart −0.100
  (−0.222 to −0.053), liver −0.356 (−0.608 to −0.193), lung −0.290 (−0.472 to −0.184); dispersion −0.066, −0.117,
  −0.324, −0.193; core essential entropy −0.137, −0.075, −0.288, −0.245. Non-essential entropy +0.093 (n.s.),
  +0.229, +0.397, +0.326. Shield share F per 10 months (housekeeping): cortex −0.013 (p 0.16), heart −0.003
  (0.68), liver −0.114 (0.014), lung −0.049 (0.010).

## R62 — one OSKM cycle, RRBS (RM17), 2026-11-13

- Commit dc71c35; `figures reprogramming_rrbs --scale paper --seed 42`; 4 cores, 32 GB; job 6469554, 43 s.
  GSE156557 (Chondronasiou 2022): female OSKM mice, young controls 13 weeks, old controls 52 weeks, old OSKM 54
  weeks (one cycle, 2 weeks off), 5 per arm; promoter CpGs covered 10+ times in all 15 mice of a tissue: pancreas
  60,109, spleen 22,903, liver 1,508. Old mice scored against the young; 5,000 label permutations.
- Result, as measured (housekeeping; effect, p): age (old − young) B pancreas −0.001 (0.77), spleen −0.001 (0.96),
  liver +0.021 (0.17); F −0.041 (0.06), −0.033 (0.31), −0.064 (0.20). OSKM (OSKM − old) B +0.004 (0.55), −0.013
  (0.19), −0.023 (0.20); F −0.005 (0.81), +0.024 (0.45), +0.067 (0.07).

## R59 — the species shield off the 0 % floor (RM13 (4)), 2026-11-13

- Commit 676bff5; `figures species_floor --scale paper --seed 42`; 4 cores, 48 GB; job 6469130, 55 min, 25.5 GB.
  GSE223748, 139 species-tissue groups, 4,485 promoter CpGs; per group only CpGs whose young-quarter mean lies at
  10-90 % (median share kept: housekeeping 24 %, core essential 29 %, non-essential 50 %).
- Result, as measured (entropy gap, housekeeping): mean −0.295 SD → −0.236; negative in 96 % → 93 % of groups,
  interval below 0 in 89 % → 66 %. Core essential −0.225 → −0.175; non-essential +0.368 → +0.212. Lifespan
  coefficient (log maximum lifespan, weight and tissue held): housekeeping entropy +0.029 (p 0.14) → −0.014
  (p 0.48); dispersion +0.036 (p 0.11) → −0.002 (p 0.92). Human marrow and lymph node stay positive (+0.28, +0.37).

## R60 — the per-person scores against Seeboth 2020's confounders (RM13 (5)), 2026-11-13

- Commit 9a75a70's code (its docstring wrapped after; the figure redrawn with a chance marker); `figures
  score_confounders --scale paper --seed 42`; 4 cores, 16 GB; job 6469169, 7.7 min. An earlier run on the login
  node (diagnosis not held, survivor instead) is superseded. SATSA 1,410 draws with smoking known; SJLIFE 4,518
  people on 604 slides. 500 resamples of people.
- Result, as measured (housekeeping; partial R² beyond age and sex, chance in brackets): SJLIFE slide B 36 %, F
  23 %, S 20 % (13 %); row on slide B 1.7 %, F 0.6 %; diagnosis B 7.7 %, F 1.6 %, S 2.8 %; SATSA smoking B 0.0 %,
  F 0.2 %, S 2.1 %; array type F 1.5 %. Terms (SD): current smoker S −0.41 (−0.64 to −0.17), F +0.15 (n.s.);
  EPIC array F +0.28; alkylating agents within diagnosis B +0.15, F −0.16 (−0.24 to −0.08), S −0.22;
  anthracyclines B −0.22, F +0.21 (0.09 to 0.33). Core essential the same.

## R58 — essentiality gradient beyond chromatin and expression (RM13 (1)), 2026-11-13

- Commit 864bda4 (decile panel redrawn by decile number after); `figures essentiality_gradient --scale paper
  --seed 42`; 4 cores, 48 GB; job 6469111, 3.7 min, 4.6 GB. DepMap 24Q4 CRISPRGeneEffect (Figshare 51064667),
  mean over lines; GTEx v8 gene median TPM; genes at 1+ TPM in 50+ of 54 tissues: 9,241 genes, about 100,500
  single-gene promoter CpGs in Hannum blood and Jaffe cortex. 200 gene resamples.
- Result, as measured (drift in SD per unit of gene effect; positive = more essential drifts less). Blood entropy
  0.130 → 0.046 (+ own chromatin state) → 0.022 (+ blood expression, breadth; 95 % 0.011-0.033); blood dispersion
  0.084 → 0.012 → −0.009 (n.s.). Cortex entropy 0.044 → 0.034 → 0.018 (p 0.06); cortex dispersion 0.044 →
  0.039 → 0.027 (p 0.01). After all held, the most essential decile sits about 0.04 SD below 0, the least 0.03
  above.

## R57 — divergence across a stem-cell transplant (RM12, bottleneck), 2026-11-13

- Commit 51d3559; `figures hsct_bottleneck --scale paper --seed 42`; 4 cores, 48 GB; job 6469087, 47 s, 1.9 GB.
  Holland et al. 2024 (Mendeley j5krtjfj6y): 334 blood samples, 92 transplants, six batches; 75,436 CpGs measured
  in all (promoter set CpGs and look-alikes). Comparisons within batch: 103 donor → recipient (median 1.0 year
  after transplant), 196 consecutive recipient draws (median 0.5 year), 4 donor duplicates.
- Result, as measured (F = set minus matched squared difference over matched; donor → recipient against recipient
  → recipient, log years and batch held, 1,000 transplant resamples): core essential −0.058 (−0.091 to −0.028),
  housekeeping −0.039 (−0.064 to −0.017), non-essential −0.057 (−0.098 to −0.018). Divergence B +50-70 %.

## R56 — a second CR blood series, GSE282661 (RM16 step 4), 2026-11-13

- Commit 7204317 (the figure redrawn after, labels only); `figures cr_blood_replication --scale paper --seed 42`;
  4 cores, 32 GB; job 6469050, 23 s. C57BL/6J whole blood, mammalian array, one run: 8 CR, 8 TLS+CR, 8 controls,
  8 at the start; long-term TLS 7 against 6. 77,888 autosomal promoter probes. Each animal scored against the
  controls (leave one out for controls); 0 and 3 composition components; 2,000 label permutations.
- Result, as measured (3 components): CR against controls, F core essential −0.090 (p_perm 0.006), housekeeping
  −0.027 (0.026); B +0.002 (n.s.). TLS+CR −0.043, −0.014 (n.s.). Controls after 12 weeks against the start:
  B −0.003, F +0.02 (n.s.): no measurable drift in 12 weeks.

## R55 — CR's blood effect with chromatin openness matched (RM12), 2026-11-13

- Commit 58ddfee; `figures cr_chromatin --scale paper --seed 42`; 4 cores, 48 GB; job 6469026, 2 min, 11.5 GB.
  GSE80672's 542,008 promoter CpGs; openness = mean ImmGen ATAC signal (GSE100738: B, naive CD4 and CD8 T, NK,
  Ly6C-high monocytes, granulocytes) within ±250 bp, in quartiles; strata island × young level decile × quartile.
- Openness: 76 % of housekeeping CpGs sit in the two most open quartiles against 41 % of background; 80 % of
  non-essential ones in the least open.
- Result, as measured. Aging in standard C57BL/6 (set minus matched, SD): housekeeping entropy −0.35 → −0.14,
  dispersion −0.42 → −0.20; core essential −0.30 → −0.14, −0.37 → −0.20; non-essential +0.36 → +0.02.
  CR in B6D2F1 (12 against 10, one series; F with 5 components removed): housekeeping +0.025 (p 0.002) → +0.008
  (p 0.066), relative to the controls' F −23 % → −18 %; core essential +0.017 (p 0.09) → +0.004 (p 0.62).
  GHRKO and Snell unchanged (n.s. both ways).

## R54 — the human marrow and lymph-node anomaly (RM16 step 2), 2026-11-13

- Commit 95c24e7 (after cb18908's run 6468984 failed in 3 min on a tissue of one age: zero spread, a division by
  zero; tissues now need a 15-year span); `figures marrow_anomaly --scale paper --seed 42`; 4 cores, 32 GB; job
  6468993, 3.5 min, 4.3 GB. GSE223748's human samples: 13 tissues of 20-282 donors, 4,485 promoter CpGs measured in
  every marrow and lymph-node sample. 300 donor resamples per tissue; 400 draws of 20 donors (ages 38-74) from each
  of 9 larger tissues.
- Result, as measured (housekeeping entropy gap, SD): marrow +0.175 (donors resampled −0.27 to +0.29), lymph node
  +0.198 (−0.08 to +0.24); chip held +0.20 and +0.20. 20-donor draws of other tissues reach marrow's gap in 6.8 %,
  lymph node's in 5.1 %. Seven tissues refitted on marrow's donors (13-25 shared) move from −0.15 to +0.004 on
  average (+0.16). Dispersion gaps of both lie inside the other tissues' range (+0.04, +0.06). Marrow's per-CpG age
  slopes follow the adipose-minus-blood direction, the most positive correlation of all tissues (r 0.28; next 0.13).

## R53 — read-level shield in sorted atlas cell types (RM22), 2026-11-13

- Commit ada6425; `figures read_level --scale paper --seed 42`; 4 cores, 48 GB; job 6460007 (6459198 at 8 cores
  cancelled after waiting on the group CPU cap), 51 min. 101 sorted samples of the atlas (GSE186458) in 27 cell-type
  and tissue groups with at least 3 samples and a 15-year age spread (ages 19-83); 69,373 windows with 10+ reads in
  76+ samples; age slope with each group a fixed effect (sex held). Some groups share donors (pancreas, blood T cells).
- Result, as measured (set minus look-alikes, SD; p): entropy slope housekeeping -0.084 (0.001), core essential
  -0.044 (0.001), non-essential +0.19 (0.001); discordance slope -0.060 (0.001), -0.029 (0.002), +0.23 (0.001).
  Levels at matched methylation near 0 for the two sets. WBC and cfDNA rows reproduced from cache, unchanged.
  Figure redrawn on the login node from the table (sorted panel labelled samples, not donors).

 over time, SATSA (RM28), 2026-11-13

- Commit 95240c5; `figures within_person --scale paper --seed 42`; 4 cores, 48 GB; job 6457453, 2 min 58 s. SATSA
  (S-BSST1206; 6.7 GB, preprocessed by the depositors: noob, functional normalisation, cell counts and batch
  regressed out): 1,469 draws of 596 people, 392 with 2+ draws; 65,386 promoter CpGs without gaps; 873 consecutive
  pairs, median 4.9 years apart.
- Result, as measured (set minus look-alikes, SD; p): divergence per year (change one person makes beyond the shared
  trend, measurement floor removed) housekeeping -0.017 (0.001), core essential +0.013 (0.10), non-essential +0.11
  (0.001); |directed change| about 0 for both sets; floor -0.03 (0.001), -0.01 (0.14). Divergence > 0 for 81 % of
  CpGs. Within-person slope per decade (SD; person-level bootstrap 95 %): F housekeeping -0.14 (-0.23 to -0.05),
  core essential -0.07 (-0.18 to 0.02); S -0.19 (-0.27 to -0.11), -0.16 (-0.24 to -0.07); B +0.03, +0.02 (both
  span 0). Between people at first draws: F -0.13, -0.06; B +0.17, +0.15; S -0.06, -0.01.
- R52b (positive control; commit 384ae50; job 6457917, 4 min 17 s, after a wait on the group CPU cap): the usual
  cross-sectional comparison on first draws, same data. Dispersion slope housekeeping -0.12 (0.001), core essential
  -0.04 (0.001); entropy slope +0.00 (0.92), +0.02 (0.19); |age slope| +0.00, +0.00. The preprocessed SATSA betas
  carry a third of Hannum's dispersion gap and none of its entropy gap: SATSA is a weak test of the CpG-level
  shield. Figure redrawn on the login node from the tables (layout only).

## R45 — when does drift begin? change point against sexual maturity (RM29), 2026-11-12

- Commit cc4bc08; `figures drift_onset --scale paper --seed 42`; 4 cores, 64 GB; job 6439005, 24 min. 68
  species-tissue groups with animals before maturity and after twice maturity, 8,671 postnatal animals.
- Result, as measured: noise per animal (development removed) rises with age in 84 % of groups for look-alike
  CpGs (42 of 68 significantly, 6 significantly falling; median Spearman 0.28) and in 68 % for housekeeping CpGs
  (median 0.14). A supported flat-then-rising onset in 6 groups (background) and 4 (housekeeping); the 6 sit at
  2.4-18 times maturity (beluga skin, marmot blood, minipig blood, mouse liver, muscle, cortex). Too few for the
  across-species regression (needs 8). Tables in `data/output/drift_onset/`. Panel titles a and b overlap in the
  figure; shortened in the code for the next run.

## R39 — health states and the shield, cross-sectional (RM24), 2026-11-12

- Commit 9996400; `figures health --scale paper --seed 42`; 4 cores, 64 GB; job 6436141, 1 h 42 min, 2.7 GB.
  Breast left out of the Hardy contrast (too few).
- Result, as measured: 522 effects (B, F, S x 3 sets x raw/5 components x 29 contrasts); 50 at p < 0.05, 19 at
  FDR q < 0.05: schizophrenia in GSE41037 blood 8, in Jaffe cortex 3; smoking in GTEx 4 (lung S housekeeping
  -0.67 SD raw, -0.50 adjusted; colon B -0.53 adjusted); slow against sudden death 4 (skeletal muscle F +1.3 to
  +1.45 SD, 12 against 12 donors). Alzheimer's and Braak stage (Lunnon, 5 tissues): none; Braak effects per stage
  |0.12| SD or less. Housekeeping F, 5 components: Jaffe cortex schizophrenia +0.57, GSE41037 blood +0.45.
  Tables in `data/output/health_states/`.

## R51 — read-level replication in cfDNA of the same donors (RM22), 2026-11-13

- Commit ba26233; `figures read_level --scale paper --seed 42`; 8 cores, 32 GB; job 6447493, 12 min. cfDNA of the
  23 donors (GSM6810003-25); 51,902 windows with coverage (WBC counts reused from cache).
- Result, as measured (set minus look-alikes, SD; p): entropy slope housekeeping -0.063 (0.001), core essential
  -0.058 (0.002); discordance slope -0.054 (0.001), -0.033 (0.019). White blood cells: -0.075, -0.049, -0.059,
  -0.034. Non-essential: cfDNA -0.02 (0.48) and +0.06 (0.06), against WBC +0.15 and +0.22.

## R50 — molecule-level disorder at essential promoters, Loyfer WBC (RM22), 2026-11-13

- Commit 7d95064; `figures read_level --scale paper --seed 42`; 8 cores, 32 GB; job 6445425, 11 min. 23 WBC donors
  aged 21-75; 109,068 4-CpG windows at 450k promoter probes, 73,885 with 10+ reads in 18+ donors (median 41 reads).
- Check: window mean methylation against Hannum's array mean for the same probes, r 0.963 (73,804 probes): the CpG
  index and the window placement are right.
- Result, as measured (set minus look-alikes matched on island and level, SD; permutation p): age slope of epiallele
  entropy, housekeeping -0.075 (0.001), core essential -0.049 (0.001), non-essential +0.148 (0.001); of discordant
  reads, -0.059 (0.001), -0.034 (0.002), +0.222 (0.001). Levels at matched methylation: entropy -0.007 / +0.013,
  discordance +0.016 / +0.017. Tables in `data/output/read_level/`.

## R48 — SJLIFE deaths (RM28), 2026-11-13, first pass

- Commit 60b685b; `figures sjlife --scale paper --seed 42`; 4 cores, 96 GB; job 6444474, 9 min. 4,518 first draws,
  277 deaths (222 over 2 years after the draw), 92,030 CpGs loaded from the 28 GB table.
- Result, as measured (odds ratio of death per SD; age, sex, 4 exposures, diagnosis held; no control has died, so
  "everyone" equals "survivors"): PhenoAge 1.62 (p 1e-10), DunedinPACE 1.78 (p 1e-17), Zhang 1.28 (p 8e-4) -- the
  positive control holds. Drift B 1.06-1.09 (p 0.12-0.41). Shield share F 0.84-0.87 (p 7e-4 to 0.005) and signed
  shift S 0.82-0.84 (p 0.002-0.011): a stronger shield goes with more deaths here, unexpected.
- Rerun with mean methylation, 11 cell fractions and chip held (R48b; module changed in the commit below).
- R48b: commit ebafd8e; job 6445940, 9 min 48 s; survivors only (no control died); 573 EpiDISH reference CpGs
  (PBMC: neutrophils 0.03). Adjusted (mean beta of loaded CpGs, 11 cell fractions, measure centred within chip):
  DunedinPACE still 1.70 (p 2e-11), PhenoAge 1.45 (p 7e-5); drift B housekeeping 1.24 (p 8e-4); shield share F 0.76
  (p 6e-5); signed shift S 0.71 (p 4e-5).
- Joint models on the run's per-person table (login node, seconds; chip-centred, age, sex, exposures, diagnosis):
  with DunedinPACE beside them, B 1.01 (p 0.83) and F 1.00 (p 0.95); all four together B 1.02, F 1.01,
  DunedinPACE 1.63 (p 2e-9), PhenoAge 1.22 (p 0.02). Correlations among survivors: B with DunedinPACE 0.31, F with
  B -0.45. Same check on EPIC-Italy (women, breast over 2 years, chip-centred): B 1.46 alone, 1.44 (p 9e-4) with
  DunedinPACE (0.99) and PhenoAge (1.12) beside it; B with DunedinPACE r 0.02.

## R49 — the shield by blood lineage, in silico (RM31), 2026-11-13

- Commit 1a753f2; `figures lineage_shield --scale paper --seed 42`; 4 cores, 32 GB; job 6444485, 3 min 21 s. Hannum
  656 donors; lineage fractions from R38 (myeloid 0.66, T cells 0.20, B and NK 0.14 on average).
- Result, as measured (set minus look-alikes on |age slope within the lineage|, SD; permutation p): T cells
  housekeeping -0.103 (0.001), core essential -0.092 (0.001); B and NK -0.043 (0.001), -0.055 (0.001); myeloid
  -0.020 (0.001), -0.018 (0.008). Non-essential also below its look-alikes here (-0.09 myeloid, -0.03 T, -0.04 B/NK).
  Table `data/output/lineage_shield/gaps`.

## Index — hg38 CpG index for read-level WGBS (RM22), 2026-11-13

- Commit 45814a7; job 6445117 (`entropy-cpgindex`), 46 s. 29,401,795 CpGs over chr1-22, X, Y, M from the GRCh38
  primary assembly: exactly the length of Loyfer's `.beta` files. Chromosome starts agree with the first reads of
  each chromosome in a PAT file (chr3 4,567,830; chr4 6,241,123; chr5 7,744,552). Written to
  `data/interim/annotation/cpg_index_hg38.parquet` (not in git).

## R47 — directed against noisy aging under exenatide and rapamycin (RM27), 2026-11-12

- Commit 54e0d43; `figures intervention_split --scale paper --seed 42`; 4 cores, 32 GB; job 6441201, 36 s; the
  figure redrawn with stacked panels (job 6442268). 11 tissues; rapamycin in 8 (5 mice each), exenatide in 11.
- Result, as measured (reversal minus permutation null, one-sided p): exenatide, no tissue p < 0.05 for either
  part. Rapamycin, noise part: colon +0.64 (p 0.04), adipose +1.30 (p 0.02), hippocampus +0.51 (p 0.015), blood
  +0.55 (p 0.08), frontal cortex +0.54 (p 0.07), muscle +0.27 (p 0.07), heart -0.60 (p 0.94); 7 of 8 positive.
  Directed part: adipose +0.71 (p 0.05), muscle +0.11 (p 0.03), others p > 0.17. Table
  `data/output/intervention_split/reversals`.

## R46 — incident cancer in EPIC-Italy (RM28), 2026-11-12

- Commit fe9e9a1; `figures cancer --scale paper --seed 42`; 4 cores, 64 GB; job 6440396, 6 min 32 s. 845 people,
  421 later cancers (235 breast, 166 colorectal), 31,564 promoter CpGs without gaps.
- Result, as measured (odds ratio per SD, age and sex held): drift B of look-alike CpGs predicts breast cancer, 1.43
  (1.19-1.71, p 1e-4) for all cases and 1.59 (1.31-1.93, p 2e-6) for cases diagnosed over 2 years after the draw;
  any cancer 1.17 (p 0.03); colorectal 0.75 (p 0.04). Shield share F and signed shift S: none (0.97-1.00). Clocks:
  Horvath colorectal 1.23-1.27 (p 0.02-0.03), DunedinPACE 1.13-1.21 (p 0.04-0.08), PhenoAge 1.13 (p 0.08).
- Checks on the login node from the run's per-person table (minutes): cases are balanced over the 72 chips (median
  case share 0.5; chip mean B unrelated to case share, rho -0.13); B centred within chip still predicts breast
  cancer, 1.36 (p 0.001); women only, B 1.43 / 1.62 / 1.66 for diagnoses over 0 / 2 / 5 years after the draw
  (p 3e-4, 1e-5, 2e-5), PhenoAge 1.09-1.14 and DunedinPACE 1.13-1.25 beside it.
- Open: van Veldhoven 2015 reported genome-wide hypomethylation before breast cancer in this cohort (not replicated
  in NOWAC). B may carry it, or blood-cell mix. Rerun with mean beta, 11 cell fractions and chip held (R46b).
- R46b (commit 4e5c3b4; job 6441619) failed after 9 min: 72 chip dummies in a logistic fit went singular where a
  chip held only controls. Now the measure is centred within its chip instead, and the Newton step uses least
  squares (test added); resubmitted.
- R46b rerun: commit 7656dd1; job 6442500, 10 min. Mean beta over all CpGs 0.533; 600 EpiDISH reference CpGs;
  neutrophils 0.51. Breast cancer, diagnosed over 2 years after the draw, adjusted (mean beta, 11 cell fractions,
  measure centred within chip): B housekeeping 1.44 (1.17-1.76, p 6e-4), B core essential 1.46 (p 3e-4); all cases
  1.37-1.40 (p 0.001). DunedinPACE breast 1.21 -> 1.04, Horvath colorectal 1.27 -> 1.01, B colorectal 0.78 -> 1.00
  once adjusted. Over all 168 effects, the 8 breast-B variants alone pass FDR 5 % (q 0.0004-0.027).

## R42 — survival in Vitality 90+ (RM28), 2026-11-12, partial

- Commit cc4bc08; `figures survival --scale paper --seed 42`; 4 cores, 32 GB; job 6439004, 1 min 27 s. 115
  nonagenarians, 39 deaths in about 2.5 years.
- Result, as measured (hazard ratio per SD, sex and batch held; then with CD8+CD28- share): shield share F 1.08
  (0.78-1.48) core essential, 1.08 housekeeping; drift B 0.84 (0.58-1.19), 0.84; signed shift S 1.19 (0.87-1.63),
  1.02; all p > 0.26, unchanged with CD8+CD28-. MMSE: no measure p < 0.2.
- Every clock was left out: a clock CpG absent from this series' table is all NaN, and the linear score turned NaN
  for everyone. `linear_score` now skips CpGs no sample has (test added); rerun for the clocks.
- Rerun with clocks: commit 6cf2348; job 6440151, 1 min 48 s. Hazard ratio of death per SD (sex, batch held):
  PhenoAge 0.98 (0.72-1.34), DunedinPACE 1.03 (0.75-1.43), Horvath 1.16, Hannum 0.86, Zhang 0.92, StocH 0.92,
  StocP 0.72 (0.50-1.03, p 0.07); MMSE: none p < 0.18. Tables in `data/output/survival_shield/`.

## R44 — the shield, drift and clocks in sorted blood cells (RM31), 2026-11-12

- Commit d891500; `figures sorted_cells --scale paper --seed 42`; 4 cores, 64 GB; job 6436739, 31 min.
  MESA CD14+ monocytes 1,202 people (GSE56046), CD4+ T cells 214 (GSE56581), ages 44-83.
- Result, as measured (set minus matched look-alikes, SD; permutation p 0.001 unless given): monocytes
  housekeeping entropy -0.14, dispersion -0.14, |mean change| -0.07; core essential -0.08, -0.06, -0.06. T cells
  housekeeping -0.22, -0.04, -0.09; core essential -0.09, -0.01 (p 0.39), -0.06.
- Clocks, r with age in monocytes / T cells: Zhang 2019 0.92 / 0.92, Hannum 0.84 / 0.80, Horvath 0.78 / 0.82,
  PhenoAge 0.77 / 0.71, StocZ 0.81 / 0.82, DunedinPACE 0.09 / 0.27; residual contamination scores explain 0.4-1.6 %
  of acceleration in monocytes, 3.6-10 % in T cells.
- Health, monocytes (per SD of the measure; age and race-sex-site held): DunedinPACE plaque +0.37 (p 4e-8), log
  coronary calcium +0.37 (p 1e-7); PhenoAge +0.14 for both (p 0.04, 0.05); Zhang 2019 log CAC +0.16 (p 0.02);
  Horvath, Hannum, Stoc clocks n.s.; shield F housekeeping plaque -0.14 (p 0.054), S -0.12 (p 0.09), the rest
  p > 0.09. Tables in `data/output/sorted_cells/`.

## R41 — organs of one donor aging together (RM26), 2026-11-12, rerun

- Commit 522d4e6; `figures organ_coupling --scale paper --seed 42`; 4 cores, 64 GB; job 6438960, 10 min 45 s.
  Lunnon 5 tissues (49-84 shared donors per pair), GTEx 9 tissues (22 pairs with at least 20 shared donors,
  median 33).
- Result, as measured (Pearson r of residual acceleration between two tissues of the same donors; median over
  pairs, GTEx / Lunnon): Zhang 2019 0.31 / 0.17 (GTEx all 22 pairs positive), Horvath 0.29 / 0.22, Cortical
  0.17 / 0.33, DunedinPACE 0.16 / 0.14, Hannum 0.15 / 0.17, PhenoAge 0.13 / 0.07; clock-free drift B -0.03 / 0.07
  and shield share F 0.04 / 0.02, no GTEx pair significant. Lunnon cortical regions couple strongly with each
  other (Cortical clock 0.75-0.76); blood with brain weakly (-0.15 to 0.34, intervals include 0 for most).
  Tables in `data/output/organ_coupling/`.

## R38 — the stochastic component of clocks, information stripped (RM23), 2026-11-12

- Commit 9099a66; `figures stochastic_component --scale paper --seed 42`; 4 cores (resubmitted from 8; job
  6435307 never started), 32 GB; job 6438042, 3 min 41 s. Hannum blood, 631 test donors; sklearn reported
  non-convergence for some noise and shield fits (those clocks read no age anyway).
- Result, as measured, RR2 = R²(stochastic) / R²(original) for Horvath / PhenoAge / Zhang 2019: published Stoc
  clocks 0.67 / 0.65 / 0.89; Tong's model 0.64 / 0.85 / 0.86; direction only 0.71 / 0.70 / 0.80; undirected
  noise 0.02 / 0.10 / 0.00; shield only 0.07 / 0.03 / 0.00. Correlation of the stochastic clock's acceleration with
  the original's, before and after 11 of 12 cell fractions are removed (EAA -> IAA), PhenoAge: published StocP
  0.44 -> 0.45, Tong model 0.51 -> 0.49.
- Added on the login node from the run's cell fractions and Hannum betas (minutes, 2 threads): the share of each
  clock's acceleration (residual on age) that the 11 fractions explain: PhenoAge 0.155, Horvath 0.049, Zhang 2019
  0.043, Hannum 0.138, DunedinPACE 0.170; StocH 0.158, StocP 0.049, StocZ 0.072. Not yet in the module.
- Tables in `data/output/stochastic_component/`. The figure's title runs past the right edge; to shorten.

## R37 — Dnmt3a gain of function (RM12), 2026-11-12

- Commit 75469a1; `figures dnmt3a --scale paper --seed 42`; 4 cores (resubmitted from 8; job 6435278 never
  started), 32 GB; job 6438041, 19 s.
- Result, as measured (12 W326R/+ against 12 wild types per tissue, age group held, CpGs centred on wild types of
  each age): shield share F -0.53 to -0.61 for core essential and housekeeping sets in bone marrow, intestine and
  liver (all p < 1e-16), +0.56 to +0.67 for non-essential; drift B +0.06 to +0.10 everywhere; signed shift S
  -0.06 to -0.09 for essential and housekeeping, +0.07 to +0.09 for non-essential. Tables in
  `data/output/dnmt3a_gain/`.

## R36 — Tet3 signed shift with array compression held (RM18), 2026-11-12

- Commit 5b769d5; `figures mmc_interventions --scale paper --seed 42`; 4 cores (resubmitted from 8; job 6435202
  never started), 16 GB; job 6438040, 1 min 26 s.
- Result, as measured (core essential, S raw -> S with compression mean |beta - 0.5| as a covariate): Tet3 cortex
  against same-chip littermates +0.0120 -> +0.0101 (p 0.003), striatum +0.0142 -> +0.0128 (p 0.001); pooled
  controls cortex +0.0116 -> +0.0107. Other arms shift S too, in both directions: GHRKO kidney +0.023, GHRKO
  cortex -0.010, high-fat liver -0.005. Table `data/output/mmc_interventions/effects`.

## R42 and R45 — failed once each, 2026-11-12

- R42 (survival, Vitality 90+; commit c4691ff; job 6436543) failed after 1 min 30 s: a score with missing values
  reached the MMSE regression ("SVD did not converge"). Such a measure is now left out with a line in the log.
- R45 (drift onset; commit 8ae0bf6; job 6437491) failed after 1 min 45 s: the consortium dates 267 fetal samples
  with negative ages, and the trend's log age became NaN. Only postnatal animals are used now.
- Both fixed in the commit below and resubmitted.

## R43 — methylation spread against expression spread, GTEx same donors (RM30), 2026-11-12

- Commit 3b73b5d; `figures meth_expr_spread --scale paper --seed 42`; 4 cores, 64 GB; job 6436659, 49 s, 1.3 GB.
- Donors with both layers: lung 161, ovary 133, colon 98, prostate 62, blood 45 (testis, kidney, breast, muscle
  under 40, left out); 8,500-11,200 genes each.
- Result, as measured: signed coupling (control, expected < 0): lung -0.032, prostate -0.035, colon -0.027, ovary
  -0.011 (54-61 % of genes negative), blood +0.058. Spread coupling per donor, mean |r| minus permutation null:
  lung +0.017 (p 0.025), prostate +0.022 (p 0.06), ovary +0.015 (p 0.06), colon +0.002 (p 0.45), blood +0.026
  (p 0.15). Spread-slope coupling across genes (Spearman, level held): prostate +0.052 (p 4e-8), colon +0.028
  (p 0.004), ovary +0.020 (p 0.046), lung -0.019 (p 0.04), blood -0.014. Tables in `data/output/meth_expr_spread/`.

## R41 — organ coupling (RM26), 2026-11-12, failed once

- Commit 9599c31; `figures organ_coupling --scale paper --seed 42`; 4 cores, 64 GB; job 6436308 failed after 31 s:
  the cohort's display name ("Lunnon") was passed where the loader expects its key ("lunnon"). Fixed in the
  commit below; resubmitted.

## R40 — tissue aging pace across species and layers (RM25), 2026-11-12

- Uncommitted code, then the commit below; `figures tissue_pace --scale paper --seed 42` on the login node
  (reads result tables only, seconds).
- Result, as measured (Spearman over shared tissues, median slope of features outside the gene sets): directed
  change (|mean slope|) mouse against rat RNA rho 0.56 (20 tissues, p 0.011); noise (dispersion slope) -0.27
  (p 0.25). Human methylation against human RNA, 9 tissues: 0.03 and 0.35 (both p > 0.3). Human against mouse
  methylation: 4 tissues, 0.80 (p 0.20, too few). Tables in `data/output/tissue_pace/`.

## R33 — clocks of our own (RM11), 2026-11-12

- Commit 2bc70dc; `figures own_clocks --scale paper --seed 42`; 8 cores, 64 GB; job 6433984, 3 h 12 min, 4.3 GB.
- Result, as measured (out-of-fold, 2 draws; pools of 38,548 blood and 38,356 cortex CpGs): Pearson r with age
  blood 0.884 shielded, 0.882 look-alike, 0.885 random; cortex 0.929, 0.943, 0.938. Acceleration overlap: blood
  look-alike with random 0.80, shielded with either 0.63; cortex 0.62 and 0.53-0.57.

## R35 — Tet1 knockout and catalytic-dead cortex (RM12), 2026-11-12

- Commit 3d9c573; `figures tet1 --scale paper --seed 42`; 8 cores, 32 GB; job 6434508. 24 animals, 65,779
  promoter CpGs.
- Result, as measured (against 8 wild types; raw / 3 components): KO housekeeping B +0.033/+0.040, F
  −0.130/−0.105, S −0.018/−0.018; core-essential F −0.170/−0.135; Hxd housekeeping F −0.130/−0.125, B +0.029/+0.031;
  all p ≤ 0.014.

## R34 — Tet1 knockout and catalytic-dead cortex (RM12), 2026-11-12

- Commit 4cb8a47; failed after 8 s: the betas file writes missing values as `NA`. Read as missing now; rerun R35.

## R32 — partial reprogramming and the shield (RM17), 2026-11-12

- Commit 803d1da; `figures reprogramming --scale paper --seed 42`; 8 cores, 16 GB; job 6433970, 8 s.
- Result, as measured (treated against untreated old mice, 8-11 against 8-14 per tissue): F housekeeping −0.115
  (muscle, p 0.02), +0.104 (spleen, p 0.046), core-essential +0.211 (spleen, p 0.001), other tissues n.s.; B
  unchanged in every tissue (p ≥ 0.14). Young against old (4 each): B higher in skin (p 0.001) and liver
  core-essential (p 0.011), S higher in all six tissues.

## R31 — clocks of our own (RM11), 2026-11-12

- Commit 1637fe8; job 6432527 cancelled after 2 h 32 min: one cohort-draw took 2.5 h on one core (the
  elastic-net search ran single-threaded), six would have taken 15 h. Now 2 draws, search on every core. R33.

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

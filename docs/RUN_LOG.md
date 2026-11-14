# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

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

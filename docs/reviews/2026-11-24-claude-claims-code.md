run claude opus medium (the engine's own answer). Read: copy at $S/src (paper.md whole, REASSESSMENT s2, CLAIMS C1-C2 head,
protection.py, stochasticity.py, chromatin_match.py, common.py, series.py + config/series, brain_replication.py,
data/reprocessed.py, tissue_contrast.py helpers, floor_check.py). Tables opened in the repository at check time.
Ran: read-only uv/polars reads of data/output tables (species_shield parquet, lnlifespan range). No writes.

FINDINGS, by severity
1 CLAIM  paper.md:32,440,468-470,576 "not tuned to species lifespan": species_floor/lifespan_fits.csv off-floor kl_slope
  housekeeping +0.004/ln-unit, se 0.020 -> 95% up to +0.044; 5-95% of species span 2.9 ln-units, so the interval admits a
  change of ~0.13 SD, the size of the mean gap (-0.13). Say "not resolved; cannot exclude a lifespan relation of the
  shield's own size".
2 CODE/CLAIM  band SD: brain_replication.group_gaps -> chromatin_match.matched_gaps:124 and series.py:119, and
  tissue_contrast._matched_sd:195 via floor_check, scale each band's gap by the SD of that band's CpGs. floor_check row
  blood hk: SD 1.43e-4 (all) vs 1.90e-4 (10-90 %). Raw gaps -2.02e-5 vs -1.48e-5: in all-CpG SD the off-floor gap is
  -0.104, not -0.078; cortex -0.080 not -0.066. paper.md:90-93 "much of the blood gap lies near the floor" compares
  numbers in different units. Report band gaps on one SD, or say per band.
3 CLAIM  paper.md:92 "in cortex the gap kept its size off the floor (-0.07, -0.11 to +0.00)": floor_check row is
  -0.066 (-0.113, +0.003), interval includes 0; 20-80 % cortex -0.048 (-0.111, +0.054) unreported.
4 CLAIM  paper.md:319-321 Tet3: "Essential promoters gained methylation ... S +0.012 and +0.014" are core_essential
  rows; mmc_interventions/effects.csv same_chip housekeeping S cortex +0.0023 p 0.14 (only with compression held
  +0.0044 p 0.006). Say which set; housekeeping S in cortex is not resolved. Discussion 608 "in brain it depends in part
  on Tet3" is causal language from one study of 8 het mice per tissue, no raw files.
5 CLAIM  internal contradiction: paper.md:245-248 consortium mouse array "separated brain ... though not by turnover
  class" vs Limits paper.md:588-589 "on the consortium's array it [ordering by turnover] is seen with chip intercepts"
  and CLAIMS C1 "ordered by turnover within chips, R118". Pick one.
6 HYPOTHESIS  "goes with cell division" (title, abstract 22, Discussion 553): rests on one 2-tissue mouse cohort (gaps
  in each tissue's own SD: colon SD 0.0032 vs hippocampus 0.0006, jena_rrbs/gaps.csv), one glial fraction p 0.025 at
  the BH edge (nuclei_paired: SOX10+ -0.003; second series -0.044, -0.108..+0.015), and culture series with no arrested
  control. Against: Huang 11 tissues reverse within chip, consortium not by turnover class, GTEx testis ~0, Endicott
  shield gone on rebuild, BJ mortal cells -0.026/+0.087. Discussion 553-556 drops the caveats ("larger in glia than in the
  same donors' neurons"). Verdict: supported in part (single contrasts), not across tissues. Unconsidered alternative:
  tissue-specific promoter activity (colon vs hippocampus differ in which promoters are active).
7 CLAIM  "cultures kept from dividing show no measurable drift" (abstract 25, Discussion 555): arrest_timecourse 2
  donors x 6 arrays, permutation of 720 orders; third donor (drift ratio 1.96, p 0.007) excluded for array order. An
  absence from a low-power test; say "no drift detected in two donors".
8 CLAIM  abstract 23-24 "accrues per population doubling in fibroblasts": Endicott's 7 lines lose the shield on rebuilt
  arrays (paper 448-449), BJ mortal cells -0.026 / +0.087. Restrict to "in one fibroblast series (Sturm)".
9 CODE  p-values: stochasticity._permute_within_strata:336 shuffles CpG labels; CpGs of one promoter are correlated, so
  permutation p (e.g. "all p <= 0.001", paper 83) is anti-conservative and feeds the global BH (13,955 tests).
  Permute genes, or report gene-resampling intervals only.
10 CODE  resampling not identical across pipelines: series.donor_intervals fixes strata and band from the full sample;
  tissue_contrast._boot_lifespan:212 and floor_check.gaps:47-50 recompute young reference, deciles and band membership
  in every donor resample. Point estimates unaffected; intervals differ in kind. Pick one and say which.
11 CODE/PROVENANCE  common.write_table:248-253 writes CSV only for <=1000 rows and never removes an old CSV:
  species_shield/species_gaps.csv (Nov 12) is stale and has no kl_slope rows; the paper's 89 % / -0.13 come from the
  gitignored parquet (Nov 14: 88.5 %, -0.128). Lifespan fits (Nov 14) predate the relative_entropy fix 77aa712 (bound
  0.005 SD measured on headline gaps only, RUN_LOG R81).
12 PROVENANCE  paper 282-283 culture-resampling intervals for BJ hTERT and astrocytes come from untracked files
  (bj_telomerase/culture_bootstrap.csv, astrocyte_passage/culture_bootstrap.csv: git ??). Commit them.
13 CLAIM  astrocyte_htert culture interval (-0.089, -0.080) is narrower than the gene interval (-0.092, -0.074) for 3
  replicate lines of one immortalisation: resampling within lines cannot carry between-line variance; n is ~1.
14 CLAIM  chromatin, cortex: paper 18, 348, 539 "a quarter to a half" and 358-359 "25 % / 56 %" are deposited Jaffe
  (chromatin_match/matched_gaps.csv: 24.8 %, 55.6 %); on the rebuilt arrays state alone removes 43.9 % (38-52 %)
  and expression alone 45.6 %, all three 77.2 % (joint_cortex_rsesame/gaps.csv). "A quarter to a half" still covers 44 %
  but the text quotes deposited numbers where a rebuild exists; quote 44 % rebuilt. Blood 83-91 % holds (82.9, 91.2).
  Monocytes: E062 state alone removes 101 % (mesa_expression), more than expression's 82 %, unreported.
15 CLAIM  breast cancer: holds as worded (cancer_matched rows: 1.76 deposited; rebuilt 1.11, boot 0.69-2.18). The loss
  is from adjustment on rebuilt arrays (age-only rebuilt 1.84), not processing alone; F "protective" 0.69 (p 0.015)
  is one of many uncorrected tests.
16 CODE (latent)  reprocessed.read_qc:62-63 passes every sample when frac_dt is absent; all 32 rebuilt QC files have
  frac_dt today. cohort_betas:121 returns a cached parquet forever, so a changed QC or rebuild is not picked up.
Holds as checked: expression 71-82 % (mesa 82.4, sturm 72.2, epidermis 71.0/70.7); jena colon -0.189 vs hippocampus
  -0.078, diff -0.111 (-0.222,-0.085); blood off-floor -0.078 (-0.153,-0.027), 20-80 % +0.0001; neuron series; SD
  for resamples = full-sample SD everywhere read (scaled_intervals, rescale, series, chromatin_match).
Not read: scripts/sesame_preprocess.R (not in copy), most experiment modules, supplement S3/S4 in full, RUN_LOG beyond
  R77-R81, senescence_modes/donors, mouse_atlas_rrbs tables, SATSA, GTEx, placebo tables.

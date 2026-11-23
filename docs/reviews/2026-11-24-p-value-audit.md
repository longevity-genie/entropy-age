# P-value audit: manuscript/paper.md lines 1-680

Repository: /data/uq2347/sources/entropy-age. Read-only audit; nothing in the repository was edited.

Rule applied. A p-value is SET-LABEL (invalid as evidence for an age association) when it is computed with the samples
held fixed and CpGs or genes as the units:
- the stratified permutation of gene-set labels among CpGs (`compare_groups` → column `p_permutation`);
- the per-gene residual test (`p_gene`);
- a test or bootstrap across genes of per-CpG age slopes, where the samples are again fixed (Spearman over genes,
  `p_spearman`; gene-resampling `p_boot`). These share the defect: the same fixed samples' noise is spread across
  every gene, so the test counts genes, not donors.

A p-value is VALID when donors, animals or cultures are resampled or permuted, treatment labels are permuted across
samples, or it comes from a model fitted across samples or groups (OLS on per-animal scores, logistic regression
across people, OLS across species-tissue groups).

Counts: 33 p-values (counting each number quoted); 25 VALID, 8 SET-LABEL (at lines 142, 156, 250 ×3, 412, 443 and
449). There are 6 replacement edits; line 250 takes 2.

## Table

| Line | Fragment | Source table | Method | Verdict | Old string → new string |
|---|---|---|---|---|---|
| 24 | "(p 0.025; not resolved in a second series)" | nuclei_paired/paired.csv (glia NeuN-/SOX10-, kl_slope, housekeeping, p_bootstrap 0.025) | nuclei_paired.py: donor bootstrap of both fractions together, 400 resamples | VALID | keep |
| 142 | "(−0.015 SD, p 0.005" | sorted_cells/shield.csv (CD4+ T, kl_slope, housekeeping: matched_diff_sd −0.0146, p_permutation 0.006, p_gene 0.0022) | sorted_cells.py → tissue_contrast compare_groups: set labels permuted among CpGs | SET-LABEL | `(−0.015 SD, p 0.005; Supplementary Fig. S10)` → `(−0.015 SD, −0.024 to −0.006 over genes; no donor resampling interval computed; Supplementary Fig. S10)` (gene interval from mesa_expression/matched_gaps.csv, CD4+ island x level) |
| 156 | "(−0.060, −0.029; all p ≤ 0.002" | read_level/gaps.csv (pdr_slope; white blood cells −0.059/−0.034, sorted cells −0.060/−0.029; p_permutation 0.001-0.002) | read_level.py: compare_groups-style stratified permutation of set labels | SET-LABEL (all 4) | `(−0.060, −0.029; all p ≤ 0.002; Supplementary Fig. S11)` → `(−0.060, −0.029; no resampling interval computed; Supplementary Fig. S11)` |
| 234 | "p 0.04" | species_floor/lifespan_fits.csv (all CpGs, kl_slope, housekeeping, gap_sd ~ log max lifespan: +0.031, p 0.042) | species_floor.py: OLS across 139 species-tissue groups (weight, tissue held) | VALID | keep |
| 234 | "0.05 with groups clustered by species" | not found in the current tables or code (no clustering term in the species_*.py modules) | regression across groups, species-clustered SE | VALID (method). The number is not traceable in the current tables. | keep; source to be confirmed |
| 235 | "(+0.004, p 0.83" | species_floor/lifespan_fits.csv (off floor and ceiling, kl_slope, housekeeping: +0.0042, p 0.834) | OLS across groups | VALID | keep |
| 250 | "cortex −0.03, −0.11; p ≤ 0.02)" (headroom, "closer to where the site started") | equilibrium_drift/per_cpg_gaps.csv (headroom, island x level; p_permutation 0.001; cortex core 0.020) | equilibrium_drift.py: per-CpG fits, set labels permuted among CpGs | SET-LABEL | `−0.42; cortex −0.03, −0.11; p ≤ 0.02)` → `−0.42; cortex −0.03, −0.11; 95 % donor intervals blood −0.31 to −0.12 and −0.47 to −0.17, cortex −0.11 to −0.02 and −0.16 to −0.07)` (per_cpg_interval/per_cpg_gaps.csv, 100 donor resamples, R-entry at RUN_LOG line 1007) |
| 250 | "p 0.24" (cortex core k_loss) | equilibrium_drift/per_cpg_gaps.csv (k_loss, cortex, core_essential, p 0.241) | as above | SET-LABEL | see next row |
| 250 | "and 0.02" (cortex housekeeping k_loss) | same (k_loss, cortex, housekeeping, p 0.019) | as above | SET-LABEL | `(−0.02 and −0.02; p 0.24 and 0.02)` → `(−0.02 and −0.02; 95 % donor intervals −0.05 to +0.05 and −0.03 to +0.04)` (per_cpg_interval) |
| 293 | "together, p 0.025)" | nuclei_paired/paired.csv (same row as line 24) | donor bootstrap | VALID | keep |
| 365 | "(p 0.006)" Tet3 F cortex | mmc_interventions/effects.csv (same_chip, CerebralCortex, housekeeping, F +0.168, p 0.0056) | consortium_interventions.py: OLS across animals of per-animal F, chip covariate | VALID | keep |
| 365 | "(p 9·10⁻⁶)" Tet3 F striatum | same (Striatum F +0.245, p 8.6e-6) | OLS across animals | VALID | keep |
| 367 | "B +0.009 in cortex, p 0.0008" | same (CerebralCortex B +0.0093, p 0.00077) | OLS across animals | VALID | keep |
| 367 | "+0.005 in striatum, p 2·10⁻⁵" | same (Striatum B +0.0051, p 2.2e-5) | OLS across animals | VALID | keep |
| 372 | "(F −0.11 to −0.17, p ≤ 0.014)" | tet1_cortex/effects.csv (F −0.105 to −0.170) | tet1_cortex.py: arm_scores + OLS across animals | VALID. The current table's largest F p is 3·10⁻⁵, so "p ≤ 0.014" is loose. 0.014 matches the Hxd core B row with 3 components. | keep (could tighten to p ≤ 3·10⁻⁵) |
| 372 | "(S −0.021 and −0.018, p < 10⁻⁷)" | tet1_cortex/effects.csv (S, max p 3.6e-8) | OLS across animals | VALID | keep |
| 379 | "(F −0.09 to −0.17, p ≤ 0.002" | tet1_cortex_rsesame/effects.csv (F −0.093 to −0.167, max p 2.9e-4) | OLS across animals | VALID | keep |
| 382 | "(p ≤ 3·10⁻⁵)" Tet1 shared genes | tet_shared_rsesame/effects.csv (shared genes, housekeeping F: KO −0.190 p 5.7e-8, Hxd −0.193 p 3.3e-5) | tet_shared.py: OLS across animals | VALID | keep |
| 382 | "(p ≤ 0.009)" Tet3 shared genes | same (Tet3 cortex +0.169 p 0.0089, striatum +0.246 p 5e-6) | OLS across animals | VALID | keep |
| 412 | "(−0.007 SD per unit, −0.021 to +0.008, p 0.33)" | essentiality_gradient/gradient.csv (blood, kl_slope, + expression and breadth: −0.0074, −0.0218 to +0.0071, p_boot 0.33) | essentiality_gradient.py: regression over CpGs with samples fixed; interval and p from resampling genes | SET-LABEL (genes as units) | `(−0.007 SD per unit, −0.021 to +0.008, p 0.33)` → `(−0.007 SD per unit, −0.021 to +0.008 over genes; no donor resampling interval computed)`. The cortex "+0.022, +0.009 to +0.040" two lines later is also a gene interval: add "over genes". |
| 443 | "correlation −0.12, p 5·10⁻³⁷)" | mesa_expression/gradient.csv (CD14+, kl_slope, spearman −0.117, p_spearman 4.6e-37) | mesa_expression.gradient: Spearman over 11,678 genes of per-gene mean relative slope, samples fixed | SET-LABEL (genes as units) | `correlation −0.12, p 5·10⁻³⁷)` → `correlation −0.12 over genes; no resampling interval computed)` |
| 449 | "(rank correlation −0.06, p 1·10⁻⁹)" | sturm_expression_rsesame/gradient.csv (kl_slope, spearman −0.059, p_spearman 1.2e-9) | sturm_expression.py → mesa_expression.gradient (same) | SET-LABEL (genes as units) | `(rank correlation −0.06, p 1·10⁻⁹)` → `(rank correlation −0.06 over genes; no resampling interval computed)` |
| 510 | "|F| ≤ 0.022, p ≥ 0.15" | chip_methylation/effects.csv | chip_methylation.py: OLS across patients (carriers vs controls) | VALID | keep |
| 510 | "|F| ≤ 0.031, p ≥ 0.19" | chip_methylation_rsesame/effects.csv | OLS across patients | VALID | keep |
| 517 | "−0.09, p 0.16" | cr_blood_replication/effects.csv (RUN_LOG R98) | within-chip randomisation of treatment with each animal's reference refitted (720 exact assignments) | VALID | keep |
| 517 | "−0.08, p 0.36" | cr_blood_replication_rsesame/effects.csv (core F −0.078, p 0.36) | same | VALID | keep |
| 525 | "+0.004 per log-unit, 95 % interval −0.036 to +0.044, p 0.84" | species_floor/lifespan_fits.csv (off floor, kl_slope, housekeeping, +0.0042, SE 0.020, p 0.834) | OLS across groups | VALID | keep |
| 526 | "(+1.43, p 0.23)" share of drift | RUN_LOG R68 (line 1273) only. The current species_shield/lifespan_fits.csv (kl_slope, housekeeping, relative gap) gives +0.59, p 0.29 (131 groups); species_shield_chip gives +0.005, p 0.96. | OLS across groups | VALID (method). The number does not match the current table. | keep the method; update the number from the table (+0.59, p 0.29) after checking |
| 531 | "+0.045, p 0.0002" | gill_reprogramming/effects.csv (B_housekeeping +0.0445, p_perm 0.0002) | gill_reprogramming.donor_permutation_p: labels permuted within donors | VALID | keep |
| 531 | "(F −0.016, p 0.25" | same (F_housekeeping −0.0158, p_perm 0.25) | within-donor label permutation | VALID | keep |
| 531 | "core-essential −0.015, p 0.35" | same (F_core_essential −0.0153, p_perm 0.347) | within-donor label permutation | VALID | keep |
| 562 | "1.19-1.71, p 10⁻⁴" | incident_cancer/effects.csv | incident_cancer.logit: logistic regression across people, age and sex held | VALID | keep |
| 564 | "1.17-1.76, p 6·10⁻⁴" | incident_cancer/effects.csv (chip, mean methylation and cell fractions held) | logistic regression across people | VALID | keep |

## Replacements (exact strings)

1. Line 142
   - old: `(−0.015 SD, p 0.005; Supplementary Fig. S10)`
   - new: `(−0.015 SD, −0.024 to −0.006 over genes; no donor resampling interval computed; Supplementary Fig. S10)`
2. Line 156
   - old: `(−0.060, −0.029; all p ≤ 0.002; Supplementary Fig. S11)`
   - new: `(−0.060, −0.029; no resampling interval computed; Supplementary Fig. S11)`
3. Line 250 (headroom)
   - old: `−0.42; cortex −0.03, −0.11; p ≤ 0.02)`
   - new: `−0.42; cortex −0.03, −0.11; 95 % donor intervals blood −0.31 to −0.12 and −0.47 to −0.17, cortex −0.11 to −0.02 and −0.16 to −0.07)`
4. Line 250 (per-CpG loss rate)
   - old: `(−0.02 and −0.02; p 0.24 and 0.02)`
   - new: `(−0.02 and −0.02; 95 % donor intervals −0.05 to +0.05 and −0.03 to +0.04)`
   - Optional: the blood loss rates in the same sentence have donor intervals too, so `(+0.07 core-essential, +0.11 housekeeping SD)` → `(+0.07 core-essential, +0.11 housekeeping SD; +0.02 to +0.10 and +0.05 to +0.15)`. The resting levels (blood −0.24, −0.33; cortex −0.08, −0.15) have intervals −0.26 to −0.13, −0.36 to −0.17, −0.11 to −0.05 and −0.17 to −0.10.
5. Line 412
   - old: `(−0.007 SD per unit, −0.021 to +0.008, p 0.33)`
   - new: `(−0.007 SD per unit, −0.021 to +0.008 over genes; no donor resampling interval computed)`
6. Line 443
   - old: `correlation −0.12, p 5·10⁻³⁷)`
   - new: `correlation −0.12 over genes; no resampling interval computed)`
7. Line 449
   - old: `(rank correlation −0.06, p 1·10⁻⁹)`
   - new: `(rank correlation −0.06 over genes; no resampling interval computed)`

## Side findings (not p-value verdicts)

- Line 142: the paper says p 0.005; the table says 0.006. This does not matter once the p is replaced.
- Line 372: "p ≤ 0.014" is looser than the current tet1_cortex table, where the largest F p is 3·10⁻⁵.
- Line 526: "+1.43, p 0.23" is not in any current species table. RUN_LOG R68 records it, but the regenerated
  species_shield/lifespan_fits.csv (dated Nov 23 10:42) gives +0.59, p 0.29.
- Line 234: the "0.05 with groups clustered by species" fit has no table or code I could find.
- Several unchanged sentences quote gene-resampling intervals without saying so:
  - line 412, the cortex value "+0.009 to +0.040";
  - line 441, the monocyte "75 to 92 %";
  - mesa/sturm matched_gaps, which resample genes.
  They may want "over genes".
- Methods line 814 (outside the audited range) states that the within-person gaps' p-values take genes as units. Any
  Results sentence that quotes them falls under the same rule. None was found in lines 1-680.

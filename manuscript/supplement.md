# Supplementary information

Companion to `paper.md`. Sizes are set minus matched CpGs (same CpG-island status and young-level decile), in SD of
the drift measure, unless stated otherwise. Run identifiers (R-numbers) refer to the project's run log
(`docs/RUN_LOG.md`); output paths are relative to the root of the repository at https://github.com/longevity-genie/entropy-age.

## Supplementary Note 1. The metric audit: binary entropy against relative entropy

Binary entropy H(m) measures how far a CpG is from 50 %. At a low-methylated site, a loss of methylation therefore
counts as a gain in order. Essential promoters lose methylation faster than their look-alikes (main text, gain and loss
rates), so every entropy gap mixed how far a site moved with which way it moved. The corrected measure is relative
entropy from the young state, KL(x ‖ m_young) in bits; its age slope counts moves up and down alike. It was added beside
the old measures, and every result whose headline rested on entropy was rerun with all four measures on the same seeds
(the old numbers reproduce). Per-person scores (P, B, F, S), clocks, survival, accessibility, and results on |mean
change| or dispersion already measure distance from a reference and were not affected.

Table S2. Results that rested on binary entropy, old against corrected (housekeeping unless stated).

| result | entropy (old) | relative entropy (corrected) | verdict |
|---|---|---|---|
| blood (Hannum), housekeeping / core-essential | −0.35 / −0.23 | −0.15 / −0.09 | kept, about 2.3× smaller |
| cortex (Jaffe), housekeeping / core-essential | −0.21 / −0.13 | −0.06 / −0.07 | kept, 1.7-3.6× smaller |
| non-essential genes drift more | blood +0.31, cortex +0.08 | blood +0.13, cortex −0.00 (p 0.94) | blood kept; cortex overturned |
| Jaffe deposited → reprocessed | −0.22 → −0.15 | −0.07 → −0.09 | steadier across preprocessing |
| chromatin share removed, blood | 77-84 % | 83-91 % | kept |
| chromatin share removed, cortex | 43-69 % | 25 % own state, 56 % with sorted ATAC | resized |
| essentiality beyond chromatin, blood | +0.022 (0.008-0.033) | −0.007 (−0.021, 0.008), p 0.33 | overturned |
| essentiality beyond chromatin, cortex | +0.018, p 0.08 | +0.022 (0.009-0.040) | kept |
| sorted CD4+ T cells | −0.22 | −0.015 (p 0.005) | resized; monocytes −0.10 |
| GTEx whole blood / lung / muscle | −0.29 / −0.27 / −0.07 | −0.14 / −0.09 / −0.00 | kept, about halved |
| mouse array, renewing tissues | −0.12 to −0.31 | −0.05 to −0.19 | kept |
| mouse array, post-mitotic and adipose | −0.22 to +0.10 | +0.03 to +0.09, all p ≤ 0.002 | reframed (hypothalamus −0.08) |
| species, share of 139 groups negative | 96 % (−0.30) | 89 % (−0.13) | kept, resized |
| species, with chip intercepts (122 groups that can be fitted) | n.a. | 86 % (−0.11) | kept |
| species off the floor (10-90 %) | 93 % (−0.24) | 73 % (−0.10) | kept, resized |
| relative shield grows with lifespan | −2.42, p 0.04 | +1.43, p 0.23 | overturned |
| no lifespan relation off the floor | −0.014, p 0.48 | +0.004, p 0.84 | kept |
| mouse RRBS and array series (GSE80672, BXD, six RRBS tissues, Meer) | −0.09 to −0.52 | −0.03 to −0.32 | all keep sign |

Source: `docs/METRIC_AUDIT.md`; runs R68, R73.

## Supplementary Note 2. Statements withdrawn after the audit

Three earlier statements are not made in the paper: (1) essentiality protects beyond chromatin in blood (R58; the
remainder is entropy-only, and on relative entropy it is cortex's); (2) the shield's share of drift grows with species
lifespan (R22; +1.43, p 0.23 on relative entropy); (3) a large shield in sorted CD4+ T cells (R44; −0.015 SD on
relative entropy). Read-level epiallele entropy (R50-R53) is not yet recomputed; the paper quotes its partner measure,
the share of discordant reads, instead.

## Supplementary Note 3. Further results

**Clocks.** Retraining clocks on simulated change (after Tong et al. 2024) in Hannum blood kept, as a share of the
original R² with age: published stochastic clocks 0.67 / 0.65 / 0.89 (Horvath / PhenoAge / Zhang 2019); each CpG's own
direction at a common rate 0.71 / 0.70 / 0.80; undirected noise 0.02 / 0.10 / 0.00; noise damped at shielded promoters
0.07 / 0.03 / 0.00 (R38). Clocks trained on equal-size pools of shielded, look-alike and random CpGs read age alike (blood
r 0.884 / 0.882 / 0.885; cortex 0.929 / 0.943 / 0.938), while the shielded clock's acceleration overlapped less with the
others (blood r 0.63 against 0.80; R33). At equal look-alike drift, faster chronological and stochastic clocks went with a
stronger shield (partial r −0.15 to −0.29), and the outcome clocks were unrelated to it (PhenoAge −0.07, DunedinPACE 0.00,
Zhang mortality −0.04; R29). In 11 of 12 published clocks, protected CpGs carry the weight of any CpG with the same slope
and starting level; DamAge (Ying et al. 2024) alone weights them up (1.66, 1.14-2.36; R64). Eleven blood-cell fractions explain 15.5 % of
PhenoAge's acceleration and 4.9 % of StocP's (R38).

**Chromatin writers and DNMT3A gain of function.** Blood of carriers against controls (5 composition components;
R18; on arrays rebuilt from raw IDATs, R95, every drift effect within 0.002 and the F effects 0.01-0.02 stronger): Sotos drift B +0.021 (p 10⁻⁴²), F −0.043; other NSD1 variants F −0.039; Weaver F −0.037; Kabuki B +0.011, F
+0.003 (n.s.); CHARGE F −0.012; KMT2D and CHD7 variant carriers ≈ 0. GSE74432 has no ages. Dnmt3a W326R/+ against
littermates (R37): F −0.53 to −0.61 for essential and housekeeping sets in bone marrow, intestine and liver, +0.56 to
+0.67 for non-essential; S −0.06 to −0.09 and +0.07 to +0.09; B +0.06 to +0.10.

**Tet3 and Tet1 detail.** Tet3 heterozygous against same-chip, same-age wild types (litters not recorded in GEO; R15, R16, R36): housekeeping F +0.17
(cortex, p 0.006), +0.25 (striatum, p 9·10⁻⁶); core-essential F +0.31, +0.49; core-essential S +0.012 (p 3·10⁻⁴),
+0.014 (p 2·10⁻⁴), with compression held +0.010 (p 0.003), +0.013 (p 0.001). Signed shifts also occur in arms without a
shield change (GHRKO kidney +0.023), so S alone is not specific. Tet1 knockout / catalytic-dead (R35): housekeeping B
+0.033 / +0.029, F −0.130 / −0.125; core-essential F −0.170; S −0.018. The global correction (R71) dropped some weaker
Tet3 and GHRKO contrasts (p 0.027-0.036); the strongest Tet3 contrasts survive.

**Mouse interventions and selection tests.** CR in B6D2F1 blood (R8; 12 against 10 mice): housekeeping F +0.025 (p 0.002, 5 components),
with blood-cell openness matched +0.008 (p 0.07; R55). CR in C57BL/6 crossed sequencing series and is not read. CR
liver (59 against 36 mice): F −0.008 (p 0.77; R13). Second CR series, young mice, 12 weeks (R56): controls at 12 weeks against the start B −0.003, but all eight start controls sat on one chip and the end controls on two others, so this aging contrast is confounded with chip and unresolved; CR
against controls, core-essential F −0.090, housekeeping −0.027; under a within-chip randomisation that refits the
controls' reference per assignment (720 assignments, R91) p 0.16 and 0.23 (the earlier label shuffle on fixed
scores gave 0.006 and 0.026). GHRKO liver strengthening was batch (same-chip +0.01, p 0.67). Clonal haematopoiesis (R25): 64
patients, |F| ≤ 0.022, p ≥ 0.15. Transplant (R57): F core-essential −0.058 (−0.091 to −0.028), housekeeping −0.039,
non-essential −0.057. Per-mouse shield share grows with age: GSE80672 standard C57BL/6 F −0.052 per 10 months
(−0.041 with components); BXD liver −0.28 (−0.055 with components) (R8, R9); re-processed from the IDATs with SeSAMe's mouse masks (R106) −0.25 (−0.109 with components),
and the high-fat diet −0.086 (−0.035 with components, p 0.07; deposited −0.017, p 0.26). Chemotherapy in SJLIFE moved B and F in
opposite directions (alkylating agents B +0.15, F −0.16; anthracyclines B −0.22, F +0.21; R60).

**Reprogramming.** Gill 2022 from IDATs (R65): reprogrammed − failed B housekeeping +0.045, core +0.036,
non-essential +0.125 (p 0.0002 each); F −0.016 (p 0.25), −0.015 (p 0.35), non-essential −0.060 (p 0.014); Horvath
skin & blood (Horvath et al. 2018) −0.36 transformed units (p 0.085). Cyclic OSKM in mouse tissues (Browder et al. 2022) changed the shield
inconsistently, and one cycle in one-year-old mice (Chondronasiou et al. 2022) left too little aging to test.
Browder cyclic OSKM (R32): muscle F −0.115 (p 0.02), spleen core-essential +0.211 (p 0.001), other tissues n.s.;
re-processed from the IDATs with SeSAMe's mouse masks (R106) muscle −0.107 (p 0.07), spleen core-essential +0.171 (p 0.007). One
OSKM cycle (R62): no measurable aging between 13 and 52 weeks.

**Species.** Phylogenetic GLS (taxonomy covariance, Pagel's λ) on entropy: housekeeping +0.025 per log lifespan
(p 0.40). Look-alike drift rate falls with lifespan (exponent −1.42 on entropy, −1.15 on dispersion; R22). Human
marrow and lymph node reversals are within donor noise (R54). Drift onset (R45): noise
also builds from early life rather than switching on at maturity: in 68 species-tissue groups, look-alike noise rose
with age in 84 %, with a late change point in only 6 (at 2.4-18 times the age of maturity), and housekeeping CpGs
gained noise more slowly in the same animals (median rank correlation with age 0.14 against 0.28).

**Directed loss outside promoters.** Solo-WCGW open-sea CpGs lose methylation 2.7× faster than other open-sea CpGs in
blood (−0.21 against −0.08 percentage points per decade, p 4·10⁻⁵), none in cortex (R14); not part of the shield claim.

**Accessibility (R23, R24, R27, ATAC variance).** Blood bulk ATAC (157 donors): housekeeping −0.048 SD (p 0.02) with
components, −0.031 raw (p 0.14). Cortical snATAC: core-essential +0.07 to +0.23 SD in 6 of 7 cell types. Hippocampus:
signs depend on adjustment.

**RNA.** Housekeeping minus matched expression spread gain (mean over class): renewing GTEx −0.08, mouse −0.29, rat
−0.10; brain +0.06, +0.01, −0.04. Tabula Muris Senis (R75): 15 of 17 organ intervals include 0; renewing − brain +0.20
(p 0.21). Methylation-expression spread coupling (R43): mean |r| above null by 0.015-0.022.

**Probe design and preprocessing (R20, R61).** The two probe designs of the 450k array differ in range, and
essential-promoter CpGs are type I more often (42-45 % against 35 %). Adding design to the matching strata reduced
blood's gap by about a sixth, moved cortex's by under 0.01 SD and changed no sign. Masking unreliable probes (Zhou et
al. 2017) moved no essential-set gap by more than 0.012 SD and the non-essential gaps by up to 0.03 (Supplementary Fig.
S13). We rebuilt Jaffe's cortex series from its 675 raw arrays with a SeSAMe-like pipeline (R61), and the 213 adult controls
again with R SeSAMe (R93): on 208 identical donors and 89,446 complete CpGs the two raw routes agree within 0.004 SD,
and both strengthen the gap over the deposited betas (housekeeping relative entropy −0.064 → −0.081, spread −0.072
→ −0.119 SD; paired difference on spread −0.047, −0.081 to −0.003). In the first rebuild the housekeeping
relative-entropy gap moved from −0.07 (deposited) to −0.09 SD (reprocessed), although single CpGs' slopes agreed only
at r 0.56-0.82 between the two versions (Supplementary Fig. S12). Matching each set to look-alikes cancels such shifts.
In this one rebuild the housekeeping gap moved by 0.07 SD on binary entropy and 0.02 on relative entropy; judging from
it, a gap's size carries roughly ±0.05 SD of preprocessing uncertainty.

**Mouse array tissues (GSE279437; R68a).** Housekeeping gap on relative entropy: white blood cells −0.08 SD, colon
−0.18, spleen −0.19, kidney −0.18 and liver −0.05; skeletal muscle +0.09, heart +0.08, frontal cortex +0.07,
hippocampus +0.03 and adipose tissue +0.03; the hypothalamus, with four aged controls, −0.08. Re-processed from the raw
IDATs (R97) the deposited betas proved to be R SeSAMe's own output, identical to float precision; three young
controls fail detection there; relative-entropy contrasts move by at most 0.02 (entropy slopes by up to 0.04). Within chip (R92; animals
on chips shared by both ages, chip intercepts; identifiable in 7 tissues): white blood cells +0.04, liver +0.01,
skeletal muscle −0.04, heart +0.09, frontal cortex −0.01, hippocampus −0.05, hypothalamus −0.14; colon, spleen, kidney
and adipose tissue drop out. Renewing minus post-mitotic +0.02 (−0.06 to +0.19), slowly renewing minus post-mitotic
+0.03 (−0.04 to +0.18). Within-chip dispersion slopes are not reported: in the hypothalamus each shared chip holds one
aged and one young animal, so with chip intercepts its dispersion slope is zero by construction.

**Health states and outcomes.** Cross-sectional (R39): 19 of 522 effects pass FDR; schizophrenia cortex F +0.57 SD
(5 components); smoking lung S −0.67 raw; Alzheimer's none. Vitality 90+ (R42): positive-control clocks fail, so
uninformative. MESA monocytes (R44): plaque and calcium track DunedinPACE (+0.37 per SD), not the shield. Per-person
scores against confounders (R60): slide explains 36 % of B, 23 % of F, 20 % of S (13 % by chance); smoking moves S only.

**Within one person, deposited betas (R52).** On the 450k draws alone (R119; 647 pairs of 447 people, no array
change inside any pair) the shield share F still falls within a person, housekeeping −0.126 per decade (−0.226 to
−0.030), while the divergence gap below does not hold (−0.003, p 0.41). Pooled over both arrays: divergence per year housekeeping −0.017 (p 0.001), core-essential
+0.013 (p 0.10); within-person F slope −0.14 per decade (−0.23 to −0.05); B +0.03 (spans 0). Positive control on first
draws kept a third of Hannum's dispersion gap (−0.12) and none of its entropy gap.

## Supplementary Note 4. Array processing and quality-control practice

The practice the analyses follow, item by item, with where it stands (R numbers are runs; RM numbers are roadmap
items still open).

- **Raw data, one recipe.** Every array series with public raw IDATs is rebuilt with R SeSAMe 1.24.0 (sesameData
  1.24.0, R 4.4.3), following SeSAMe's recommendation per platform: QCDPB for human 450k, EPIC and EPIC v2; TQCDPB for
  the mouse MM285 design (T masks probes hit by the inferred strain's variants); SHCDPB for the 40k mammalian array in
  mouse (S the species mask, H cg/ch probes only). The steps applied are recorded for each array. The mouse series were
  first rebuilt with the human recipe and are being rebuilt again (RM43). SATSA and Gill were rebuilt with methylprep
  1.7.1, labelled as such; on Jaffe's cortex the two implementations agreed within 0.004 SD (R93). Cohorts without
  raw IDATs keep the authors' processing and are marked wherever reported.
- **Sample flow.** Per series: deposited samples, samples with both IDAT channels, samples passing QC, samples in each
  model, and donors (Supplementary Table S4 names the analyses; per-series counts are in each run's log).
- **Sample QC.** Detection by pOOBAH: an array passes when it detects at least 90 % of its probes (for the
  multi-species 40k array, of the probes the species mask keeps). Recorded beside it: the bisulfite-conversion index
  (GCT, near 1 when conversion is complete), mean intensities, dye-bias statistics and each sample's beta quantile
  curve with its distance from the series' median curve (R103). Three human series have a high conversion index
  (clonal haematopoiesis 1.58, the two writer-syndrome series 1.39-1.47, against about 1.07 elsewhere); their results
  are rerun with high-index arrays excluded (RM43).
- **Probe QC and one CpG universe.** The CpGs analysed are those measured in the samples of a contrast, fixed before
  any group comparison, and set and background CpGs come from that same universe. CpGs missing in more than a small
  share of samples (5 % for most array analyses; 10 % for RRBS, the two mouse blood RRBS analyses and Gill) are dropped; the remaining
  sporadic missing values are filled with that CpG's median across the contrast's samples, alike for set and
  background. Rerun on complete CpGs only (R117), every headline gap agreed to the third decimal; only the rebuilt
  cortex arrays had missing values (4,577 of 94,513 CpGs). No clock training means or values pooled across studies are imputed.
- **Batch.** Within each series, the share of the leading principal components explained by chip, row and position,
  and how far each design variable (age, group, tissue, sex) is confounded with chip (R103). Where chip and the
  contrast overlap, the contrast is estimated within chip (Huang mouse tissues, R92; CR blood, R91; EPIC-Italy with
  same-chip matching, R84b; Tet3 against same-chip wild types). Where chip fully determines the contrast (CR blood's
  aging controls; Lunnon blood against brain) it is reported as unidentifiable; batch correction does not create the
  missing comparison.
- **Platforms and studies.** Effects are estimated within one study and platform and reported side by side, with two
  pooled exceptions named in Methods (SATSA, the transplant series). Comparisons across arrays are to be redone on
  the CpGs both measure (RM41: Tet3 against Tet1, GTEx EPIC against 450k cohorts); RRBS and array estimates are not
  ranked against each other.
- **Second tracks.** For each analysis Supplementary Table S4 says whether chromatin, expression or outcome comes
  from the same samples, the same study or another study, and from the same cell or tissue type.
- **Published clocks.** Scoring with each clock's own inputs, reference and missing-data rule, and an audit of how
  much coefficient weight is lost per person, are not yet done (RM41); results that use published clocks are
  reported as indices on our matrices, not as the clocks' published outputs.

## Supplementary figures

| figure | file | content |
|---|---|---|
| S1 | `docs/figures/global_fdr.png` | one BH correction over 7,626 tests (R71) |
| S2 | `docs/figures/drift_hotspots.png` | region classes of the fastest-rising spread (R72) |
| S3 | `docs/figures/mmc_interventions.png` | consortium interventions incl. Tet3 (R13, R15, R16) |
| S4 | `docs/figures/tet1_cortex.png` | Tet1 knockout and catalytic-dead cortex (R35) |
| S5 | `docs/figures/composition_gaps.png` | 11 blood-cell fractions held (R70) |
| S6 | `docs/figures/mouse_tissues.png` | 11 mouse tissues, turnover classes (R68a) |
| S7 | `docs/figures/species_floor.png` | species shield off the 0 % floor (R59, R68c) |
| S8 | `docs/figures/gtex_tissues.png` | GTEx nine tissues (R68d) |
| S9 | `docs/figures/meer_tissues.png` | four mouse tissues by RRBS (R63, R73) |
| S10 | `docs/figures/sorted_cells.png` | sorted monocytes and CD4+ T cells (R44, R68b) |
| S11 | `docs/figures/read_level.png` | read-level WGBS, WBC, cfDNA, sorted cell types (R50-R53) |
| S12 | `docs/figures/reprocessed_jaffe.png` | Jaffe deposited vs rebuilt from IDATs (R61) |
| S13 | `docs/figures/probe_robustness.png` | probe design and probe mask (R20) |
| S14 | `docs/figures/passaged_cells.png`, `passaged_cells_sturm.png` | culture over doublings (R1 and earlier) |
| S15 | `docs/figures/essentiality_gradient.png` | essentiality beyond chromatin (R58, R68b) |
| S16 | `docs/figures/gill_reprogramming.png` | transient reprogramming (R65) |
| S17 | `docs/figures/hsct_bottleneck.png`, `chip_methylation.png`, `cr_chromatin.png` | selection tests (R57, R25, R55) |
| S18 | `docs/figures/stochastic_component.png`, `own_clocks.png` | what clocks read (R38, R33) |
| S19 | `docs/figures/rna_tms.png`, `rna_gtex.png`, `rna_rodents.png` | RNA tests (R75 and earlier) |
| S20 | `docs/figures/sjlife_deaths.png` | SJLIFE deaths (R48) |
| S21 | `docs/figures/syndromes.png`, `dnmt3a_gain.png` | chromatin writers, DNMT3A W326R (R18, R37) |

## Supplementary Table S1. Datasets

| dataset | accession | platform | values |
|---|---|---|---|
| Hannum blood, 656 people | GSE40279 | 450k | deposited; **no raw IDATs on GEO, not rebuilt** |
| Jaffe prefrontal cortex | GSE74193 | 450k | deposited; rebuilt from IDATs (R61 methylprep, R93 R SeSAMe) |
| Lunnon blood and four brain regions | GSE59685 | 450k | deposited; **no raw IDATs on GEO, not rebuilt** |
| GSE41037 blood | GSE41037 | 27k | deposited; **no raw IDATs on GEO, not rebuilt** |
| GTEx nine tissues | GSE213478 | EPIC | deposited (noob + BMIQ); **raw IDATs controlled access, not rebuilt** |
| Muscle fibres | GSE304045 | EPIC v2 | deposited; rebuilt from IDATs (R101) |
| SATSA twins, longitudinal | E-MTAB-7309 / S-BSST1206 | 450k | deposited (R52); rebuilt from IDATs (R66) |
| Gill 2022 transient reprogramming | GSE165179 | EPIC | rebuilt from IDATs (R65) |
| Cellular Lifespan Study fibroblasts | GSE179847 | EPIC | deposited; rebuilt from IDATs (R94) |
| Endicott primary cells | GSE197512 | EPIC | deposited; rebuilt from IDATs (R94) |
| EPIC-Italy incident cancer | GSE51032 | 450k | deposited; rebuilt from IDATs (R96) |
| SJLIFE survivors | GSE314261 | EPIC | deposited |
| MESA sorted monocytes, CD4+ T cells | GSE56046, GSE56581 | 450k | deposited |
| Vitality 90+ survival | GSE68194 | 450k | deposited |
| Chromatin-writer syndromes | GSE74432, GSE97362 | 450k | deposited; rebuilt from IDATs (R95) |
| Stem-cell transplant series | Mendeley Data, https://doi.org/10.17632/j5krtjfj6y.1 | 450k/EPIC | deposited |
| Loyfer 2023 WGBS | GSE186458 | WGBS | authors' read files, our per-read summaries |
| Mammalian Methylation Consortium (65 species; interventions incl. Tet3) | GSE223748 | Mammal40 | deposited |
| Huang 11 mouse tissues | GSE279437 | Mammal320 | rebuilt (R97; deposited = R SeSAMe output) |
| Petkovich mouse blood | GSE80672 | RRBS | deposited |
| CR blood replication | GSE282661 | Mammal320 | deposited; rebuilt from IDATs (R98) |
| Mouse RRBS six tissues | GSE120132 | RRBS | deposited |
| Meer 2018 four tissues | GSE121141 | RRBS | deposited |
| BXD liver | GSE199979 | mammalian methylation array (40k) | rebuilt (R100) |
| Browder 2022 cyclic OSKM | GSE190665 | Mammal40 | deposited |
| Chondronasiou 2022 one OSKM cycle | GSE156557 | RRBS | deposited |
| Tet1 cortex | GSE271375 | MM285 | rebuilt (R99; 19 of 24 arrays pass) |
| Dnmt3a gain of function | GSE324160 | mouse array | rebuilt (R99) |
| Clonal haematopoiesis | GSE315366 | EPIC | rebuilt (R99) |
| ImmGen ATAC (mouse blood cells) | GSE100738 | ATAC-seq | deposited peaks |
| Sorted neuron/glia chromatin | GSE96614 | ATAC-seq | deposited |
| Reference epigenomes | Roadmap E062, E030, E073 | ChromHMM | deposited |
| Blood ATAC (Morandini 2024) | GSE193140 | ATAC-seq | deposited |
| Cortical snATAC (Catching 2026) | Zenodo, https://doi.org/10.5281/zenodo.18394349 | snATAC | deposited |
| Hippocampal snATAC | GSE278576 | snATAC | deposited |
| Rodent RNA atlases | GSE314817, GSE314970 | RNA-seq | deposited counts |
| Tabula Muris Senis bulk | GSE132040 | RNA-seq | deposited |
| GTEx RNA | GTEx v10 | RNA-seq | deposited |
| DepMap CRISPR gene effect | DepMap 24Q4 Public, https://doi.org/10.25452/figshare.plus.27993248.v1 | CRISPR screens | deposited |
| GTEx expression breadth | GTEx v8 gene median TPM | RNA-seq | deposited |

Source: `docs/PROVENANCE.md`, `docs/extension.md`.

## Supplementary Table S3. Traceability: claim → run → output table

| claim in paper.md | claim id | run(s) | output table(s) |
|---|---|---|---|
| Blood/cortex gaps on relative entropy (−0.15, −0.09; −0.06, −0.07; non-essential +0.13, −0.00) | C1 | R68a | `data/output/tissue_contrast/point_estimates.csv` |
| Blood minus cortex, relative entropy (−0.09; −0.02) | C1 | R68a | `data/output/tissue_contrast/tissue_differences.csv` |
| Same-donor dispersion gaps; per-person 86-91 % | C1 | R3/tissues, R2 | `data/output/tissue_contrast/`, `data/output/per_person/same_donor_turnover.csv` |
| Placebo sets | C1 | R69 | `data/output/placebo_sets/summary.csv` |
| Cell fractions held | C1 | R70 | `data/output/composition_gaps/gaps.csv` |
| Sorted monocytes −0.10, CD4+ −0.015 | C1 | R44, R68b | `data/output/sorted_cells/shield.csv` |
| Probe design and mask | C1 | R20, R68b | `data/output/probe_robustness/gaps.csv` |
| Jaffe from IDATs | C1 | R61, R68b | `data/output/reprocessed_jaffe/gaps.csv`, `per_cpg_agreement.csv` |
| Read-level discordance | C1 | R50, R51, R53 | `data/output/read_level/` |
| Global FDR, 7,626 tests | all | R71 | `data/output/global_fdr/by_folder.csv` |
| GTEx tissues | C1 | R68d | `data/output/gtex_tissues/point_estimates.csv` |
| Mouse array 11 tissues, renewing − post-mitotic −0.19 | C1, C2 | R68a | `data/output/mouse_tissues/point_estimates.csv`, `contrasts.csv` |
| RRBS and BXD gaps | C1 | R73 | `data/output/meer_tissues/aging_gaps.csv`, `mouse_rrbs_tissues/aging_gaps.csv`, `interventions/aging_gaps.csv`, `bxd_liver/aging_gaps.csv` |
| Culture, share against doublings | C3 | Endicott and Sturm runs (pre-R-numbering), R1 | `data/output/passaged_cells/coupling_across_lines`, `data/output/passaged_cells_sturm/` |
| Species 89 % / 73 % off floor; lifespan +0.004 | C6 | R12, R59, R68c | `data/output/species_floor/species_gaps.parquet`, `data/output/species_floor/lifespan_fits.csv`, `share_off_floor.csv` |
| Phylogenetic GLS | C6 | species_phylo | `data/output/species_phylo/lifespan_fits` |
| Gain and loss rates; resting level | C2, C4 | R67 | `data/output/equilibrium_drift/pooled_rates.csv`, `per_cpg_gaps.csv`, `identified.csv` |
| Faster loss within active TSS | C2 | R74 | `data/output/equilibrium_drift/pooled_rates.csv`, `per_cpg_gaps.csv` |
| Tet3 F and S | C4 | R13, R15, R16, R36 | `data/output/mmc_interventions/effects.csv` |
| Tet1 | C4 | R35 | `data/output/tet1_cortex/effects.csv` |
| DNMT3A W326R | C4 | R37 | `data/output/dnmt3a_gain/` |
| Chromatin-writer syndromes | C5 | R18 | `data/output/syndromes/effects` |
| Chromatin share 83-91 % blood, 25 %/56 % cortex | C2 | R3, R4, R68b | `data/output/chromatin_match/matched_gaps.csv` |
| Essentiality beyond chromatin | C2 | R58, R68b | `data/output/essentiality_gradient/gradient.csv` |
| Drift hotspots, Polycomb 2.1× | C1 | R72 | `data/output/drift_hotspots/enrichment.csv` |
| SATSA from IDATs | C1 | R66 | `data/output/satsa_idats/gaps.csv`, `cross_sectional.csv`, `person_slopes.csv` |
| SATSA deposited, within-person F slope | C1 | R52 | `data/output/within_person/` |
| CHIP null | C3 | R25 | `data/output/chip_methylation/effects` |
| Transplant bottleneck | C3 | R57 | `data/output/hsct_bottleneck/effects` |
| CR blood, openness-matched, liver, second series | C3 | R8, R55, R13, R56 | `data/output/interventions/intervention_effects`, `cr_chromatin/intervention_effects`, `mmc_interventions/effects.csv`, `cr_blood_replication/effects` |
| Lifespan scaling of look-alike drift | C6 | R22 | `data/output/species_shield/lifespan_fits.csv` |
| Drift onset | — | R45 | `data/output/drift_onset/` |
| Transient reprogramming | C10 | R65 | `data/output/gill_reprogramming/effects.csv` |
| Cyclic / one-cycle OSKM | C10 | R32, R62 | `data/output/reprogramming/effects`, `data/output/reprogramming_rrbs/effects` |
| RNA (GTEx, rodents, TMS) | — | earlier rna runs, R75 | `data/output/rna_gtex/`, `rna_rodents/`, `rna_tms/contrasts` |
| Methylation-expression spread | — | R43 | `data/output/meth_expr_spread/` |
| Accessibility | C11 | R23, R24, R27, atac_variance | `data/output/atac_variance/gaps`, `atac_cortex/gaps`, `atac_hippocampus/` |
| Clocks read direction | C8, C9 | R38, R33, R29, R64 | `data/output/stochastic_component/`, `own_clocks/`, `shield_clocks/correlations`, `clock_weights/summary` |
| Breast-cancer lead (1.43; 1.62/1.66; adjusted 1.44) | C12 | R46, R46b | `data/output/incident_cancer/effects.csv`, `per_person.csv` |
| B beside DunedinPACE/PhenoAge (1.44), SJLIFE (1.01) | C12 | R48, R48b | `data/output/sjlife_deaths/`, `data/output/incident_cancer/per_person.csv` (joint models fitted ad hoc from the per-person tables; not in a module) |

## Supplementary Table S4. Provenance of every comparison

For each analysis: the methylation data and its platform and processing; any second track (chromatin, expression,
outcome, cell fractions) and its source; whether the comparison stays within one batch or study; whether the second
track comes from the same samples, the same study or another study; and whether it comes from the same cell or tissue
type. Generated from `docs/provenance.csv` by `scripts/provenance_table.py`.

<!-- provenance:start -->

| analysis | runs | claim | tissue / species | methylation (study, accession) | platform | processing | second track | its source | batch | second track from | cell or tissue type of second track | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Blood headline gap | R68a | C1 | human whole blood | Hannum, GSE40279 | 450k | deposited (no IDATs) | none |  | one study | none | n.a. |  |
| Cortex headline gap | R68a R93 | C1 | human prefrontal cortex | Jaffe, GSE74193 | 450k | deposited; rebuilt in the paired pilot (R93) and the full rerun R102 | neuron fraction (covariate) | same arrays' metadata | one study | same samples | same | cross-check: rebuilt on 208 of the 213 donors (R93): gap stronger |
| Blood minus cortex | R68a | C1 | blood vs cortex | Hannum vs Jaffe | 450k | deposited (both); two pipelines | none |  | across studies | none | n.a. | two studies, different pipelines and age ranges; neuron fraction held in cortex only; descriptive |
| Same-donor blood and brain | R68a R2 | C1 | human blood and 4 brain regions | Lunnon, GSE59685 | 450k | deposited (no IDATs) | none |  | one study; tissue confounded with slide and per-tissue normalisation | none | n.a. | blood drawn in life, brain at death |
| Cell fractions held | R70 | C1 | human whole blood | Hannum, GSE40279 | 450k | deposited | 12 cell fractions | EpiDISH reference panel (Luo 2023) | one study | other study | related | reference panel of sorted cells |
| Sorted monocytes and T cells | R44 R68b | C1 | human CD14+ monocytes / CD4+ T cells | MESA, GSE56046 GSE56581 | 450k | deposited | none |  | one study | none | n.a. | text sets them beside Hannum (other study) |
| Read-level blood | R50 | C1 | human white blood cells | Loyfer, GSE186458 | WGBS | authors' reads | Hannum array means (check) | GSE40279 | across studies | other study | same | cross-platform used only as a validation (r 0.963) |
| GTEx nine tissues | R68d | C1 | human 9 tissues | GTEx, GSE213478 | EPIC | deposited (IDATs controlled) | none |  | one study | none | n.a. | same donors across tissues |
| Huang mouse tissues | R68a R92 R97 | C1 C2 | mouse 11 tissues | Huang, GSE279437 | Mammal320 | rebuilt (= deposited) | none |  | one study; age confounded with chip | none | n.a. | within chip (R92) the turnover ordering vanishes |
| Mouse RRBS tissues | R63 R10 R11 R7 R73 | C1 | mouse tissues and blood | Meer GSE121141; Thompson GSE120132; Petkovich GSE80672 | RRBS | deposited | none |  | one study each | none | n.a. | reported per study; not ranked against arrays |
| BXD liver | R9 R100 | C1 | mouse liver | BXD, GSE199979 | Mammal40 | rebuilt (R SeSAMe) | cell-mix components | same arrays | one study | same samples | same | HFD effect partly survives components on rebuild; Mammal40 rebuilt without SeSAMe's species masks; sample QC weak |
| Mouse turnover across platforms | paper text | C1 C2 | mouse tissues | Huang vs Meer/Thompson/Petkovich vs BXD | Mammal320 vs RRBS vs Mammal40 | mixed | none |  | across studies and platforms | none | n.a. | cross-platform: side by side only (RM41 9b) |
| Culture Endicott | R94 | C3 | human primary cell lines | Endicott, GSE197512 | EPIC | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | shield-in-culture vanishes on rebuild |
| Culture Sturm | R1 R94 | C3 | human fibroblasts | Sturm, GSE179847 | EPIC | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | unchanged on rebuild |
| Species shield | R12 R68c R59 R22 R118 | C6 | 65 mammals | MMC, GSE223748 | Mammal40 | deposited (no IDATs) | lifespan and weight; NCBI taxonomy | same metadata | one study; chip partly confounded with age (56 of 139 groups); with chip intercepts 86 % of 122 fittable groups stay negative (R118) | same samples | n.a. | mouse mm10 promoter and island annotation applied to every species |
| Gain and loss rates | R67 R83 | C2 C4 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | none |  | one study each | none | n.a. |  |
| Within active promoters | R74 | C2 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | chromatin state (TssA) | Roadmap E062 (PBMC) / E073 (DLPFC) | one study each | other study | blood: related (PBMC vs whole blood); cortex: same | reference epigenomes from other people |
| Chromatin match blood | R3 R68b | C2 | human whole blood | Hannum, GSE40279 | 450k | deposited | chromatin state | Roadmap E062 PBMC, E030 neutrophils; E073 control | one study | other study | related | E062 lacks granulocytes; E073 a deliberate other-tissue control |
| Chromatin match cortex | R3 R4 R68b | C2 | human prefrontal cortex | Jaffe, GSE74193 | 450k | deposited | chromatin state and sorted ATAC | Roadmap E073 DLPFC; GSE96614 NeuN+/- ATAC BA9; E062 control | one study | other study | same |  |
| Chromatin match Lunnon brain | R3 R68b R104 | C2 | human FC STG EC (cerebellum excluded) | Lunnon, GSE59685 | 450k | deposited | chromatin state | Roadmap E073 (FC), E072 (STG), E071 (EC); BA9 ATAC for FC only | one study | other study | FC same; STG EC nearest region map (related) | R104: region maps remove 5-10 %, blood's map 25-27 % |
| Essentiality gradient | R58 R68b | C2 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | DepMap gene effect; GTEx expression; chromatin state | DepMap 24Q4 (cancer lines); GTEx v8; E062/E073 | one study each | other study | expression same tissue; essentiality other (cancer cell lines); chromatin blood related, cortex same | three external tracks |
| Drift hotspots | R72 | C1 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | chromatin state | E062 / E073 | one study each | other study | blood related; cortex same | Polycomb enrichment rests on E062 labels |
| Expression-matched placebo | R82 | C1 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | expression | GTEx v8 whole blood / frontal cortex BA9 | one study each | other study | same tissue |  |
| Within-person SATSA deposited | R52 R60 | C1 | human blood | SATSA, S-BSST1206 | 450k + EPIC pooled | deposited | none |  | batches pooled with a chip-change term | none | n.a. | pools two platforms within people |
| Within-person SATSA rebuilt | R66 | C1 | human blood | SATSA, E-MTAB-7309 | 450k only | rebuilt (methylprep) | none |  | one platform | none | n.a. |  |
| Tet3 heterozygous | R13 R15 R16 R36 | C4 | mouse cortex / striatum | MMC, GSE223748 | Mammal40 | deposited (no IDATs) | none |  | one study; same-chip wild types | none | n.a. |  |
| Tet1 knockout | R35 R99 | C4 | mouse cortex | Foong, GSE271375 | MM285 | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | 19 of 24 arrays pass |
| Tet3 vs Tet1 | paper text | C4 | mouse brain | GSE223748 vs GSE271375 | Mammal40 vs MM285 | mixed | none |  | across studies and platforms | none | n.a. | cross-platform: stated as different arrays |
| Dnmt3a W326R | R37 R99 | C4 | mouse marrow / intestine / liver | Sarni, GSE324160 | MM285 | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | littermate controls |
| Writer syndromes | R18 R95 | C5 | human blood | GSE74432, GSE97362 | 450k | rebuilt (R SeSAMe) | none |  | one study each vs own controls | none | n.a. | concordance pools the two series' effects; bisulfite conversion index high (median GCT 1.4-1.6 vs ~1.07) |
| Clonal haematopoiesis | R25 R99 | C3 | human blood | Basrai, GSE315366 | EPIC | rebuilt (R SeSAMe) | CH calls | same patients | one study | same samples | same | bisulfite conversion index high (median GCT 1.4-1.6 vs ~1.07) |
| Stem-cell transplant | R57 | C3 | human blood | Holland, Mendeley j5krtjfj6y | 450k + EPIC batches | deposited (ENmix) | none |  | batches pooled with batch terms; pairs within batch | none | n.a. | platform per batch not recorded |
| CR mouse blood RRBS | R7 R8 | C3 | mouse whole blood | Petkovich, GSE80672 | RRBS | deposited | none |  | one study | none | n.a. |  |
| CR openness-matched | R55 | C3 | mouse whole blood | Petkovich, GSE80672 | RRBS | deposited | chromatin accessibility | ImmGen ATAC GSE100738 (6 sorted blood cell types) | one study | other study | related (sorted lineages vs whole blood) | strain of ImmGen mice vs B6D2F1 not recorded |
| CR blood second series | R56 R91 R98 | C3 | mouse whole blood | Jayne, GSE282661 | Mammal320 | rebuilt run pending (R98) | none |  | one study; aging contrast confounded with chip | none | n.a. | R91: CR p 0.006 -> 0.16 under within-chip randomisation |
| Cyclic OSKM | R32 R100 | C10 | mouse 6 tissues | Browder, GSE190665 | Mammal40 | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | Mammal40 rebuilt without SeSAMe's species masks; sample QC weak |
| One OSKM cycle | R62 | C10 | mouse pancreas / spleen / liver | Chondronasiou, GSE156557 | RRBS | deposited | none |  | one study | none | n.a. |  |
| Transient reprogramming | R65 | C10 | human fibroblasts | Gill, GSE165179 | EPIC | rebuilt (methylprep) | epigenetic clock | same arrays | one study | same samples | same |  |
| Breast cancer | R46 R84b R96 | C12 | human blood | EPIC-Italy, GSE51032 | 450k | rebuilt (R SeSAMe) | registry outcomes; cell fractions; clocks | same participants | one study | same samples | same | adjusted OR 1.76 -> 1.11 on rebuild |
| SJLIFE deaths | R48 | C12 | human PBMC | SJLIFE, GSE314261 | EPIC | deposited | deaths; cell fractions | same participants | one study | same samples | same |  |
| Health states | R39 | – | human blood / brain / GTEx | Lunnon, Jaffe, GSE41037, GTEx | 450k / 27k / EPIC | mixed | disease and smoking labels | same samples | one study per contrast | same samples | same | GSE41037 is 27k |
| Methylation-expression spread | R43 | – | human 5 tissues | GTEx, GSE213478 | EPIC | deposited | expression | GTEx v10 RNA (same donors) | one study | same donors | same |  |
| Blood ATAC spread | atac_variance | C11 | human PBMC | none | – | – | bulk ATAC | Morandini, GSE193140 | one study | the measure itself | own tissue |  |
| Cortex snATAC | R23 R24 | C11 | human prefrontal cortex | none | – | – | single-nucleus ATAC | Catching 2026, Zenodo 18394349 | one study | the measure itself | own tissue |  |
| Hippocampus snATAC | R27 | C11 | human hippocampus | none | – | – | single-nucleus ATAC | Zemke, GSE278576 | one study | the measure itself | own tissue |  |
| RNA spread GTEx and rodents | rna_gtex rna_rodents R75 | – | human 47 tissues; mouse; rat | none | – | – | RNA-seq | GTEx v10; GSE314817; GSE314970; GSE132040 | one study each | the measure itself | own tissue |  |
| Joint test in cortex | R108 | C2 | human prefrontal cortex | Jaffe, GSE74193 | 450k | rebuilt (R SeSAMe) | expression; chromatin state; sorted ATAC | GTEx v8 frontal cortex BA9; Roadmap E073; GSE96614 BA9 | one study | other study | same | all three tracks from other people; gap mostly removed |
| Sorted cortical nuclei | R109 R110 | C1 C3 | human cortex: neurons, oligodendrocyte lineage, other glia, unsorted | Walker, GSE279509 | EPIC | rebuilt (R SeSAMe) | none |  | one study; cell types share chips; age mostly between chips | none | n.a. | 524 of 749 arrays pass detection; neurons against glia paired within donors |
| Second cortex series | R109b | C1 | human prefrontal cortex and nucleus accumbens | White, GSE252501 | EPIC | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | same brain bank as Jaffe: shared donors excluded; depression diagnosis held |
| Second sorted series | R115 | C3 | human cortex neurons and glia | Gasparoni, GSE66351 | 450k | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | 31 paired donors of whom 2 under 60: low power |
| Senescence modes in culture | R111 | C3 | human fibroblasts | GSE197723 | EPIC | rebuilt (R SeSAMe) | none |  | one study; every arm on every chip | none | n.a. | donors and doublings not in the record; noise floor from untreated pairs |
| Arrested fibroblasts over time | R112 | C3 | human fibroblasts (confluent) | E-MTAB-15859 | EPIC | rebuilt (R SeSAMe) | none |  | one study; each donor within one chip | none | n.a. | 3 donors with 6 arrays each; month tied to array position in 2 donors; RNA-seq of the same cultures not yet used |
| Colon and hippocampus RRBS | R114 | C1 C3 | mouse colon and hippocampus | GSE233734, GSE304754 | RRBS | deposited per-CpG tables | none |  | one cohort; both tissues on the same flowcells; age balanced | none | n.a. | per-mouse link between tissues unverified |
| Breast case-control | R113 | C12 | human blood | GSE283951 | EPIC | rebuilt (R SeSAMe) | case status | same participants | one study; 9 of 16 chips hold both groups | same samples | same | 32 cases; draw not shown to precede diagnosis |

<!-- provenance:end -->

## Open items for the authors

- Reference metadata marked [verify] in `paper.md` (authors, journals or years not given in the project docs).
- Read-level epiallele entropy and the per-person `K` measure are not yet recomputed on relative entropy.
- Joint models with DunedinPACE/PhenoAge (EPIC-Italy, SJLIFE) and the clock-composition shares (R38 addendum) were run
  ad hoc from per-person tables and are not yet in the pipeline.
- The culture runs (GSE197512) predate the R-numbered run log; their run identifiers are not recorded there.

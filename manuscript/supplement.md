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
| species, 13 human groups, human annotation instead of mouse | n.a. | 9 of 13 negative under both (−0.05 against −0.02) | kept |
| same donors, blood minus brain region (Lunnon), tissues normalised together | n.a. | −0.08 to −0.10 against −0.07 to −0.10 on the deposited values | kept; still confounded with slide |
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

**Expression in the same monocytes (R122).** MESA measured methylation (GSE56046) and expression (GSE56045) in
the same 1,202 monocyte samples. Join check: per-gene correlation of promoter methylation and expression across
people, mean −0.007 (52.5 % of 15,890 genes negative; permutation p 0.005). Housekeeping gap on relative entropy:
−0.099 SD with island x level, −0.017 with the gene's expression decile added, +0.001 with the mononuclear
chromatin state, +0.017 with both. Core-essential: −0.035, +0.020, +0.028, +0.039. Non-essential: +0.218, +0.124,
+0.080, +0.073. Outside the gene sets, the slope relative to matched CpGs by expression decile, lowest to highest:
+0.07, +0.04, +0.03, −0.02, −0.10, −0.10, −0.11, −0.11, −0.12, −0.20 SD. Of housekeeping CpGs 75 % and of
core-essential CpGs 82 % lie in expression deciles 6 to 10, against 27 % of CpGs outside the sets. In CD4+ T cells of
214 of the same people (GSE56581 with GSE56580; R123) the housekeeping gap on relative entropy is −0.015 SD and
the expression decile leaves −0.013; outside the gene sets there is no gradient with expression (rank correlation
+0.02). On the entropy slope the T-cell gap, −0.219, falls to −0.069 with the expression decile, and outside the gene
sets that slope falls with expression (rank correlation −0.20).

**Expression in the same fibroblast cultures (R124).** Sturm et al. 2022: 89 untreated cultures of six healthy
lines at ambient oxygen (GSE179847, rebuilt), slopes on population doublings with a line intercept; expression is
the mean over 77 untreated RNA-seq samples of the same lines (GSE179848). Relative entropy, island x level against
the expression decile added: housekeeping −0.102 and −0.028 SD; core-essential −0.055 and −0.006; non-essential
+0.113 and +0.023. Outside the gene sets, by expression decile: +0.11, +0.09, +0.01, −0.06, −0.08, −0.03, −0.08,
−0.08, −0.09, −0.10. Two lines give 60 of the 89 cultures. With chip intercepts (R125; 50 chips, 27 with two or
more cultures, 0.33 of the variance of doublings left within line and chip): housekeeping −0.102 and −0.024,
core-essential −0.053 and −0.002, non-essential +0.164 and +0.056.

**Expression in the same epidermis (R129).** Holzscheck et al. 2020: suction-blister epidermis of healthy women
aged 21 to 79, on the EPIC array (E-MTAB-8993; 65 women, 61 after the detection cut) or the 450k array (E-MTAB-8992;
56 women on 450k; up to 31 of the EPIC women are women of the 450k group measured again three years later), rebuilt from IDATs, with RNA-seq of every woman (E-MTAB-8995; the depositors' TPM table,
batch-corrected by them, values below zero set to zero). The two array groups are analysed apart. Age slopes hold the
depositors' processing batch. Relative entropy, EPIC then 450k, island x level against the expression decile added:
housekeeping −0.106 and −0.031, −0.105 and −0.031; core-essential −0.058 and −0.003, −0.054 and −0.001;
non-essential +0.218 and +0.088, +0.203 and +0.081. With the chromatin state of keratinocytes (Roadmap E057) in
place of expression: housekeeping −0.052 and −0.055; with both: −0.019 and −0.021. Outside the gene sets, by
expression decile (EPIC): +0.17, +0.06, +0.01, −0.06, −0.08, −0.07, −0.11, −0.09, −0.13, −0.14 (rank correlation
−0.07; 450k −0.09). The depositors renamed the IDAT files; the chip was read from each file's header (9 EPIC chips,
5 450k chips, every chip spanning at least 27 years of age; chip explains 0.24 and 0.18 of the variance of age).
With chip intercepts: housekeeping −0.102 and −0.028 (EPIC), −0.110 and −0.034 (450k). Across women, promoter
methylation and expression of the same gene correlate weakly (mean r −0.003 and −0.028; with chip held −0.022 and
−0.024), so the join is checked by the gradient across genes more than by this correlation. Limits: one
laboratory, women only, cell composition not held.

**One line with and without telomerase (R130).** GSE268145: BJ fibroblasts passaged as mortal cells and with hTERT,
EPIC, rebuilt; slopes on doublings with a trajectory intercept. hTERT cells (20 arrays, doublings 2 to 38):
housekeeping −0.136 SD on relative entropy (−0.143 to −0.129), −0.138 with chip intercepts; core-essential −0.064.
Mortal cells on the same chips (12 arrays): −0.026 (−0.036 to −0.018); core-essential −0.025. Mortal cells of an
earlier batch (23 arrays, other chips): +0.087 on relative entropy and −0.164 on the entropy slope, so that group
gives no reading. One cell line; no arrested control; hTERT and mortal cells are not compared.

**Tissues of the same adults on shared chips (R131).** GSE280465 (Apsley et al. 2025): buccal swab, saliva, dried
blood spot and blood mononuclear cells of adults aged 19 to 70 on EPIC v2, rebuilt; 17 of 27 chips carry all four
tissues. Saliva is left out: 16 of 37 arrays fall below 90 % detection. Per pair of tissues, donors with both;
housekeeping gap on between-person dispersion, first tissue, second tissue and their difference with a donor
bootstrap: buccal −0.046, dried blood spot −0.071, +0.026 (−0.004 to +0.049; 33 donors); buccal −0.042, mononuclear
cells −0.107, +0.064 (+0.012 to +0.095; 43 donors); dried blood spot −0.071, mononuclear cells −0.079, +0.008
(−0.032 to +0.039; 32 donors). With chip intercepts no pair differs (+0.012, +0.042, +0.039; every interval
includes zero). Age slopes on relative entropy: −0.07 to −0.13 in every tissue, no pair differing (p 0.34 to 0.64).

**A third cortex series (R132).** GSE254293: bulk orbitofrontal cortex of 90 people aged 24 to 84, EPIC, rebuilt;
diagnosis held. Housekeeping gap on relative entropy −0.017 SD (−0.030 to −0.008), −0.010 (−0.022 to −0.001) with
chip intercepts; on the dispersion slope −0.037 and −0.022. Core-essential −0.024 and −0.016 (the second interval
includes zero).

**Passaged astrocytes (R133).** GSE202554: primary fetal astrocytes of two donors, 45 unsorted cultures over up to
15 population doublings, EPIC, rebuilt; slopes on doublings with a donor intercept. Housekeeping gap per doubling
−0.102 SD on relative entropy (−0.111 to −0.094), −0.101 with chip intercepts; core-essential −0.072; non-essential
+0.342. Most arrays (37) come from one donor.

**Astrocytes immortalised with telomerase (R135).** GSE226079: fetal astrocytes immortalised with hTERT, 43 cultures
in three replicate lines over 29 to 72 cumulative population doublings, EPIC, rebuilt; slopes on doublings with a
line intercept. The 12 chips each hold cultures of several doublings, and 81 % of the variance of doublings is left
within line and chip. Housekeeping gap per doubling −0.083 SD on relative entropy (−0.092 to −0.074), −0.084 with
chip intercepts; core-essential −0.041 (−0.060 to −0.023); non-essential +0.141 (+0.089 to +0.197). The series'
primary astrocytes (five passages) are not used.

**Three brain regions, fourth brain series (R137).** GSE349373: prefrontal cortex, amygdala and hippocampus of about
300 donors aged 16 to 80 (per region about 100 controls, 100 with major depression and 100 with post-traumatic
stress disorder), EPIC, rebuilt; one array per donor and region; donors among the adult controls of Jaffe's cohort
(same brain bank) left out; age slopes with sex and diagnosis held. Housekeeping gap on relative entropy, without
and with chip intercepts: prefrontal cortex −0.086 and −0.081 (−0.089 to −0.069; 294 donors, 62 chips), amygdala
−0.065 and −0.070 (−0.081 to −0.058; 297, 56 chips), hippocampus −0.079 and −0.069 (−0.081 to −0.056; 300, 57
chips). Core-essential with chip: −0.059, −0.043, −0.053. Controls alone with chip: housekeeping −0.043 (−0.052 to
−0.034; 100), −0.036 (−0.046 to −0.025; 99), −0.078 (−0.089 to −0.069; 98). Off the 0 % floor (10-90 %), every
donor, with chip: −0.061, −0.055, −0.045, each interval excluding zero; in the amygdala controls alone the off-floor
gap is not resolved (−0.012, −0.040 to +0.016). With diagnosis and chip held, 73 to 79 % of the variance of age is
left (36 to 49 % in the controls alone).

**Sorted cortical neurons, third series (R134).** GSE112179: NeuN+ nuclei of frontal cortex, 100 arrays of 83
donors, EPIC, rebuilt; one array per donor is used (the first by accession among those passing detection), giving 82
donors aged 23 to 77 (28 controls, the others with schizophrenia or bipolar disorder); age slopes with sex and
diagnosis held, on the 88,190 promoter CpGs measured in every array. Housekeeping gap −0.051 SD on relative entropy
(−0.060 to −0.042), −0.030 (−0.040 to −0.020) with chip intercepts; core-essential −0.043 (−0.060 to −0.027), −0.009
(−0.028 to +0.014) with chip. Off the 0 % floor (10-90 %): housekeeping −0.078, −0.043 with chip. On the dispersion
slope the housekeeping gap is -0.048 without and +0.032 with chip intercepts. With diagnosis and chip held, 62 %
of the variance of age is left. In the 28 controls alone 18 % is left within chips; without chip their housekeeping
gap is −0.030 (−0.041 to −0.018), and with chip it is not estimable (+0.024).

**Senescence by division and by damage in six donors (R136).** GSE227160: skin fibroblasts of two donors and
mesenchymal stromal cells of four (bone marrow and adipose), each as a young culture, a middle-passage culture, two
cultures in replicative senescence and two each made senescent by irradiation and by doxorubicin; EPIC, rebuilt, 47
of 48 arrays passing detection, 89,960 promoter CpGs measured in all. Every donor has one array on each of the 8
chips. Per CpG, methylation was regressed on donor, chip and treatment (young as reference); for each treatment the
squared coefficient minus its squared standard error estimates the squared change without array noise, scaled by
m(1 − m) of the young level. Mean over CpGs matched to housekeeping promoters: 0.0141 in replicative senescence,
0.0006 after irradiation, −0.0015 after doxorubicin and −0.0003 at middle passage, so only replicative senescence
shows drift beyond noise. There the share of squared drift is −0.39 for housekeeping promoters (−0.36 to −0.46
across the six fits that leave one donor out) and −0.35 for core-essential ones; non-essential promoters carry more
than their matched CpGs (+0.51). In the other three treatments the share is a ratio of two values near zero and is
not interpretable. The gap of the absolute coefficient is negative in every treatment (housekeeping −0.06 to
−0.12 SD), including those without drift; it reflects lower array noise at these promoters and is not read as a
shield.

**Arrays passing the 90 % detection cut (passing of rebuilt).** GSE279509 524 of 749; GSE252501 232 of 232; GSE66351
189 of 190; GSE197723 313 of 316; E-MTAB-15859 79 of 79; GSE283951 85 of 87; GSE112179 99 of 100 (82 donors used); GSE227160 47 of
48; E-MTAB-8993 61 of 65; GSE280465 saliva 21 of 37; Jaffe adult controls 208 of 213; EPIC-Italy 839 of 845; BXD
liver 288 of 338; Browder mouse tissues 155 of 166; Tet1 19 of 24; Dnmt3a 72 of 72; Huang
mouse tissues 581 of 583.

**Blood with the neutrophil chromatin map (R121).** Whole blood has no reference map of its own; the analyses use
mononuclear cells (Roadmap E062). With neutrophils (E030) in its place: loss per decade at active-TSS island
promoters, set minus matched, core-essential +0.0053 (+0.0026 to +0.0081; E062 +0.0063), housekeeping +0.0041
(+0.0022 to +0.0061; E062 +0.0015); drift hotspots 2.02 times expected in Polycomb-repressed chromatin (E062 2.11)
and 0.32 in active TSS (0.59); essentiality per unit with the state held +0.054 on entropy (E062 +0.046) and +0.007
on relative entropy (p 0.28; E062 +0.004). The housekeeping gap on relative entropy, −0.147 unmatched, is −0.013
with E062, −0.025 with E030 and −0.006 with both maps held at once.

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
  first rebuilt with the human recipe and were then rebuilt with the recipe of their platform. SATSA and Gill were rebuilt with methylprep
  1.7.1, labelled as such; on Jaffe's cortex the two implementations agreed within 0.004 SD (R93). Cohorts without
  raw IDATs keep the authors' processing and are marked wherever reported.
- **Sample flow.** Per series: deposited samples, samples with both IDAT channels, samples passing QC, samples in each
  model, and donors (Supplementary Table S4 names the analyses; per-series counts are in each run's log).
- **Sample QC.** Detection by pOOBAH: an array passes when it detects at least 90 % of its probes (for the
  multi-species 40k array, of the probes the species mask keeps). Recorded beside it: the bisulfite-conversion index
  (GCT, near 1 when conversion is complete), mean intensities, dye-bias statistics and each sample's beta quantile
  curve with its distance from the series' median curve (R103). Three human series have a high conversion index
  (clonal haematopoiesis 1.58, the two writer-syndrome series 1.39-1.47, against about 1.07 elsewhere). Every array of these
  series is above 1.2, so none can be excluded; the index was instead held as a covariate of the group contrast
  (R126). In clonal haematopoiesis it does not differ by group (p ≥ 0.77) and no effect changes. In the syndrome
  series it differs between cases and controls (Sotos +0.93 SD, Kabuki −2.30 SD), which shows cases and controls
  were not processed alike; the effects are unchanged with the index held (Sotos housekeeping F −0.066 → −0.068,
  Weaver −0.054 → −0.054, CHARGE −0.019 → −0.019), and are reported with that caveat.
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
  the CpGs both measure (GTEx EPIC against 450k cohorts; the Tet3 and Tet1 arrays share only 150 CpG positions, so
  those two models are compared on the genes both arrays cover, R120); RRBS and array estimates are not
  ranked against each other.
- **Second tracks.** For each analysis Supplementary Table S4 says whether chromatin, expression or outcome comes
  from the same samples, the same study or another study, and from the same cell or tissue type.
- **Published clocks.** Scoring with each clock's own inputs, reference and missing-data rule, and an audit of how
  much coefficient weight is lost per person, are not yet done (RM41); results that use published clocks are
  reported as indices on our matrices, not as the clocks' published outputs.

## Supplementary figures

| figure | file | content |
|---|---|---|
| S1 | `docs/figures/global_fdr.png` | one BH correction over 13,955 tests |
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
| Lunnon blood and four brain regions | GSE59685 | 450k | deposited; **no raw IDATs on GEO, not rebuilt**; all tissues normalised together from the deposited raw intensities as a check (R127) |
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
| Global FDR, 13,955 tests | all | R71 | `data/output/global_fdr/by_folder.csv` |
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
| Same-donor blood and brain | R68a R2 R127 | C1 | human blood and 4 brain regions | Lunnon, GSE59685 | 450k | deposited (no IDATs); joint normalisation from raw intensities as a check | none |  | one study; tissue confounded with slide | none | n.a. | blood drawn in life, brain at death; R127: joint normalisation, blood minus frontal cortex -0.105 against -0.092 |
| Cell fractions held | R70 | C1 | human whole blood | Hannum, GSE40279 | 450k | deposited | 12 cell fractions | EpiDISH reference panel (Luo 2023) | one study | other study | related | reference panel of sorted cells |
| Fibroblast cultures with own expression | R124 R125 | C2 | human dermal fibroblasts in culture | Sturm, GSE179847 | EPIC | rebuilt (R SeSAMe) | expression (RNA-seq) | Sturm, GSE179848 | one study; holds with chip intercepts (R125) | same cultures (mean over cultures) | same | 89 untreated cultures of 6 lines; 2 lines give 60 |
| Epidermis with own expression | R129 | C2 | human epidermis (women) | Holzscheck, E-MTAB-8993 / E-MTAB-8992 | EPIC / 450k (apart) | rebuilt (R SeSAMe) | expression (RNA-seq); chromatin state | Holzscheck, E-MTAB-8995 (depositors' TPM); Roadmap E057 keratinocytes | one study; holds with chip intercepts (chip read from IDAT headers) | same women; chromatin other study | same tissue; cultured keratinocytes for chromatin (related) | 61 + 56 women; cell composition not held |
| Fibroblast line with and without telomerase | R130 | C2 | human BJ fibroblasts in culture | GSE268145 | EPIC | rebuilt (R SeSAMe) | none |  | one study; hTERT and 12 mortal arrays share chips; holds with chip intercepts in hTERT | none | n.a. | one cell line; mortal groups small or inconsistent across measures; no arrested control |
| Tissues of the same adults on shared chips | R131 | C1 | human buccal, dried blood spot, blood mononuclear cells | GSE280465 | EPIC v2 | rebuilt (R SeSAMe) | none |  | one study; tissues share chips (17 of 27 hold all four) | none | n.a. | 32 to 43 donors per pair; saliva out (detection); no tissue difference with chip held; all tissues renew |
| Third cortex series | R132 | C1 | human orbitofrontal cortex (bulk) | GSE254293 | EPIC | rebuilt (R SeSAMe) | none |  | one study; age spread over 12 chips; holds with chip intercepts | none | n.a. | 90 people; diagnosis held; cell mix not held |
| Passaged astrocytes | R133 | C2 | human fetal astrocytes in culture | GSE202554 | EPIC | rebuilt (R SeSAMe) | none |  | one study; doublings spread over 8 chips; holds with chip intercepts | none | n.a. | 45 cultures of 2 donors (37 from one); 15 doublings; no arrested control |
| Astrocytes with telomerase | R135 | C2 | human fetal astrocytes immortalised with hTERT | GSE226079 | EPIC | rebuilt (R SeSAMe) | none |  | one study; doublings spread over 12 chips; holds with chip intercepts | none | n.a. | 43 cultures in 3 replicate lines; 29 to 72 doublings; no arrested control |
| Three brain regions | R137 | C1 | human prefrontal cortex, amygdala, hippocampus (bulk) | GSE349373 | EPIC | rebuilt (R SeSAMe) | none |  | one study; age spread over 56 to 62 chips per region; holds with chip intercepts | none | n.a. | about 300 donors per region (100 controls); diagnosis held; Jaffe donors left out; bulk tissue |
| Sorted neurons third series | R134 | C1 | human frontal cortex NeuN+ nuclei | GSE112179 | EPIC | rebuilt (R SeSAMe) | none |  | one study; age spread over 13 chips; smaller with chip intercepts | none | n.a. | 82 donors, one array each (28 controls); diagnosis held; neurons only |
| Senescence modes in six donors | R136 | C2 | human skin fibroblasts and mesenchymal stromal cells | GSE227160 | EPIC | rebuilt (R SeSAMe) | none |  | one study; every donor on each of 8 chips; donor and chip held | none | n.a. | 6 donors, 47 cultures; one young culture per donor; doublings not deposited |
| Monocytes and T cells with own expression | R122 R123 | C2 | human CD14+ monocytes / CD4+ T cells | MESA, GSE56046 GSE56581 | 450k | deposited | expression (HT-12 array) | MESA, GSE56045 GSE56580 | one study | same samples | same | 1,202 and 214 samples joined one to one; expression deposited normalised values |
| Sorted monocytes and T cells | R44 R68b | C1 | human CD14+ monocytes / CD4+ T cells | MESA, GSE56046 GSE56581 | 450k | deposited | none |  | one study | none | n.a. | text sets them beside Hannum (other study) |
| Read-level blood | R50 | C1 | human white blood cells | Loyfer, GSE186458 | WGBS | authors' reads | Hannum array means (check) | GSE40279 | across studies | other study | same | cross-platform used only as a validation (r 0.963) |
| GTEx nine tissues | R68d | C1 | human 9 tissues | GTEx, GSE213478 | EPIC | deposited (IDATs controlled) | none |  | one study | none | n.a. | same donors across tissues |
| Huang mouse tissues | R68a R92 R97 | C1 C2 | mouse 11 tissues | Huang, GSE279437 | Mammal320 | rebuilt (= deposited) | none |  | one study; age confounded with chip | none | n.a. | within chip (R92) the turnover ordering vanishes |
| Mouse RRBS tissues | R63 R10 R11 R7 R73 | C1 | mouse tissues and blood | Meer GSE121141; Thompson GSE120132; Petkovich GSE80672 | RRBS | deposited | none |  | one study each | none | n.a. | reported per study; not ranked against arrays |
| BXD liver | R9 R100 | C1 | mouse liver | BXD, GSE199979 | Mammal40 | rebuilt (R SeSAMe) | cell-mix components | same arrays | one study | same samples | same | HFD effect partly survives components on rebuild; Mammal40 rebuilt without SeSAMe's species masks; sample QC weak |
| Mouse turnover across platforms | paper text | C1 C2 | mouse tissues | Huang vs Meer/Thompson/Petkovich vs BXD | Mammal320 vs RRBS vs Mammal40 | mixed | none |  | across studies and platforms | none | n.a. | cross-platform: side by side only (RM41 9b) |
| Culture Endicott | R94 | C3 | human primary cell lines | Endicott, GSE197512 | EPIC | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | shield-in-culture vanishes on rebuild |
| Culture Sturm | R1 R94 | C3 | human fibroblasts | Sturm, GSE179847 | EPIC | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | unchanged on rebuild |
| Species shield | R12 R68c R59 R22 R118 R128 | C6 | 65 mammals | MMC, GSE223748 | Mammal40 | deposited (no IDATs) | lifespan and weight; NCBI taxonomy | same metadata | one study; chip partly confounded with age (56 of 139 groups); with chip intercepts 86 % of 122 fittable groups stay negative (R118) | same samples | n.a. | mouse mm10 promoter and island annotation applied to every species; the 13 human groups agree under the human annotation (R128) |
| Gain and loss rates | R67 R83 | C2 C4 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | none |  | one study each | none | n.a. |  |
| Within active promoters | R74 | C2 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | chromatin state (TssA) | Roadmap E062 (PBMC) / E073 (DLPFC) | one study each | other study | blood: related (PBMC vs whole blood); cortex: same | reference epigenomes from other people |
| Chromatin match blood | R3 R68b | C2 | human whole blood | Hannum, GSE40279 | 450k | deposited | chromatin state | Roadmap E062 PBMC, E030 neutrophils; E073 control | one study | other study | related | E062 lacks granulocytes; E073 a deliberate other-tissue control |
| Chromatin match cortex | R3 R4 R68b | C2 | human prefrontal cortex | Jaffe, GSE74193 | 450k | deposited | chromatin state and sorted ATAC | Roadmap E073 DLPFC; GSE96614 NeuN+/- ATAC BA9; E062 control | one study | other study | same |  |
| Chromatin match Lunnon brain | R3 R68b R104 | C2 | human FC STG EC (cerebellum excluded) | Lunnon, GSE59685 | 450k | deposited | chromatin state | Roadmap E073 (FC), E072 (STG), E071 (EC); BA9 ATAC for FC only | one study | other study | FC same; STG EC nearest region map (related) | R104: region maps remove 5-10 %, blood's map 25-27 %; R127 joint normalisation: 12-18 % and 35-40 % |
| Essentiality gradient | R58 R68b | C2 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | DepMap gene effect; GTEx expression; chromatin state | DepMap 24Q4 (cancer lines); GTEx v8; E062/E073 | one study each | other study | expression same tissue; essentiality other (cancer cell lines); chromatin blood related, cortex same | three external tracks |
| Drift hotspots | R72 | C1 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | chromatin state | E062 / E073 | one study each | other study | blood related; cortex same | Polycomb enrichment rests on E062 labels |
| Expression-matched placebo | R82 | C1 | human blood / cortex | Hannum, Jaffe | 450k | deposited (both) | expression | GTEx v8 whole blood / frontal cortex BA9 | one study each | other study | same tissue |  |
| Within-person SATSA deposited | R52 R60 | C1 | human blood | SATSA, S-BSST1206 | 450k + EPIC pooled | deposited | none |  | batches pooled with a chip-change term | none | n.a. | pools two platforms within people |
| Within-person SATSA rebuilt | R66 | C1 | human blood | SATSA, E-MTAB-7309 | 450k only | rebuilt (methylprep) | none |  | one platform | none | n.a. |  |
| Tet3 heterozygous | R13 R15 R16 R36 | C4 | mouse cortex / striatum | MMC, GSE223748 | Mammal40 | deposited (no IDATs) | none |  | one study; same-chip wild types | none | n.a. |  |
| Tet1 knockout | R35 R99 | C4 | mouse cortex | Foong, GSE271375 | MM285 | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | 19 of 24 arrays pass |
| Tet3 vs Tet1 | R120 | C4 | mouse brain | GSE223748 vs GSE271375 | Mammal40 vs MM285 | mixed | none |  | across studies and platforms | none | n.a. | arrays share 150 CpG positions (42 promoter), so no shared-CpG test; opposite signs hold on the 314 genes both cover |
| Dnmt3a W326R | R37 R99 | C4 | mouse marrow / intestine / liver | Sarni, GSE324160 | MM285 | rebuilt (R SeSAMe) | none |  | one study | none | n.a. | littermate controls |
| Writer syndromes | R18 R95 R126 | C5 | human blood | GSE74432, GSE97362 | 450k | rebuilt (R SeSAMe) | none |  | one study each vs own controls; conversion index differs between cases and controls (R126) | none | n.a. | concordance pools the two series' effects; bisulfite conversion index high (median GCT 1.4-1.6 vs ~1.07) |
| Clonal haematopoiesis | R25 R99 R126 | C3 | human blood | Basrai, GSE315366 | EPIC | rebuilt (R SeSAMe) | CH calls | same patients | one study | same samples | same | bisulfite conversion index high (median GCT 1.4-1.6 vs ~1.07) |
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

**Supplementary Table S5. Which analysis was run in which series.** One row per methylation series; each cell names
the run that did the analysis, or says it was not run. A statement about what explains the gap (chromatin state,
expression, kinetics) rests only on the rows whose cell for that analysis is filled; most tissues have only the gap
itself. Runs are described in the run log of the repository.

<!-- coverage:start -->

| series | species | tissue or cells | accession | platform | processing | samples | gap | off the 0 % floor | chip or batch held | chromatin state matched | expression from the same samples | accessibility | gain and loss rates | interval resamples |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Hannum blood | human | whole blood | GSE40279 | 450k | deposited | 656 donors | R68a | R80 | plate (R68a) | R3 R68b R121 | not run | not run | R67 R83 | donors |
| Jaffe cortex | human | prefrontal cortex | GSE74193 | 450k | rebuilt | 208 donors | R68a R93 | R80 | plate (R68a); chip R93 | R3 R4 R108 | not run (GTEx expression; R108) | sorted reference (R4 R108) | R67 R83 | donors |
| Lunnon blood and brain | human | blood and four brain regions | GSE59685 | 450k | deposited and jointly normalised | 71-80 donors per tissue | R2 R127 | not run | not separable (slide = tissue) | R104 | not run | not run | not run | donors |
| Three brain regions | human | prefrontal cortex, amygdala, hippocampus | GSE349373 | EPIC | rebuilt | about 300 donors per region | R137 R142 | R137 | chip (R137) | not run | not run | not run | not run | genes; donors (R142) |
| Orbitofrontal cortex | human | orbitofrontal cortex | GSE254293 | EPIC | rebuilt | 90 donors | R132 | not recorded | chip (R132) | not run | not run | not run | not run | genes |
| Second DLPFC series | human | dorsolateral prefrontal cortex | GSE252501 | EPIC | rebuilt | 52 donors | R109b | R109b | chip (R109b) | not run | not run | not run | not run | genes |
| Sorted cortical nuclei | human | neurons, oligodendrocyte lineage, other glia, unsorted | GSE279509 | EPIC | rebuilt | 89-148 donors per fraction | R109 R110 | R109 | chip (R109) | not run | not run | not run | not run | genes; donors for paired contrast (R110) |
| Second sorted series | human | neurons and glia | GSE66351 | 450k | rebuilt | 31 donors | R115 | not run | not run | not run | not run | not run | not run | donors (paired) |
| Sorted neurons | human | NeuN+ nuclei | GSE112179 | EPIC | rebuilt | 82 donors | R134b R142 | R134b | chip (R134b) | not run | not run | not run | not run | genes; donors (R142) |
| MESA monocytes and T cells | human | CD14+ monocytes and CD4+ T cells | GSE56046 GSE56581 | 450k | deposited | 1202 and 214 people | R44 R122 R123 | not run | not run | E062 (R122) | R122 R123 | not run | not run | genes |
| Epidermis | human | suction-blister epidermis | E-MTAB-8993 E-MTAB-8992 | EPIC and 450k (not pooled) | rebuilt | 117 arrays (up to 31 women twice) | R129 | not run | chip (R129) | E057 (R129) | R129 | not run | not run | genes |
| GTEx tissues | human | nine tissues incl. lung, colon, kidney, muscle, testis | GSE213478 | EPIC | deposited | 47-224 donors per tissue | R68d | not run | not run | not run | not run | not run | not run | donors (32 per resample) |
| Fibroblast cultures | human | primary fibroblasts | GSE179847 | EPIC | rebuilt | 89 cultures of 6 lines | R124 R125 | not run | chip (R125) | not run | R124 R125 | not run | not run | genes |
| Telomerase fibroblasts | human | BJ fibroblasts with and without hTERT | GSE268145 | EPIC | rebuilt | 20 hTERT cultures | R130 R138 | not run | chip (R130) | not run | not run | not run | not run | genes; cultures (R138) |
| Passaged astrocytes | human | fetal astrocytes | GSE202554 | EPIC | rebuilt | 45 cultures of 2 donors | R133 R138 | not run | chip (R133) | not run | not run | not run | not run | genes; cultures (R138) |
| Astrocytes with telomerase | human | fetal astrocytes with hTERT | GSE226079 | EPIC | rebuilt | 43 cultures of 3 lines | R135 R138 | not run | chip (R135) | not run | not run | not run | not run | genes; cultures (R138) |
| Senescence modes | human | fibroblasts | GSE197723 | EPIC | rebuilt | 12 chips x 4 arms | R111 | not run | within chip (R111) | not run | not run | not run | not run | chips (sign-flip) |
| Senescence in six donors | human | fibroblasts and stromal cells | GSE227160 | EPIC | rebuilt | 6 donors x 8 cultures | R136 | not run | donor and chip (R136) | not run | not run | not run | not run | leave one donor out |
| Arrested fibroblasts | human | confluent fibroblasts | E-MTAB-15859 | EPIC | rebuilt | 3 donors (2 usable) | R112 | not run | within chip (R112) | not run | not run | not run | not run | month reorderings |
| Four tissues of the same adults | human | buccal, saliva, dried blood spot, PBMC | GSE280465 | EPIC v2 | rebuilt | 32-43 donors per pair | R131 | not run | chip (R131) | not run | not run | not run | not run | donors |
| Read-level blood | human | white blood cells and 27 sorted cell types | GSE186458 | WGBS | deposited reads | 23 donors | R50 | not run | not run | not run | not run | not run | not run | genes |
| Within one person | human | blood over 20 years | E-MTAB-7309 | 450k | rebuilt | 444 people | R66 | not run | batch (R66) | not run | not run | not run | not run | people |
| Single muscle stem cells | mouse | muscle stem cells | GSE121437 | scBS-seq | deposited coverage | 275 cells of 6 mice | R140 | not run | not run | not run | R140 (same cells) | not run | not run | cells within mouse |
| Colon and hippocampus | mouse | colon and hippocampus | GSE233734 GSE304754 | RRBS | deposited coverage | 82 and 83 mice | R114 | R114c | shared flowcells (not modelled) | not run | not run | not run | not run | mice |
| Four tissues of the same mice | mouse | lung, heart, cortex, liver | GSE225166 | RRBS | deposited coverage | 30-32 mice per tissue | R141b | R141c | lane (R141b) | not run | not run | not run | not run | mice |
| Huang mouse tissues | mouse | 11 tissues | GSE279437 | mouse array | rebuilt | young and aged controls | R68a R92 R97 | not run | chip (R92 R97; ordering not identifiable) | not run | not run | not run | not run | mice |
| Mouse RRBS series | mouse | blood; six tissues; four tissues | GSE80672 GSE120132 GSE121141 | RRBS | deposited coverage | about 20 mice per tissue | R63 R10 R11 R7 R73 | not run | not run | not run | not run | not run | not run | mice |
| Mammals | 65 species | 139 species-tissue groups | GSE223748 | mammalian array | deposited | 15 or more animals per group | R12 R68c | R59 | chip (R118) | not run | not run | not run | not run | groups |
| Colon mucosa next to cancer | human | normal colon mucosa | GSE101764 | 450k | rebuilt | 146 patients | R143 | R143 | chip (R143) | not run | not run | not run | not run | genes; donors |
| Lung next to cancer | human | normal lung of non-smokers | GSE314841 | 450k | rebuilt | 114 people | R143 | R143 | chip and batch (R143) | not run | not run | not run | not run | genes; donors |
| Human tissues on the mammalian array | human | heart, muscle, liver, kidney, lung, spleen, fat, blood | GSE184211 | mammalian array | rebuilt | 32-112 per tissue (about 1400 promoter CpGs) | R143 | R143 | chip (R143) | not run | not run | not run | not run | genes; donors |

<!-- coverage:end -->

## Open items for the authors

- Reference metadata marked [verify] in `paper.md` (authors, journals or years not given in the project docs).
- Read-level epiallele entropy and the per-person `K` measure are not yet recomputed on relative entropy.
- Joint models with DunedinPACE/PhenoAge (EPIC-Italy, SJLIFE) and the clock-composition shares (R38 addendum) were run
  ad hoc from per-person tables and are not yet in the pipeline.
- The culture runs (GSE197512) predate the R-numbered run log; their run identifiers are not recorded there.

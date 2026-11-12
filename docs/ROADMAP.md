# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.
The order and the rule for each item are argued in `docs/RESEARCH_PLAN.md`.

## RM15 — Directed demethylation outside promoters: the culture and cortex follow-ups

**State:** after RM19 · **Owner:** lead

Progress: first part done (`figures intergenic`, R14, `docs/extension.md`): in blood, solo-WCGW CpGs lose
methylation 2.7 times as fast as other open-sea CpGs (directed, not noise); none in cortex; the PCA does not single
it out. Left: the same test in the culture series (EPIC manifest needed for the tags; the loss should grow with
doublings), and what cortex's two strongly age-tracking components are.
## RM12 — Selection pressure under interventions

**State:** unattended · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): do interventions change the shield through selection, apart from
chromatin? What we have: in B6D2F1 blood CR weakens the shield by a quarter (within one sequencing series, robust
to composition), without slowing drift; the shield's share also grows with age in liver, which rarely renews. So
neither proof nor disproof. Tests that would separate selection from chromatin, in order of cost:
- Done (R13, `docs/extension.md`): CR in liver does not weaken the shield; Tet3 heterozygous mice lose it in
  cortex and striatum (controls pooled across studies, a batch could do it); GHRKO strengthens it in liver.
- Done (R15): against same-chip littermates the Tet3 loss holds and grows (housekeeping F +0.17, +0.25); the GHRKO
  liver effect was batch.
- Done (R16): under Tet3 loss essential CpGs gain methylation, but S rises in other arms too (array compression
  could do it); not specific yet.
- Next: RM14's syndromes in human blood (built: `figures syndromes`).
- Second model found (`docs/REPROG_TET_SURVEY.md`): GSE271375, Tet1 knockout and catalytic-dead cortex, 8 per
  group, mouse array; then Dnmt3a gain of function (GSE324160). Next after the quota resets.
- Done (R35): Tet1 knockout and catalytic-dead cortex keep or strengthen the shield while drifting more;
  the opposite of Tet3. Next: Dnmt3a gain of function (GSE324160).
- **The same intervention in a renewing and a non-renewing tissue of the same animals.** Selection predicts CR
  or rapamycin changes the shield in blood or gut and not in brain or muscle; chromatin predicts the same change
  (or none) everywhere. Search for a multi-tissue CR or rapamycin series with at least 8 animals per arm.
- **Interventions that raise selection pressure.** Chemotherapy, irradiation or transplantation shrink the stem
  cell pool and let clones expand. Selection predicts a stronger shield afterwards, chromatin no change. Search
  for blood methylation before and after such bottlenecks (for example bone-marrow transplant donors and
  recipients, or patients after chemotherapy).
- **Chromatin held fixed.** Match on chromatin first (RM7's method), then test the intervention on the remainder.
  If the CR effect survives the matching, it is not a chromatin change.

## RM11 — Clocks of our own

**State:** building · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): can we build clocks from the entropy-age ideas, and do they capture
more than conventional clocks? Candidates from the results: a stochastic clock from Polycomb CpGs with each
tissue's offset calibrated, and the per-person shield score beside it. "More" needs a readout beyond age: which
clock tracks interventions (CR mice, GSE80672) and, later, outcomes in a cohort with follow-up.
- antonkulaga (2026-11-12): among people with the same shield, are those with a higher methylation age less
  protected? Plan, first tomorrow (no downloads): (1) per-person P against age acceleration of Horvath, Hannum,
  PhenoAge, GrimAge in Hannum and Jaffe; (2) stochastic clocks from look-alike CpGs only and from shielded CpGs
  only; (3) stratify acceleration on P.
- Done (R29): at equal drift, faster chronological and stochastic clocks go with a stronger shield; outcome
  clocks (PhenoAge, DunedinPACE, mortality) are blind to it. Next: train our own clocks, (2).

## RM13 — Answer the prior work (novelty review)

**State:** after RM12 · **Owner:** lead

Progress: `docs/NOVELTY_REVIEW.md`. To make the paper new, in order: (1) does essentiality still predict protection
after matching on chromatin, expression level and breadth? First answer from R4's tables (`docs/extension.md`):
in blood, core-essential genes keep almost nothing beyond chromatin, housekeeping genes a small clear remainder,
brain both; left: a DepMap gradient within broadly expressed genes, and expression level as a matching variable.
(2) absolute drift slopes beside the CR share; (3) clonality in the culture series where reads allow (epiallele
structure), against Hollmann 2020; (4) read the species result against Crofts 2026 (stem-cell dynamics) and the
floor of near-0 % essential promoters; (5) the shield score against Seeboth 2020's confounders.

## RM18 — Normalization and probe robustness

**State:** building · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12): are results comparable across arrays and array versions, and do we
follow best practice? By design every comparison stays within one dataset and platform, on the deposited betas,
matched on island status and starting level. Gaps:
1. Infinium probe type (I/II) is not in the strata; island CpGs at essential promoters are more often type I.
   Add it, rerun the headline gaps (Hannum, Jaffe, Lunnon).
2. Cross-reactive and SNP-affected probes are not masked: apply Zhou et al. 2017's masks, rerun.
3. One cohort reprocessed from raw IDATs (SeSAMe, or noob + BMIQ) against its deposited betas.
4. Each sample's beta compression as a covariate for the signed shift S.
- Done (R20): items 1-2. Design explains about a sixth of blood's gaps; the mask nothing; no sign changes.
5. A short review of current practice, recorded in `docs/`.

## RM19 — Does chromatin accessibility become stochastic with age?

**State:** after RM18 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12): can we collect enough comparable ATAC-seq data to test stochastic
changes of chromatin opening? Test: per region, the age slope of between-individual (or between-cell)
accessibility variability, essential promoters against matched regions, as for methylation. First a survey of
public data (one tissue, one protocol, many individuals across ages, processed counts).
- Survey done (`docs/ATAC_SURVEY.md`): no published test of per-region accessibility variance against age, none
  against essential promoters. First test: Morandini 2023 PBMC bulk ATAC (GSE193140, 157 libraries, 20-74 y,
  10 MB counts). Then Catching 2026 prefrontal cortex snATAC (362 donors, 15-100 y, Zenodo 18394349, 14.5 GB)
  and Zemke 2024 hippocampus (GSE278576, 40 donors). Controlled-access cohorts (Márquez 2020, CIMA) later.
- First test done (`figures atac_variance`, `docs/extension.md`): with cell mix removed, housekeeping and
  core-essential promoter peaks spread apart more slowly with age (−0.05 SD, p 0.02; −0.08 SD, p 0.05); raw not
  significant. Next: per cell type in the 362-donor cortex set (a 14.5 GB download; a Slurm job).
- Done (R23, R24): in all seven cortical cell types essential promoters spread apart with age *more* than
  matched ones (core-essential +0.07 to +0.23 SD after quality and composition adjustment). No general
  accessibility shield. Next: hippocampus (Zemke 2024, GSE278576) to replicate in brain; a sorted blood set.
- Done (R27): hippocampus does not replicate cleanly: signs flip with the adjustment in 20-40 donors.
  Accessibility spread is too sensitive to sample quality to read the shield. RM19 parks; revisit with a
  cohort that has per-sample TSS enrichment and fragment QC, or sorted bulk ATAC.

## RM20 — Why would long life or CR weaken the shield?

**State:** after RM16 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12). The pattern: no lifespan effect overall (p 0.40), a lean toward a
weaker shield in blood of long-lived species (p 0.054), CR weakening it in one blood study and not in liver,
Tet3 loss weakening it in brain. Hypotheses, each with its test:
1. Catch-up: better fidelity everywhere lowers drift at unprotected CpGs toward the protected floor. Test: species
   look-alike drift B per year and per lifespan fraction against lifespan; CR blood already shows B unchanged.
2. Purging: part of the blood shield is selection in fast-renewing marrow; long life and CR damp proliferation.
   Test: blood against skin and liver within species; rapamycin and GHRKO blood; RM16 clonal haematopoiesis.
3. Investment: long-lived species favour genome-wide fidelity (DNMT1, UHRF1, repair) over targeted promoter
   protection (TET, CFP1, H3K4me3). Test: cross-species expression of that machinery against lifespan (public
   multi-species RNA, e.g. Fushan 2015, Ma 2016).
4. Tet3 removes the protection directly, in non-dividing cells: a separate mechanism (done, R15-R16).
- Done (R22): hypothesis 1 not supported: look-alike drift falls as lifespan^−1.4, but the relative gap holds or
  grows.
- The blood lean rests largely on human (jackknife: +0.085 → +0.047 without it, p 0.30); 12 and 7 species
  with paired tissues are too few for hypothesis 2. RM20 parks until the pattern is firmer (Anton told).

## RM22 — Stochasticity at the level of single molecules (long reads)

**State:** after RM11 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12): can nanopore data add insight on the stochasticity of clocks? Long
reads give many CpGs on one molecule from one cell: within one sample, do essential promoters stay uniform
molecule to molecule while look-alike regions become mixtures of epialleles with age (epiallele entropy,
methylation haplotype load, discordant reads)? Read-level WGBS does the same on shorter reads.
- Survey done (`docs/LONGREAD_SURVEY.md`): no open long-read cohort with ages (CARD brain ONT, 351 donors,
  controlled). First test: Loyfer 2023 WGBS read-level PAT files (GSE186458), 23 white-blood-cell donors aged
  21-75, ~13 GB, `wgbstools homog`/entropy, essential promoters against matched regions. No published
  single-molecule test of epiallele heterogeneity against age.

## RM23 — The stochastic component of clocks, with the shield (extending Tong et al. 2024)

**State:** building · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12): extend Tong et al. 2024 (Nature Aging; stochastic clocks StocH, StocP,
StocZ) and recompute the stochastic part of the same clocks. They simulate uniform stochastic DNAm change from a
young ground state and retrain each clock's CpGs on it; the share of the clock's accuracy reproduced is the
"stochastic component" (66-75 % Horvath, 90 % Zhang, 63 % PhenoAge). Our extension:
1. Reproduce their fractions on Hannum blood with the pinned Stoc and original clocks.
2. Shield-aware simulation: CpG-specific rates (shielded CpGs slower by the measured gap; or each CpG's observed
   drift rate) instead of uniform; retrain; recompute the fractions. Does shaped noise explain more?
3. Split each clock's prediction into shielded and unprotected CpGs: which carries age, which acceleration.
4. Cell composition (antonkulaga, 2026-11-12: PhenoAge, the least stochastic clock, may mostly reflect blood-cell
   mix): estimate 12 blood-cell fractions by reference-based deconvolution (EpiDISH cent12CT or Salas 2022 IDOL),
   then recompute each clock's stochastic fraction and acceleration with and without them. Is PhenoAge's
   non-stochastic part composition?
Methods read (scratch notes; key points): rates per CpG are set by its observed old-minus-young effect and
changes go in that direction (probit-scale half-normal steps, P = 1 - exp(-gamma |effect|)); ground state from the
youngest MESA monocytes; elastic net alpha 0.5 on simulated data; fraction RR2 = R2(Stoc) / R2(original);
intrinsic acceleration removes 11 of 12 EpiDISH fractions. Variants to compare: their model; direction only;
pure noise; rates by gene class (shield) only.

## RM24 — Does the shield track health? Cross-sectional tests on data at hand

**State:** building · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12): with the data already downloaded, can we estimate whether the shield
predicts health? No cohort on disk follows people forward; cross-sectional health states can be tested:
Alzheimer's pathology (Lunnon, AD against control and Braak stage, four brain regions and blood), schizophrenia
(Jaffe cortex; GSE41037 blood), smoking (GTEx). Per cohort, cases against controls scored as the syndromes are
(controls leave themselves out), age and technical covariates; shield share F, drift B, signed shift S.
Cross-sectional: a weaker shield may be a consequence of disease or treatment, not a predictor.

## RM25 — Which tissues age fastest, and is the order shared across species and layers?

**State:** done (R40) · **Owner:** lead

Progress: from docs/BEYOND_SHIELD.md idea 5 (antonkulaga asked for questions beyond the shield, 2026-11-12). Rank
tissues by drift and clock pace from result tables on disk (human GTEx, mouse array and RRBS, RNA in GTEx and
rodents) and compare the orders.

## RM26 — Do all organs of one person age at the same rate?

**State:** unattended (next) · **Owner:** lead

Progress: idea 1 of docs/BEYOND_SHIELD.md. Per donor, acceleration (clock residuals, drift B, entropy) in each
tissue; correlations between tissues of the same donors (Lunnon blood and 4 brain regions, GTEx, Huang mice).
Lunnon blood was drawn before death: each sample keeps its own age. Literature (docs/LITERATURE_COVERAGE.md): clock
coupling across tissues is published (Horvath 2022 GeroScience; Richardson 2025 on this GTEx set), so the clock arm
is a replication; the new part is drift and entropy coupling without clocks.

## RM27 — Do interventions reverse disorder or shift a programme?

**State:** waits until RM26 · **Owner:** lead

Progress: idea 3 of docs/BEYOND_SHIELD.md. Split aging CpGs into directed (consistent mean change) and noisy
(variance rise only); ask which part each intervention moves.

## RM28 — Does the shield (and drift) predict death and disease? Open follow-up cohorts

**State:** unattended · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12) for longitudinal cohorts with age and health; survey in
docs/LONGITUDINAL_SURVEY.md. Open: SJLIFE GSE314261 (4,520 people, 294 deaths with age at death, 378 with 2-5
draws), Vitality 90+ GSE68194 (115 nonagenarians, 39 deaths, follow-up in days), EPIC-Italy GSE51032 (845, 421
incident cancers, time to diagnosis), SATSA S-BSST1206 (596 people, up to 8 draws over 20 years; no mortality).
Order: Vitality 90+ (small, a true survival test), EPIC-Italy, SJLIFE, SATSA (within-person change). Controlled
cohorts with mortality (Generation Scotland, Lothian, Framingham, Normative Aging Study) need antonkulaga's
application; links in the survey.

## RM29 — When does drift begin? A change point at sexual maturity across species

**State:** unattended, after RM28's first cohorts · **Owner:** lead

Progress: idea 6 of docs/BEYOND_SHIELD.md; open per docs/LITERATURE_COVERAGE.md (drift rate scales with lifespan,
Bertucci-Richter 2023; no change-point test at maturity). Fit per species a drift curve with a breakpoint and
compare it with age at sexual maturity.

## RM30 — Does methylation spread reach expression spread? Same GTEx donors

**State:** unattended · **Owner:** lead

Progress: idea 7 of docs/BEYOND_SHIELD.md; the dispersion coupling is open per docs/LITERATURE_COVERAGE.md (mean-
level coupling is published: Oliva 2023). Per gene, the age rise in between-donor spread of promoter methylation
against that of expression, in lung, colon and ovary donors with both.

## RM16 — Selection by cellular fitness, in vivo

**State:** after RM22 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12, #entropy-age > general chat): what about selection by cellular
fitness? The evidence so far is mostly against selection (culture, turnover, chromatin matching, Tet3, liver),
with two leads for it (CR in one blood study; human bone marrow and lymph node reversing the shield). Steps:
1. Clonal haematopoiesis: blood methylation of CHIP carriers (DNMT3A, TET2 clones) against non-carriers, on GEO.
2. The bone-marrow and lymph-node anomaly of R12: donors' ages and any clonality data.
3. Sorted stem cells against mature blood, young against old, where public; else a DFG work package.
4. A second within-batch CR blood series.
- Survey done (`docs/CHIP_SURVEY.md`): open per-person CH status only in GSE315366 (Basrai 2026; 29 CH, 35 not,
  EPIC blood, cardiogenic-shock patients, no VAF); within-person mutant against wild-type cells in GSE158067 (Nam
  2022, scRRBS, 2 donors). Large cohorts with VAF are controlled access (TOPMed, Lothian, Danish twins).
  First test: GSE315366 with the `syndromes` scoring (carriers against non-carriers, age covariate).
- Done (R25): CH carriers do not differ (any, DNMT3A, TET2; p ≥ 0.15); diluted clones, no VAF.
- GSE158067 (2 donors) has no age range within a person, so it cannot test drift with age; it could only
  compare cell-to-cell deviation of mutant and wild-type cells. Lower priority than RM19's replication.

## RM17 — Partial reprogramming: does a reset restore the shield?

**State:** after RM12 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12). The consortium's OSKM arms are 2-5 mice (3.4F+7mDox kidney F −0.17,
p 0.10; skin n.s.): too small. Public array series to find: maturation-phase transient reprogramming of human
fibroblasts (Gill et al. 2022), in-vivo cyclic OSKM mice (Browder et al. 2022; Ocampo et al. 2016). Chromatin
maintenance predicts the shield restored with the clock; selection predicts no change (no cells removed).
- Survey done (`docs/REPROG_TET_SURVEY.md`). First: GSE190665 (Browder 2022, cyclic OSKM in vivo, 6 mouse
  tissues, Mammal40 array as in the consortium, 2-7 mice per group, 49 MB); then GSE156557 (RRBS, 5 per group).
- Done (R32): no consistent change (muscle stronger, spleen weaker, four tissues unchanged); the regimen barely
  moves promoter CpGs toward youth. Next if wanted: GSE156557 (one OSKM cycle, RRBS).

## RM14 — Natural experiments in chromatin: broken writers, knockouts, edits

**State:** after RM12 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): cell lines tell little about clocks; can public data show what
happens when chromatin is changed? Three sources, in order:
1. Blood methylation of people with a constitutional variant in a chromatin writer: Sotos (NSD1), Weaver (EZH2),
   Tatton-Brown-Rahman (DNMT3A), Kabuki (KMT2D). Martin-Herranz et al. 2019 found clocks accelerated in Sotos.
   Test: is the essential-gene shield weaker than in matched controls? Search GEO for the cohorts.
2. Mouse knockouts: Tet3 heterozygous brain (GSE223748, at hand after R12); Tet2 and Dnmt3a mutant blood.
3. Edited cells over time: CRISPRoff (Nunez et al. 2021) tracked imposed promoter methylation for months.
Rule: a cohort with at least 10 carriers and matched controls on one platform.
- Done: Tet3 heterozygous brain loses the shield against littermates (R15, R16; see RM12).
- Done (R18): Sotos, NSD1 variants, Weaver, Kabuki and CHARGE carriers drift more in blood but keep the shield
  (F unchanged or more negative); KMT2D and CHD7 variants show nothing within the controls' ages.
- Next: the DNMT3A loss-against-gain pair below; a TET2 or TET3 human cohort if one exists.
- Found, below the rule: DNMT3A loss of function (TBRS) GSE128801, 5 carriers against 7 controls; DNMT3A gain
  of function (microcephalic dwarfism) GSE120428, 4 against 5 in blood. Together a loss-against-gain contrast of
  one writer; read only as a pair, after R17.

## RM3 — Interventions: replicate the CR result

**State:** after RM2 · **Owner:** lead

Progress: GSE80672 done (`figures interventions`, `docs/extension.md`): in living mouse blood the shield's share
grows with age, and CR weakens it by a quarter in the one within-batch comparison (B6D2F1, 12 against 10 mice,
robust to composition), without slowing drift. Left, to replicate it:
- **GSE89275** (liver RRBS: Ames dwarf, rapamycin, CR) and its WGBS companions GSE89273/89274: about 4 mice
  per arm, below the plan's 8; not used.
- Another blood series with CR and controls sequenced together (to find).
- GSE199979 (BXD liver) done (`figures bxd`): the shield's share grows with age in liver too, so the growth needs
  no renewing tissue; HFD's effect is the cell mix. It has no CR arm.
- Left: a second CR (or rapamycin) series with treated and control animals sequenced together, at least 8 per
  arm, in blood or another tissue (to find).
## RM2 — A turnover gradient across many tissues (methylation)

**State:** building · **Owner:** lead

Progress: mouse (GSE279437), human GTEx (GSE213478) and muscle fibres (GSE304045) done; see
`docs/extension.md`. In mouse the protection is missing in cortex and hippocampus but kept in heart;
in GTEx muscle is weakly protected, as are breast, ovary and testis, so the split does not follow
turnover. Left:

- **GSE184211** (human, 11 tissues, mammalian array, exact ages 26–91, 190 MB). Not paired by donor, but
  it shares its CpGs with the mouse array, so the comparison across species is direct.
- **GSE50498** (450k muscle, ages 18–89): a muscle cohort with exact ages and more donors than GTEx's 47.
  Like GTEx it says nothing of exercise, and satellite-cell divisions follow loading, injury and
  wasting, so muscle tests selection only in a cohort that records physical activity or training
  history beside age (to find; the fibre cohort GSE304045 records training but not age).
- **Mouse replication:** Thompson 2018 (GSE120132) done (`figures rrbs_tissues`): the shield in all six tissues,
  ordered by neither size nor growth by turnover. Left: Meer 2018 (GSE121141; liver, lung, heart, cortex, 6-30 months).
- The fibre contrast (GSE304045) was uninformative: a dozen adults of similar age show too little
  between-person divergence. A within-tissue turnover contrast needs aged and young people.

## RM4 — Stochasticity in RNA

**State:** after RM2 · **Owner:** lead

Progress: the turnover test is done in GTEx v10 (`figures rna_gtex`, 47 tissues) and in the mouse and rat
atlases (`figures rna_rodents`); see `docs/extension.md`. Mouse repeats the methylation result (brain
without protection, heart with it), rat leans the same way, GTEx shows no turnover effect in spread.
Left:
- **Health at death (antonkulaga, 2026-11-11):** some GTEx donors died fast, others after long illness. Hardy is a
  covariate already; rerun GTEx RNA and methylation on fast deaths only (Hardy 1-2) and say what moves.

- **A stochastic RNA clock**, trained on simulated transcriptional noise as the paper does for
  methylation, tested on the rodent atlases (exact ages).
- **GTEx with more power:** the same-donor design (every tissue of one donor against that donor's
  blood), and subsamples larger than 70 donors in the tissues that allow them.
- Tabula Muris Senis bulk **GSE132040**, a third rodent cohort.
- Interventions: **GSE280382** (rapamycin and exenatide) and **GSE131754** (eight liver interventions).

## RM5 — Direct selection test: protection against division count in culture

**State:** after RM4 · **Owner:** lead

Progress: GSE197512 and its replication GSE179847 done (`figures passage`, `figures passage_sturm`,
`docs/extension.md`): in both the protection is a roughly constant fraction of drift and does not grow
with doublings; the MMC arms (25 days) and contact inhibition (120 days) carry too little drift to decide
whether protection needs divisions. Left:

- An arrest with more drift than GSE179847's contact inhibition gives (quiescence over a year, or senescent
  cells kept for months), so that a gap without divisions can be measured above the noise floor.
- GSE268145 (BJ, 175 MB, downloaded), GSE202554 (astrocytes), GSE131280 (glucose): smaller replications.

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
## RM13 — Answer the prior work (novelty review)

**State:** building · **Owner:** lead

Progress: `docs/NOVELTY_REVIEW.md`. To make the paper new, in order: (1) does essentiality still predict protection
after matching on chromatin, expression level and breadth? First answer from R4's tables (`docs/extension.md`):
in blood, core-essential genes keep almost nothing beyond chromatin, housekeeping genes a small clear remainder,
brain both; left: a DepMap gradient within broadly expressed genes, and expression level as a matching variable.
(2) absolute drift slopes beside the CR share; (3) clonality in the culture series where reads allow (epiallele
structure), against Hollmann 2020; (4) read the species result against Crofts 2026 (stem-cell dynamics) and the
floor of near-0 % essential promoters; (5) the shield score against Seeboth 2020's confounders.
- (1) done (R58): within genes expressed in 50+ GTEx tissues, DepMap essentiality still predicts less drift with
  chromatin state, expression and breadth held (blood entropy 0.022 SD per unit, cortex 0.018-0.027), a tenth or
  less of the raw contrast, in cortex too. Next: (2) absolute drift slopes beside the CR share.
- (2) answered: the CR tables give drift B beside F; CR does not slow absolute drift in the within-batch series
  (B6D2F1 B −0.004, p 0.37; −0.008 with 5 components, p 0.12; R55), so its effect is on the share, not the rate.
- (4) done (R59): off the 0 % floor (10-90 % CpGs) the species shield holds (−0.24 SD, negative in 93 % of
  groups) and the lifespan lean vanishes (−0.014, p 0.48): the lean was the floor; fits Crofts 2026.
- (3) no read-level culture series at hand (the passage series are arrays); waits for one.
- (5) done (R60): slide batch explains 36 % of B against 13 % by chance, 20-23 % of F and S; smoking moves only S;
  chemotherapy moves F against B (alkylating F −0.16, anthracyclines +0.21). Batch must be held in every use.

## RM19 — Does chromatin accessibility become stochastic with age?

**State:** unattended · **Owner:** lead

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

**State:** blocked: desk #28 (CALERIE), parked by antonkulaga until the pattern is firmer · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12). The pattern: no lifespan effect overall (p 0.40), a lean toward a
weaker shield in blood of long-lived species (p 0.054), CR weakening it in one blood study and not in liver,
Tet3 loss weakening it in brain. Hypotheses, each with its test:
1. Catch-up: better fidelity everywhere lowers drift at unprotected CpGs toward the protected floor. Test: species
   look-alike drift B per year and per lifespan fraction against lifespan; CR blood already shows B unchanged.
2. Purging: part of the blood shield is selection in fast-renewing marrow; long life and CR damp proliferation.
   Test: blood against skin and liver within species; rapamycin and GHRKO blood; clonal haematopoiesis (R25), bottleneck (R57).
3. Investment: long-lived species favour genome-wide fidelity (DNMT1, UHRF1, repair) over targeted promoter
   protection (TET, CFP1, H3K4me3). Test: cross-species expression of that machinery against lifespan (public
   multi-species RNA, e.g. Fushan 2015, Ma 2016).
4. Tet3 removes the protection directly, in non-dividing cells: a separate mechanism (done, R15-R16).
- Done (R22): hypothesis 1 not supported: look-alike drift falls as lifespan^−1.4, but the relative gap holds or
  grows.
- The blood lean rests largely on human (jackknife: +0.085 → +0.047 without it, p 0.30); 12 and 7 species
  with paired tissues are too few for hypothesis 2. RM20 parks until the pattern is firmer (Anton told).
- Since (2026-11-13): hypothesis 2 tested three ways without support: the CR effect goes mostly with chromatin
  openness (R55), a second CR blood series does not repeat it (R56), a transplant bottleneck changes nothing
  specific to essential genes (R57). The pattern got weaker, not firmer; hypothesis 3 is the one left untested.

## RM28 — Does the shield (and drift) predict death and disease? Open follow-up cohorts

**State:** unattended (the SATSA raw-IDAT rerun, with the IDAT pipeline R61 validated) · **Owner:** lead

Progress: asked by antonkulaga (2026-11-12) for longitudinal cohorts with age and health; survey in
docs/LONGITUDINAL_SURVEY.md. Open: SJLIFE GSE314261 (4,520 people, 294 deaths with age at death, 378 with 2-5
draws), Vitality 90+ GSE68194 (115 nonagenarians, 39 deaths, follow-up in days), EPIC-Italy GSE51032 (845, 421
incident cancers, time to diagnosis), SATSA S-BSST1206 (596 people, up to 8 draws over 20 years; no mortality).
Order: Vitality 90+ (small, a true survival test), EPIC-Italy, SJLIFE, SATSA (within-person change). Controlled
cohorts with mortality (Generation Scotland, Lothian, Framingham, Normative Aging Study) need antonkulaga's
application; links in the survey.
- Done (R42): Vitality 90+ uninformative (mortality clocks fail too). Done (R46, R46b): EPIC-Italy, drift B
  predicts breast cancer years ahead, robust to chip, cell mix and mean beta; the shield share does not. Next:
  a second breast-cancer cohort (NOWAC or Breakthrough Generations, controlled: antonkulaga's application).
- Done (R48, R48b): SJLIFE deaths; drift and shield add nothing beyond DunedinPACE.
- Done (R52, R52b): SATSA. A person's own shield share grows with age (F -0.14 SD per decade, as between
  people); CpG-level divergence gap small (-0.017), but the deposited betas keep little of the shield even
  cross-sectionally. Open data left: none with outcomes; a raw-IDAT SATSA rerun (E-MTAB-7309) is optional.

## RM33 — Metric audit: correct what rests on binary entropy

**State:** after RM34 · **Owner:** lead

Asked by antonkulaga (2026-11-14): go through previous discoveries, evidence and hypotheses, correct the
computations, update the artifact. Binary entropy measures distance from 0.5; R67 showed essential promoters lose
methylation faster, which lowers entropy at low levels, so entropy gaps mix direction with stability. Relative entropy
from the young state (`kl_slope`, per-person `K`) added (e3ad29c); R68a-d rerun every figure whose headline rests on
entropy. Then `docs/METRIC_AUDIT.md` (each result kept, resized or overturned), CLAIMS, state, artifact.
- Done (R68): the shield keeps its sign everywhere at about half the size; three statements withdrawn (R58 blood,
  R22, CD4+ T size); cortex chromatin share cut; mouse array now orders the shield by turnover.
- Done (R73): the mouse RRBS, BXD, intervention and Meer gaps keep their sign, smaller. Left: per-person `K` in
  effect tables, read-level entropy.

## RM34 — Bias checks from the omics re-analysis practice

**State:** unattended (left: download checksums in `docs/PROVENANCE.md`) · **Owner:** lead

Asked by antonkulaga (2026-11-14), against `longevity-skills/omics-reanalysis/SKILL.md`, cheapest first: (1)
placebo gene sets (random promoter sets with the same size, island and level mix must give 0); (2) Hannum's gaps
with 12 Houseman/EpiDISH cell fractions as covariates; (3) one FDR over every p reported; (4) stochastic or
specific: the CpGs whose spread rises most, enrichment by region class against a mean-matched background; (5) a
provenance table per dataset.
- (1) done (R69): placebo sets sit at or above 0; none reaches an essential gap.
- (4) done (R72): in blood the fastest-spreading CpGs are a Polycomb programme (2×) kept out of active chromatin;
  in cortex the spread is broad (0.7-1.3× across classes). "Stochastic" holds for cortex, not blood.
- (5) done: `docs/PROVENANCE.md`, one row per dataset (deposited or IDATs, platform, figures); checksums next.
- (3) done (R71): one BH over 7,626 reported tests; p up to 0.024 survives. Lost: CR housekeeping share in the
  second blood series, a few knockout contrasts, any-cancer and colorectal drift, GTEx smoking and death type.
- (2) done (R70): cell fractions leave the blood gaps unchanged; dispersion's shrinks by a third.

## RM32 — Drift toward each site's own equilibrium, not toward 0.5

**State:** unattended · **Owner:** lead

Asked by antonkulaga (2026-11-13): methylation repair is not symmetric (maintenance by DNMT1, de novo by DNMT3A/B,
active loss by TET), so should the Markov model's entropy be adjusted? The model already has per-site maintenance
and de novo rates and so a per-site equilibrium z (`docs/model.md` §2); the binary-entropy slope assumes drift
toward 0.5 and misreads sites whose z is far from it, which essential promoters are. Fit z and λ per CpG
(`fit_relaxation`, z bounded by island status where the age range cannot pin it) in Hannum blood and Jaffe's adult
cortex; compare essential promoters with look-alikes on the fraction of the way to z, on the divergence from z,
and on z itself. Answers whether the shield is a slower approach to the same equilibrium or a different
equilibrium (repair biased toward unmethylated).
- Done (R67): a different resting level. Essential island promoters gain no faster than look-alikes and lose
  faster (blood +0.0025, cortex +0.0054 per decade), so they rest lower with less room to drift.
- Done (R74): within active-TSS chromatin the faster removal remains (blood core +0.006, cortex +0.006-0.007 per
  decade); the lower resting level is mostly chromatin (four fifths in blood). Next: Tet-deficient cells.

## RM14 — Natural experiments in chromatin: broken writers, knockouts, edits

**State:** waits until 2026-12-13 (search again for a cohort meeting the rule) · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): cell lines tell little about clocks; can public data show what
happens when chromatin is changed? Three sources, in order:
1. Blood methylation of people with a constitutional variant in a chromatin writer: Sotos (NSD1), Weaver (EZH2),
   Tatton-Brown-Rahman (DNMT3A), Kabuki (KMT2D). Martin-Herranz et al. 2019 found clocks accelerated in Sotos.
   Test: is the essential-gene shield weaker than in matched controls? Search GEO for the cohorts.
2. Mouse knockouts: Tet3 heterozygous brain (GSE223748, at hand after R12); Tet2 and Dnmt3a mutant blood.
3. Edited cells over time: CRISPRoff (Nunez et al. 2021) tracked imposed promoter methylation for months.
Rule: a cohort with at least 10 carriers and matched controls on one platform.
- Done: Tet3 heterozygous brain loses the shield against littermates (R15, R16).
- Done (R18): Sotos, NSD1 variants, Weaver, Kabuki and CHARGE carriers drift more in blood but keep the shield
  (F unchanged or more negative); KMT2D and CHD7 variants show nothing within the controls' ages.
- Next: the DNMT3A loss-against-gain pair below; a TET2 or TET3 human cohort if one exists.
- Checked and not run (2026-11-13): the DNMT3A pair is smaller than listed. GSE120428's blood holds 2 TBRS
  carriers, 1 person with the gain-of-function dwarfism and 2 controls (its fibroblasts are replicates of one
  patient); GSE128801 is one family, 5 TBRS carriers against 7 relatives. A gain arm of one cannot be read.
- Waits until a TET2, TET3 or DNMT3A cohort meeting the rule (10+ carriers, matched controls) is found.

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

- **GSE184211**: done within GSE223748 (all 661 of its samples are the consortium's human set; R12, R54).
- **GSE50498** (450k muscle, ages 18–89): a muscle cohort with exact ages and more donors than GTEx's 47.
  Like GTEx it says nothing of exercise, and satellite-cell divisions follow loading, injury and
  wasting, so muscle tests selection only in a cohort that records physical activity or training
  history beside age (to find; the fibre cohort GSE304045 records training but not age).
- **Mouse replication:** Thompson 2018 (GSE120132) done (`figures rrbs_tissues`): the shield in all six tissues,
  ordered by neither size nor growth by turnover. Done (R63, Meer 2018): the shield in all four tissues, cortex
  too (−0.09 SD), smaller in cortex and heart than in liver and lung, and growing with age only in the latter.
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

# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.
The order and the rule for each item are argued in `docs/RESEARCH_PLAN.md`.

## RM37 — Manuscript round v0.2: the authors' comments and the referees' concerns

**State:** waits until 2026-11-23 (on 2026-11-18 the v0.1 Doc had no comments and desk #48 no answer; a v0.2 release before he has read v0.1 would cross his reading) · **Owner:** lead

Round v0.1 (RM36, done 2026-11-16) is on Drive. v0.2 narrows the claims the simulated referees questioned (the 0 %
floor, faster loss resting on a pooled fit, placebo sets not matched on expression, the cancer cohort's design),
adds R77 (held loci protected in blood), R78 (no sex difference), R79 (no universal hotspots), R80 (the blood shield is mostly the floor) and R82 (in blood, mostly expression), and fills the authors' decisions as they come.

## RM42 — A provenance table for every comparison: same samples or not, platform, chromatin source

**State:** urgent (asked by antonkulaga 2026-11-18: "have a table for comparisons to clearly separate: where several
tracks (chromatin and methylation) is from same sample or study and where not, what are the arrays. Also, how
chromatin is decided by a tissue-specific atac assay, or from general ensembl track (and if general - from same cell
or tissue type or general) these things must be clearly shown in paper and taken into consideration in paper,
manuscript and comparisons") · **Owner:** lead

1. Inventory every reported analysis: methylation study and platform; second track (chromatin, expression, outcome)
   and its source; same samples, same study or other study; chromatin as own-tissue assay, reference epigenome of
   the same tissue, of another tissue, sorted cells, or none.
2. The table goes into the supplement (and a short version in Methods); every Results sentence that joins tracks
   from different studies says so; cross-platform and cross-tissue pairings are flagged in claims and wording.
   The artifact gets the same table, and each result there carries icons: same or different batch (chip);
   second modality (chromatin etc.) from the same samples or not; same cell or tissue type or not (antonkulaga,
   2026-11-18).
3. Comparisons that join another tissue's chromatin reference, or two platforms, are redone on matched sources
   where data exist, or reported as such.

## RM43 — The independent methods review of 2026-11-18: reruns on comparable data

**State:** urgent (antonkulaga 2026-11-18 asked for "a separate reviewer agent to check if we compare comparable stuff
and if our selection for chromatin tracks, cell/tissue type definitions, array processing methods are proper"; its
findings are defects against what the paper says, so decided) · **Owner:** lead

Text was corrected at once (abstract, headings, Lunnon and blood-minus-cortex caveats, Tet1/Tet3 heading, provenance
labels; claims C1, C2, C12). Reruns, in order:
1. **Processing guard (M1):** a sesame-mode run that mixes rebuilt and deposited cohorts in one contrast writes the
   processing of each side into its tables and refuses a cross-cohort difference unless both match.
2. **Blood minus cortex (H2):** both cohorts on the shared age range (21-85), EpiDISH fractions held in Hannum as
   neuron fraction is in Jaffe; reported side by side.
3. **Chromatin references (M2, M3):** Lunnon entorhinal cortex with E071, superior temporal gyrus with E072, cerebellum
   out of chromatin matching, BA9 accessibility for frontal cortex only; blood rerun with E030 alone and with CpGs
   whose state agrees between E062 and E030 (R74, R72, R58 too). Done: R104 (brain maps), R121 (blood: no conclusion
   depends on the map).
4. **Array recipe (M6, M7):** mouse series with SeSAMe's mouse steps (MM285: T strain mask; Mammal40 in mouse: S and H
   masks); sample QC with expected probes from the species mask, not the data, plus a bisulfite-conversion (GCT) and
   intensity cut; rerun BXD, Browder, Tet1, Dnmt3a, CR, Huang; check CH and the syndromes against their high GCT.
5. **Species shield (M5):** group gaps with chip intercepts: done, R118 (holds). Left: human groups with the human
   annotation.
6. **SATSA (M8):** done, R119.
7. **Lunnon (H1):** renormalise all tissues jointly from GEO's raw intensities; slide stays confounded with tissue, so
   the same-donor contrast is reported as batch-confounded whatever the result.
8. **Tet3 against Tet1 (M9):** done, R120: the arrays share 150 CpG positions (42 at promoters), so shared CpGs cannot
   carry the comparison; on the 314 genes both cover, the opposite signs hold.
9. **Missing values (best-practice check):** `FeatureMatrix` median-fills sporadic missing values (up to 5-10 % per
   CpG); rerun the headline gaps on complete CpGs only and report the change. Done, R117 (no change).

## RM41 — One preprocessing from raw IDATs for every array series, with QC and batch checks

**State:** urgent (in progress: 16 series rebuilt, reruns R93-R109 written up; asked by antonkulaga 2026-11-18: "use raw idat, preprocess data same way, check for batch
effects ... apply it to all samples and validate if things are done properly for existing one") · **Owner:** lead

1. **Environment:** R with Bioconductor `sesame` in `/data/uq2347/envs/sesame` (micromamba), its platform data cached.
2. **One recipe:** `scripts/sesame_preprocess.R`, `prepSesame("QCDPB")` (quality mask, channel inference,
   non-linear dye bias, pOOBAH detection masking, noob) for HM450, EPIC, EPICv2, MM285 and Mammal40; per-sample
   `sesameQC_calcStats`; chip and position from the IDAT names. Replaces methylprep's SeSAMe imitation (R61).
3. **Samples with raw IDATs on GEO:** Jaffe GSE74193, EPIC-Italy GSE51032, writer syndromes GSE74432 and GSE97362,
   clonal haematopoiesis GSE315366, Gill GSE165179, Endicott GSE197512, Sturm GSE179847, fibres GSE304045, Huang
   GSE279437, CR blood GSE282661, Sarni GSE324160, Tet1 GSE271375, mammal subseries GSE190665 and GSE199979, SJLIFE
   GSE314261 (5,013 EPIC, 70 GB, last).
4. **No raw IDATs, kept on deposited betas and labelled so:** Hannum GSE40279, Lunnon GSE59685, MESA GSE56046 and
   GSE56581, GSE68194, GSE41037 (27k), GTEx GSE213478 (IDATs controlled access, dbGaP), the mammal superseries
   GSE223748 (only some subseries carry IDATs). Human lifespan blood from IDATs: EPIC-Italy and SJLIFE.
5. **Comparable arrays:** within a comparison, one platform; across 450k, EPIC and EPICv2 only their shared probes
   (EPICv2 replicates collapsed to the cg id); mouse MM285 against Mammal40 by mm10 CpG position.
6. **QC and batch:** samples failing detection (< 90 % of probes) or bisulfite conversion dropped and listed; per
   series, principal components of betas and of control-probe intensities against chip, position, plate and the
   design (age, group): the variance each explains, and age-on-chip confounding as in the 2026-11-18 audit.
10. **Progress (2026-11-18 17:00):** R SeSAMe rebuilds done for GSE74193 adult controls (pilot R93), GSE197512,
    GSE179847, GSE282661, GSE324160, GSE271375, GSE315366, GSE279437, GSE190665, GSE199979, GSE304045, GSE51032,
    GSE74432, GSE97362; full GSE74193 running. QC: Mammal40 detection read over the probes a series measures (all
    mouse samples failed the whole-array cut); Tet1 keeps 19 of 24 standard arrays. `ENTROPY_AGE_BETAS=sesame` makes
    every cohort table read its rebuild (no-IDAT cohorts stay deposited); the GEO-series readers (cancer,
    syndromes, CR blood, clonal haematopoiesis, Tet1, Sarni) still need the same switch. Then each figure reruns
    on rebuilt betas and the manuscript and claims follow (item 8).
8. **Decided by antonkulaga (2026-11-18):** results from cohorts without raw IDATs (Hannum, Lunnon, MESA, GSE68194,
   GSE41037, GTEx, the mammal superseries) are kept, each carrying a warning that it rests on the authors'
   processing and was not rebuilt; every series with IDATs is rebuilt. After each rebuilt series lands, the
   manuscript (`manuscript/paper.md`, supplement) and the hypotheses (`docs/CLAIMS.md`,
   `docs/evidence_matrix.csv`) are updated with it, not batched.
9. **Comparisons across platforms (antonkulaga, 2026-11-18: "you also had mistakes where you compared samples on
   different arrays"):** every claim that sets one array's result against another's is redone on the CpGs both
   measure, or reported only as separate within-study estimates with their heterogeneity. Known cases: (a) Tet3
   (Mammal40, GSE223748) against Tet1 (MM285, GSE271375): both rescored on the mouse CpGs the two arrays share
   (by mm10 position); (b) mouse turnover, Huang (Mammal320, MM285 probes) against RRBS (Meer, Thompson, Petkovich)
   and BXD liver (Mammal40): within-platform estimates side by side, no pooled ranking across platforms; (c) the
   species array (R59) against human 450k blood (R80): stated as different arrays and CpG universes; (d) GTEx
   (EPIC) gap sizes beside the 450k cohorts': on the 450k-EPIC shared probes. Status: (a)-(d) open; (b) within
   Huang the turnover ordering itself vanishes within chip (R92), so no cross-platform ranking is drawn; the paper's
   RRBS-against-array ranking sentence removed (R100). (b) done. Within-chip
   dispersion-slope intervals were degenerate (±1e12): hypothalamus has one aged-young pair per shared chip, so its
   within-chip dispersion slope is 0 by construction; dropped there from the next run on (paired_only).
7. **Validation of what exists:** per cohort with both, deposited against reprocessed betas (per-probe and
   per-sample correlation) and each headline result rerun on the reprocessed betas, the change reported.

## RM40 — Replacement data for the chip-confounded designs

**State:** urgent (asked by antonkulaga 2026-11-18: "can you search for replacement samples?"; on 2026-11-19 he pointed to the external search's final list: fetch and analyse it now) · **Owner:** lead

1. **Mouse turnover within chips, from data on disk:** GSE223748's untreated mice have numeric ages and young
   (<= 6 months) and old (>= 18 months) animals on shared chips in liver (72/100), muscle (20/54), cortex (21/46),
   blood (17/18), striatum (6/18), cerebellum (6/11), lung (6/9); spleen (4/2), kidney (6/2), heart (3/3) thin; no
   colon. A per-tissue age slope with chip intercepts, then the turnover-class contrast.
2. **Renewing tissues from other series** (GEO search 2026-11-18; batch read from IDAT chip ids or FASTQ flowcells):
   GSE233734 colon RRBS, 3-28 months, n 16-18 per age, every flowcell holds every age (best); GSE238172 mammal
   array, sorted marrow (each cell type on one chip with all ages) and blood (one chip), check overlap with
   GSE223748 first; GSE200527 MM285 blood, both B6 chips hold young and old; GSE198638 intestinal stem-cell WGBS,
   3 lanes hold both ages; GSE161141 rat blood RRBS, every lane holds both. Unknown or weak: GSE213628, GSE129712,
   GSE236617, GSE225166. Arrays enter through RM41's SeSAMe recipe; RRBS and WGBS need a sequencing pipeline.
4. **Search widened (antonkulaga 2026-11-18: "start searching for the samples that actually will be good enough in
   terms of criteria and will test selection and other things that we had to reject"; "do we have any scmethylation
   data"):** candidates must have public raw IDATs (or per-CpG sequencing files), one array type per contrast, and
   the contrast spread across chips or flowcells, tabulated before anything is fetched. Needs, in order: a second
   brain cohort and sorted neurons; culture passages not confounded with batch; clonal expansions with methylation;
   tissue turnover within animals; CR over time; a second prospective cancer cohort; single-cell methylation with
   ages (none in the project yet: `fig3` is a simulation, R50-R53 are read-level, not cell-level).
5. **Matched modalities from the same samples (antonkulaga 2026-11-18: "datasets where people got several modalities
   together like methylation+rna-seq, or methylation+atac-seq ... extremely important"):** methylation with RNA-seq,
   ATAC-seq or histone marks from the same animals, donors, cells or cultures across ages or passages; open-access
   animal and cell data first; single-cell multi-omics too. These replace the reference maps from other people
   (Table S4's weakest column) and rank with the second brain cohort.
6. **Scope widened and a gap recorded (external search, handoff revision 2, 2026-11-19):** for matched modalities
   the study need not be about aging: any dataset with methylation plus RNA or chromatin from the same samples and an
   age per sample qualifies, IDAT-backed arrays first. Sequencing-based methylation (bisulfite, enzymatic, nanopore,
   PacBio) with a second modality from the same samples is a recorded gap: 0 of the 46 analyses in
   `docs/provenance.csv` have it (arrays: 1 of 36, GTEx). Every fetched series is audited first with the external
   agent's `geo_audit.py` (and `sra_flowcell.py` for sequencing), with the filter that will be analysed, and the
   output kept in `data/input/geo/audit/<GSE>/`. Fetched so far: GSE252501, GSE279509 (IDATs), GSE233734 (per-CpG
   files). Corrections to the first search are at the top of `docs/DATASET_SEARCH_2026-11-18.md`.
7. **Final external results, 2026-11-19:** fetch order from its manifest: GSE279509 (analysed, R109); GSE304754 with
   GSE233734 (colon and hippocampus RRBS of one mouse cohort); GSE197723 and E-MTAB-15859 (division against arrest in
   culture); expression for cohorts in hand after donor linkage is established (GSE179848, GSE197471, GSE56045,
   GSE56580) and the epidermis set; GSE283951 (breast, pre-diagnostic). Rows marked [P] are held until audited here.
   Sequencing-based methylation with a same-sample second modality remains a gap: 0 of 46 analyses, and nothing
   adequate found in the open. The external search's unfinished tracks are to be completed or scheduled, not read
   as negative results.
3. **CR blood over time:** no public series has CR and ad-lib blood at more than one age in shared batches
   (closest: GSE75495, CR at one old age, batch unknown; liver only: GSE92486, GSE89274). The CR 12-week aging
   contrast stays unresolved.

## RM15 — Directed demethylation outside promoters: the culture and cortex follow-ups

**State:** after RM19 · **Owner:** lead

Progress: first part done (`figures intergenic`, R14, `docs/extension.md`): in blood, solo-WCGW CpGs lose
methylation 2.7 times as fast as other open-sea CpGs (directed, not noise); none in cortex; the PCA does not single
it out. Left: the same test in the culture series (EPIC manifest needed for the tags; the loss should grow with
doublings), and what cortex's two strongly age-tracking components are.
## RM13 — Answer the prior work (novelty review)

**State:** waits until 2026-12-13 (left: item 3 needs a read-level culture series; search again) · **Owner:** lead

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

**State:** waits until 2026-12-14 (search again for a cohort with per-sample ATAC QC, or sorted bulk ATAC) · **Owner:** lead

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

**State:** blocked: desk #40 · **Owner:** lead

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

- Done (R66): SATSA rebuilt from 1,072 IDATs: within-person directed change −0.14 SD (deposited −0.002); the
  cross-sectional shield reappears. Next: the mortality cohorts (waits on Anton).
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
- Done: Tet3 heterozygous brain loses the shield against same-chip wild types (R15, R16).
- Done (R18): Sotos, NSD1 variants, Weaver, Kabuki and CHARGE carriers drift more in blood but keep the shield
  (F unchanged or more negative); KMT2D and CHD7 variants show nothing within the controls' ages.
- Next: the DNMT3A loss-against-gain pair below; a TET2 or TET3 human cohort if one exists.
- Checked and not run (2026-11-13): the DNMT3A pair is smaller than listed. GSE120428's blood holds 2 TBRS
  carriers, 1 person with the gain-of-function dwarfism and 2 controls (its fibroblasts are replicates of one
  patient); GSE128801 is one family, 5 TBRS carriers against 7 relatives. A gain arm of one cannot be read.
- Waits until a TET2, TET3 or DNMT3A cohort meeting the rule (10+ carriers, matched controls) is found.

## RM3 — Interventions: replicate the CR result

**State:** blocked: desk #28 · **Owner:** lead

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
## RM4 — Stochasticity in RNA

**State:** waits until 2026-12-14 (no RNA shield replicated: GTEx null, rat weak, Tabula Muris Senis null (R75);
the intervention and RNA-clock items test a shield not yet shown, so they wait for a reason to expect one) · **Owner:** lead

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
- Done (R75): Tabula Muris Senis bulk (GSE132040): no RNA shield and no turnover order; the mouse-atlas pattern
  does not replicate.
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

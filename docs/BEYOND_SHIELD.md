# Beyond the shield: what else the data at hand can answer

Asked by antonkulaga (2026-11-12): other than the chromatin shield, what could this project investigate with what
it has already built and downloaded, to learn about health and aging? A fresh read of `docs/`, the module list
and the data on disk (2026-11-12). Every idea below reuses on-disk data and existing code; none needs a download
unless said. Ideas the data cannot support are listed at the end with the reason.

The reusable machinery: per-sample drift `B`, signed shift `S` and entropy `E` against a reference
(`entropy_age.protection`, `per_person`, `consortium_interventions.arm_scores`, `health_states`); 17 published
clocks at a pinned commit (`data/clocks.py`: Horvath 1/2, Hannum, PhenoAge, DunedinPACE, Zhang, Lin, Cortical,
DNAmTL, EpiTOC1, StocH/P/Z, Ying CausAge/DamAge/AdaptAge); own elastic-net clocks (`own_clocks`); 12-type blood
deconvolution (`data/input/epidish`); Roadmap chromHMM for blood, cortex and others (`data/chromatin`); tissue
contrasts with bootstraps (`tissue_contrast`, `human_tissues`, `mouse_tissues`); GTEx v10 RNA read from rna-age
(`data/gtex_rna.py`) and the rodent RNA atlases (`rna_rodents`).

Novelty notes are from memory unless a PMID is given in `docs/NOVELTY_REVIEW.md`; each is marked **unverified**
where not checked this session. A literature check comes before any of these is written up.

## Ranked list (insight × feasibility)

### 1. One body, one aging rate? Coupling of age acceleration across tissues of the same individual

- **Question.** If a person's (or mouse's) blood is epigenetically "old for its age", are their other organs too?
  Or does each organ age at its own pace?
- **Why it matters.** Blood is the only tissue sampled in living people. Every blood clock used as a health
  marker assumes that blood speaks for the body. Strong coupling supports that; weak coupling means organ-specific
  aging and limits what a blood test can say about the brain.
- **Data.** Lunnon GSE59685 (`lunnon_*`): whole blood plus entorhinal cortex, superior temporal gyrus, frontal
  cortex and cerebellum, 67 donors with all five, 44 with four. GTEx GSE213478 (`gtex_*`): 424 donors, 272 with
  2-7 tissues (colon-lung 153 pairs, colon-prostate 85, lung-blood 37). Huang GSE279437 (`huang_*`): 11 tissues of
  the same 72 mice on one array (about 42 aged, 30 young).
- **Code.** Published clocks (`data/clocks.py`, applied as in `shield_clocks` and `clock_cpgs`, which already
  computed Lunnon's brain-minus-blood offset); per-sample `B` and `E` (`per_person.score_against`); own per-tissue
  mouse clocks (`own_clocks`) since no mouse clock is on disk.
- **Analysis.** (1) Per tissue, residual of each clock on age (and of `B`, `E`), plus covariates (sex, plate,
  composition components). (2) Within-donor correlation matrix of residuals across tissues, with the bootstrap
  over donors. (3) A one-factor model: share of each tissue's acceleration explained by a donor-wide factor.
  (4) In mice, the same within the aged group only (age is fixed, so residuals are pure individual variation).
- **Effort.** 10-14 h (Lunnon and GTEx 6 h; mouse clocks 4-6 h; one Slurm job for clocks on the large betas).
- **Informative either way.** r ≥ 0.3 between blood and brain: blood clocks carry brain information. r ≈ 0 with
  tight intervals: organ-specific aging, a direct caution for blood biomarkers. A shared factor in mice but not in
  people would point at genetic/environmental heterogeneity versus a systemic signal.
- **Caveats.** Lunnon blood was drawn before death (`age.blood` and `age.brain` differ in the series matrix; the
  samples table keeps one age): residuals must use each sample's own age. Lunnon is AD-enriched (adjust for
  status and Braak). GTEx ages are decade bins (residual noise up to ±5 y). Clocks are blood- or pan-tissue
  trained; use `B`/`E` (clock-free) beside them.
- **Novelty risk.** Moderate. Hannon et al. 2015 (Epigenetics) compared blood and brain methylation per CpG in
  the same London cohort; organ-specific aging rates are published from proteomics (Oh et al. 2023, Nature) and
  imaging (Tian et al. 2023, Nat Med); within-donor correlation of epigenetic age across GTEx tissues may be in
  the GTEx methylation paper (Oliva et al. 2023) or later work (**unverified**). A three-dataset answer with
  clock-free drift and a mouse arm at fixed age is likely new.

### 2. Does illness before death age the methylome? GTEx's Hardy scale as a health contrast

- **Question.** Donors who died after a long chronic illness (Hardy 4, 68 donors with methylation) against
  those who died suddenly (Hardy 1-2, 132): at equal age, are their clocks, drift `B` and entropy higher, and in
  which of eight tissues?
- **Why it matters.** It is the only health readout on disk that spans many human tissues of the same people.
  Which tissues show systemic illness, and whether stochastic drift or directed clocks register it, says which
  measure reads health rather than time.
- **Data.** GTEx methylation (`gtex_samples.parquet`, 424 donors) joined to `DTHHRDY` from
  `rna-age/data/download/gtex/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt` (0: 200, 1: 22, 2: 110,
  3: 24, 4: 68). Smoking (168 current, 128 never) is a second exposure in the same run.
- **Code.** `health_states.py` already scores GTEx current against never smokers with leave-one-out references;
  add Hardy as a contrast. Clocks as in idea 1.
- **Analysis.** (1) Per tissue, `B`, `E`, `S` and clock residuals against Hardy 4 vs 1-2, with age bin, sex,
  autolysis, plate, ischemic time (sample attributes) and composition components. (2) Ventilator cases (Hardy 0)
  as a separate arm, not pooled. (3) Same for smoking; ask whether smoking moves directed sites only (AHRR-type)
  or raises drift everywhere.
- **Effort.** 4-6 h (one short job).
- **Informative either way.** A chronic-illness effect in drift but not in clocks (or the reverse) separates
  noise from programme as health readouts; no effect at all says terminal illness does not reach the methylome
  within months, which bounds how fast these markers can respond.
- **Caveats.** Hardy correlates with ischemic time and autolysis; both must be covariates. Cross-sectional: an
  effect may be treatment or agonal.
- **Novelty risk.** Low to moderate. Hardy is commonly a nuisance covariate in GTEx expression work; as an
  exposure for methylation age and drift **unverified**, likely little done.

### 3. Noise or programme: what do interventions actually reverse?

- **Question.** When an intervention makes a clock younger, does it reduce stochastic drift (between-animal
  spread, distance from the young state at drifting CpGs) or shift directed CpGs (mean moves toward young)?
- **Why it matters.** Rejuvenation claims rest on clocks. If interventions move only the directed part, they are
  resetting a regulated programme and not repairing accumulated disorder; this changes how to read every
  intervention trial.
- **Data.** Huang GSE279437: exenatide (13 aged vs 14 aged controls, 8 young treated) and rapamycin (5) in 11
  tissues; Petkovich GSE80672: CR, GHRKO, Snell blood; Browder GSE190665: cyclic OSKM in 6 tissues; BXD GSE199979:
  high-fat diet and body weight in liver; mammal consortium mouse arms (GSE223748).
- **Code.** `mouse_tissues`, `interventions`, `reprogramming`, `consortium_interventions.arm_scores`,
  `bxd_liver.ols`; the directed/stochastic split of `stochastic_component` (Tong-style rates) and `intergenic`.
- **Analysis.** (1) Per tissue, split age-associated CpGs into directed (consistent mean change, low dispersion
  growth) and stochastic (dispersion growth, mean small), from young vs aged controls. (2) For each arm, effect
  on the mean of directed CpGs (toward young?) and on dispersion/`B` at stochastic CpGs. (3) Compare across
  tissues: does exenatide act in the same tissues in both components?
- **Effort.** 8-10 h (all data parsed; one job).
- **Informative either way.** Directed only: interventions retune, not repair. Both: clocks' rejuvenation is
  partly real repair of disorder. Already seen in pieces: CR leaves `B` unchanged (R8), OSKM leaves `B` unchanged
  (R32), which leans toward "directed only".
- **Caveats.** Rapamycin and KD arms are 5 animals; read exenatide and CR as the main arms.
- **Novelty risk.** Moderate. Tong et al. 2024 (PMID 38724732) report stochastic clocks under reprogramming and
  disease; Huang et al. (data owners) presumably report clocks per tissue (**unverified**); the explicit
  directed-vs-stochastic split across 11 tissues of the same animals is likely new.

### 4. Alzheimer's: does disorder follow the pathology's regional order?

- **Question.** In the same AD brains, Braak staging says entorhinal cortex is hit first, then temporal, then
  frontal; cerebellum is spared. Do drift `B`, entropy and clock acceleration rise with Braak stage in that
  order, and does blood carry any of it?
- **Why it matters.** It separates "the disease ages tissue" (effect graded by local pathology) from "old-for-age
  people get AD" (a donor-wide effect, the same in every region and in blood). That is a causal-direction hint
  from a cross-sectional cohort.
- **Data.** Lunnon GSE59685: about 60 AD and 21-26 controls per region, Braak 0-VI; blood 48 AD, 9 controls.
- **Code.** `health_states.py` (built for RM24; contrasts AD and Braak per region); clocks including the Cortical
  clock (`DNAmClockCortical.csv`); idea 1's donor factor.
- **Analysis.** (1) Per region, slope of `B`, `E` and clock residuals on Braak, with age, sex, slide and neuron
  proportion (the cortex composition must be estimated or taken as components). (2) Test the regional gradient
  of the slopes (EC > STG > FC > CER) by a donor-level mixed model. (3) Add the donor-wide factor of idea 1: does
  Braak predict the region-specific part or the shared part?
- **Effort.** 4-6 h after RM24's run.
- **Informative either way.** Graded by region: pathology drives local disorder. Flat across regions and present
  in cerebellum: a donor trait. Neuron loss is the main confounder; a graded effect that vanishes with neuron
  fraction is composition.
- **Novelty risk.** High for clocks: Levine et al. 2015 (Aging) linked prefrontal epigenetic age to neuritic
  plaques in ROSMAP; Shireby et al. 2020 built the Cortical clock partly on these cohorts (**unverified** details).
  The regional-gradient and drift/entropy angle is less covered.

### 5. Which tissues age fastest, and is the order the same in human, mouse and rat, in methylation and RNA?

- **Question.** Rank tissues by drift rate per year (and per lifespan fraction) in methylation and by
  transcriptional dispersion growth in RNA, and test whether the ranking is conserved across species and between
  the two layers.
- **Why it matters.** A conserved order (e.g. colon and blood fast, brain slow) points at a shared cause
  (division, metabolism) and says which organs to watch first; a non-conserved order says tissue aging rate is
  species-specific.
- **Data.** Tables already in `data/output/`: `gtex_tissues`, `mouse_tissues`, `mouse_rrbs_tissues`,
  `rna_gtex` (47 tissues), `rna_rodents` (24 mouse, 28 rat tissues). Look-alike drift `B` per tissue is already
  estimated for the shield's matched null.
- **Code.** None new beyond reading the point-estimate tables; rank correlations with bootstrap intervals from the
  `bootstrap.parquet` files.
- **Analysis.** (1) Extract each tissue's background (non-gene-set) drift slope from each table. (2) Map tissues
  across species (about 8-12 shared). (3) Spearman of the orders, methylation vs RNA within species, and mouse vs
  rat vs human within layer, with bootstrap intervals.
- **Effort.** 2-4 h (mostly tables at hand).
- **Informative either way.** A conserved order is a cheap, citable result; discordance between methylation and
  RNA says the two layers age by different rules.
- **Novelty risk.** Moderate. Schaum et al. 2020 (Nature, Tabula Muris Senis bulk) timed tissue transcriptome
  aging in mouse; Vaidya et al. 2023 (PMID 36797759) tied methylation entropy to division rate across mouse
  tissues. The cross-species, cross-layer comparison is likely new (**unverified**).

### 6. When does disorder begin? The drift trajectory from fetus to old age, and its timing across species

- **Question.** Does between-person drift start at birth, at maturity, or later, and does it grow linearly with
  age? Across species, does drift onset scale with age at sexual maturity?
- **Why it matters.** If disorder accumulates from maturity at a constant rate, aging is a post-developmental
  process with a measurable start; that sets when prevention could begin.
- **Data.** Jaffe GSE74193: 36 fetal, 94 under 20, 404 adults, prefrontal cortex. Petkovich: mouse blood from
  0.67 to 35 months. Mammal consortium GSE223748: age, max lifespan and **age at sexual maturity** per species
  (`ageatsexualmaturityyears` field).
- **Code.** `atlas_real_data` (per-CpG drift, entropy and dispersion slopes), `per_person`, `species_shield`.
- **Analysis.** (1) In Jaffe, per-sample distance from a young-adult reference at background CpGs, fitted with a
  segmented or spline model in age; locate the change point. (2) Same in mouse blood. (3) Across species,
  dispersion slope before vs after sexual maturity, against log lifespan, with `species_phylo`'s covariance.
- **Effort.** 6-10 h.
- **Informative either way.** A kink at maturity supports a developmental-to-aging transition; a continuous rise
  from fetal life says disorder begins with development.
- **Novelty risk.** Moderate. Horvath 2013's log-linear age before 20; Kerepesi et al. 2021 (Sci Adv) "ground
  zero" of rejuvenation in early embryogenesis; Lu et al. 2023 universal clocks use relative age with maturity.
  Drift/entropy (not clock) onset is less studied (**unverified**).

### 7. Does epigenetic disorder reach transcription? Methylation-expression coupling with age in the same donors

- **Question.** In GTEx donors with methylation and RNA from the same tissue, do genes whose promoter CpGs drift
  most show the largest rise in expression variance with age, and does a donor's drift `B` predict their
  expression noise?
- **Why it matters.** Drift matters for health only if it changes what cells do. A link from methylation noise to
  transcriptional noise makes drift a mechanism; no link makes it a passive record (a clock, not a cause).
- **Data.** Same-tissue methylation and RNA donors: lung 161, colon 139, ovary 133, prostate 62, whole blood 45,
  muscle 37, breast 36, testis 31 (counted from the sample attributes, RNASEQ freeze). GTEx RNA via
  `data/gtex_rna.py` (rna-age folder).
- **Code.** `rna_gtex` (residualised log CPM with Hardy, RIN, ischemic time), `atlas_real_data` per-CpG slopes,
  `human_tissues` matching.
- **Analysis.** (1) Per gene, promoter dispersion slope (methylation) and residual-variance slope with age (RNA)
  in lung, colon, ovary. (2) Correlate across genes, against matched background. (3) Per donor, `B` against a
  donor-level expression-noise score (mean absolute residual of tissue-typical genes).
- **Effort.** 10-16 h (joining two pipelines; memory-heavy, a Slurm job).
- **Informative either way.** Coupling: drift has transcriptional consequences. None: methylation disorder is
  mostly silent, which itself argues clocks read a by-product.
- **Novelty risk.** Moderate. Oliva et al. 2023 (GTEx methylation, Cell Genomics) mapped eQTMs; age-related
  dispersion coupling between the layers in the same donors is **unverified**, possibly new.

### 8. Time, divisions or metabolism: what sets the rate of drift in cells?

- **Question.** In cultured human cells, how much of drift and of each clock's ticking is per doubling, per day,
  and how much is raised by mitochondrial dysfunction, glucocorticoids or hypoxia?
- **Why it matters.** It tells whether stress and metabolism accelerate aging markers directly or only by
  changing division rate, a question behind every lifestyle-and-clock association.
- **Data.** Sturm GSE179847 (`sturm_*`): 9 donors, SURF1 (Leigh syndrome) lines, 19 agents (oligomycin,
  dexamethasone, contact inhibition, hypoxia...), doublings and days per sample. Endicott GSE197512: MMC arrest,
  serum withdrawal, 3 % oxygen, TERT. BJ GSE268145.
- **Code.** `passaged_cells`, `stochastic_component`, clocks (EpiTOC1 is a mitotic clock; DNAmTL).
- **Analysis.** (1) Per culture, `B`, `E` and each clock regressed on doublings and days jointly (lines as
  random effects). (2) Treatment interactions: does oligomycin or SURF1 raise drift per doubling or per day?
  (3) Rank clocks by time-vs-division loading.
- **Effort.** 6-8 h.
- **Informative either way.** Per-day drift raised by mitochondrial stress at equal divisions: metabolism ages
  the methylome directly. Purely per-doubling: stress acts through proliferation.
- **Novelty risk.** Moderate to high. Sturm et al. 2023 (Commun Biol) reported accelerated clocks under OxPhos
  defects in this series; Endicott et al. 2022 (PMID 36347867) separated time and division for late-replicating
  loss (**unverified** beyond the review's notes). The drift/entropy decomposition across all agents is the new part.

### 9. Progeroid and machinery mutants: accelerated aging or a different pattern?

- **Question.** Do the methylation changes of Dnmt3a W326R/+ progeria (bone marrow, intestine, liver), Tet3+/-,
  Tet1 KO and the human writer syndromes point in the same direction as natural aging at the same CpGs?
- **Why it matters.** Progeroid models are used as aging proxies; a weak concordance says they are not.
- **Data.** GSE324160, GSE223748 Tet3 arm, GSE271375, GSE74432 and GSE97362; natural age slopes from Thompson
  (6 tissues), Huang (11 tissues) and Petkovich.
- **Code.** `dnmt3a_gain`, `tet1_cortex`, `syndromes`, `mouse_tissues` slopes.
- **Analysis.** (1) Per tissue, mutant-minus-control effect per CpG. (2) Correlate with the natural age slope in
  the matched tissue; split directed and stochastic CpGs as in idea 3. (3) Express as "years of aging
  equivalent" by projection.
- **Effort.** 4-6 h.
- **Novelty risk.** Moderate to high: Sarni et al. 2026 link Dnmt3a hypermethylation to age-related pathology
  and probably test concordance (**unverified**); Martin-Herranz et al. 2019 found Sotos clock acceleration.

### 10. Sex differences in drift, across tissues and species

- **Question.** At equal age, do females drift more slowly than males, in which tissues, and is the sign the
  same in mouse and across mammals?
- **Why it matters.** Women live longer; whether their methylome is less disordered, or only their clocks
  slower, is open.
- **Data.** Hannum, Jaffe, GTEx (non-sex-specific tissues), Lunnon; Thompson mice (both sexes, 6 tissues); rodent
  RNA atlases (both sexes); mammal consortium (`female` field). Autosomal CpGs only.
- **Code.** `per_person` (`B`, `E`), `species_shield` grouping.
- **Analysis.** (1) Sex term on `B`, `E` and clock residuals per cohort and tissue, with composition. (2)
  Meta-analyse; (3) across species, the sex difference in drift slope.
- **Effort.** 4-6 h.
- **Novelty risk.** Moderate: men's faster clocks are known (Horvath et al. 2016, Genome Biol); drift and the
  cross-species comparison are less covered (**unverified**).

### 11. Universal noise hotspots: CpGs that become disordered in every mammal

- **Question.** Which conserved CpGs gain between-individual variance with age in most species and tissues, and
  what genes and functions do they mark?
- **Why it matters.** Universally fragile loci are candidate weak points of mammalian aging and targets for
  biomarkers that transfer across species (from mouse trials to people).
- **Data.** GSE223748 (15,043 samples, 349 species); human (Hannum, GTEx) and mouse sets for the overlap.
- **Code.** `species_shield` (per-group slopes), `atlas_real_data` (per-CpG dispersion), `data/mouse.py`
  annotation.
- **Analysis.** (1) Per species-tissue group, per-CpG dispersion slope on relative age. (2) Count groups with a
  positive slope per CpG; rank, with phylogenetic weighting. (3) Gene and chromatin-state enrichment.
- **Effort.** 8-12 h (large matrix; a Slurm job).
- **Novelty risk.** Moderate to high: Haghani et al. 2023 (Science) and Lu et al. 2023 (universal clocks, PRC2
  targets) did mean changes; variance across species is less done (**unverified**).

### 12. Schizophrenia: blood and brain

- **Question.** Clock, drift and entropy in schizophrenia: blood (GSE41037) and prefrontal cortex (Jaffe).
- **Why it matters.** Schizophrenia shortens life by 10-20 years.
- **Data/code.** `health_states.py` (RM24) already covers both; add clocks and `E`.
- **Effort.** 2-3 h after RM24.
- **Novelty risk.** High: McKinney et al. 2017 (Jaffe cortex, no acceleration) and Higgins-Chen et al. 2020
  (Biol Psychiatry, blood) (**unverified** details). Low ranking because different people provide blood and brain.

## Dropped: what the data cannot support

- **Mortality or incident disease prediction.** No cohort on disk follows people forward (RM24 says so).
- **Strain lifespan in BXD.** GSE199979 has weight and diet but no strain lifespan; Mozhui et al. 2022 (the data's
  own paper) already report weight and diet effects on clocks.
- **Blood cell composition with age or disease as a result.** Well published (immune aging); keep it as a
  covariate.
- **Cell-type aging from the 362-donor cortex snATAC.** RM19 found accessibility spread dominated by sample
  quality; cell proportions with age are the data owners' own result.
- **Ethnicity in Hannum.** Two groups by site; Horvath et al. 2016 covered it.
- **Clonal haematopoiesis beyond R25.** No allele fractions; cardiogenic-shock patients only.
- **Long-read epialleles.** Not downloaded (RM22's WGBS set is 13 GB).

## Top 8 at a glance

| idea | data | effort | novelty risk |
|---|---|---|---|
| 1 Cross-tissue coupling of aging rate | Lunnon 5 tissues; GTEx 2-7 tissues; Huang 11 mouse tissues | 10-14 h | moderate |
| 2 Hardy (illness before death) and smoking | GTEx methylation + DTHHRDY | 4-6 h | low-moderate |
| 3 Interventions: noise or programme | Huang, Petkovich, Browder, BXD | 8-10 h | moderate |
| 4 AD: disorder by Braak region | Lunnon 4 regions + blood | 4-6 h | high (clocks), moderate (drift) |
| 5 Tissue aging order across species and layers | existing output tables | 2-4 h | moderate |
| 6 When drift begins (fetal to old; maturity) | Jaffe, Petkovich, mammals | 6-10 h | moderate |
| 7 Methylation noise to expression noise | GTEx same-donor meth + RNA | 10-16 h | moderate |
| 8 Time vs divisions vs metabolism in culture | Sturm, Endicott, BJ | 6-8 h | moderate-high |

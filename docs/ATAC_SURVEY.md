# ATAC-seq data for stochastic chromatin opening with age (roadmap RM19)

Survey of public data, 2026-11-12, by a delegated agent for the lead; accessions checked against GEO and Zenodo where stated, unverified items marked. Not yet a result.

# Public ATAC-seq data for testing age-related accessibility drift and a housekeeping "shield"

Survey date: 2026-11-12. Report only; nothing was changed in any repository.

Question: within one tissue/cell type and one protocol, does per-region variability of accessibility
(between individuals, or between cells) grow with age, and does it grow less at promoters of
essential/housekeeping genes than at matched regions?

How things were verified (labels used below):
- **[GEO]** read directly from the GEO FTP listing and series-matrix header through E-utilities/FTP on 2026-11-12
  (file names, sizes, sample characteristics, ages counted from the matrix).
- **[file]** a small file was read (Zenodo demographics table; first lines of a GEO count matrix).
- **[paper]** read from the paper's full text (PMC) via fetch.
- **[secondary]** from a search summary, abstract, or news item only; treat as unverified until the paper is read.
- **[unverified]** stated from memory or partial evidence; check before relying on it.

---

## Ranked candidates

### 1. Human prefrontal cortex snMultiome aging cohort (Catching et al. 2026)

- **Publication:** Catching A, Weller CA, Hu F, et al. "Single-nucleus multiome analysis in the human prefrontal
  cortex identifies gene expression and cis-regulatory elements associated with aging." *Cell Reports* 2026,
  PMID 41832957 (PubMed search hit), DOI 10.1016/j.celrep.2026.117110 [secondary].
- **Accession:** Zenodo record **18394349** (open access) [file]. Title on Zenodo: "Single nucleus multiome analysis
  in the human pre-frontal cortex identifies cell-type restricted gene expression and cis-regulatory elements
  associated with aging"; creators Reed X, Weller C, Catching B; published 2026-01-28.
- **Files:** `final_atac_data.h5ad` **14.5 GB**, `final_rna_data.h5ad` 18.4 GB, `TableS1_CohortDemographics.xlsx`,
  `TableS2_NucleusCountandPercentagebyCohort.xlsx`, `TableS3.Significant_genes_and_peaks_age.xlsx` [file].
  Whether the ATAC h5ad has raw cell-by-peak counts (rather than only normalized/LSI values) is **not verified**.
  The paper reports >521,000 peaks, ~1.5 M nuclei [secondary, via the 2026 bioRxiv clocks preprint].
- **Organism / tissue:** human, postmortem dorsolateral/prefrontal cortex; neurologically unaffected donors.
  Cell types (annotation table): astrocytes, excitatory neurons, inhibitory neurons, microglia,
  oligodendrocytes, OPCs, vascular cells [file].
- **Individuals:** TableS1 lists **362** donors (the paper text says 357 [secondary]) [file].
- **Ages:** exact (one decimal for some), **15–100 years** [file]. Per decade: ≤20: 47, 21–30: 48, 31–40: 48,
  41–50: 80, 51–60: 56, 61–70: 20, 71–80: 16, 81–90: 31, >90: 16 [file].
- **Sex:** 230 male, 132 female [file].
- **Other covariates in TableS1:** PMI, ethnicity, race, brain bank (HBCC 158, UMARY 148, SH 29, UKY 27) [file].
  Brain bank is partly collinear with ancestry (HBCC contributes the Black donors) and should be a covariate.
- **Protocol:** 10x Genomics Single Cell Multiome ATAC + Gene Expression [secondary].
- **Why rank 1:** the only cohort found that is open, large (>350), spans adolescence to centenarians in one tissue
  and protocol with exact ages, and resolves cell types. It supports both designs: (a) between-individual
  variance of per-donor pseudobulk accessibility per cell type vs age; (b) between-cell variability within donor
  vs age. The matched RNA in the same nuclei lets the "shield" be checked against expression level.
- **Caveats:** postmortem tissue (PMI, agonal state); per-donor nuclei counts per cell type vary; one tissue only;
  12 old-age donors per decade above 60 is thinner than below; h5ad content (raw counts? per-donor labels?)
  to be confirmed by reading the `.obs`/`.var` of the file.

### 2. Human PBMC bulk ATAC-seq aging clock cohort (Morandini et al. 2023)

- **Publication:** Morandini F, Rechsteiner C, Perez K, ... Ocampo A. "ATAC-clock: An aging clock based on chromatin
  accessibility." *GeroScience* 2023 (PMID 37924441, linked from the GEO series) [GEO][secondary].
- **Accession:** **GSE193140** (ATAC-seq subseries; SuperSeries GSE193142, RNA-seq GSE193141; BioProject
  PRJNA795127) [GEO].
- **Organism / tissue:** human, bulk PBMC [GEO].
- **Samples / individuals:** 157 ATAC samples [GEO]; column names like `CR_027_1`, `CR_027_2` show some donors have two
  libraries, so the donor count is somewhat below 157 (the abstract says blood from 159 donors overall)
  [file][secondary]. Exact donor count for ATAC **not verified**.
- **Ages:** exact integer ages in the series matrix, **20–74 years**; per decade: 20s: 20, 30s: 39, 40s: 29, 50s: 36,
  60s: 26, 70s: 7 [GEO].
- **Sex:** 115 M, 42 F [GEO].
- **Protocol / platform:** ATAC-seq, Illumina NovaSeq 6000 [GEO].
- **Processed matrix:** **yes** — `GSE193140_atac_counts_qval_0.05.txt.gz`, **9.8 MB**, space-delimited,
  ~42,557 peaks × 157 samples of raw integer counts (columns `chr start end` + samples) [GEO][file].
  It covers only peaks passing q < 0.05, so very weak/rare peaks are absent; that suits a promoter-level test.
- **Cell composition:** the abstract says cell-composition data were generated [secondary]; whether it is on GEO
  is **not verified** (GSE193141 has a CSV supplementary file; not opened).
- **Why rank 2:** smallest and cheapest first test (one 10 MB file), continuous exact ages, >150 libraries.
  Bulk PBMC confounds variance with cell-type composition (CD8 memory expansion, naive T loss with age), so
  composition (or a deconvolution estimate) must be regressed out before variance is compared.

### 3. Human hippocampus snMultiome across adult life (Zemke et al. 2024)

- **Publication:** "Epigenetic and 3D genome reprogramming during the aging of human hippocampus", PMID 39463924
  (GEO-linked); contributors Zemke NR, Lee S, ... Ren B [GEO]. Journal/year: 2024 bioRxiv per PMID date
  [unverified whether peer-reviewed].
- **Accession:** **GSE278576** [GEO].
- **Organism / tissue:** human postmortem hippocampus; cell types include DG, CA1, CA2–3, astrocytes, chandelier,
  endothelial and others (from bigWig names) [GEO].
- **Individuals:** **40** donors, each with ATAC and RNA libraries (80 GSM) [GEO].
- **Ages:** exact, **20–95 years**, read from the cell metadata file: 20, 20, 25, 26, 28, 28, 31, 33, 38, 38, 41, 41,
  43, 44, 46, 46, 48, 50, 54, 55, 65, 66, 68, 69, 69, 71, 75, 75, 75, 79, 82, 83, 86, 86, 87, 91, 91, 92, 94, 95
  [GEO file]. Age groups used by the authors: 20–40, 40–60, 60–80, 80–100 (~10 per group).
- **Sex:** 19 M, 21 F (from the same metadata) [GEO file].
- **Protocol:** 10x Single Cell Multiome ATAC + Gene Expression [GEO].
- **Processed data:** per-donor `*_raw_feature_bc_matrix.h5` (10x output, includes the Peaks feature type;
  49–83 MB each, ~2.5–3 GB for 40), RNA Seurat object 3.1 GB, cell metadata 12 MB (with Age, Gender,
  subclass), bigWigs per cell type and per age bin [GEO]. Per-donor peaks from Cell Ranger ARC are called per
  library, so a common peak set must be built (or counts re-quantified from fragments, which are in the 52 GB
  RAW tar only if included — **not verified**).
- **Why rank 3:** clean, evenly spread exact ages 20–95 with cell-type resolution; moderate n (40).

### 4. Mouse organism-wide scATAC aging atlas (Lu, Zhou, ... Cao 2025)

- **Accession:** **GSE288730** [GEO]; PMID 40463164 (GEO-linked). Title: "Organism-wide cellular dynamics and
  epigenomic remodeling in mammalian aging".
- **Organism / tissues:** C57BL/6 mouse, 21 tissues (BAT, bone marrow, cecum, colon, esophagus, eye, heart, intestine,
  kidney, liver, lung, muscle, ovary/uterus, pancreas, skin, spleen, stomach, thymus, three WAT depots) [GEO].
- **Individuals / ages:** 32 mice — 1 month (n=10, 5M/5F), 5 months (n=10, 5M/5F), 21 months (n=12, 6M/6F) [GEO].
  Ages exact by design (three time points).
- **Protocol:** EasySci-ATAC-seq (combinatorial indexing), ~7 M cells, ~1.3 M cCREs [GEO].
- **Processed data:** per-tissue `*_peak_count.h5ad` cell-by-peak matrices, **1.3–37 GB each** (liver 9.1 GB,
  kidney 37 GB, heart 2.9 GB, muscle 10 GB), plus gene-count h5ad; cell metadata 121 MB [GEO].
- **Fit:** between-cell variability within a cell type and per mouse, per tissue, with ~10 mice per age; the
  1-month group is juvenile, so "young adult vs old" is 5 vs 21 months only. Sparse sci-ATAC per-cell counts make
  cell-level variance estimates noisy; per-mouse pseudobulk gives n=10–12 per age.
- A related Cao-lab series **GSE331488/GSE331489** ("A global view of cell-type-specific temporal dynamics in
  mammalian aging", 614 + ~725 ATAC GSM, per-tissue peak h5ad 2–55 GB) has ages "1, 5 and 21-month-old" and sex
  mostly "pooled male and female" in the metadata [GEO]; individual-level structure **not verified**.

### 5. Mouse liver lifespan snMultiome (5 ages, both sexes)

- **Accession:** **GSE328206** (SuperSeries GSE341533, BioProject PRJNA1454349) [GEO]. Title: "A lifespan-resolved
  multi-omic study reveals dynamic transcription factor networks and the regulatory logic of liver aging".
  Publication **not found/verified**.
- **Design:** C57BL/6JN liver, 40 multiome libraries (ATAC + RNA), age groups Y, M, O, P, G (8 each: 4M/4F),
  two batches [GEO]. Exact ages per group **not verified** (labels only).
- **Processed data:** `GSE328206_objects.tar.gz` 26 GB (analysis objects), RAW tar 144 GB [GEO].
- **Fit:** one tissue, one protocol, 5 ages × 8 animals; inbred mice remove genetic variance, so between-individual
  variance with age is cleaner than in humans. Small n per age.

### 6. Human PBMC ATAC, sex-dimorphic immune aging (Márquez et al. 2020)

- **Publication:** Márquez EJ, ... Ucar D. "Sexual-dimorphism in human immune system aging." *Nat Commun* 11:751
  (2020) [paper].
- **Accession:** dbGaP **phs001934.v1.p1** (all PBMC ATAC/RNA, controlled) and EGA **EGAS00001002605** (subset)
  [paper]. No GEO series.
- **Individuals:** 120 profiled by ATAC (54 M, 66 F), 100 after QC [paper]; ages 22–93 (cohort); age groups in the
  full cohort: young 22–40 (54), middle 41–64 (59), older ≥65 (59) [paper]. Exact ages are in dbGaP phenotype
  files [unverified].
- **Protocol:** Omni-like ATAC on 50,000 nuclei, Nextera Tn5, HiSeq 2500, ≥30 M reads, ≥2 technical replicates
  merged; hg19, MACS2 [paper].
- **Processed matrix:** **not public** as far as found; Source Data and the Shiny app (immune-aging.jax.org) give
  derived results only [paper]. Requires a dbGaP application and realignment-free counts only if dbGaP holds them
  [unverified].
- **Fit:** ideal design (n=100, 22–93, one protocol), blocked by controlled access.

### 7. Human PBMC and sorted immune ATAC with age (Ucar et al. 2017)

- **Publication:** Ucar D, et al. "The chromatin accessibility signature of human immune aging stems from CD8+ T
  cells." *J Exp Med* 214:3123 (2017), PMID 28904110 [paper].
- **Accession:** EGA **EGAS00001002605** (controlled) [paper]. (OmicsDI also lists GSE98983/GSE98984 — these are
  unrelated HepG2 NF-κB series [GEO]; do not use.)
- **Individuals:** PBMC ATAC 49 (28 young 22–40, 21 old ≥65), 44 after QC (25/19); sorted naive/memory CD4 and
  CD8 T cells from 8 donors (4/4); monocytes n=20, naive B n=7 (age split not given) [paper]. Ages two-group only.
- **Processed matrix:** not described [paper].
- **Fit:** extreme-group design; small sorted sets. Mainly valuable as prior work (below).

### 8. Chinese Immune Multi-Omics Atlas (CIMA) PBMC scATAC

- **Publication:** CIMA consortium, *Science* 2026 391(6781) [secondary]; DOI 10.1126/science.adt3130
  [unverified].
- **Accession:** NGDC BioProject **PRJCA030049**; sample-level data controlled via GSA-Human **HRA009968**
  [secondary]. Portal db.cngb.org/trueblood/cima (open processed resources; contents not checked).
- **Individuals:** 428 adults, 20–77 years; scATAC on ~401 [secondary]. DNBelab C4 scATAC.
- **Fit:** large, wide age range, cell-type resolved; access is controlled for per-donor data. A 2026 paper used it
  for a single-cell accessibility clock (Wei et al., *npj Aging*, see prior work).

### 9. TenK10K phase 1 PBMC scATAC caQTL atlas

- medRxiv 2025 (10.1101/2025.08.27.25334533): ~4.5 M nuclei scATAC from **952 donors**, multiome for 90, with WGS
  [secondary]. Age range, sex and data access **not verified** (preprint page returned 403). Probably controlled.
  The largest between-individual scATAC resource found; worth a direct enquiry.

### 10. Human brain FANS-sorted neuron/glia bulk ATAC caQTL (Zeng, Bendl et al.)

- bioRxiv 2023 (10.1101/2023.03.02.530826): **1,932 samples, 616 donors**, NeuN+ and NeuN− nuclei, four regions
  (CMC: DLPFC, ACC; AMP-AD: STG, PHG) [paper/preprint]. Data: PsychENCODE Synapse **syn25955362** [paper];
  access tier, age range and diagnosis mix **not verified** (cohort includes SCZ/BD/AD cases). Peer-reviewed
  version not checked.
- **Fit:** sorted cells, very large n, adult–elderly range; disease and controlled access are the obstacles.

### Lower suitability (checked, kept for completeness)

| Dataset | Accession | n / ages | Why lower |
|---|---|---|---|
| Liver caQTL bulk ATAC (138 donors) | GSE264684 | 477 libraries; ages not in GEO | only consensus peak BED (3.6 MB), narrowPeaks, bigWigs (151 GB); no count matrix; ages not given [GEO] |
| SEA-AD MTG snMultiome | AD Knowledge Portal / brain-map.org | 84 donors, 65–102 | no young donors; AD pathology confounds; processed snATAC open on AWS [secondary] |
| Human skeletal muscle atlas (Lai et al. 2024 Nature) | CNGB CNP0004394, CNP0004395, CNP0004494, CNP0004495 | single-nucleus part 22 individuals total, 15–99 | snATAC donor count below 22 (not verified); processed data in a portal [paper][secondary] |
| Human HSPC sorted bulk ATAC across life | GSE243641 | 31 libraries, 7 ages 21–79 (6 cell types) | ~7 donors; QuantileNormalizedSignalOverBackground.tsv 68 MB [GEO] |
| PsychENCODE single-cell 388 brains (Emani 2024) | GSE261983 | GEO has 20 ATAC, ages 17–90+ | rest on PsychENCODE, controlled [GEO] |
| Mouse brain snATAC, 8 regions (Amaral et al. 2025/2026) | GSE294772 | 2, 9, 18 mo; replicates are pools of 3–14 mice | pooled nuclei: no between-individual variance; between-cell only; h5ad per sample 0.4–0.9 GB [GEO] |
| Mouse hippocampus bulk ATAC, M/F young vs aged | GSE244506 (sub GSE244504/5) | 29 ATAC libraries, 2 ages | small, two ages; processed matrix not seen [GEO] |
| Mouse HSC scATAC/bulk ATAC aging (Itokawa 2022) | GSE190424 (SuperSeries GSE162662) | young vs aged, 4 bulk replicates | two ages, few mice [paper] |
| Calderon 2019 immune stimulation ATAC | GSE118189 | 175 libraries, 25 cell types, few donors | not an aging design [GEO][unverified donor count] |
| Benayoun 2019 Genome Res mouse aging | — | 3/12/29 mo | **histone ChIP-seq (H3K4me3, H3K27ac) + RNA, not ATAC** [secondary]; drop from the ATAC list |
| ENCODE / ENTEx / Roadmap | ENCODE portal | ENTEx = 4 adult donors | too few donors per tissue for variance-with-age [unverified in detail] |
| Allen Human Immune Health Atlas / Sound Life | apps.allenimmunology.org | 108 / 234 donors | released modality is scRNA; scATAC/TEA-seq release not found [secondary] |
| Gate 2018 stimulated CD4 T cells | GEO GSE115980 holds Hi-C only | up to 105 donors | ATAC accession not found (likely dbGaP) [GEO][unverified] |
| iPSC-derived neurons/microglia ATAC | GSE335939 | 137 lines | iPSC-derived: donor age erased by reprogramming; unsuitable |

---

## Prior work on accessibility variability ("chromatin noise") with age

- **Ucar et al. 2017 (JEM), human PBMC bulk ATAC:** aging combines consistent closing at promoters/enhancers
  (T-cell signalling, IL7R) with "a potentially stochastic chromatin opening mostly found at quiescent and repressed
  sites"; opening peaks are "low-frequency and often subject-specific" [paper]. The closest direct precedent:
  it implies between-individual divergence at repressed chromatin but did not test variance per region vs age.
- **Cheung et al. 2018 (Cell), EpiTOF single-cell histone-mark profiling of human immune cells:** increased
  cell-to-cell variability of chromatin marks with age, partly environmental (twin comparisons) [secondary].
  Histone marks, not accessibility.
- **Wei et al. 2026 (npj Aging), sc-ChromAging clock on CIMA scATAC (401 donors):** immune aging spans "highly
  predictable patterns to largely stochastic changes"; lineages with weak age prediction (naive B, myeloid) are
  attributed to extrinsic heterogeneity; no variance test [paper].
- **Amaral et al. 2025/2026 (Cell Reports), mouse brain snATAC:** heterochromatin instability with age — gains of
  accessibility at heterochromatin domains in specific neuron types [GEO][secondary]; earlier Cell Research 2022
  (s41422-022-00719-6) reported age-associated decay of heterochromatin domains in excitatory neurons
  [secondary; first author not verified].
- **Bioarxiv 2026 "Cell type-resolved chromatin accessibility clocks for brain aging":** cites the information
  theory of aging (noise accumulation) as motivation but builds clocks on donor pseudobulk means; no variance
  result [paper/preprint].
- **Mouse liver (Nature Aging 2023, s43587-023-00513-y):** increased cell-to-cell transcriptional variability with
  age in diploid hepatocytes; accessibility changes reported locus by locus [secondary].
- **Itokawa et al. 2022 (Nat Commun), mouse HSC scATAC:** aged-associated open/closed regions shift in the same
  direction across all HSC clusters; heterogeneity not the question [paper].

**Gap:** no study found that measures, per region, how between-individual or between-cell variance of ATAC
accessibility changes with age, nor one comparing housekeeping/essential-gene promoters against matched regions.
The test would be new as far as this search reached.

---

## What a first analysis would need (top three)

### Catching 2026 PFC (Zenodo 18394349)
- **Download:** `final_atac_data.h5ad` (14.5 GB) + `TableS1_CohortDemographics.xlsx`; optionally the RNA h5ad
  (18.4 GB) for expression-matched controls. Store under /data, not $HOME. Open backed (`anndata.read_h5ad(...,
  backed="r")`) first to check `.X` dtype (raw counts?), `.obs` (donor ID, cell type, QC) and `.var` (peak
  coordinates, genome build).
- **Design:** per donor × cell type pseudobulk (sum counts, CPM/TMM); per-region variance across donors in
  sliding age windows or a model of the residual variance on age (e.g. a GAMLSS/Levene-type test or a
  mean-variance-corrected dispersion vs age); separately within-donor cell-to-cell dispersion per region
  (accounting for per-cell depth and sparsity).
- **Covariates:** sex, PMI, brain bank (batch/ancestry), race, nuclei count per donor and cell type, TSS
  enrichment/FRiP, sequencing depth; regress the mean-age trend before comparing variance.
- **Shield test:** promoters (±1 kb TSS) of essential genes (DepMap common-essential, Hart/Blomen sets) or
  housekeeping genes (Eisenberg & Levanon / HRT Atlas) vs promoters of non-essential genes matched on mean
  accessibility, CpG density/GC, peak width and expression level (from the matched RNA).

### Morandini 2023 PBMC (GSE193140)
- **Download:** `GSE193140_atac_counts_qval_0.05.txt.gz` (9.8 MB) and the series matrix (age, sex); check
  GSE193141/GSE193142 for the cell-composition table.
- **Design:** one library per donor (collapse or pick `_1` replicates; check the replicate structure); variance
  across donors in age bins (20s–70s; 70s has 7 only, so bin 60–74) or modelled continuously.
- **Covariates:** sex (unbalanced 115/42), cell composition (flow or a deconvolution of ATAC using sorted-cell
  references), library size/FRiP, batch if recorded (not in the GEO characteristics).
- **Limitation:** 42.5k peaks only (q < 0.05 set); bulk tissue.

### Zemke 2024 hippocampus (GSE278576)
- **Download:** 40 `*_raw_feature_bc_matrix.h5` (49–83 MB each, ~2.5–3 GB) + cell metadata (12 MB; Age, Gender,
  subclass per barcode). Build a common peak set (merge per-donor peaks or use cell-type bigWig peak calls), or
  requantify from fragment files if they are in the 52 GB RAW tar (not verified).
- **Design:** per donor × subclass pseudobulk, variance vs age across 40 donors (~10 per 20-year bin); cell-level
  dispersion per donor.
- **Covariates:** sex, nuclei per subclass, TSS enrichment/depth; PMI and cause of death if obtainable from the
  paper's supplementary table (not in GEO).

---

## Sources consulted
- GEO FTP and series matrices: GSE193140, GSE278576, GSE288730, GSE331488, GSE328206/GSE341533, GSE264684,
  GSE294772, GSE244506, GSE243641, GSE261983, GSE338484, GSE335939, GSE247455, GSE118189, GSE115980,
  GSE98983/4 (checked via E-utilities and https://ftp.ncbi.nlm.nih.gov/geo/series/).
- Zenodo API record 18394349 and its TableS1/TableS2.
- PMC full texts: PMC7005316 (Márquez 2020), PMC5626401 (Ucar 2017), PMC10120699 (Zeng/Bendl preprint),
  PMC11062927 (Lai 2024), PMC9110722 (Itokawa 2022), PMC13376647 (Wei 2026).
- bioRxiv 10.64898/2026.08.22.746439 (brain accessibility clocks).
- PubMed: PMID 37924441 (ATAC-clock), 41832957 (Catching 2026).
- Search summaries (secondary): CIMA (NGDC PRJCA030049), TenK10K preprint, SEA-AD portal, Allen Immune Health
  Atlas, Benayoun 2019, Cheung 2018.

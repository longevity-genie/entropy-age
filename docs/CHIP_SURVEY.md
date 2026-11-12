# Public data for selection in vivo: clonal haematopoiesis and methylation (roadmap RM16)

Survey of public data, 2026-11-12, by a delegated agent for the lead; unverified items marked. Not yet a result.

# Public methylation data with clonal-haematopoiesis (CHIP) status: survey

Date: 2026-11-12. Question: do CHIP clones change how methylation at essential/housekeeping-gene promoters drifts with
age (the "shield")? Carriers vs non-carriers, clone size (VAF), and mutant vs wild-type cells within a person.

How this was checked: GEO E-utilities sweeps (esearch/esummary on `gds`, ~1,100 series screened), GEO SOFT
sample annotations, the GEO FTP listings for file sizes, PubMed/Europe PMC full texts for data-availability
statements, and 60 KB range-reads of two processed files to see their layout. Nothing over 50 MB was
downloaded. "Verified" means seen in GEO metadata, the FTP listing or the paper's own text. "[unverified]"
means I did not see it.

Corrections to the brief's candidate list (verified):
- Nachun et al. 2021 is in **Aging Cell** (20:e13366), not Genome Biology. Its cohorts are FHS, JHS, WHI and MESA (TOPMed), n=5,522. ARIC is not among them.
- Robertson et al. 2019 is in **Current Biology** (29:R786), on the Lothian Birth Cohorts, not Genome Biology.
- Nam et al. 2022 is in **Nature Genetics** (54:1514): single-cell RRBS plus genotyping of DNMT3A R882 CH. Nam 2019 in Nature is GoT (RNA only), and GoT-ChA measures chromatin, not methylation.
- Tovy 2022 (Haematologica) covers germline DNMT3A (Tatton-Brown-Rahman syndrome, TBRS). The DNMT3A-mutant clone study is Tovy et al. 2020, Cell Stem Cell: a mosaic individual with "epigenetic erosion".
- Izzo et al. 2020 (Nat Genet, GSE124822): all methylation in it is **mouse**. The only human sample is one RNA-seq/GoT sample.
- Jakobsen et al. 2024 (Cell Stem Cell) used TARGET-seq (RNA and genotype). It has no methylation.
- Mitchell 2022 (Nature), Fabre 2022 (Nature) and Spencer Chapman 2024 (Nature) are phylogeny or longitudinal-VAF studies with no methylation (checked from titles and abstracts).
- "Uddin/Liu": Uddin et al. 2022 (Nat Commun) is verified. I found no CHIP-methylation study first-authored by Liu. [unverified]

---

## Ranked candidates

### 1. GSE315366: Basrai et al. 2026, Nat Commun 17:2754 (PMID 41690920), "Blood-based epigenetic instability linked to human aging and disease" (cardiogenic-shock cohort)
- **Design:** 64 cardiogenic-shock patients. 29 have CH (DNMT3A 21, TET2 12, both 4) and 35 have no CH. Verified by parsing the GEO sample annotation.
- **Per person:** age, sex, `dnmt3a mutation` TRUE/FALSE, `tet2 mutation` TRUE/FALSE, `clonal hematopoiesis` TRUE/FALSE.
- **VAF:** not in GEO. CH was called in Scolari et al. 2022, Eur J Heart Fail (PMID 35729851). VAFs might be in that paper's supplement, but its patient IDs are not mapped to GEO's `Patient_n`. [unverified]
- **Age:** carriers 33-89 (median 61), non-carriers 18-93 (median 61).
- **Platform and cells:** EPIC (GPL21145), peripheral blood (buffy coat according to the paper).
- **Access:** open.
- **Files:**
  - `GSE315366_cardiogenic_shock_beta_values.txt.gz`, **462 MB**: tab-separated, `ID_REF` plus 64 `Patient_n` columns, one beta per column. Checked by range-read.
  - `GSE315366_RAW.tar`, 1.0 GB (IDATs).
  - Series matrix, 4.5 KB (annotations only).
- **Why first:** it is the only open array set I found with per-person CHIP status and betas. The paper also measures nearly the quantity in question. Its "DMI" is the per-sample SD of beta across 31,744 CpGs that are stably unmethylated in young blood: 92.6% lie within 1.5 kb of a TSS, CGI OR=15. DMI rises with age, tracks clonal burden in AML, and predicts mortality in this cohort independently of CH status.
- **Caveats:**
  - Patients are acutely ill, which raises inflammation and shifts cell composition.
  - Only DNMT3A and TET2 were reported.
  - n is small, and there is no VAF.
- **Sister series** (myeloid neoplasms, the fallback set):
  - GSE315451: 57 MDS + 69 CML at diagnosis, EPIC, mostly bone marrow. Beta file 908 MB, RAW 1.9 GB.
  - GSE315367: 14 AML patients, longitudinal, 46 samples, EPIC. Beta file 334 MB. The paper reports VAFs for 4 patients.
  - The other 29 cohorts are listed in the paper's Supplementary Data 1, which I did not open. [unverified] It may name more public sets.

### 2. GSE158067: Nam et al. 2022, Nat Genet 54:1514 (PMID 36138229), single-cell multi-omics of DNMT3A R882 CH
- **Design:** CD34+ HSPCs from G-CSF-mobilised blood or bone marrow of people with DNMT3A R882 CH (5 individuals in the study).
- **Single-cell methylation (MspI scRRBS) with per-cell genotype**, from the GEO metadata xlsx (cells passing QC, >25k CpGs):
  - CH02: 155 mutant, 65 WT, 64 not genotyped.
  - CH04: 102 mutant, 137 WT, 40 not genotyped.
  So there are about 257 mutant and 202 WT cells from 2 people.
- **Other assays:** scRNA-seq with genotyping (10x GoT and Smart-seq2).
- **VAF:** per cell rather than per person. Bulk VAFs are in the paper. [unverified detail]
- **Age:** older adults. Exact ages are not in the GEO metadata. [unverified]
- **Platform:** scRRBS (GPL24676/GPL11154). RRBS is enriched for CpG-island promoters, which suits the housekeeping-promoter question.
- **Access:** processed data open; raw reads in EGA (controlled).
- **Files:**
  - `GSE158067_all_CpG_meth.csv.gz`, **282 MB**: all CpGs by cell (not downloaded).
  - `GSE158067_scDNAme_cell_metadata.xlsx`, 41 KB: downloaded; it holds the genotype per cell.
  - `GSE158067_RAW.tar`, 10 GB.
- **Why second:** this is the cleanest natural experiment on selection. Mutant and WT cells come from the same person, the same age and the same environment. Published finding: R882 causes *selective* hypomethylation that is focal, at PRC2/early-progenitor regulatory elements, not genome-wide.
- **Caveats:**
  - Only 2 donors have methylation.
  - Coverage per cell is sparse, so promoter-level pooling is needed.
  - The cells are CD34+, not mature blood.

### 3. GSE210430 (in SuperSeries GSE210435): Mohammed Ismail et al. 2025, Blood Adv (PMID 39631069), and Ben-Crentsil et al. 2024, Cancer Discov (PMID 39189614); TET2/DNMT3A CH and COVID-19
- **Design:** 18 PBMC samples: 5 CH only, 9 CH + COVID-19, 4 COVID-19 only. EPIC after TrueMethyl oxBS, so the readout is 5mC proper.
- **Mutations and VAFs per patient** are drawn in Fig. 2A of the Blood Adv paper, not given in GEO. They would need transcribing.
- **Age:** the cohort ranged 19-99 (median 60). Per-sample ages are not in GEO. [unverified]
- **Access:** open.
- **Files:**
  - `GSE210430_GEO_metadata_Methylation_Matrix_signal_intensities.csv.gz`, **85 MB**: methylated and unmethylated signal per CpG and sample (checked by range-read).
  - RAW tar, 412 MB.
  - Series matrix, 80 MB.
- **Caveats:**
  - Tiny.
  - Only 4 non-CH samples, all COVID-positive. A clean non-CH control is missing.
  - PBMC excludes granulocytes, where myeloid CH clones are largest.
  - The SuperSeries also holds scRNA, multiome and GoT-ChA from the same patients.

### 4. EPI-Clone: Scherer et al. 2025, Nature 643:478 (PMID 40399669)
- **Design:** targeted single-cell methylation (scTAM-seq, HhaI-based) plus CITE-seq plus myeloid-panel genotyping (DNMT3A, TET2, ASXL1, JAK2, PPM1D and others). Human cohorts:
  - A.1-A.7: total bone marrow, ageing donors.
  - B.1-B.5: CD34+, donors mostly aged 50-60.
  - X.1-X.2.
  CH clones are identified as epimutation clones; for example, DNMT3A C666Y and R659H in donor A.4.
- **Access:** processed Seurat objects open on Figshare (A.1-A.7: 10.6084/m9.figshare.25526899.v2; B.1-B.5: 10.6084/m9.figshare.28082048.v1). Count matrices are in GSE282971. Human raw reads are in EGAS00001008056 (controlled; use restricted to haematology and methods research).
- **File sizes:** [unverified]
- **Limitation:** the panel is a few hundred CpGs chosen *because* they vary (clonal-epimutation and differentiation CpGs). Housekeeping promoters that stay unmethylated are mostly absent by design, so this set does not test the shield directly. It is useful to show clonal stability of epimutations: the paper says they arise randomly, stay stable, and do not increase in number across divisions.

### 5. Controlled cohorts with CHIP calls from WGS/WES and 450k/EPIC: the best-powered, but dbGaP
| Study | Cohorts | n | Genes / VAF | Access |
|---|---|---|---|---|
| Nachun 2021, Aging Cell | FHS, JHS, WHI, MESA (TOPMed) | 5,522 | all CHIP genes; VAF from WGS | dbGaP phs000974/phs000007, phs000964, phs001237/phs000200, phs001416 |
| Uddin 2022, Nat Commun | CHS (discovery), ARIC (replication) | 582 + 2,655 | DNMT3A, TET2 and others; VAF | CHS phs001368 (CHIP calls), ARIC phs000668; methylation via ancillary-study proposals (4-9 weeks) |
| Kirmani 2025, Nat Commun | FHS, JHS, CHS, ARIC | 8,196 | EWAS of CHIP | phs000007, phs000964, phs001368, phs000668 |
| Mack 2024, Nat Aging | TOPMed | 4,370 | clonal expansion rate (PACER) | dbGaP TOPMed |

- All whole blood (FHS and others), 450k/EPIC, ages about 40-90.
- These are the right sets for the full test (thousands of people, VAF, many genes), but they need a dbGaP application.

### 6. Lothian Birth Cohorts 1921/1936: Robertson et al. 2019, Curr Biol (PMID 31430471); Robertson et al. 2022, Nat Med (PMID 35788175)
- **Design:** about 6% of participants carry CHIP. The cohorts have repeated 450k/EPIC whole-blood methylation and, in the 2022 paper, longitudinal VAF trajectories over about 12 years.
- **Value:** within-person drift could be compared against within-person clone growth, which no other set allows.
- **Access:** controlled, through LBC data access. Methylation is deposited in EGA. Accession numbers and CHIP-call release were not checked. [unverified]

### 7. Danish twins (LSADT): Tulstrup et al. 2021, Nat Commun (PMID 34663818)
- **Twins:** 305 aged ≥70; 116 carry CHIP (55 DNMT3A, 44 TET2), VAF ≥2%. Whole-blood 450k.
- **CCUS:** granulocytes (5 TET2 with VAF 33-50% vs 8 controls) and BM-MNC (20; 10 TET2), EPIC.
- **Access:** restricted, via the Danish Twin Registry (about 6 weeks plus legal review) and the authors (about 6 months). There is no GEO or EGA accession.

### 8. Other CHIP-plus-methylation studies, less suitable
- **GSE253126:** Singapore Chinese Health Study, MeDIP-seq on blood about 9.9 years before diagnosis. 66 people later developed AML/MDS (38 pre-MDS, 28 pre-AML) and 167 controls. Region-level only (MeDIP), no CpG betas. Processed xlsx 14 MB. Mutations from targeted sequencing are in the paper (PMID 38437679). Open.
- **Feldkamp 2022, Haematologica (BASE-II):** CHIP-related age acceleration rises with clone size. Data "on reasonable request".
- **Parker 2025, Blood (PMID 39565983):** a methylation-sequencing assay that infers cell-type proportions to read CH VAF trajectories. Data availability not found. [unverified]
- **Hasan 2026, Geroscience:** 56 CH and 12 controls with VAF, but PBMC RNA only. No methylation.

### 9. Sorted HSPC vs mature cells, young vs old (no CHIP status)
- **GSE104405:** Adelman et al. 2019, Cancer Discov. ERRBS of bone-marrow HSCs (Lin-CD34+CD38-), young (19-27, n≈5-6) vs aged (63-71, n≈5-6), no known haematological malignancy. CH was not assessed. RAW 356 MB. Open. This is the HSC-level baseline: does the shield exist before any selection in the blood?
- **GSE87196 (BLUEPRINT):** WGBS of HSC, MPP, CMP, GMP, MEP, CLP, MLPs and mature T, B, NK, neutrophils and monocytes. Processed open (Uddin 2022 used it); raw in EGA. Few donors.

### 10. Germline or engineered DNMT3A/TET2 models (no selection, gene dose only)
- **GSE128801:** Jeffries 2019, Genome Res. 450k blood of an Amish TBRS sibship (DNMT3A R771Q) and their *mosaic* father. 12 samples. Beta matrix open. Reports accelerated epigenetic ageing.
- **Smith 2021, Nat Commun:** TBRS WGBS, dbGaP phs000159 (controlled).
- **Kaasinen 2019, Nat Commun:** germline TET2 haploinsufficiency family. EGA EGAS00001003454 (controlled).
- **Tovy 2020, Cell Stem Cell:** mosaic DNMT3A, mutant vs WT LCL clones. Accession not found. [unverified]
- **GSE168581:** Lim 2021, JEM. EPIC of hESC-derived DNMT3A+/- and -/- myeloid cells, 40 samples, RAW 727 MB. Open.

### 11. Myeloid-neoplasm fallback (open)
- GSE315451 and GSE315367 (above).
- GSE31600: CMML, TET2.
- GSE105420: CMML.
- GSE129828: MDS.
- GSE319461 / GSE319444: AML WGBS methylation stochasticity, with normal CD34+, CMP, GMP and MEP. RAW 14-32 GB; no processed matrix listed.

---

## What a first analysis needs (top three)

**GSE315366**
1. Download the beta file (462 MB, too large to pull here) and parse the 64-sample annotation (done once already).
2. Map the shield CpG set and its matched controls onto EPIC probes.
3. Per sample, compute the dispersion of shield CpGs vs matched CpGs: SD, or |beta − young-reference mean| using the young-blood reference medians.
4. Model it as dispersion ~ age × CH + sex + estimated cell fractions (EpiDISH/FlowSorted.Blood.EPIC). Then split by DNMT3A and TET2.
5. Optionally, request VAF (Scolari 2022 / Abelson lab) to fit dispersion against VAF.
6. Watch for acute-illness confounding: adjust for neutrophil fraction.

**GSE158067**
1. Download `all_CpG_meth.csv.gz` (282 MB) and join it to the per-cell genotype in the metadata xlsx.
2. Pool the sparse per-cell calls into promoter bins of shield genes vs matched genes.
3. Compare mutant vs WT cells within each donor, for mean and for the cell-to-cell variance or epipolymorphism.
4. Control for the progenitor-state annotation (from Smart-seq2 RNA), because R882 shifts states.
5. Read the result as a direct test: if mutant clones lose the shield, their promoters should show more variance or loss than WT cells from the same person.

**GSE210430**
1. Compute betas from the signal-intensity matrix (85 MB) or the IDATs. Note that the input is oxBS, so it is 5mC only.
2. Transcribe the per-patient genes and VAFs from Fig. 2A of Mohammed Ismail 2025.
3. Use it as a VAF-gradient replication of the GSE315366 test. With n=18 this is descriptive only. Separate the COVID effect (4 COVID-only vs 9 CH+COVID).

---

## Published findings on CHIP and epigenetic drift, entropy or age (one line each)
- **Basrai 2026, Nat Commun:** variability at CpGs stably unmethylated in young blood (promoter/CGI) rises with age, tracks clonal burden in AML, and predicts mortality independently of CH.
- **Scherer 2025, Nature:** stochastic epimutations at specific CpGs act as stable clonal barcodes. Expanded clones with and without CH drivers form one age-related spectrum, becoming detectable at ~50 and inevitable after ~60.
- **Nam 2022, Nat Genet:** DNMT3A R882 causes selective, not global, hypomethylation in early progenitors.
- **Tovy 2020, Cell Stem Cell:** DNMT3A-mutant clones show conserved focal hypomethylation ("epigenetic erosion") compared with WT clones from the same person.
- **Uddin 2022, Nat Commun:** DNMT3A-CHIP and TET2-CHIP have opposite genome-wide associations (hypo vs hyper).
- **Tulstrup 2021, Nat Commun:** TET2-CHIP brings enhancer hypermethylation in whole blood, with no change in global mean methylation.
- **Kirmani 2025, Nat Commun:** an EWAS of CHIP in 8,196 people finds CHIP-associated CpGs, 261 of them linked to cardiovascular traits and mortality.
- **Nachun 2021 and Robertson 2019:** CHIP carriers are epigenetically older by about 1.3-4.5 years, depending on the clock.
- **Simonson 2026, Ageing Res Rev (meta-analysis):** larger clones go with more age acceleration, and TET2 has a larger effect than DNMT3A.
- **Mack 2024, Nat Aging:** faster clonal expansion goes with measured and genetically predicted epigenetic age acceleration.
- **Du/Teschendorff 2026, Nat Commun (OmniAge):** mitotic clocks correlate with DNAm-based CH proxies. Note that CH was not sequenced.
- **Feldkamp 2022, Haematologica:** CHIP-related age acceleration scales with clone size.
- **Jeffries 2019, Genome Res:** germline DNMT3A (TBRS) gives strongly accelerated epigenetic ageing.
- **Adelman 2019, Cancer Discov:** aged human HSCs lose methylation and H3K27ac at enhancers. CH was not assessed.

No paper I found tests drift at essential or housekeeping promoters by CHIP status. Basrai 2026 is the closest, and their stable-unmethylated CpG set probably overlaps the shield CpGs heavily (worth intersecting).

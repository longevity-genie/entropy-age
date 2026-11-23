<!-- Search by a delegated agent, 2026-11-23. Verified from GEO sample records (brief text view), GEO download pages (archive sizes) and the GDC files API only; papers not read; nothing downloaded. Not committed. -->

# Raw-IDAT methylation series for non-tumour solid tissues (GTEx validation)

**Goal.** Find public 450k, EPIC or EPIC v2 series with raw IDATs for the solid tissues where our age results come only
from GTEx deposited betas (GSE213478; GTEx IDATs are dbGaP-only). Priority tissues: lung, kidney, testis, breast, ovary,
prostate and skeletal muscle, then any other GTEx tissue.

**Criteria.** A series meets when it has all of the following:

- raw IDATs on GEO or ArrayExpress;
- a numeric age per sample;
- adult donors spanning 30 or more years;
- 30 or more donors of one tissue;
- chip (Sentrix) IDs recoverable;
- one sample per donor, or donor IDs given.

Healthy tissue is preferred to tumour-adjacent tissue.

**How it was checked.** I ran GEO ESearch/ESummary (`db=gds`) on the 450k, EPIC and EPIC v2 platforms (GPL13534,
GPL16304, GPL21145, GPL23976, GPL33022, GPL34394) for 19 tissue term sets, then ran free-text searches. 1,662 series
were found and screened by title, IDAT flag and size. Every candidate named below was then read sample by sample from
the GEO brief text view. That reading gave the n, the numeric age range, the subset labels and the chip IDs from IDAT
file names. Each archive size comes from the GEO download page. ArrayExpress was searched through the BioStudies API.

**Age-by-chip column.** "age-by-chip η²" is the share of age variance explained by chip within the subset. The chance
level is roughly k/n, where k is the number of chips and n the number of samples; both are given. "Chips <45 & ≥65"
counts the chips that hold both a donor under 45 and a donor aged 65 or over.

**Excluded as already decided.** GSE314841, GSE101764, GSE296589, GSE184211, GSE50498, GSE213478 and GSE61454 are not
re-reported. Neither is any series with a verdict in `DATASET_SEARCH_2026-11-18.md`, `DATASET_SEARCH_2026-11-24_more_samples.md` or
`entropy-age-audits/2026-11-23-muscle-search/CANDIDATES.md`. That covers GSE244996 (20 donors × 9 tissues), GSE235414
(lung, adjacent), GSE131013, E-MTAB-8993/8992, GSE180474, GSE61258, GSE90124, GSE51954, GSE197670 and GSE164197.

## Lung

| Accession | PMID | Platform | n donors (tissue) | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE175458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175458) | [35816432](https://pubmed.ncbi.nlm.nih.gov/35816432/) | EPIC | 145 `dx: normal` (112 Pittsburgh Lung Donor Core, 25 UCSF, 7 National Jewish, 1 LTRC) + 6 LTRC `control`; one array per title | 14–88 (normals; median 53; 48 <45, 35 ≥65) | Non-diseased control lung; the series also has 345 IPF and 51 "other" | Yes; RAW.tar 7.7 GB for 547 arrays, normals ≈2.0 GB | Yes (`<chip>_R01C01_noid_Grn.idat`, also `methylationid` field); 60 chips for normals, age-by-chip η² 0.37 (k/n 0.41), 16 chips <45 & ≥65 | **Meets.** Largest open non-tumour lung series found. Donor-lung origin of the Pittsburgh controls is inferred from the centre name; read the paper's methods before relying on "not resection". RNA-seq of the same lungs is in super-series GSE175459. Smoking ever/never is recorded. |
| [GSE173356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173356) | [34038697](https://pubmed.ncbi.nlm.nih.gov/34038697/) | EPIC | 13 normal (24 IPF) | 52–75 | Normal control lung | Yes, 0.52 GB | Yes; 27 chips for 37 arrays, η² 0.82 (k/n 0.73) | **Fails**: n 13, span 23 y. |
| [GSE268362](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268362) | [39970931](https://pubmed.ncbi.nlm.nih.gov/39970931/) | EPIC | 6 control donors × 6 regions (35 arrays); 10 IPF | 43–72 | Control lung | Yes, 1.5 GB | Yes, 14 chips | **Fails**: 6 control donors. |
| [GSE52401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52401) | [24572595](https://pubmed.ncbi.nlm.nih.gov/24572595/) | 450k | 244 "normal lung" | No age field | Normal lung (resection, not stated) | **No IDATs** (processed only) | No | **Fails**: no IDATs, no ages. |
| [GSE63704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63704) | [26425700](https://pubmed.ncbi.nlm.nih.gov/26425700/) | 450k | 43 healthy | No age field | Healthy, fibrosis, COPD, cancer | **No IDATs** | No | **Fails**. |
| [GSE66836](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66836) | [26601720](https://pubmed.ncbi.nlm.nih.gov/26601720/) | 450k | 19 normal | No age field | Tumour-adjacent | Yes, 1.6 GB (whole series) | Yes, 16 chips | **Fails**: n 19, no ages. |
| [GSE138348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138348) | [31619441](https://pubmed.ncbi.nlm.nih.gov/31619441/) | EPIC | 29 healthy volunteers × 2 weeks (airway, e-cig trial) | 21–30 | Healthy | Yes, 0.97 GB | Yes, 10 chips | **Fails**: age span 9 y. |

TCGA fallback, outside GEO and ArrayExpress (see the end note): LUSC has 42 normal-adjacent arrays with open IDATs and
LUAD has 32.

## Kidney

| Accession | PMID | Platform | n donors (tissue) | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE79100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79100) | none in GEO | 450k | 31 "normal kidney" (16 F, 15 M; 16 EA, 15 AA) | 18–78 (median 56; 8 <45, 8 ≥65) | "Normal human kidney"; source (nephrectomy margin or donor) **not stated in GEO**; likely tumour-nephrectomy normal | Yes, 0.95 GB | Yes; 5 chips, η² 0.04 (k/n 0.16), 4 of 5 chips <45 & ≥65 | **Meets with a caveat**: n just above 30; check tissue origin in the paper (title "Sex, Age, and Ethnic Differences"). |
| [GSE61441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61441) | [26515236](https://pubmed.ncbi.nlm.nih.gov/26515236/) | 450k | 46 matched normal | No age field | Tumour-adjacent (ccRCC) | Yes, 0.94 GB | **No** (files renamed `2014ccRCC01N_Grn.idat`) | **Fails**: no ages, no chip IDs. |
| [GSE253101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253101) | [40939817](https://pubmed.ncbi.nlm.nih.gov/40939817/) | EPIC | 43 normal (25 cortex, 18 medulla) | No age field | Tumour-adjacent | Yes, 5.6 GB (411 arrays) | Yes | **Fails**: no ages. |
| [GSE105288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105288) | [30131446](https://pubmed.ncbi.nlm.nih.gov/30131446/) | 450k | 17 normal renal | No age field | Tumour-adjacent | Yes, 0.66 GB | Yes, 4 chips | **Fails**. |
| [GSE233758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233758) | [38961842](https://pubmed.ncbi.nlm.nih.gov/38961842/) | EPIC | 33 microdissected proximal tubules (8 normal, 4 DN, 21 diabetic) | No age field | Diabetic / normal | Yes, 0.61 GB | Yes, 8 chips | **Fails**: no ages. |
| [GSE50874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50874) | [24098934](https://pubmed.ncbi.nlm.nih.gov/24098934/) | 450k | 85 tubule samples (CKD/control) | No age field | CKD and control | **No IDATs** | No | **Fails**. |
| [GSE292896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292896) | [39665570](https://pubmed.ncbi.nlm.nih.gov/39665570/) | EPIC | 330 kidney arrays | Children | Predisposed children, Wilms | Yes, 6.3 GB | Yes | **Fails**: paediatric. |
| [GSE269180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269180) | none | EPIC | 121 | 23–84 | **All ccRCC tumour** (no normals despite the title) | Yes, 1.8 GB | Yes | **Fails**: tumour only. |

TCGA fallback: KIRC has 160 normal-adjacent kidneys with open IDATs, ages 34–89, about 1.4 GB; KIRP adds 45.

## Testis

No 450k, EPIC or EPIC v2 series of adult testis tissue was found on GEO or ArrayExpress, with or without IDATs. The 30
GEO hits on testis terms are germ-cell tumours, sperm, cell lines or "cancer-testis antigen" studies. TCGA-TGCT has no
normal-adjacent arrays. **No usable series.**

## Breast

| Accession | PMID | Platform | n donors (tissue) | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE88883](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88883) | [28693600](https://pubmed.ncbi.nlm.nih.gov/28693600/) | 450k | 100 healthy breast biopsies (Komen Tissue Bank IDs K10xxxx) | 18–82 (median 37; 72 <45, 5 ≥65) | **Healthy donors**, cancer-free; family history recorded | Yes, 1.0 GB | Yes; 9 chips, η² 0.09 (k/n 0.09), 4 of 9 chips <45 & ≥65 | **Meets.** Age skews young, few over 65. |
| [GSE101961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101961) | [29383109](https://pubmed.ncbi.nlm.nih.gov/29383109/) | 450k | 121 cancer-free women (81 white, 40 black) | 17–76 (median 38; 80 <45, **1 ≥65**) | Healthy, cancer-free | Yes, 1.1 GB | Yes; 28 chips, η² 0.25 (k/n 0.23), 1 chip <45 & ≥65 | **Meets with a caveat**: the oldest decade is nearly empty, so slopes rest on 17–64 y. No chip or position overlap with GSE88883 (checked). |
| [GSE287331](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287331) | [40500734](https://pubmed.ncbi.nlm.nih.gov/40500734/) | EPIC | 182 healthy donated breast (HDB, Komen IDs); also 60 AN, 67 OQ, 68 CUB and 69 tumour from 69 cases | **No age field in GEO** | Healthy (HDB); case tissues are tumour-adjacent or contralateral | Yes, 6.3 GB (446 arrays), HDB ≈2.6 GB | Yes, 63 chips | **Fails as deposited**: ages absent. Worth a request or a check of the paper's supplement, because 182 healthy EPIC breasts would be the largest set. 5 Komen IDs are also in GSE88883. |
| [GSE60185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60185) | [25146004](https://pubmed.ncbi.nlm.nih.gov/25146004/) | 450k | 46 normal (17 reduction mammoplasty, 29 needle biopsy) | Age **NA for all normals** | Healthy / reduction | Yes, 2.5 GB | Yes | **Fails**: no ages for the normals. |
| [GSE69914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69914) | [26664652](https://pubmed.ncbi.nlm.nih.gov/26664652/) | 450k | 50 normal (healthy women), 42 normal-adjacent | No age field | Healthy and adjacent | **No IDATs** | Chip in characteristics | **Fails**: no IDATs. |
| [GSE67919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67919) | [25953686](https://pubmed.ncbi.nlm.nih.gov/25953686/) | 450k | 96 (71 mastectomy-adjacent, 18 reduction, 7 prophylactic) | 19–84 | Mostly tumour-adjacent | **No IDATs** | No | **Fails**: no IDATs. |
| [GSE141338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141338) | [32019179](https://pubmed.ncbi.nlm.nih.gov/32019179/) | 450k | 6 normal | — | Non-tumoural | Yes, 0.58 GB | Yes | **Fails**: n 6. |
| [GSE74214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74214) | none | 450k | 18 | 13–80 | Healthy | Yes, 0.33 GB | Yes, 3 chips | **Fails**: n 18. |

TCGA fallback: BRCA has 97 normal-adjacent arrays with open IDATs.

## Ovary (and fallopian tube)

| Accession | PMID | Platform | n donors (tissue) | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE293125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293125) | none | GPL34510 (Infinium, new chip) | 44 granulosa cells (IVF) | 27–41 | Fertility patients | Yes, 0.65 GB | Yes, 7 chips | **Fails**: cells, not tissue; 14-y span. |
| [GSE155761](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155761) | [32817081](https://pubmed.ncbi.nlm.nih.gov/32817081/) | EPIC | 12 adjacent normal fallopian tube; 13 cancer-free endometrium, 12 cervix | 34–79 (part) | Tumour-adjacent / cancer-free | Yes, 1.9 GB | Yes | **Fails**: n ≤13 per tissue. |
| [GSE80468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80468) | [27192117](https://pubmed.ncbi.nlm.nih.gov/27192117/) | 450k | — | — | Whole blood (PCOS) | Yes | Yes | **Fails**: blood. |

No normal ovary tissue series with arrays was found. GEO and ArrayExpress hold only ovarian tumours, granulosa or
cumulus cells, and PCOS blood. TCGA-OV has no 450k normals. **No usable series.**

## Prostate

| Accession | PMID | Platform | n donors (tissue) | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE183015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183015) | [38112617](https://pubmed.ncbi.nlm.nih.gov/38112617/) | EPIC | 58 `_BENIGN` prostate, one per patient (patient IDs given; tumour and buffy coat of the same men also present) | 52–75 (median 64; 23 ≥65) | **Tumour-adjacent** benign tissue from prostatectomy | Yes, 2.5 GB (174 arrays), benign ≈0.83 GB | Yes; 16 chips, η² 0.48 (k/n 0.28), no chip spans <45 & ≥65 | **Meets with a caveat**: tumour-adjacent; 23-y span, below the 30-y aim; age is somewhat chip-bound. Best prostate option. Also in super-series GSE183040 with RNA. |
| [GSE323385](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE323385) | none | EPIC | 53 normal (African American men) | Age `#N/A` for all normals | Tumour-adjacent | Yes, 1.8 GB | Yes, 15 chips | **Fails**: normals lack ages. |
| [GSE269244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269244) | [39162297](https://pubmed.ncbi.nlm.nih.gov/39162297/) | EPIC | 118 normal | No age field | Tumour-adjacent | Yes, 3.4 GB | Yes, 32 chips | **Fails**: no ages (largest prostate normal set; ages may be in the paper). |
| [GSE304293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304293) | [41057544](https://pubmed.ncbi.nlm.nih.gov/41057544/) | EPIC + EPIC v2 | 70 non-tumour | No age field | Tumour-adjacent | Yes, 4.2 GB | Yes | **Fails**: no ages. |
| [GSE323379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE323379) | none | EPIC | 18 normal, 15 BPH | BPH 55–78; normals no age | Adjacent / BPH | Yes, 1.4 GB | Yes | **Fails**: n, ages. |
| [GSE116338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116338) | none | 450k | 10 benign epithelium | 46–81 | Adjacent | Yes | Yes | **Fails**: n 10. |
| [GSE112047](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112047), [GSE73549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73549), [GSE101908](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101908), [GSE245609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245609), [GSE262522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262522) | 29740534, 28081708, –, –, 38566104 | 450k / EPIC | 13–32 normal each | No ages | Adjacent | Yes | Mostly yes (GSE262522 renamed, no) | **Fail**: n and no ages. |
| [GSE76938](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76938) | [28412973](https://pubmed.ncbi.nlm.nih.gov/28412973/) | 450k | 63 benign | 43–73 | Adjacent | **No IDATs** | No | **Fails**: no IDATs. |

TCGA fallback: PRAD has 50 normal-adjacent arrays with open IDATs.

## Skeletal muscle

Nothing new meets. Beyond the series already judged in the muscle search:

- [GSE135066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135066) (PMID 31519890): 57 EPIC arrays,
  metabolic-surgery time course. No age field. Fails.
- [GSE343068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343068): the super-series of GSE341192, already
  judged. No ages.
- [GSE173613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173613) (PMID 34802257): 16 CABG patients
  pre- and post-surgery. No ages. Fails.
- [GSE303899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303899) (PMID 41578880): sarcopenia study, but
  the samples are peripheral blood. Fails.

**No muscle series meets the full criteria.** The muscle audit's conclusion stands. Its nearest options are GSE61454
(fails QC here) and GSE244996 (20 donors).

## Other GTEx tissues

| Tissue | Accession | PMID | Platform | n donors | Age range | Status | IDATs, size | Chip IDs | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| Colon (transverse/sigmoid analogue) | [GSE151732](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151732) | [33377907](https://pubmed.ncbi.nlm.nih.gov/33377907/) | EPIC | 128 donors × right and left colon (pair IDs given) | 34–80 (median 55; 5 <45, 25 ≥65) | **Healthy** screening colonoscopy, not tumour-adjacent; 48 donors with adenoma history; 89 AA | Yes, 3.5 GB (one side ≈1.75 GB) | Yes; 36 chips; left η² 0.24 (k/n 0.28), right 0.43 (0.28) | **Meets.** A healthy-colon replicate for the tumour-adjacent GSE101764 result. |
| Colon (rectum) | [GSE216024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216024) | [36612042](https://pubmed.ncbi.nlm.nih.gov/36612042/) | EPIC | 140 normal rectum | 34–80 | Healthy; 58 with adenoma | Yes, 2.1 GB | Yes; 19 chips, η² 0.20 (k/n 0.14) | **Meets with a caveat**: **123 of its 140 donor IDs are GSE151732 donors**, so it is the same cohort, not an independent replicate. |
| Thyroid | [GSE97466](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97466) | [28938489](https://pubmed.ncbi.nlm.nih.gov/28938489/) | 450k | 50 non-neoplastic | 20–77 (median 40; 31 <45, 5 ≥65) | **Tumour-adjacent** (PTC and other nodules) | Yes, 1.3 GB (141 arrays), normals ≈0.46 GB | Yes; 9 chips, η² 0.06 (k/n 0.18), 4 chips <45 & ≥65 | **Meets with a caveat**: tumour-adjacent; few over 65. [GSE86961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86961) is a subset (all 82 arrays identical by chip and position), so do not pool them. TCGA-THCA adds 56. |
| Liver | [GSE107038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107038) | [29554203](https://pubmed.ncbi.nlm.nih.gov/29554203/) | 450k | 40 healthy liver donors | 13–90 (11 <45, 17 ≥65; 1 aged 13) | **Healthy** (transplant donors) | Yes, 0.51 GB | Yes; 5 chips, η² 0.16 (k/n 0.12), all 5 chips <45 & ≥65 | **Meets.** Not in the earlier liver lists. Small but clean, with matched expression in super-series GSE107039. |
| Liver | [GSE61278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61278) | [25282492](https://pubmed.ncbi.nlm.nih.gov/25282492/) | 450k | 52 accident victims, 44 cancer | Only "Adult/Fetal" | Mixed | Yes, 1.0 GB | Yes | **Fails**: no numeric age. |
| Liver | [GSE123995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123995) | [31798965](https://pubmed.ncbi.nlm.nih.gov/31798965/) | EPIC | 56 primary hepatocyte donors (AA) | 0.25–67 | Cultured hepatocytes | Yes, 0.79 GB | Yes | **Fails**: cells; children included. |
| Pancreas | [GSE49149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49149) | [24500968](https://pubmed.ncbi.nlm.nih.gov/24500968/) | 450k | 29 adjacent | 33–87 (18 ≥65) | Tumour-adjacent (PDAC) | Yes, 1.7 GB, adjacent ≈0.25 GB | Yes, 13 chips | **Meets with a caveat**: n 29, tumour-adjacent. |
| Pancreas | [GSE155353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155353) | [33635431](https://pubmed.ncbi.nlm.nih.gov/33635431/) | 450k | 34 non-cancerous | Decade bins 40–49 … 70–79 | Tumour-adjacent / non-cancer pancreas of PDAC patients | Yes, 1.1 GB | Yes, 9 chips | **Meets with a caveat**: age in 10-y bins; adjacent. Possible donor overlap with GSE49149 not checked. |
| Pancreas | [GSE143209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143209) | [32736653](https://pubmed.ncbi.nlm.nih.gov/32736653/) | 450k | 64 islet preparations | No age field | Donor islets | Yes, 0.70 GB | Yes | **Fails**: no ages. |
| Uterus (endometrium) | [GSE223817](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223817) | [37587191](https://pubmed.ncbi.nlm.nih.gov/37587191/) | EPIC | 347 controls (no endometriosis); 637 cases | 18–55 (median 37) | Premenopausal, cycle phase recorded | Yes, 13.6 GB, controls ≈4.8 GB | Yes; 133 chips for controls, η² 0.43 (k/n 0.39) | **Meets with a caveat**: endometrium, not GTEx uterus body; 37-y span, all premenopausal; cycle phase must be a covariate. |
| Uterus (endometrium) | [GSE73950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73950) | [26759613](https://pubmed.ncbi.nlm.nih.gov/26759613/) | 450k | 17 without endometriosis | No age field | — | Yes | Yes | **Fails**. |
| Stomach | [GSE99553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99553) | [29574700](https://pubmed.ncbi.nlm.nih.gov/29574700/) | 450k | 42 healthy controls + 42 cases (normal mucosa) | No age field | Healthy control mucosa; H. pylori status recorded | Yes, 0.87 GB | Yes, 7 chips | **Fails as deposited**: no ages. Ages may be in the paper (case-control matched). |
| Stomach | [GSE103186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103186), [GSE127857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127857) | 29290541, 34199386 | 450k | ~55 normal biopsies; 24 adjacent | No age field | Normal / adjacent | Yes | Yes | **Fail**: no ages. |
| Oesophagus | [GSE72872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72872) | [26905591](https://pubmed.ncbi.nlm.nih.gov/26905591/) | 450k | 64 normal squamous | No age field | Adjacent (EAC patients) | Yes, 2.1 GB | Yes | **Fails**: no ages. |
| Oesophagus | [GSE89181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89181) / [GSE121932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121932) | 27795744 / 31962291 | 450k | 11 squamous / 15 normal | 40–93 / – | Barrett's / ESCC-adjacent | Yes | Yes / no (renamed) | **Fail**: n. |
| Adipose (cells) | [GSE222595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222595) | [37188674](https://pubmed.ncbi.nlm.nih.gov/37188674/) | 450k (discovery) + EPIC (replication) | 48 control subcutaneous adipocyte preps (24 per cohort and platform) | 22.6–70.3 for controls | Healthy controls vs extreme obesity; **isolated adipocytes** | Yes, 2.4 GB | Yes; controls 20 chips, η² 0.37 (k/n 0.21) | **Meets with a caveat**: isolated cells, two platforms of 24 each; visceral and subcutaneous preparations may share donors (no donor ID). |
| Skin | [GSE151617](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151617), [GSE268703](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268703), [GSE52980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52980) | 32660606, –, 25886480 | EPIC / EPIC / 450k | 16 biopsies / 53 post-mortem / 13 with IDATs | 6–84 / none / – | — | Yes / yes / mostly no | Yes | **Fail**: n, ages or IDATs. |
| Heart | [GSE197672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197672) | [36656640](https://pubmed.ncbi.nlm.nih.gov/36656640/) | 450k + EPIC | 7 non-failing LV | 28–76 | Heart failure series | Yes | Yes | **Fails** (super-series of GSE197670, already judged). |
| Tibial nerve, artery, spleen, adrenal, pituitary, salivary gland, vagina, bladder | — | — | — | — | — | — | — | — | **No series found.** E-MTAB-13433 has normal adrenal medulla only as controls for phaeochromocytoma; not opened. |

ArrayExpress turned up nothing new that meets. E-MTAB-13564 (blood plus subcutaneous and omental fat from bariatric
surgery, EPIC) and E-MTAB-11037/11841 (adipose, type 2 diabetes in obesity) were seen by title only and are not
assessed.

### End note: TCGA normal-adjacent IDATs (outside GEO and ArrayExpress)

The GDC files API (queried 2026-11-23) lists 450k "Masked Intensities" IDATs for "Solid Tissue Normal" samples as
**open access**, two files per array. Normal-adjacent arrays by project:

- KIRC 160 (ages 34–89 by `demographic.age_at_index`, read over the first 400 KIRC normal IDAT files; about 1.4 GB)
- BRCA 97
- THCA 56
- LIHC 50
- PRAD 50
- UCEC 46
- KIRP 45
- LUSC 42
- COAD 38
- LUAD 32
- ESCA 16
- PAAD 10

File names are UUIDs (`<uuid>_noid_Grn.idat`), so **chip IDs are not recoverable from names**; the TCGA aliquot
barcode's plate is the usable batch proxy. "Masked" means the probe intensities overlapping germline variants are
blanked. All of this is tumour-adjacent. It is the only route to arrays of kidney or prostate with more than 50 donors.

## Ranked shortlist: download first

1. **[GSE175458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175458), lung**: 145 non-diseased lungs,
   112 of them from a lung donor core, aged 14–88. EPIC, chip IDs in names, ≈2.0 GB for the normals. The only
   non-tumour lung series with IDATs and wide ages, and it gives an independent check of GSE314841's adjacent lung.
2. **[GSE88883](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88883), breast**: 100 healthy Komen donors,
   aged 18–82. 9 chips with age balanced across them; 1.0 GB.
3. **[GSE101961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101961), breast**: 121 cancer-free women,
   aged 17–76, 1.1 GB. It does not overlap GSE88883, so the two make a second, independent breast cohort, though
   with almost no donors over 65.
4. **[GSE151732](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151732), colon**: 128 healthy colonoscopy
   donors aged 34–80, right and left colon each. A healthy, non-adjacent counterpart to GSE101764. One side is
   ≈1.75 GB.
5. **[GSE79100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79100), kidney**: 31 normal kidneys aged
   18–78, with age spread over 5 chips; 0.95 GB. It is the only kidney series on GEO or ArrayExpress with IDATs and
   ages, but its tissue origin is unverified.
6. **[GSE97466](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97466), thyroid**: 50 tumour-adjacent
   thyroids aged 20–77, ≈0.46 GB for the normal subset. Do not add GSE86961 to it.

**Next in line:**

- [GSE107038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107038): 40 healthy livers aged 13–90, 0.51 GB.
- [GSE183015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183015): 58 benign, tumour-adjacent prostates
  aged 52–75. The only prostate option, short of the 30-y span.
- [GSE223817](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223817): 347 control endometria aged 18–55.

**Tissues with no usable public IDAT series:**

- **Testis** and **ovary**: none of any kind.
- **Skeletal muscle**: nothing meets.
- **Prostate**: only tumour-adjacent tissue with a 23-y span.
- **Stomach, oesophagus and skin**: IDATs exist but without ages.
- **Nerve, artery, spleen, adrenal, pituitary, salivary gland, vagina and bladder**: nothing.

Kidney and prostate gain size only through TCGA normal-adjacent IDATs, which are open but tumour-adjacent and carry no
chip IDs.

## What I could not check

- Papers were not read. Unverified points:
  - the tissue origin of GSE79100;
  - that the Pittsburgh controls in GSE175458 are donor lungs rather than resections;
  - whether ages for GSE287331, GSE269244, GSE99553 and GSE323385 are in paper supplements.
- Archive sizes are from the GEO download pages. Subset sizes (≈) are pro-rated by array count, not summed per file.
- Donor overlap was checked only for GSE88883 against GSE101961 (none), GSE88883 against GSE287331 (5 Komen IDs),
  GSE86961 against GSE97466 (full) and GSE216024 against GSE151732 (123 IDs). It was not checked for GSE49149 against
  GSE155353 or for GSE222595's depots.
- TCGA ages were read for KIRC only, from `demographic.age_at_index`.

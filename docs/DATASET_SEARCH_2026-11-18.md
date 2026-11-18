<!-- Search by a delegated agent, 2026-11-18 (RM40 items 4-5). Verified from GEO records and file lists only; papers not read; nothing downloaded. Statements marked unverified are to be checked before use. -->

# Dataset search for the shield tests: verified candidates per need

Six needs have at least one usable candidate; clonal haematopoiesis with VAF, a breast-specific prospective cohort, and calorie restriction at several ages with n ≥ 8 have nothing adequate in the open.

**What I checked.** For every accession below I read the GEO series record, `suppl/filelist.txt` and the per-sample records (characteristics and supplementary file names). Chip ids come from the IDAT file names; citations come from PubMed via the PMID in the GEO record. I did not read the papers, did not open SRA run tables, and downloaded no data.

**How to read the overlap numbers.**
- "Both tertiles" is the number of chips that hold at least one sample from the lowest and one from the highest age (or passage) tertile.
- "eta²" is the share of age variance that lies between chips; the value in brackets is what random allocation would give.
- Every candidate was found by GEO search, so none repeats your known or used lists.

## Need 1: second human cortex cohort, and sorted neurons versus glia

Sorted nuclei, ranked:

| Accession | Citation | Tissue, assay | n and ages | Raw data | Verdict |
|---|---|---|---|---|---|
| GSE279509 | Walker 2025, Brief Bioinform | Sorted nuclei, EPIC | NeuN+ 182, Sox10+ 168, double-neg 164, total 218; 18–108 y | 751 IDAT pairs, 11.2 GB | Meets with a caveat |
| GSE107729 | Mendizabal 2019, Genome Biol | BA46 NeuN+ and OLIG2+, WGBS | 53 NeuN, 42 OLIG2; 25–85 y | Per-CpG txt, 11.7 GB | Meets with a caveat |
| GSE66351 | Gasparoni 2018, Epigenetics Chromatin | Sorted and bulk cortex, 450k | 31 neuron, 31 glia, 128 bulk; 18–97 y | 190 IDAT pairs, 1.8 GB | Meets with a caveat |
| GSE112179 | Pai 2019, Nat Commun | Frontal cortex NeuN+, EPIC | 100 (33 controls); 23–77 y | 100 IDAT pairs, 1.6 GB | Meets with a caveat |
| GSE330111 | Xu 2026, bioRxiv | PFC Glu and MGE-GABA neurons, WGBS + oxBS | 310 libraries, 103 donors | Per-CpG tsv, 1.38 TB | Meets with a caveat |
| GSE234520 | Hannon 2024, BMC Biol | PFC sorted nuclei, EPIC | 28 NeuN+; 32–95 y | 112 IDAT pairs, 1.65 GB | Fails for age |

- **GSE279509**
  - The neuron-versus-glia contrast is within chip: 132 of 143 chips hold two or more cell types (737 of 751 samples).
  - Age in NeuN+ sits mostly between chips, close to chance: 115 chips, 66 with two or more NeuN+ samples, 10 with both tertiles, eta² 0.73 (0.63), median within-chip age range 12 y.
  - GEO lists 237 individual ids; the design text says 287 donors.
  - Brain region and diagnosis are not in the sample fields (unverified).
- **GSE107729**
  - 54 subject ids; 45 control and 50 schizophrenia libraries.
  - 89 libraries ran on HiSeq X Ten and 6 on HiSeq 2500 (all NeuN, 3 control and 3 schizophrenia).
  - Flowcell and lane are not in GEO; GEO links SRA SRP126677, which I did not open.
- **GSE66351**
  - Neurons are 16 control and 15 AD, on 9 chips, 4 with both tertiles, eta² 0.08 (0.27); most are old (lower quartile 73 y).
  - Bulk controls are 52, aged 18–88, on 10 chips, 7 with both tertiles.
  - Frontal and temporal bulk samples sit on separate chips (0 mixed), so regions cannot be compared.
- **GSE112179**
  - Diagnosis is mixed on 12 of 13 chips.
  - Age across all 100 samples: 8 of 13 chips with both tertiles, eta² 0.39 (0.12).
  - Controls alone span only 31–68 y (interquartile 42–51), and there are no glia.
- **GSE330111**
  - Ages are not in the GEO fields (unverified).
  - FASTQ is controlled access at NDA; the processed per-CpG files are open.
- **GSE234520**
  - NeuN+ samples are on 23 chips with 0 holding both tertiles, eta² 0.96.
  - Donor overlap with GSE279509 is unverified.

Bulk cortex, ranked:

| Accession | Citation | Tissue, assay | Controls and ages | Raw data | Verdict |
|---|---|---|---|---|---|
| GSE252501 | White 2024, Biol Psychiatry GOS | DLPFC and nucleus accumbens, EPIC | 59 DLPFC controls; 24.8–78.0 y | 232 IDAT pairs, 3.6 GB | Meets criteria |
| GSE349373 | No PMID in GEO | PFC, amygdala, hippocampus, EPIC | 101 PFC controls; 17.5–80.0 y | 904 IDAT pairs, 13.5 GB | Meets with a caveat |
| GSE203332 | PMID 40162278 (citation not fetched) | Frontal cortex, EPIC | 77 controls; 49–99 y | 492 IDAT pairs, 7.4 GB | Fails: old donors only |
| GSE341563 | No PMID in GEO | PFC and anterior cingulate, EPIC | 87 Lewy-stage-0 PFC; 47–97 y | 779 IDAT pairs, 11.6 GB | Fails: old donors only |

- **GSE252501**
  - Region is mixed on 30 of 30 chips and case/control on 28 of 30.
  - Age in DLPFC controls: eta² 0.48 (0.50), 7 chips with both tertiles.
  - Age in all 117 DLPFC samples: eta² 0.26 (0.25), 18 of 30 chips with both tertiles.
- **GSE349373**
  - Age in PFC controls: 43 chips, 7 with both tertiles, eta² 0.58 (0.42); diagnosis is mixed on 47 of 62 PFC chips.
  - The GEO record is defective: the tissue field says "Plasma cfDNA" and the design text describes a breast cancer study.
  - The contributors include LIBD investigators, so donors may overlap the already used GSE74193 (unverified).
- **Others that fail**
  - GSE197305 (Brains for Dementia Research, 1,221 arrays) has no age in GEO.
  - GSE41826, GSE80970 and GSE61431 have no IDATs in the file list.
- **Not tabulated**
  - GSE125895 (LIBD, 450k, four regions) and GSE134379 (middle temporal gyrus and cerebellum, 450k, 54–103 y) have IDATs, but I did not compute chip by age.

## Matched modalities from the same samples

Same cell or nucleus:

| Accession | Citation | System | Modalities | Ages and n | Files |
|---|---|---|---|---|---|
| GSE247988 (20 sub-series, e.g. GSE193274) | Chien 2024, Neuron | Human DLPFC BA46 neurons | Methylation + RNA, snmCT-seq | 11 donors: 5 young (23–30 y), 6 aged (70–74 y) | Per-nucleus allc; 365 GB for one sub-series |
| GSE121437 (GSE121436 + GSE121364) | Hernando-Herraez 2019, Nat Commun | Mouse muscle stem cells | scBS + scRNA | 3 mice at 2 mo, 3 at 24 mo; 275 methylomes | Bismark cov, 1.73 GB |
| GSE154370 (GSE154271 + GSE154368) | Castillo-Fernandez 2020, Aging Cell | Mouse GV oocytes | scBS + scRNA | 39 young, 41 old methylomes; 88 RNA | cov, 1.15 GB |
| GSE225173 (GSE225171 + GSE225172 + GSE225166) | No PMID in GEO | Mouse blood cells, plus bulk tissues | scM&T-seq, plus bulk RRBS | 16 single-cell samples; ages not in GEO | 8 tar archives, 3.15 GB |

- **GSE247988**
  - About 3,050 nuclei per sub-series and two sub-series per donor.
  - I extrapolate 6–7 TB for the whole set from the one sub-series I sized.
  - Plate ids are in the sample titles; coverage per nucleus is unverified.
- **GSE121437**
  - 40–48 cells per mouse.
  - The id mapping between the methylome and RNA of each cell is unverified.
- **GSE154370**
  - The number of mice and the exact ages are not in the GEO fields.

Same study, linkage across modalities by replicate label only (same animal unverified unless stated):

| Accession | Citation | System | Modalities | Ages and n |
|---|---|---|---|---|
| GSE233879 | Wang 2025, Cell Res | Mouse HSC, young and old | RRBS 28, RNA 28, ATAC 12, H3K4me3/H3K27ac 12, scRNA | Per-age n not tabulated |
| GSE47819 | Sun 2014, Cell Stem Cell | Mouse HSC | WGBS 13, ChIP 16 (H3K4me3, H3K27me3, H3K36me3), RNA 8 | 4 and 24 mo |
| GSE92486 | Hahn 2017, Genome Biol | Mouse liver, AL versus DR | WGBS 18, RNA 12 | 5 and 26 mo, n = 3 per group |
| GSE156558 (GSE156557 + GSE155740) | Chondronasiou 2022, Aging Cell | Mouse pancreas, spleen, liver | RRBS 45, plus one more modality | Young, old, old + OSKM; n = 5 each |
| GSE242660 (GSE129712 + GSE129767) | Krepelova 2025, Nat Aging | Mouse small-intestinal crypts | RRBS 24, plus one more modality | 2, 18, 26 mo; n = 8 each |
| GSE213723 | Vaidya 2023, Genome Biol | Mouse intestinal epithelium and 7 other tissues | RRBS 120, RNA 24 (intestine only) | 4–25 mo |
| GSE179880 | Wang 2022, Protein Cell | Mouse skin fibroblasts, P1–P7 | RNA, RRBS, ATAC | 3–4 replicates per passage |
| GSE165180 | Gill 2022, eLife | Human fibroblasts, transient reprogramming | EPIC 118 IDAT pairs, RNA 219 | Donors 0–53 y |
| GSE156409 | McCauley 2021, Nat Aging | Human MSC, PD12 versus PD32 | Histone ChIP 24, RNA 18, WGBS 2 | One WGBS library per state |
| GSE104408 | Adelman 2019, Cancer Discov | Human HSC | ERRBS 12, plus 5 more sub-series | Young 19–27 y, aged 63–71 y |
| GSE89275 | Cole 2017, Genome Biol | Mouse hepatocytes, dwarf, CR, rapamycin | Bisulfite 32, RNA 16 | Ages not in GEO fields |
| GSE171998 | Sahinyan 2022, eLife | Mouse muscle stem cells and myoblasts | WGBS 6, ATAC 4, RNA 28 | Young versus 21–26 mo |

- **GSE47819:** the mice were pooled, so the modalities come from the same study but different pools.
- **GSE92486:** the 26-month group was sequenced as separate paired-end and single-end sets.
- **GSE156558:** each RRBS sample is one mouse; I did not open GSE155740, so the second modality (RNA by the paper title) is unverified.
- **GSE242660:** I did not open GSE129767, so the second modality is unverified.
- **GSE165180:** chip by group was not tabulated.
- **GSE104408:** the other five sub-series were not opened, so their modalities are unverified.
- **GSE295059** (mouse intestinal organoids, passages P0–P27): the design names RRBS and RNA-seq, but no RNA series is linked in the record.
- **GSE152798** (Benjamin 2023, Cell Metab; muscle stem cells at 4 versus 22–23 mo; RNA 14, bisulfite 8): no processed files in the super-series listing.
- **GSE181797** (mouse kidney and liver aging) holds only RNA and ATAC, so it does not qualify.

## Need 2: selection test in culture

| Accession | Citation | Cells, assay | Passages and n | Raw data | Verdict |
|---|---|---|---|---|---|
| GSE268145 | No PMID in GEO | BJ fibroblasts with and without hTERT, EPIC | PDL 2–104, two trajectories | 66 IDAT pairs, 1.18 GB | Meets with a caveat |
| GSE202554 | Higgins-Chen 2022, Nat Aging | Fetal astrocytes, EPIC | P1–10, 47 bulk samples, 2 donors | 63 IDAT pairs, 1.09 GB | Meets with a caveat |
| GSE226079 | Minteer 2023, Sci Adv | Astrocytes, EPIC | hTERT arm P12–27 (45); primary arm P2–10 (25) | 70 IDAT pairs, 0.98 GB | hTERT arm meets; primary arm fails |
| GSE227160 | Kwiatkowska 2023, Cells | Fibroblasts and MSC, EPIC | 6 donors × young, middle, senescent, irradiated, doxorubicin | 48 IDAT pairs, 0.88 GB | Meets with a caveat |
| GSE295059 | No PMID in GEO | Mouse intestinal organoids, RRBS | P0–P27 from 4- and 24-month mice, 70 libraries | Bismark cov, 5.1 GB | Meets with a caveat |
| GSE116375 | Franzen 2021, Commun Biol | MSC and HUVEC, 450k | 16 arrays on 2 chips | 16 IDAT pairs | Fails: too small |

- **GSE268145:** 17 chips, and nearly every chip mixes early and late doublings (for example PDL 2, 3, 20, 30 and 104 on one chip). It is one cell line with no arrested control.
- **GSE202554:** 8 of 8 chips hold both passage tertiles, eta² 0.07. The span is short (15 cumulative doublings or fewer). It also has 16 sorted fractions (senescent, galactosidase high and low).
- **GSE226079:** the hTERT arm has 9 of 12 chips with both tertiles, eta² 0.16 (0.25). The primary arm has 0 of 12, eta² 0.89 (0.46), so passage is confounded with chip there. Each arm is one donor.
- **GSE227160:** all 8 chips hold several treatments, and every chip holds several donors. There are only three passage levels and 6 donors.
- **GSE295059:** library batch is not recorded, and the number of animals is small (unverified count).
- I found nothing new with barcoded lineages or parallel arrested controls; the already used GSE197512 remains the only series of that kind.

## Need 3: clonal expansions in vivo

| Accession | Citation | System, assay | n and ages | Verdict |
|---|---|---|---|---|
| GSE141254 | Lewis 2020, Cell Mol Gastroenterol Hepatol | Human intestinal crypts, mucosa, organoids; 450k and EPIC | 36 crypt, 27 mucosa, 21 organoid; 0.5–82.5 y | Meets with a caveat |
| GSE290999 | Masika 2026, Nat Commun | Single mouse crypts and single human hair follicles, WGBS | 80 crypts from 2 mice (22 and 28 wk); 50 follicles from one 53-year-old | Fails: no young/old contrast |
| GSE282971 | Scherer 2025, Nature | Mouse bone marrow, scTAM-seq | 2 young and 2 old native-aging samples, plus barcoded transplants | Fails for promoters |
| GSE210435 | PMID 39189614 (citation not fetched) | Clonal haematopoiesis with COVID-19, EPIC | 18 arrays | Fails: too small |

- **GSE141254**
  - 84 IDAT pairs, 1.28 GB.
  - Age is spread within chips: 5 of 5 chips (450k, 23 crypt samples) and 3 of 4 chips (EPIC, 13 crypt samples) hold both tertiles.
  - The two platforms must be analysed separately.
  - Whether a crypt sample is a single crypt or a pooled prep is unverified.
- **GSE290999:** the deposited files are aggregated bedgraphs totalling 3 MB, not per-CpG calls.
- **GSE282971:** the assay is a targeted CpG panel, so there is no promoter-wide coverage; it is useful only for clonal structure.
- **GSE222486** (Li 2023, Cell; barcoded mouse line with joint single-cell profiling): I did not find an age contrast in the record.

## Need 4: tissue turnover

| Accession | Citation | Tissues, assay | Ages and n | Batch overlap | Verdict |
|---|---|---|---|---|---|
| GSE225166 | No PMID in GEO | 14 mouse tissues and cell types, RRBS | 3, 15, 24 mo (plus 6, 9 mo); 43 mice | Each of 16 lanes holds 3–5 ages | Meets for four tissues |
| GSE213628 | Vaidya 2023, Genome Biol | Intestinal epithelium, blood, spleen, heart, kidney, liver, lung, muscle; RRBS | 4, 12, 18, 24–25 mo; 120 libraries | Not recorded | Meets with a caveat |
| GSE184410 | Zhou 2022, Cell Genom | Mouse methylation array, 44 tissues | 1,239 IDAT pairs; old (≥ 540 d) n ≤ 11 per tissue | Young and old share a chip on 4 of 50 (small intestine), 0 of 21 (colon) | Fails |
| GSE290585 | Lee 2025, Genome Biol | Mouse methylation array, 29 tissues | 8 wk versus 76 wk; old n ≤ 6 per tissue | 0–1 chips per tissue hold both | Fails: age confounded with chip |
| GSE93957 | Stubbs 2017, Genome Biol | Cortex, heart, liver, lung; RRBS | 1–41 wk | Not checked | Fails: no old animals |

- **GSE225166**
  - Liver, heart, cortex and lung have 10–12 animals at each of 3, 15 and 24 months.
  - Kidney, muscle, skin, cerebellum, intestinal stem cells, lung macrophages and CD4, CD8 and B cells have only 2 per age, which fails the n criterion for the renewing tissues.
  - The same mouse ids recur across tissues: 25 mice give 4 tissues and 9 mice give 8–13.
  - Lanes are readable from the file names (lane 1 holds 7, 7 and 5 samples of the three ages).
  - Bismark cov files, 5.0 GB.
- **GSE213628**
  - n per tissue and age is small.
  - cov files, 3.3 GB.
  - Whether the animal letters in the titles mark the same animal across tissues is unverified.
- **Single tissue, several ages**
  - GSE304754: mouse hippocampus RRBS, 3–28 mo, n = 16–18 per age, but no processed files on GEO.
  - GSE129712: crypts at 2, 18 and 26 mo.

## Need 5: lifespan intervention over time

| Accession | Citation | System, assay | Arms and time points | Raw data | Verdict |
|---|---|---|---|---|---|
| GSE316651 | No PMID in GEO | Mouse PBMC, mouse methylation array | 89 mice, up to 5 time points; control 153, NMN 162 | 315 IDAT pairs, 1.7 GB | Meets with a caveat |
| GSE92486 | Hahn 2017, Genome Biol | Mouse liver, WGBS | AL versus DR at 5 and 26 mo, n = 3 | cov, 3.1 GB | Fails: n = 3 per group |
| GSE89275 | Cole 2017, Genome Biol | Mouse hepatocytes, bisulfite | Dwarf, CR, rapamycin | txt, 3.9 GB | Fails: small n, ages not in GEO |
| GSE75499 | Maegawa 2017, Nat Commun | Mouse, monkey, human blood | 8 CR samples | 142 MB | Fails: not bisulfite per-CpG, small |

- **GSE316651**
  - Arms are mixed on 29 of 31 chips.
  - 24 chips hold three or more time points, and 23 hold baseline together with the fourth or fifth.
  - Samples per time point are 83, 76, 79, 59 and 18.
  - NMN is not an established lifespan intervention, and ages at each time point are not in GEO.
  - The record is inconsistent: titles say "cell-free DNA" while descriptions say PBMC.
- A reprogramming series (GSE156557; pancreas, spleen, liver; young, old and old + OSKM, n = 5 each) is listed under matched modalities.

## Need 6: prospective blood cohort with incident cancer

| Accession | Citation | Cohort, assay | Cases and non-cases | Raw data | Verdict |
|---|---|---|---|---|---|
| GSE280206 | Suárez-Pérez 2024, Genome Med | Stroke patients, blood at stroke, 15-year follow-up; 450k and EPIC | 83 incident cancers, 565 none | 648 IDAT pairs, 7.2 GB | Meets with a caveat |

- 80 of 83 cases share a chip with a non-case, and 59 of 157 chips hold both.
- The two platforms never share a chip: 450k has 410 samples (53 cases) and EPIC has 238 (30 cases), so analyse them separately.
- Age, cancer type and time to diagnosis are not in the GEO fields, and every participant is a stroke patient.
- I found no breast-specific prospective cohort with open IDATs.

## Need 7: single-cell methylation with age

| Accession | Citation | Species, cells | Ages and donors | Cells | Files | Verdict |
|---|---|---|---|---|---|---|
| GSE247988 | Chien 2024, Neuron | Human DLPFC neurons, snmCT-seq | 5 young (23–30 y), 6 aged (70–74 y) | ~6,100 per donor | Per-nucleus allc, ~6–7 TB (extrapolated) | Meets with a caveat: size |
| GSE121436 | Hernando-Herraez 2019, Nat Commun | Mouse muscle stem cells, scBS | 3 mice at 2 mo, 3 at 24 mo | 275 | cov, 1.73 GB | Meets with a caveat: 3 mice per age |
| GSE154271 | Castillo-Fernandez 2020, Aging Cell | Mouse GV oocytes, scBS | Young and old groups | 80 | cov, 1.15 GB | Meets with a caveat |
| GSE292803, GSE313770 | Zeng 2026, Cell | Mouse brain, 36 cell types, female and male | 8 wk, 9 mo, 18 mo | Pseudobulk only | allc, 858 GB and 585 GB | Fails for per-cell work |
| GSE225171 | No PMID in GEO | Mouse blood, scM&T-seq | Ages not in GEO | Unknown | 8 tar archives, 3.15 GB | Unverified ages |
| GSE282971 | Scherer 2025, Nature | Mouse bone marrow, scTAM-seq | 2 young, 2 old | Count matrices | 130 MB | Fails for promoters |

- Coverage per cell, and whether promoters and CpG islands are covered well enough for per-promoter heterogeneity, is unverified for every row; it needs the files or the papers.
- Batch: plate ids are in the titles for GSE247988, and mouse ids are in the fields for GSE121436; nothing else is recorded.
- **GSE154271:** the number of mice and the exact ages are not in the GEO fields.
- **GSE292803 and GSE313770:** these hold one merged file per cell type and age, so they suit neuron-type age contrasts but not cell-to-cell heterogeneity; the number of animals is unverified.
- **Leads not checked**
  - GSE211786 with GSE209656: scNMT-seq of the adult neural stem cell lineage, 2,084 cells; ages not checked.
  - GSE314211: mouse brain 5mC, 5hmC and RNA; titles show one age only, and no files are listed.
- **Starting points I could not confirm**
  - I found no GEO series for the Gravina 2016 hepatocyte data.
  - I did not verify which datasets Trapp 2021 used.

## The five to fetch first

1. **GSE279509** (11.2 GB). It is the only large sorted-neuron and glia series with IDATs across 18–108 y, and the cell-type contrast is within chip.
2. **GSE252501** (3.6 GB). It is the cleanest independent DLPFC replicate: 59 controls aged 25–78, with age no more chip-bound than chance.
3. **GSE225166 with GSE225171** (5.0 + 3.15 GB). It gives the same mice across tissues at 3, 15 and 24 months with ages mixed on every lane, plus single-cell methylation with RNA from the same study.
4. **GSE121437** (about 2 GB). It gives single-cell methylation and RNA from the same muscle stem cells, young versus old, and is small enough to test the per-promoter heterogeneity method.
5. **GSE268145 with GSE202554** (1.2 + 1.1 GB). Both are culture series where passage is spread across chips, unlike the withdrawn ones.

Next in line are GSE107729 (11.7 GB, WGBS neurons versus oligodendrocytes, a platform-independent check of item 1), GSE233879 (2.1 GB, the richest bulk multi-omic HSC aging set), GSE280206 and GSE316651.

## Needs with nothing adequate in the open

- **Clonal haematopoiesis with VAF and methylation in healthy people:** only an 18-array series exists on GEO. My recollection that the cohort data sit behind dbGaP or EGA is unverified.
- **Single-cell-derived colonies or organoids, young versus old, genome-wide:** the only crypt series with a wide age range (GSE141254) is small and split over two platforms; the single-crypt WGBS has no age contrast and only aggregated files.
- **Breast-specific prospective cohort with IDATs:** none found; GSE280206 is all cancers in stroke patients.
- **Calorie restriction at several ages with n ≥ 8 per arm:** nothing beyond the already used GSE80672; the liver DR series has n = 3.
- **Culture with barcoded lineages or arrested controls:** nothing beyond the already used GSE197512.
- **Renewing tissues with n ≥ 8 per age and recorded batch:** GSE225166 has only 2 per age for intestine, skin and immune cells; GSE213628 has no batch record.


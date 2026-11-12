# Public data for partial reprogramming and methylation-machinery models (roadmap RM17, RM12)

Survey of public data, 2026-11-12, by a delegated agent for the lead; unverified items marked [U]. Not yet a result.

# Public methylation data for the "shield" questions: partial reprogramming (A) and methylation-machinery models (B)

Survey date: 2026-11-12. Sources: NCBI GEO (E-utilities esearch/esummary, series and sample SOFT records, FTP listings for file sizes), PubMed abstracts, PMC/Europe PMC full text for data-availability statements. Nothing was downloaded beyond metadata and directory listings.

Legend: **[V]** = verified in the GEO record or FTP listing (design, n, platform, file sizes). **[P]** = from the paper or abstract only. **[U]** = could not verify.
Sizes are the GEO supplementary files (GB/MB as listed on ftp.ncbi.nlm.nih.gov). "RAW.tar" can usually be fetched per sample (GSM supplementary files), which is how to stay under the 50 MB limit.
All GEO series listed are public (open access, no application).

Platform keys: GPL28271 = HorvathMammalMethylChip40 (same as the Mammalian Methylation Consortium, MMC); GPL32543 = HorvathMammalMethylChip320; GPL30650/GPL31950/GPL32685/GPL35806 = Illumina Infinium Mouse Methylation (MM285) BeadChip variants; GPL21145/GPL23976 = EPIC; GPL13534 = 450k.

---

## Question A. Partial / transient reprogramming

### Ranked candidates

| Rank | Accession | Study | Design (groups, n) | Tissue / cells | Platform | Processed data | Notes |
|---|---|---|---|---|---|---|---|
| 1 | GSE190665 | Browder et al. 2022 Nat Aging (PMID 37118377) | 7 groups per tissue: 4F-young (3 mo, n=4), 4F-old untreated (n=2-7), B6-old (n=6-7), B6+Dox (n=2), 4F+1 mo Dox (n=1-3), 4F+7 mo Dox (n=4-5), 4F+10 mo Dox (n=2-3); old ≈ 21-27 mo; both sexes | kidney, liver, lung, muscle, skin, spleen (168 arrays) | GPL28271 (Mammal40, same as MMC) | `datBetaNormalized.csv.gz` 49 MB; RAW.tar 126 MB | [V]. Best fit: same array as the MMC Tet3 analysis, young reference (4F-young) and untreated old 4F + WT controls in the same study. Per-tissue n is small (2-7); pool tissues or use per-tissue paired contrasts. |
| 2 | GSE156557 | Chondronasiou et al. 2022 Aging Cell (PMID 35235716) | young control (13 wk, n=5), old control (52-54 wk, n=5), old OSKM single cycle + 2 wk off-dox (54 wk, n=5) per tissue; all female | pancreas, liver, spleen | RRBS (GPL17021) | per-sample Bismark `.cov.gz` ~20-22 MB each; RAW.tar 971 MB | [V]. Clean three-arm design with an in-study young reference; per-sample files fit the size limit. Old = 1 yr only (middle-aged). |
| 3 | GSE165179 | Gill et al. 2022 eLife (PMID 35390271), MPTR | per donor (38 y, 53 y) and per reprogramming length (13/15/17 d): transiently reprogrammed (n=2), failed-to-reprogram (n=2-4), negative (mock) control (n=2-4); plus intermediates and day-0 fibroblasts; 96 arrays | dermal fibroblasts | EPIC (GPL21145) | `Matrix_processed_transient.txt.gz` 530 MB (over the 50 MB limit); RAW.tar 1.5 GB (per-sample IDATs ~7-8 MB each) | [V]. The key human dataset; no young donor in the series, so a young-fibroblast reference must come from elsewhere (or use donor-matched mock controls as reference). Needs IDAT processing per sample to stay under the limit. |
| 4 | GSE247179 | Mitchell et al. 2024 eLife (PMID 38517750), chemical (2c, 7c) | young (4 mo) and old (20 mo) male mouse fibroblasts × vehicle / 2c / 7c, n=5 each; plus iPSC n=3 per age | dermal fibroblasts (C57BL/6) | GPL32543 (Mammal320k) | `data_sesame_40k.RDS.gz` 74 MB (Mammal40 subset), `data_sesame.RDS.gz` 150 MB; RAW.tar 208 MB | [V]. Only public methylation data on chemical partial reprogramming I found. Young vehicle = in-study reference; 40k subset aligns with MMC probes. R RDS format. |
| 5 | GSE252401 | Horvath et al. 2025 GeroScience (PMID 39037528), hippocampal OSKM adenovector | young (3.5 mo, n=6), old control placebo vector (25-27 mo, n=7), old OSKM 39 d (n=7); female rats | hippocampus (rat) | GPL28271 (Mammal40) | RAW.tar 15 MB (IDATs only; no beta matrix) | [V]. Brain, same array as MMC. Batch caveat: 5 of 6 young on chip ...005, most old on chip ...010. Rat, not mouse. |
| 6 | GSE231658 (SuperSeries GSE289856) | Camacho, Koldobskiy et al. 2026 Mol Syst Biol (PMID 41667860) | skin: 4F cyclic OSKM 15→22 mo (n=5), 4F untreated 22 mo (n=4), 4F 3 mo (n=3); tail epidermis WT 3 mo vs 22 mo (n=5 each, no OSKM) | dorsal skin; tail epidermis | WGBS (GPL21103, GPL34290) | per-sample bigWigs only: MML (mean methylation), NME (normalized methylation entropy), JSD; ~100 MB each; RAW.tar 5.4 GB; raw reads in SRA | [V]. The only reprogramming dataset with published methylation entropy. No per-CpG coverage files; MML bigWig gives region-level means (informME), usable for promoter-level scoring. |
| 7 | GSE142439 | Sarkar et al. 2020 Nat Commun (PMID 32210226), mRNA OSKMLN 4 d | paired untreated vs treated: fibroblasts from 4 donors (61-69 y), endothelial cells from 4 donors (47-49 y) | skin fibroblasts, vein endothelial cells | EPIC (GPL21145) | IDATs only, RAW.tar 387 MB (~7 MB per IDAT); GEO characteristics carry Horvath and skin-blood clock ages | [V]. Small n (4 pairs per cell type), no young reference. |
| 8 | GSE147436 | Lu et al. 2020 Nature (PMID 33268865) | SH-SY5Y neurons ± OSK (AAV), ± vincristine damage: n=3 per group, 6 at 9-d recovery | differentiated SH-SY5Y (cell line) | EPIC (GPL23976) | IDATs, RAW.tar 338 MB | [V]. Cell line, injury model; weak fit. Mouse RGC methylation from this paper not found in GEO [U]. |
| 9 | GSE288182 | Local brain AAV-OSK in APP/PS1 (PMID 42210397, 2026) | WT mock, APP/PS1 mock, APP/PS1 OSK, n=3 each | brain | RRBS (GPL24247) | per-sample annotated site tables 330-390 MB each; RAW.tar 3.2 GB | [V]. Disease model, n=3; low priority. |
| 10 | GSE165178 | Gill 2022, Sendai full reprogramming time course | reprogramming vs failing-to-reprogram, days 9/11/15, 3 donors (29-53 y), n=1-2 | fibroblasts | EPIC | processed 143 MB; RAW.tar 467 MB | [V]. Full (not transient) reprogramming; context only. |
| 11 | GSE326022 | SRSF1 screen (2026), not OSK | mouse liver AAV SRSF1 vs GFP cyclic dox (n=5 vs 5); human fibroblasts ± SRSF1 / rapamycin | liver; fibroblasts | listed as EPIC GPL21145 for all samples [U: mouse liver on EPIC is odd] | processed 2.3-295 MB | [V] design. An "alternative rejuvenation" comparator, not reprogramming. |

### Requested studies with no usable public methylation data

| Study | Status |
|---|---|
| Ocampo et al. 2016 Cell (PMID 27984723) | No GEO link; no DNA methylation profiling (H3K9me3/H4K20me3 by imaging) [P]. |
| Macip et al. 2024 Cell Reprogram (PMID 38381405; PMC10909732), OSK AAV in 124-wk mice | Methylation run by the Clock Foundation (liver, heart; HEK001 keratinocytes with 2 technical repeats); **no data availability statement, no accession** [V from full text]. Ask the authors / Rejuvenate Bio. |
| Yang et al. 2023 Cell "Loss of epigenetic information" (ICE mice; PMID 36638792; PMC10166133) | RRBS only as raw reads in SRA: BioProjects PRJNA554729, PRJNA655981; no GEO processed data [V from full text]. OSK methylation only in ICE MEFs (n=8 vs 8); in vivo OSK readout was H3K9me3, not DNAm [P]. |
| Yang et al. 2023 Aging "Chemically induced reprogramming" (PMC10373966) | Transcriptomic clocks only; no DNA methylation measured [P]. If "Yang 2023 Cell chemical" meant this paper, there is nothing to take; GSE247179 is the usable chemical dataset. |
| Alle et al. 2022 Aging Cell (single early-life OSKM; PMC9649606) | "Available from the corresponding authors on reasonable request" [V]. |
| Olova et al. 2019 Aging Cell (PMID 30450724) | Reanalysis of public human reprogramming time-course arrays (likely GSE54848, 450k) [U on exact accession]; full reprogramming, not transient. |
| Parras et al. 2023 Nat Aging; Rodriguez-Matellan 2020 | No methylation datasets linked [V: no GEO links]. |

### Published findings on clocks / entropy (A), one line each
- Gill 2022: MPTR rejuvenated the DNAm (skin & blood) clock by about 30 years, alongside the transcriptome clock [P].
- Sarkar 2020: 4-day OSKMLN mRNA lowered Horvath and skin-blood clock ages modestly (GEO-listed values: 0-8 y per pair) [V values in GEO].
- Browder 2022: long-term (7-10 mo) cyclic OSKM reversed mammalian epigenetic clocks in some tissues (skin, kidney) [P].
- Chondronasiou 2022: one OSKM cycle reversed a fraction of age-associated methylation changes in pancreas, liver, spleen and blood [P].
- Camacho/Koldobskiy 2026: age- and rejuvenation-related changes converge on PRC2 targets, with **increased methylation level and entropy** there; hypomethylated, highly entropic regions bounded by H3K9me2 LOCKs [P]. Only entropy result in this set.
- Mitchell 2024: 7c chemical cocktail lowered epigenetic and transcriptomic clock age of mouse fibroblasts [P].
- Horvath 2025 (rat hippocampus): marginally significant clock rejuvenation; partial reversal of age-related hypermethylation; chromatin-state analysis suggested rejuvenation [P].
- Macip 2024: LUC clock lower in liver (p=0.014) and heart (p=0.041) of OSK mice; no public data [P].
- Yang 2023 Cell: ICE accelerates DNAm clocks ~50%; OSK reverses DNAm age in ICE MEFs up to 57% [P].
- Olova 2019: epigenetic age declines steadily early in reprogramming, before loss of somatic identity [P].

---

## Question B. Second methylation-machinery models (loss or gain of DNMT/TET/guides)

### Ranked candidates

| Rank | Accession | Model | Groups, n | Tissue, age | Platform | Processed data | Notes |
|---|---|---|---|---|---|---|---|
| 1 | GSE271375 | Tet1 KO and Tet1 catalytic-dead (HxD) vs WT; Foong et al. 2024 Epigenetics (PMID 38970823) | WT 8, HxD 8, KO 8 animals (48 arrays = BS + bACE/5hmC per animal) | adult cortex [age not in GEO record, U] | Mouse MM285 (GPL31950) | `Tet1_array_betas.csv.gz` 145 MB (over limit; per-sample IDATs in RAW.tar 245 MB, ~5 MB each) | [V]. Best Tet analogue to the MMC Tet3 brain result: brain, n=8, three genotypes, 5mC and 5hmC separable. Littermate status [U]. Promoter coverage of MM285 is good for housekeeping CGIs. |
| 2 | GSE324160 (SuperSeries GSE324236) | Dnmt3a W326R/+ **gain-of-function** (Heyn-Sproul-Jackson progeria model); Sarni et al. 2026 Nat Genet (PMID 42286141) | +/+ vs W326R/+, 3 per genotype per tissue per age (12 per genotype per tissue) | bone marrow, intestine, liver at 4 d, 23 d, 3 mo, 1 yr | Mouse Methylation BeadChip (GPL32685) | `sarni_et_al_mouse_infinium_beta.tsv.gz` 147 MB; detection p 6.4 MB; RAW.tar 398 MB | [V]. Gain-of-function "accelerated ageing" model with an age series and same-batch controls; ideal opposite-direction test. Companion EM-seq of HSCs (GSE324227, 3 pooled vs 3). |
| 3 | GSE147898 | Dnmt3a heterozygous KO vs WT; Christian et al. 2020 Cell Rep (PMID 33238114) | 8-wk cortex: WT 23 vs Het 14 (both sexes); developmental series 0-20 wk n=2-3; plus cerebellum, hippocampus, striatum, liver (n=6-8 per region) and many cultured-neuron mutants | brain regions, liver | sparse/low-coverage WGBS (GPL19057, GPL21493) | only summary tables (`mdt_table_GEO2.tsv.gz` 5.9 MB, per-sample global levels); per-CpG data only as raw SRA reads | [V]. Large n but shallow coverage: promoter-level per-sample scores would need pooling or aggregation. Main finding is non-CG loss. |
| 4 | GSE225372 | Dnmt3a P900L/+ and R878H/+ (TBRS / CH-hotspot) vs WT littermates; Gabel lab 2023 | cortex 8 wk: each line ~19 mutant vs ~19 WT (low depth, some high depth); 1-2 per brain region | cortex, cerebellum, hippocampus, striatum | WGBS low and high depth (+ ChIP, RNA) | DMR BEDs and an xlsx summary only; RAW.tar 7.6 GB (mostly ChIP) | [V]. Same caveat as rank 3; publication [U]. |
| 5 | GSE299673 | Dnmt3a overexpression in cardiomyocytes (AAV9-Tnnt2) vs GFP AAV | 5 vs 5, injected at 3 mo, harvested 3 mo later | sorted cardiomyocyte nuclei | RRBS (GPL34290) | Bismark `.cov.gz` 25-37 MB per sample; RAW.tar 255 MB | [V]. Gain model, adequate n, files fit the limit. Unpublished [U]. |
| 6 | GSE262342 + GSE244432 | Dnmt3a muscle transgenic (overexpression) and Dnmt3a muscle KO; Oyabu et al. 2025 iScience (PMID 40151644) | gastrocnemius: Tg vs WT at 3 mo (4 vs 4) and 26 mo (6 vs 6); tibialis anterior: Tg 3, KO 3, WT 6 | skeletal muscle | WGBS (PBAT) | per-sample `cpg.cover/methyl.bedGraph.gz` ~190-245 MB each (over limit); RAW.tar 8.1 GB / 4.7 GB | [V]. Has young and old in both genotypes: lets you ask whether gain-of-methylation mimics ageing at essential promoters. Files are large. |
| 7 | GSE117728 | Dnmt3a PWWP D329A (gain at bivalent/H3K27me3 domains) and Dnmt3a +/-; Sendžikaitė et al. 2019 Nat Commun (PMID 31015495) | hypothalamus ~14 wk: +/+ 8 vs Δ/D329A 7 (+ Δ/+ 3, +/D329A 2); P1, P25 n=2; liver, pituitary n=2-3 | hypothalamus, liver, pituitary, E7.5 epiblast | WGBS (GPL13112) | no per-sample processed methylation in series files; RAW.tar 9.6 GB | [V]. Good n in hypothalamus only. |
| 8 | GSE255947 | Dnmt3a R878H/+ (mouse CH-hotspot) vs +/+ | CD11b+: 6 vs 4; CD4+: 4 vs 4 | splenic myeloid and T cells | WGBS | merged metilene tables 8-11 MB | [V]. Mouse model of a CH mutation; include only if CH-like data are acceptable. |
| 9 | GSE272913 | Dnmt3a Het, KO, R832, E752 (Vav-Cre) vs control | HSC/KSL/BM, 2-4 per genotype | haematopoietic | WGBS/EM-seq + RNA | RAW.tar 7.6 GB | [V]. Small n per cell type. |
| 10 | GSE183675 | Tet2 KO vs WT, young vs old (Hong et al. 2023 Nat Aging, PMID 37884767) | WGBS only 6 samples (rest are ChIP, RNA, Hi-C, CMS-IP) | LSK cells | WGBS (n per group ~1-2) + CMS-IP bigWigs (~370 MB each) | [V]. Only Tet2 KO × age design found, but WGBS n is too small. |
| 11 | GSE200227 | Tet2/Tet3 double KO (Villin-Cre) vs fl/fl | 3 vs 3, 11 wk | small-intestine crypts | WGBS | per-sample cov/ratio BEDs ~260 MB; RAW.tar 2.9 GB | [V]. Tet2+Tet3 double loss; n=3. |
| 12 | GSE84549 | Uhrf1 cKO (Emx1-Cre) vs Uhrf1+/- controls | 3 vs 3 | E16 cortex | RRBS | RAW.tar | [V]. Embryonic, n=3; only Uhrf1 in vivo tissue methylation found. |
| 13 | GSE175380 | Dnmt1 + Dnmt3a deletion (Syn-Cre AAV) vs WT | 4 vs 4 | mPFC pyramidal neurons | WGBS (GPL21493) | DMR table + RAW.tar | [V]. Adult neurons, double KO. |
| 14 | GSE142834 | Dnmt1 hypomorphic missense (from zebrafish/medaka screen) | 6 vs 6 (listed as "dnmt1-/-") | lin- bone marrow progenitors, BALB/c | WGBS | DMR tables | [V] design; exact allele [U]. |

Human germline (non-CH) cohorts:

| Accession | Condition | Groups, n | Tissue | Platform | Processed | Notes |
|---|---|---|---|---|---|---|
| GSE128801 | DNMT3A R771Q Tatton-Brown-Rahman (Amish sibship); Jeffries et al. 2019 Genome Res | 4 TBRS + 1 mosaic vs 7 WT relatives | blood | 450k | `average_beta_matrix.txt.gz` 27 MB | [V]. Small but open and under 50 MB; paper reports strongly accelerated DNAm age. The 15 further TBRS patients in the paper are not in this series [U]. |
| GSE324230 | DNMT3A gain-of-function HESJAS (Sarni 2026) | 7 patients (0.7-31 y) vs 9 controls | blood | EPIC | beta 103 MB | [V]. Pairs with mouse GSE324160; age-mismatched (patients mostly children). |
| GSE95040 | ICF syndrome (ICF1 = DNMT3B; ICF2-4 = ZBTB24, CDCA7, HELLS) | ICF1 8, ICF2 4, ICF3 2, ICF4 1, controls 10 | whole blood | 450k | raw-data matrix 60 MB | [V]. DNMT3B loss of function in humans. |
| GSE126890 | DNMT1 ADCA-DN | blood 3 vs 3; fibroblasts 6 vs 4; iPSC and induced neurons 12 vs 8 | blood, fibroblasts, iPSC, iN | sequencing (GPL21697 NextSeq; assay [U]) | RAW.tar | [V] design. |

### Searched, nothing usable found
- **Tet3 heterozygous outside MMC**: none. Only maternal Tet3 KO zygotes/blastocysts (GSE62719, RRBS 2 vs 2), Tet3 knockdown neural progenitors (GSE123110), Purkinje-neuron Tet3 work (GSE166423, expression-centred). GEO does not index sample characteristics, so other knockout cohorts inside the MMC super-series GSE223748 could not be searched from here [U].
- **TET3 Beck-Fahrner syndrome** (human, hypermethylation episignature; Levy et al. 2021 NPJ Genom Med): no public dataset found (episignature data typically held by EpiSign) [U].
- **Germline human TET2 deficiency**: no methylation dataset found.
- **Dnmt3b** adult mouse tissue loss-of-function: none (only ESC/embryo work); human ICF1 above is the substitute.
- **Dnmt1 hypomorph** (chip/N alleles, Apc-Min intestinal studies): no tissue methylation series found besides GSE142834.
- **Cfp1/Cxxc1**: only germ cells (GSE240280 spermatocyte RRBS 2 vs 4; GSE159581 oocytes). **Kdm2b**: only ESC RRBS (GSE180558, 1-2 per genotype). No somatic tissue.
- **GSE218055** ("Dnmt1 A560V cortex", 5 het, 2 WT, 2 hom, 12 mo) is RNA bisulfite sequencing (m5C in RNA, meRanCall output), not DNA methylation. Excluded.
- GSE125008 (Uhrf1 hepatocyte KO, liver regeneration): sub-series are expression and "genome binding"; no DNA methylation sub-series identified [U].

### Published findings on clocks / entropy / noise (B), one line each
- Jeffries 2019 (TBRS, DNMT3A loss of function): highly accelerated DNAm (Horvath) age; strongest in an R882C carrier [P].
- Sarni 2026 (HESJAS, DNMT3A gain of function): recapitulates age-related methylation gains and phenocopies ageing in humans and mice [P]; clock details not checked [U].
- Christian 2020 (Dnmt3a+/- brain): global loss of neuronal non-CG methylation; no clock reported [P].
- Oyabu 2025 (Dnmt3a muscle overexpression): ageing muscle gains methylation; overexpression gives an ageing-like phenotype; no clock reported [P].
- Hong 2023 (Tet2 KO HSPC, young vs old): PMD methylation and heterochromatin relocalisation change with age; no clock/entropy [P].
- Foong 2024 (Tet1 KO/HxD cortex): KO and catalytic-dead perturb distinct CpG subsets; no clock/entropy [P].
- Sendžikaitė 2019 (Dnmt3a D329A): progressive hypermethylation of H3K27me3/bivalent domains with age (P1 → P25 → adult) [P], a Polycomb-target gain pattern that resembles ageing.
- No published entropy or epigenetic-noise analysis was found for any Tet or Dnmt mouse model in this list.

---

## Suggested order of work
1. A: GSE190665 (Mammal40 betas, 49 MB, same probe space as the MMC Tet3 result) → essential-promoter deviation from 4F-young per tissue, treated vs untreated old.
2. A: GSE156557 (per-sample RRBS cov ~20 MB) and GSE247179 (40k RDS, 74 MB; over limit, needs approval) as independent replication in RRBS and in fibroblasts.
3. A (human): GSE165179 via per-sample IDATs (~15 MB per sample pair), with mock-transfected same-donor controls as reference.
4. B: GSE271375 (Tet1 cortex, n=8/group; per-sample IDATs) as the closest analogue of Tet3 brain; GSE324160 (Dnmt3a gain, 4 ages × 3 tissues) as the opposite-direction test; GSE299673 RRBS (5 vs 5) as a small gain-of-function check.
5. Human support: GSE128801 (TBRS, 27 MB) and GSE95040 (ICF1, 60 MB).

Files over 50 MB that the top candidates need: GSE165179 processed (530 MB), GSE271375 betas (145 MB), GSE324160 betas (147 MB), GSE247179 RDS (74-150 MB), GSE95040 matrix (60 MB). Per-sample IDATs or coverage files avoid the limit for all except GSE262342/GSE244432 bedGraphs (~200 MB each).

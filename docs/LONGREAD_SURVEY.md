# Long-read and read-level methylation data with ages (roadmap RM22)

Survey by a delegated agent for the lead, 2026-11-12. [V] checked against a primary file, [P] from a paper via search, [U] unverified. Not yet a result.

# Public data for read-level (single-molecule) methylation heterogeneity with age

Survey date: 2026-11-12. Report only; nothing downloaded beyond metadata (<1 MB).

**Question.** Is age-related DNA methylation change stochastic at the level of single molecules? Do promoters of essential or housekeeping genes stay protected (each molecule uniformly methylated or uniformly unmethylated) while matched regions turn into mixtures of epialleles with age?

**Labels.** Each fact carries one of these labels:

- **[V]** I checked it myself in this session, against a GEO SOFT file, an S3 listing, PMC full text through E-utilities, or the paper's own page.
- **[P]** It comes from a paper's text as reported by a search or fetch summary. I did not open the primary file.
- **[U]** It is unverified: my recollection or an inference. Check it before relying on it.

## Ranking (suitability: open per-read data, many donors spanning ages, one tissue)

| # | Dataset | Platform | Tissue | Donors / ages | Per-read? | Access |
|---|---|---|---|---|---|---|
| 1 | Loyfer atlas GSE186458 (+WBC/cfDNA add-on) | WGBS PE150 | sorted cell types; WBC | 23 WBC donors 21–75; atlas 3–83 | yes, PAT files | open |
| 2 | Chan 2025 entropy GSE288139 | targeted BS PE150 | buccal | 100, 7–84 | yes (raw reads) | open |
| 3 | CARD NABEC + HBCC brain ONT | ONT R9/R10 ~40x | frontal cortex | 351, 15–96 | BAMs in AnVIL (MM/ML [U]) | controlled |
| 4 | Mouse RRBS clocks (GSE80672, GSE120137, GSE93957) | RRBS | blood; multi-tissue | 62–549 samples, full lifespan | yes (raw reads) | open |
| 5 | Găitănaru 2025 ONT blood | ONT adaptive sampling | blood | 7, 5–91 | modBAM (gofile link) | open (fragile host) |
| 6 | All of Us lrWGS | HiFi + ONT | blood | ~2,800 lrWGS (CDRv8), adult ages | BAM [U] | controlled (Workbench) |
| 7 | 1000 Genomes ONT (1kGP-LRSC) | ONT ~37x | LCL | 500, **no ages** | yes, phased BAM with 5mC tags | open |
| 8 | deCODE Iceland ONT | ONT ~20x | blood | 7,179, ages yes | on site only | closed |

---

## 1. Loyfer et al. 2023, human methylation atlas (Nature 613:355): GEO GSE186458. Best open option

- **[V] GEO SOFT, parsed in this session.** There are 253 samples, and 252 have a numeric age (range 3–83, median 54). Every sample ships `*.pat.gz` plus `.csi` for both hg19 and hg38, together with `.beta` and `.bigwig` files. **PAT is a read-level format.** Each read-pair keeps its CpG index and its string of methylated, unmethylated and unknown calls, so it carries per-molecule epialleles.
- **[V] Ages by cell type.** Most groups have 2–8 donors. The exceptions are:
  - Epithelium: 63 samples across many organs.
  - **White blood cells: 23 samples from 23 donors aged 21–75**, mostly female (2 male). The GSM6810026–GSM6810048 names follow the pattern `CNVS-NORM-WBC-…-WGBS`.
  - cfDNA: 23 samples from the same 23 donors.
  - Endothelium: 19 samples, ages 28–81.
  - Neurons: 10 samples, ages 19–65.
  - Macrophages: 8 samples, ages 43–79.
  - Cardiomyocytes: 6 samples, ages 54–83.
- **[P] Where the WBC samples come from.** Ochana et al. 2025 call these "23 donors subjected to WGBS of genomic DNA from blood ... from donors in the USA". They are a later addition to the series, not one of the 205 sorted-cell atlas samples.
- **[V] Size.** One hg38 WBC `pat.gz` is 555 MB, so the 23 WBC samples total about 13 GB. The whole atlas is roughly 150 GB in hg38 PAT [U, extrapolated].
- **Fit.**
  - For: open, read-level, one tissue (WBC) spanning 21–75, and sorted cell types that remove the composition confound.
  - Against:
    - WBC is a mixture whose myeloid shift with age confounds heterogeneity. The sorted types have too few donors per type for an age analysis.
    - Read pairs of 150 bp give ~3–10 CpGs per fragment in CpG islands and promoters, but few in CpG-poor regions.
- **First analysis.**
  1. Pull `*.hg38.pat.gz` + `.csi` for the 23 WBC samples (~13 GB) with `geofetch` or plain `curl` from `ftp.ncbi.nlm.nih.gov/geo/samples/GSM6810nnn/…`.
  2. Use **wgbstools** (Loyfer/Kaplan; `wgbstools homog` counts per region the reads that are mostly unmethylated, mixed or mostly methylated, i.e. the protected-versus-mixed question directly).
  3. Compute PDR, Shannon epiallele entropy over 4-CpG windows, and MHL from the PAT reads with a small script.
  4. Build the region sets: housekeeping promoters (Eisenberg & Levanon 2013) and essential genes (DepMap common essentials). Match the controls on CpG density and mean methylation.
  5. Regress each metric on age; the cfDNA samples from the same donors serve as replication.

## 2. Chan, Rubbi, Pellegrini 2025, "DNA methylation entropy is a biomarker for aging" (Aging 17:685): GEO GSE288139

- **[V] PMC text and GEO summary.** Targeted bisulfite sequencing (3,015 probes) of **buccal swabs from 100 donors aged 7–84**, run on NovaSeq PE150 and aligned with BSBolt. GEO lists 100 samples, and the raw reads are on SRA [P]. The paper computes entropy over 4-CpG windows (16 possible epialleles).
- **Fit.**
  - For: one tissue, wide age span, 100 donors, open read-level data, and a published entropy baseline to compare against.
  - Against: the probes were chosen for age association [P], so whether they cover housekeeping promoters is **[U]** (Supplementary File 1 lists the probes). Buccal samples are an epithelial/leukocyte mixture.
- **First analysis.**
  1. Download the SRA FASTQ for 100 samples (size **[U]**; a small capture is probably ~0.5–2 GB per sample).
  2. Align with Bismark or BSBolt, then convert to PAT with `wgbstools bam2pat`.
  3. Intersect the probe list with housekeeping promoters, and compute PDR, entropy and MHL against age.

## 3. CARD long-read brain: Billingsley et al. (bioRxiv 2024.12.16.628723) and Grant/Makarious/…/Billingsley 2025 long-read clocks (bioRxiv 2025.09.30.679553)

- **[P] Cohort (bioRxiv full text).**
  - Frozen post-mortem **frontal cortex**, neurologically normal donors.
  - NABEC: 205 donors of European ancestry, mean age 52.7, **range 15–96**. Sequenced on R9.4.1 with LSK110 at ~46x.
  - HBCC: 154 donors of African or African-admixed ancestry, mean age 44.2, **range 18–85**. Sequenced on R10.4.1 with LSK114 at ~39x.
  - 351 donors were used in the analyses. N50 is ~27 kb.
  - Methylation was called with `modkit pileup` (NAPU pipeline), with phased output.
  - **Results (bulk).** About 28% of 1-kb windows associate with age in NABEC; PCDHG gains methylation and FKBP5 loses it.
  - **Results (clocks).** The 2025 clock paper uses 28 M CpGs and an elastic net.
- **[V] Access (PMC text of the 2025 preprint).** "Human brain sequencing datasets are under controlled access and require a dbGaP application (phs001300.v5) (phs000979.v4) … available through the restricted AnVIL workspace." HBCC also involves NDA collection 3151 [P]. I did not verify whether the AnVIL BAMs keep their MM/ML tags or only bedMethyl is shared **[U]**. The NAPU pipeline normally keeps the tags, which makes it likely.
- **Fit.** This is the best design for the question: a single tissue, hundreds of donors from 15 to 96, and ~27 kb molecules from post-mitotic cells. It is controlled-access, and the chemistry differs between the two cohorts (R9 vs R10), so analyse each cohort separately.
- **Notes on the cohort.**
  - Neuronal/glial composition changes with age; deconvolve or stratify by read-level cell-type assignment.
  - Post-mortem interval is a covariate.
- **First analysis.**
  1. File a dbGaP Data Access Request for phs001300 and phs000979, plus NDA access for HBCC.
  2. Compute in AnVIL; the BAMs are ~100–200 GB each, so stream regions instead of downloading them.
  3. Run `modkit extract` (per-read calls) on the housekeeping-promoter and matched-control BED files. `modkit entropy` computes windowed methylation entropy directly **[U: check that the subcommand exists in the installed version]**.
  4. Then compute per-read mean methylation, the fraction of discordant reads, and haplotype-aware metrics (the data are phased).

## 4. Mouse RRBS, open and read-level (short reads, CpG-island-biased, which suits promoters)

- **[V] GEO summaries.** Sample counts in the bullets below come from GEO; the authors are my attributions [U].
  - **GSE80672.** Petkovich et al. 2017, whole-blood RRBS, 255 samples, various ages, strains and interventions. Bertucci-Richter et al. 2024 analysed "255 mice" RRBS for entropy [P], probably this set [U].
  - **GSE120137.** "A multi-tissue full lifespan epigenetic clock for mice", 549 samples [V]. I recall the paper as Thompson et al. 2018 Aging [U].
  - **GSE93957.** Stubbs et al. 2017, multi-tissue (62 samples).
  - **GSE89275.** Lifespan-intervention series (48 samples, RRBS + RNA). I recall it as Cole et al. 2017 [U].
- RRBS reads are often single-end 50–100 bp **[U per series]**, which gives 2–6 CpGs per read in CpG islands: enough for PDR and 4-CpG entropy.
- Mouse has **no public ONT or HiFi aging methylation dataset** that I could find. A BioProject search for nanopore methylation aging in mouse returned only Faulk's genome-skimming project, PRJNA927034, which has no age design.

## 5. Small or pilot long-read blood sets with ages

- **Găitănaru et al. 2025** (Epigenomes 9:25), "Exploring epigenetic ageing using direct methylome sequencing".
  - **[V] PMC text.** ONT adaptive sampling of **7 leftover blood samples, ages 5–91**, called with dorado HAC 5mC/5hmC and modkit.
  - **[P] Data.** The "download link for our modBAM input files" is in GitHub `Independent-Research-Lab/nanopore-epigenetics`, pointing to a gofile.me share (fragile host, size unknown).
  - BioProject PRJNA1227682, "Methylation profiling of human gDNA using Oxford Nanopore Adaptive Sampling", is likely the same project [U].
  - Too few donors for statistics; useful only as a pilot for the pipeline.
- **Yuen et al. bioRxiv 2023.12.10.571037** [P].
  - ONT R9.4.1 adaptive sampling over a ~113 Mb panel of aging genes, with blood from 10 donors aged 25–76.
  - Accession not found in the text **[U]**.
- **Eames, Moqri, Poganik, Gladyshev 2025, LongReadAge** (Aging Cell 24:e70084).
  - **[V] Leukocytes.** The leukocyte data are 16 PacBio Sequel II buffy-coat profiles from 11 pregnant women (CUHK, Yu et al. 2021); these are **not public**, only "on request".
  - **[V] cfDNA.** A new PacBio plasma cfDNA set is in GEO **GSE296282** (donor ages there [U]).
  - cfDNA is not one tissue, but it is open.

## 6. Controlled or closed large long-read cohorts

- **All of Us** [P].
  - The CDRv8 (Feb 2025) Controlled Tier includes **2,800 lrWGS** samples (PacBio HiFi Revio/Sequel IIe + ONT) from blood, with ages from EHR/survey data.
  - Whether the released BAMs keep their MM/ML tags is **[U]**; HiFi Revio reads normally carry them.
  - Analysis must stay inside the Researcher Workbench. The cohort deliberately oversampled under-represented groups.
- **deCODE / Sigurpalsdottir et al. 2024** (Genome Biology 25:69).
  - **[V] PMC text.** ONT whole blood from **7,179 individuals at ~20x**, with 5mC called by Nanopolish.
  - "Access to the raw Icelandic sequence data is available on request ... at the premises of deCODE genetics. The data are not publicly available because of Icelandic state law."
  - Only CpG summary statistics are public (decode.com/summarydata, Zenodo).
  - **[P]** Sigurpalsdottir's PhD thesis lists a submitted manuscript on age-associated disruption of imprinting.
- **SMaHT** (Somatic Mosaicism across Human Tissues) [P].
  - Targets 150 non-diseased donors × 19 tissue sites, with PacBio and ONT among its platforms.
  - Somatic variant calls are open (s3 `smaht-open-data-public`). **Aligned reads are controlled** (dbGaP phs004193 benchmarking, phs004194 production).
  - Donor ages and the methylation tags were not checked **[U]**. It is a multi-tissue, same-donor design worth watching.
- **CARD** is the brain set described in section 3.

## 7. LCL fallbacks (open, per-read, but no usable ages)

- **1000 Genomes ONT, 1kGP-LRSC** (Gustafson et al. 2024 Genome Res 34:2061).
  - **[V] `s3://1000g-ont` listing.** `PROCESSED_DATA/READ_ME.txt` (2025-11-04) covers **500 samples**.
  - **[V] Metadata.** `1kGP_LRSC_500_ONT_Metadata.tsv` has the columns NHGRI_ID, Sex, SubPopulation, SuperPopulation, library and pore. **There is no age column.** The 1000 Genomes Project did not collect donor ages [U, standard knowledge].
  - **[V] Files.** Phased minimap2 BAMs come from the basecallers `guppy-sup-5mC` and `dorado050_sup_5mCG_v33`, so they carry 5mC MM/ML tags. Each BAM is **134–202 GB** with a `.bai`. Modkit bedMethyl files (combined, hp1, hp2, split 5hmC) are ~0.5 GB each.
  - **[V] Access.** Open, over plain HTTPS: `samtools view https://1000g-ont.s3.amazonaws.com/...bam chr:start-end` streams regions without a full download.
  - Additional R10, UltraLong and pod5 folders exist [V].
  - **Use.** EBV transformation and culture make methylation abnormal (LCLs are clonal and hypomethylated). These BAMs are useful only for **methods development** and for a "protected promoter" baseline (housekeeping versus control heterogeneity in a single sample), not for age.
- **HPRC** (release 2, ~230 individuals) [U, recollection].
  - Open `s3://human-pangenomics`, with HiFi and ONT that carry MM/ML; LCLs, no ages.
  - Same use as 1000 Genomes ONT.
- **Kerr et al. 2023** (PLoS Genet) is the method precedent for single-molecule heterogeneity in LCL nanopore reads (highest in PMDs) [P].

## 8. Other read-level bisulfite data with ages

- **Heyn et al. 2012 PNAS** [V GEO]. GSE31263, 3 samples. WGBS of CD4+ T cells from a newborn and a centenarian, plus one middle-aged donor [U]. Open, but n=2–3.
- **BLUEPRINT WGBS** [U, recollection].
  - Purified blood cell types (monocytes, neutrophils, T and B cells) from about 50–100 donors, with ages in 5-year bins.
  - Held in EGA under controlled access (DAC application); read-level BAMs.
  - A good second blood resource if the controlled process is acceptable.
- **Ochana et al. 2025 Cell Rep 44:115958 (Dor/Kaplan/Shemer)** [V PMC preprint text].
  - Ultra-deep multiplexed bisulfite-PCR (~12,800 fragments per amplicon) of blood from **296 healthy donors aged 17–78**, plus 52 Jerusalem Perinatal Study donors sampled at ~32 and ~42.
  - The data deposit was **not found** (no GEO link from PubMed) **[U]**.
  - Read-level by design, but limited to amplicons at age-associated loci.

## Published findings on read-level heterogeneity, entropy and age (one line each)

- **Xie et al. 2011 NAR.** Defined methylation entropy over 4-CpG read patterns, the basis of later epiallele metrics [U recollection].
- **Landau et al. 2014 Cancer Cell.** PDR (proportion of discordant reads) on RRBS. Locally disordered methylation in CLL tracks clonal evolution; disorder is a stochastic process (cancer, not ageing).
- **Guo et al. 2017 Nat Genet (Shao, Zhang lab).** Methylation haplotype blocks and MHL for tissue-of-origin deconvolution, not ageing [U on the year and authors].
- **Shipony et al. 2014 Nature.** Clonal single-cell lineages show both dynamic (stochastic, heterogeneous) and static methylation memory [U].
- **Heyn et al. 2012 PNAS.** The centenarian CD4 methylome has lower, more heterogeneous methylation (hypomethylation at CpG-poor regions) than the newborn's.
- **Jenkinson et al. 2017 Nat Genet.** Entropy and potential-energy landscapes from WGBS. Entropy differs by tissue and rises in cancer; an age claim is [U].
- **Hernando-Herraez et al. 2019 Nat Commun.** scM&T-seq of mouse muscle stem cells: ageing raises promoter methylation heterogeneity (stochastic drift), linked to cell-to-cell transcriptional variability.
- **Bertucci-Richter et al. 2024 Aging.** RRBS of 255 mice: entropy and disorder change with age, but are partly a function of the mean and of cell heterogeneity, and are decoupled from the canonical clock.
- **Chan et al. 2025 Aging.** Buccal targeted BS of 100 donors aged 7–84: entropy changes reproducibly with age, at some loci up (more disorder) and at some down, independent of mean methylation; entropy clocks match beta clocks, and combining them does best.
- **Ochana et al. 2025 Cell Rep.** Age changes occur at clusters of CpGs, either stochastically (each CpG independently) or in coordinated blocks (whole molecules flip). Fifty molecules suffice to infer age, suggesting each cell carries age information.
- **Tarkhov et al. 2024 Nat Aging.** Epigenetic clocks can be built from a purely stochastic component; part of clock signal is stochastic drift [U on the details].
- **Billingsley 2024 / Grant 2025 (CARD).** Bulk long-read age associations in brain (PCDHG, FKBP5) and ancestry-aware ONT clocks. They report no read-level heterogeneity analysis yet [P]: an open gap.
- **Kerr et al. 2023 PLoS Genet.** Single long reads show large-scale heterogeneity hidden in bulk, highest in heterochromatin and PMDs, with nucleosome-scale periodicity (LCLs, no age).
- **Gap.** I found **no published study of single-molecule epiallele heterogeneity versus age on long reads** (ONT or HiFi) in any tissue. The question appears open.

## Practical notes

- Everything here was checked over the proxy. S3 listings work over HTTPS (`https://s3.amazonaws.com/1000g-ont?list-type=2&prefix=…`), and NCBI E-utilities work.
- Housekeeping and essential gene lists are not dataset-specific: Eisenberg & Levanon 2013 HK genes, and DepMap CRISPR common essentials. Matched controls: promoters with the same CpG density (CGI versus non-CGI), the same mean methylation in young donors, and the same coverage.
- Confounders for any analysis:
  - Cell-type composition (myeloid shift in blood, neuronal loss in brain). Mitigate with sorted cells (Loyfer) or per-read cell-type assignment (LongReadAge approach).
  - Coverage, which biases entropy and PDR. Downsample to a fixed number of reads per region.
  - Basecaller and chemistry (R9 vs R10), which bias per-read calls.
  - Mean methylation, which bounds entropy. Use metrics conditioned on the mean, or the Bertucci-Richter decomposition.

<!-- Search by a delegated agent, 2026-11-23. Checked against GEO series and sample records (brief text view), GEO file sizes (HTTP HEAD on the GEO FTP mirror), the figshare API, and the full texts of the papers through Europe PMC where open. Nothing was downloaded. Not committed. -->

# Lineage-resolved promoter methylation: maintenance versus selection

**Question.** Promoters of essential and housekeeping genes drift less in DNA methylation with age or with cell
division. There are two candidate explanations:

- **(1) Maintenance.** Each lineage keeps its essential promoters clean, through replication-coupled maintenance.
- **(2) Selection.** Lineages whose essential promoters drift are lost or out-competed.

To tell them apart, we need methylation measured per lineage (clone) rather than as a population average. We also need
a way to see selection: clone abundance, or clones grown without competition.

**What a discriminating dataset looks like.** It has several donors or animals. Promoter-level methylation is measured
per clone, by array, WGBS, RRBS or single-cell methylome. Lineage is given by a barcode, a somatic mutation, a
mitochondrial variant, or a physically clonal unit such as a single-cell-derived clone or a crypt. Expression, clone
size and time points (or matched doublings) are a bonus.

**Two predictions separate the explanations:**

- Under (1), essential promoters stay clean *within every clone*, small clones included, and clones grown without
  competition stay clean too.
- Under (2), clones grown without competition carry essential-promoter drift. Clone size or survival is then
  anti-correlated with that drift, and the clean population average comes from the survivors.

**How it was checked.**

- Literature: PubMed, Europe PMC (including bioRxiv), web search, and the reference list of the 2026 review
  "DNA methylation meets lineage tracing" (Quant Biol, [10.1002/qub2.70017](https://doi.org/10.1002/qub2.70017)).
- GEO: ESearch/ESummary on `db=gds`, then sample-level reads of each candidate.
- File sizes: GEO by HTTP HEAD; figshare by its API.
- EGA, NGDC and dbGaP: only the accession and the access mode were noted. File sizes there were not checked.

**Bottom line.** No public dataset has all of the following together: adult human donors, replication across donors,
per-lineage genome-wide promoter methylation, and clone abundance over time. The closest open datasets each have part
of it:

- **Single-cell-derived fibroblast clones on 450k with IDATs.** Clones are grown without competition, at matched
  doublings, from two strains. There is no clone size.
- **Barcoded single-cell genome-wide methylomes in mouse (DARLIN/Camellia-seq).** These give lineage, promoter
  methylation, RNA and clone size, but come from young animals and few of them.
- **Human clonal haematopoiesis single-cell RRBS with genotyping.** Four or five donors; an expanded clone against the
  wild-type background. The DNMT3A mutation itself confounds the methylation.
- **EPI-Clone.** It has the best clone-abundance and age design, with old against young mice and 12 human donors. Its
  ~450-CpG panel was chosen from *variable* CpGs and is not built to cover unmethylated housekeeping promoters.

## Ranked top 5 (open access first)

| Rank | Dataset | Organism | Why it ranks here | Verdict |
|---|---|---|---|---|
| 1 | Flunkert 2018 fibroblast clones, GSE112812 + GSE112873 | Human (2 fetal strains) | Per-clone 450k with IDATs; clones grown without competition at ~20 PD; parental culture arrayed | Meets with caveat |
| 2 | DARLIN Camellia-seq, GSE222328 (super GSE222486) | Mouse | Barcoded single-cell genome-wide methylome + RNA in vivo; clone size from barcodes | Meets with caveat |
| 3 | Nam 2022 DNMT3A CH, GSE158067 | Human (4 + 1 donors) | scRRBS + genotype + Smart-seq2 per cell; expanded clone vs wild-type | Meets with caveat |
| 4 | EPI-Clone, figshare + GSE282971 | Mouse + human (12 donors) | Clone sizes, old vs young mice, LARRY ground truth, a mito-lineage donor | Meets with caveat (fails as a promoter readout) |
| 5 | Shipony 2014, GSE53610 | Human (WI-38, one donor) | RRBS read-level patterns per clone vs parent population | Meets with caveat (one donor) |

## Open-access candidates

| Accession | Paper | Organism, cell type | Design | Assay | Promoter-level? | Lineage, clone size | Files (open) | Verdict |
|---|---|---|---|---|---|---|---|---|
| [GSE112812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112812), [GSE112873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112873) | Flunkert 2018, PMID 29964050, [10.1016/j.yexcr.2018.06.034](https://doi.org/10.1016/j.yexcr.2018.06.034) | Human fetal fibroblasts, 2 strains | Strain 1: 11 non-irradiated clones + 35 irradiated (2 Gy) + 2 starting-culture arrays. Strain 2: 20 + 26 + 2. Clones from single cells, ~20 PD. Karyotype stable or unstable recorded | 450k (GPL13534) | Yes, 450k promoter probes; IDATs given | Lineage = clone (physically clonal). No clone size or growth rate in GEO; survivorship of plated clones is not recorded | RAW.tar 577 + 579 MB; normalised betas 29–93 MB per file | **Meets with caveat.** Best open per-clone promoter design. Test: is essential-promoter drift per clone (vs the starting culture) as low as the population's? There is no competition between clones inside one clone. Irradiated clones add a perturbation arm. Caveats: fetal donors; 2 strains; no clone sizes; whether clones failed to grow is unknown |
| [GSE222328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222328) (super [GSE222486](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222486); scRNA [GSE221501](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221501)) | Li 2023 DARLIN, PMID 37852258, [10.1016/j.cell.2023.09.019](https://doi.org/10.1016/j.cell.2023.09.019) | Mouse HSCs: adult bone marrow, fetal liver, AGM | 1,080 single HSCs. Adult BM from 2 mice (LL731, LL901; 264 cells). Fetal liver from 1 litter, 2 embryos (480). AGM (144). No ageing time points | Camellia-seq: genome-wide single-cell DNA methylation + accessibility + RNA + DARLIN barcode | Yes after pooling cells per clone (sparse per cell) | DARLIN Cas9 barcode per cell. Clone size = cells sharing a barcode in the sample | RAW.tar 38.3 GB | **Meets with caveat.** The only open in vivo set with lineage barcodes and genome-wide methylomes plus RNA. Test: correlate per-clone essential-promoter drift with clone size. Caveats: young or developmental only; 2 adult animals; few cells per clone; drift in young HSC is small |
| [GSE158067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158067) | Nam 2022, PMID 36138229, [10.1038/s41588-022-01179-9](https://doi.org/10.1038/s41588-022-01179-9) | Human CD34+ HSPC, DNMT3A R882 clonal haematopoiesis | 4 G-CSF-mobilised donors + 1 validation marrow. Cross-sectional | Single-cell RRBS + Smart-seq2 + targeted genotyping; also 10x + genotyping | Yes; RRBS covers CpG-island promoters | Lineage = DNMT3A-mutant clone vs wild-type cells per cell. Clone size from the mutant fraction | Per-cell CpG matrix `all_CpG_meth.csv.gz` 282 MB; cell metadata; RAW.tar 10.5 GB; raw reads in EGA | **Meets with caveat.** Donor-replicated, human, with RNA. An expanded (selected) clone is compared with the polyclonal background. Caveat: the DNMT3A loss changes methylation itself, so drift and genotype are confounded. Lineage is binary, not many clones |
| figshare [24204750](https://doi.org/10.6084/m9.figshare.24204750), [25472434](https://doi.org/10.6084/m9.figshare.25472434.v1), [27917454](https://doi.org/10.6084/m9.figshare.27917454.v1), [25526899](https://doi.org/10.6084/m9.figshare.25526899.v2), [28082048](https://doi.org/10.6084/m9.figshare.28082048.v1), [28082066](https://doi.org/10.6084/m9.figshare.28082066.v1); [GSE282971](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282971); SRA PRJNA1191391 | Scherer 2025 EPI-Clone, PMID 40399669, [10.1038/s41586-025-09041-8](https://doi.org/10.1038/s41586-025-09041-8) | Mouse HSPC + mature cells; human bone marrow and CD34+ | Mouse: LARRY-barcoded transplants (M.1, M.2); native young vs old (M.7, 2 replicates); transplantation (M.8). Human: total bone marrow A.1–A.7 and CD34+ B.1–B.5 across ages (CD34+ cohort mostly 50–60), with and without CH drivers. X.2: one donor with mitochondrial variants | scTAM-seq (HhaI-digest targeted): mouse 453 CpGs from 663 amplicons; human panel of similar design; 20 surface proteins | **Mostly no.** CpGs were chosen as differentially, intermediately or within-sample heterogeneously methylated. "Static" clone-informative CpGs sit in heterochromatin. Constitutively unmethylated housekeeping CGI promoters are not the design target (inferred from the panel rules; check Supp. Table 6 for any promoter CpGs) | Yes: LARRY ground truth, EPI-Clone clone calls, relative clone size, old vs young | Seurat RDS 9–199 MB per dataset (CC BY 4.0); GSE282971 RAW 124 MB; human raw in EGA | **Meets with caveat for the selection side; fails as a promoter readout.** The best donor-replicated clone-size and age design. Use it for the expansion analysis, or test any promoter CpGs that are in the panel |
| [GSE53610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53610) (super GSE54960) | Shipony 2014, PMID 25043040, [10.1038/nature13458](https://doi.org/10.1038/nature13458) | Human WI-38 fibroblasts (one fetal donor); also hESC H1/H9 and cancer lines | WI-38 parent population + 6 single-cell clones (2–3 replicates); cell-cycle-sorted G/S fractions | RRBS, read-level methylation patterns (`.pats`) | Yes, at CpG-island promoters covered by RRBS; per-read patterns allow epipolymorphism per lineage | Lineage = clone. No clone size | Many small `.pats.txt.gz` files (sizes not read) | **Meets with caveat.** Per-clone read-level patterns measure maintenance fidelity per lineage directly. One donor; no selection readout |
| [GSE262218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262218) (+ LARRY [GSE262220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262220), RNA [GSE262219](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262219)) | Chen 2025 MethylTree, PMID 39820752, [10.1038/s41592-024-02567-1](https://doi.org/10.1038/s41592-024-02567-1) | Mouse Lin−cKit+Sca1− progenitors | One adult mouse; LARRY-labelled, 6 days of culture; 576 cells | Modified Camellia-seq: single-cell genome-wide methylation + RNA + LARRY | Yes after pooling per clone | LARRY clone and clone size per cell | RAW.tar 5.7 GB (methylation), 3 MB (RNA) | **Meets with caveat.** One animal and about 6 days, so little drift is expected. Human series [GSE259241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259241) (HEK293) and [GSE262216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262216) (H9 hESC) are cell lines and fail |
| [GSE109085](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109085) | Gaiti 2019, PMID 31092926, [10.1038/s41586-019-1198-z](https://doi.org/10.1038/s41586-019-1198-z) | Human B cells from 6 healthy donors (B01–B06, CD19+/CD27± subsets) + CLL | Cross-sectional; CLL with 2 time points | Multiplexed scRRBS (+ Smart-seq2 in the dual protocol) | Yes per cell, at CGIs | Lineage only inferred from the epimutations themselves | Processed per-cell averages only; raw on SRA (size not checked) | **Fails for this test.** Lineage is circular, because it is called from the methylation being tested. Usable as a per-cell epimutation reference in healthy B cells |
| [GSE87527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87527) | Yu 2016, PMID 27863245, [10.1016/j.cell.2016.10.045](https://doi.org/10.1016/j.cell.2016.10.045) | Mouse LKS SLAM HSCs, multicolour clonal tagging | 2 cohorts; colour-sorted clone groups (R/Y; G/P) | WGBS (6 libraries), RNA-seq, ATAC | Yes | Colour = oligoclonal group; functional output per group | RAW.tar (size not checked) | **Fails (weak).** A colour is a group of clones, n ≈ 2 per comparison, young mice |
| [GSE41048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41048) (super GSE41050) | Landan 2012, PMID 23064413 | Human immortalised fibroblasts | In vitro evolution over 300 generations | Bisulfite sequencing + MeDIP | Partly | Lineage = culture line; no clone sizes | Not checked | **Fails (weak).** Single immortalised donor line |
| [GSE129712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129712) (RRBS crypts), [GSE271646](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271646) (RNA single crypts) | Krepelova 2025, PMID 41299091, [10.1038/s43587-025-01021-x](https://doi.org/10.1038/s43587-025-01021-x) | Mouse intestinal crypts, young vs aged | 50 + 50 single crypts, 3 mice per group; laser-captured neighbouring crypts in 2 aged mice | Bisulfite pyrosequencing of single crypts: Dkk2, Sfrp1, Dkk1, Sfrp2 promoters + Actb and Dppa3 controls. RRBS on pooled crypts | Only 6 promoters per crypt (in source data). Actb is the one housekeeping control | Crypt = clone; neighbours share drift through fission | Source data with the paper | **Fails for a genome-wide test.** Shows clonal spread of drift through crypt fission. The per-crypt Actb control is a small direct look at a housekeeping promoter |
| Franzen 2021 (4C/CTCF only in [GSE144196](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144196)) | PMID 34011964, [10.1038/s42003-021-02116-y](https://doi.org/10.1038/s42003-021-02116-y) | Human umbilical-cord MSC, 2 donors | Lentiviral RGB-barcoded cultures to senescence; oligoclonal by passage 10 | Barcoded bisulfite amplicon sequencing at 10 culture-associated CpGs | No | Barcode clone frequencies (from Selich 2016) | Not in GEO | **Fails.** Only 10 CpGs, none at housekeeping promoters. Clone frequencies exist but methylation is not per clone |

## Controlled-access candidates

| Accession | Paper | Organism, cell type | Design | Assay | Promoter-level? | Lineage, clone size | Access | Verdict |
|---|---|---|---|---|---|---|---|---|
| EGAS00001005514 | Gabbutt 2022 fluctuating methylation clocks, PMID 34980912, [10.1038/s41587-021-01109-w](https://doi.org/10.1038/s41587-021-01109-w) | Human single colon crypts (31 crypts, 8 individuals), small-intestinal crypts, endometrial glands | Cross-sectional; several crypts per person | EPIC per crypt | Yes, EPIC promoters | Crypt = monoclonal unit; no clone size or fission data | EGA, controlled. Open Supp. Tables 2–3 hold betas for fCpG loci only | **Meets with caveat.** Per-clone human promoter methylation across donors. Test: between-crypt variance at essential vs other promoters within a person. There is no selection readout; crypt survival is unobserved |
| EGAS00001008056 | EPI-Clone human cohort (see above) | Human bone marrow, CD34+ | 12 donors | scTAM-seq raw | As above (panel) | As above | EGA, restricted to haematology and methods research | Processed objects are already open on figshare; request only if raw reads are needed |
| EGA (raw reads of GSE158067) | Nam 2022 | As above | As above | scRRBS raw | Yes | As above | EGA, controlled | Not needed: the per-cell CpG matrix is open in GEO |
| NGDC GSA-human HRA008624 | Chen 2025 MethylTree | Human data of the MethylTree study | Not verified | scBS | Not verified | Not verified | Controlled through a DAC | Unverified lead; open GEO human parts are cell lines |

## Other leads read and excluded

- **Krepelova 2025 single crypts.** Listed in the open table. It is the only in vivo clone-unit study that reports drift
  spreading by fission.
- **Graham 2011 and the Shibata crypt studies** (Yatabe 2001, Kim 2002/2004, Nicolas 2007). They used bisulfite clones
  of a few neutral loci and deposited no genome-wide data. Excluded.
- **Mitchell 2022 and Lee-Six 2018/2019** (HSC colony whole-genome phylogenies, clone sizes, many donors). There is no
  methylation on the colonies. They would be the ideal frame if colony methylomes were ever produced.
- **ReDeeM and other mtDNA lineage studies in human blood.** They have accessibility or RNA, not methylation. The one
  methylation and mitochondrial-lineage combination found is EPI-Clone X.2, a single donor.
- **Kraiczy 2022, 173 intestinal organoid lines.** Organoid lines are not single-cell clones, and passage is
  confounded. Excluded.
- **Meir 2020 (GSE144357).** HCT116 cancer clones. **GSE44222**: 5 adipose stem-cell clones on the 27k array, one
  donor. Both excluded.
- **"Cell-to-cell variability and gain of methylation at polycomb CGIs"** (Nat Commun 2026,
  [10.1038/s41467-026-74118-5](https://doi.org/10.1038/s41467-026-74118-5)). It reuses single-cell methylomes without
  lineage. Not a dataset for this test, but relevant background: faster-proliferating cells gain polycomb-CGI
  methylation faster.
- **Not found:** a Blood 2025 abstract on high-throughput single-cell methylation of aged marrow with TET2 CH
  (S0006497125035529; 18,268 methylomes). It may be a conference abstract from the EPI-Clone group. Neither its
  accession nor its assay could be confirmed.

## Suggested use

1. **Flunkert clones.** Compute per-clone promoter means for essential genes (DepMap common-essential) against matched
   control promoters. Compare clone-to-clone variance, and each clone's shift from its starting culture, between the
   two sets.
   - Maintenance predicts small shifts in every clone.
   - Selection predicts that clones, grown without competition, show essential-promoter shifts the population does not.
2. **DARLIN Camellia-seq adult bone marrow.** Pool cells per barcode. Regress per-clone essential-promoter drift on
   clone size.
3. **Nam 2022.** Repeat the comparison per donor for mutant against wild-type cells. Read the result against the known
   DNMT3A hypomethylation signature (PRC2 targets) so the genotype effect is not mistaken for selection.
4. **EPI-Clone (old against young mice; 12 human donors).** Use it for clone-size distributions. First check Supp.
   Table 6 for any panel CpGs in housekeeping promoters.

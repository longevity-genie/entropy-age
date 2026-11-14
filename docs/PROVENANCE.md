# Data provenance (roadmap RM34 item 5)

One row per methylation dataset: where it comes from, whether we use the authors' processed values or rebuilt them
from raw files, and which figures read it. The omics re-analysis practice asks that a result name the run it came
from and that datasets processed differently are never combined; here no test pools two datasets, every comparison
is within one dataset and platform (`docs/PREPROCESSING_PRACTICE.md`). Tool versions are pinned in `uv.lock`; the
IDAT pipeline is `scripts/reprocess_idats.py` (methylprep 1.7.1, SeSAMe-like: noob, non-linear dye bias, channel
switch, pOOBAH mask; Python 3.11, pandas < 2, numpy < 2). Runs, commits and job ids are in `docs/RUN_LOG.md`.

"Deposited" means the authors' processed betas (or M-values, counts) as GEO or the archive serves them, read as
is; "IDATs" means rebuilt by us from the raw array files. Accessions were read from the code (`src/entropy_age`).

| dataset | accession | platform | values | figures (`entropy-age figures ...`) |
|---|---|---|---|---|
| Hannum blood, 656 people | GSE40279 | 450k | deposited (no public IDATs) | tissues, per_person, placebo_sets, composition_gaps, equilibrium_drift, drift_hotspots, stochastic_component, clock_weights, chromatin, essentiality_gradient |
| Jaffe prefrontal cortex | GSE74193 | 450k | deposited; **IDATs** for R61 (`reprocessed_jaffe`) | tissues, reprocessed_jaffe, health_states, placebo_sets, equilibrium_drift, drift_hotspots |
| Lunnon blood and four brain regions | GSE59685 | 450k | deposited | tissues, health_states, organ_coupling |
| GSE41037 blood (paper replication) | GSE41037 | 27k | deposited | atlas, genesets, fig4, fig5, health_states |
| GTEx nine tissues | GSE213478 | EPIC | deposited (authors' noob + BMIQ) | gtex |
| Muscle fibres | GSE304045 | EPIC v2 | deposited (dasen) | fibres |
| SATSA twins, longitudinal | E-MTAB-7309 / S-BSST1206 | 450k | deposited (R52); **IDATs** for R66 (`satsa_idats`) | within_person, satsa_idats |
| Gill 2022 transient reprogramming | GSE165179 | EPIC | **IDATs** (R65) | gill_reprogramming |
| Cellular Lifespan Study fibroblasts | GSE179847 | EPIC | deposited | passaged_cells |
| Endicott primary cells | GSE197512 | EPIC | deposited | passaged_cells |
| EPIC-Italy incident cancer | GSE51032 | 450k | deposited | incident_cancer |
| SJLIFE survivors | GSE314261 | EPIC | deposited | sjlife_deaths |
| MESA sorted monocytes and T cells | GSE56046, GSE56581 | 450k | deposited | sorted_cells |
| Survival cohort | GSE68194 | 450k | deposited | survival_shield |
| Chromatin-writer syndromes | GSE74432, GSE97362 | 450k | deposited | syndromes |
| Loyfer 2023 sorted cells and blood (WGBS reads) | GSE186458 | WGBS | the authors' read files, our per-read summaries | read_level |
| Mammalian consortium (65 species; interventions) | GSE223748 | Mammal40 | deposited | species, species_floor, species_phylo, drift_onset, mmc_interventions, marrow_anomaly |
| Huang 11 mouse tissues | GSE279437 | Mammal40 | deposited | mouse, intervention_split, meer_tissues (reference) |
| Petkovich / Thompson mouse series | GSE80672 | RRBS | deposited | interventions, cr_chromatin |
| CR blood replication | GSE282661 | Mammal320 | deposited | cr_blood_replication |
| Mouse RRBS six tissues | GSE120132 | RRBS | deposited | rrbs_tissues |
| Meer 2018 four tissues | GSE121141 | RRBS | deposited | meer_tissues |
| BXD liver | GSE199979 | RRBS | deposited | bxd |
| Browder 2022 cyclic OSKM | GSE190665 | Mammal40 | deposited | reprogramming |
| Chondronasiou 2022 one OSKM cycle | GSE156557 | RRBS | deposited | reprogramming_rrbs |
| Tet1 cortex | GSE271375 | MM285 | deposited | tet1 |
| Dnmt3a gain of function | GSE324160 | mouse array | deposited | dnmt3a |
| CHIP methylation | GSE315366 | EPIC | deposited | chip_methylation |
| CR chromatin (ATAC, mouse blood cells) | GSE100738 | ATAC-seq | deposited peaks | cr_chromatin |
| Reference epigenomes (chromatin states) | GSE96614 / Roadmap | ChromHMM | deposited | chromatin, drift_hotspots |
| ATAC variance (Morandini 2023) | GSE193140 | ATAC-seq | deposited | atac_variance |
| ATAC hippocampus | GSE278576 | ATAC-seq | deposited | atac_hippocampus |
| Rodent RNA (methylation-expression spread) | GSE314817, GSE314970 | RNA-seq | deposited counts | rna_rodents |

## What this shows

* **Three datasets are rebuilt from IDATs** (Jaffe, Gill, SATSA); every other array dataset is the authors'
  processing, which differs from study to study (noob, BMIQ, dasen, SWAN, or unstated). That is safe for the
  within-dataset contrasts we run, and is why sizes are not compared across studies without the ±0.05 SD of R61.
* **No dataset mixes pipelines.** Jaffe and SATSA have both versions, each analysed on its own and compared
  side by side (R61; R52 against R66).
* **Not yet recorded here:** checksums of the downloaded files. `entropy-age fetch` keeps the archive's files under
  `data/input/`; adding their MD5 against the archive's is the next step of this table.

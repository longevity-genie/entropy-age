# What entropy-age code could move to other repositories

Review by a delegated agent for the lead, 2026-11-12, at antonkulaga's request. Report only; nothing moved.

# entropy-age: what is worth reusing (review 2026-11-12, report only)

Scope: `src/entropy_age` (about 16k lines with tests), 150 tests passing in 75 s. Every test uses synthetic data or tmp_path fixtures, so none needs network or real data. Two existing projects overlap with this code: **dna-seq/just-biomarkers** (polars clock scoring engine with 44 clocks, the same biolearn pin and a GEO downloader) and **rna-age** (`download.py` with resumable, checksummed fetches, and `ages.py` for GEO age parsing).

## Ranked candidates (most value for the least effort first)

1. **Matched, stratified group contrast** (b, new lib). Sources: `stochasticity.py` `baseline_strata`, `_matched_difference`, `_stratified_permutation`, `_permute_within_strata`, `compare_groups`, `cliffs_delta`, `paired_comparison`, and `protection.py` `contrast_weights`, `gene_counts`, `resampled_gaps`. It tests a gene set (or any feature group) on any per-feature metric while matching on a confounder such as level, island status, GC or length, with labels permuted within strata and genes resampled, not CpGs.
   - Dependencies: numpy, scipy, polars.
   - Tests: good. `test_gene_sets` (equals a per-stratum loop; permutation keeps strata) and `test_protection`.
   - Cleaning: none to speak of; the code is free of paths. Effort: small. Value: high, also for rna-age expression sets.
2. **GEO series-matrix and supplementary-table readers** (`data/geo.py`). (b), or better (d) to just-biomarkers, whose GEO support covers only 5 curated sets.
   - What it does: `read_series_matrix` turns `key: value` characteristics into columns; `read_supplementary_betas` handles suffix, preamble and ID-row quirks.
   - Dependencies: polars. Tests: `test_data`, `test_tissues`.
   - Cleaning: remove `ProjectPaths` from `fetch_`/`prepare_series`; use rna-age's `download.py` and `ages.py` instead of the local downloader. Effort: small. Value: high, since every GEO cohort needs it.
   - Optional (c): a `geo_series_samples(gse)` tool for holy-bio-mcp.
3. **Illumina manifest handling** (`data/annotation.py`): `_assay_block`, `read_manifest` (27k), `read_manifest_450k` (probe type, island relation, gene links exploded, promoter flag), `solo_wcgw`, `read_probe_mask` (Zhou `M_general`), `link_membership`.
   - Dependencies: polars. Tests: `test_tissues`, `test_intergenic`, `test_probe_robustness`.
   - Cleaning: split it from the gene-set and HGNC code, which pins a dated HGNC file. EPIC/EPICv2 are not handled yet. Effort: small to medium.
   - Destination: just-biomarkers (array QC before scoring). methylprep/SeSAMe do this too, but not in polars.
4. **Per-feature drift statistics** (`stochasticity.py`): `drift_metrics` (batch-adjusted slopes, partial Spearman, entropy slope), `dispersion_levels`, `expression_drift`, `fit_relaxation`, `regress_out`, `composition_components` (a reference-free proxy for the cell mix), `binary_entropy`.
   - Dependencies: numpy, scipy, polars. Tests: `test_clock_and_stochasticity` (recovers known rates and confounds), `test_rna_gtex`.
   - Cleaning: `fit_relaxation` imports `simulate.MaintenanceRates`; cut that link. Effort: small. Value: medium-high for any cross-sectional aging omics, rna-age included.
5. **Phylogenetic GLS on a taxonomy covariance**: `data/taxonomy.py` (`lookup`, `fetch_lineages`, `shared_depth`, `lineage_covariance`) and `experiments/species_phylo.py` (`gls`, `fit_lambda` for Pagel lambda, `jackknife`). (b), plus (c) a cross-species skill in longevity-skills.
   - Dependencies: numpy, scipy, stdlib urllib. Tests: `test_species_phylo` (lambda near 0 for independent noise; clade noise is found).
   - Cleaning: move `gls`/`fit_lambda` out of the experiment, make `LAMBDAS`/`TERMS` parameters, and take lineages from biothings-mcp `TaxonTools` rather than raw E-utilities.
   - Caveat: this is a taxonomy, not a dated tree, so a real tree should be an option. Effort: small. Value: medium (AnAge / synergy-age style species comparisons; Python has little for Pagel lambda).
6. **Out-of-core pseudobulk** from h5ad CSR (`atac_cortex.group_sums`, `pseudobulk`) and from 10x CSC (`atac_hippocampus.subclass_peak_sums`).
   - Dependencies: h5py, numpy. Tests: `test_atac_cortex` (checked against dense sums), `test_atac_hippocampus`.
   - Cleaning: `obs` keys `sample_id`/`cell_type`, BLOCK size and stderr prints are hard-coded; the CSR and CSC variants should be unified. Effort: medium.
   - Value: medium, for rna-pipeline's atlas only if memory bites. scanpy `sc.get.aggregate` and decoupler cover in-memory cases, so offer it (d) to scanpy as a backed mode before keeping a fork.
7. **Per-person shield score** (`protection.py` `deviations`, `fit_reference`/`Reference.adjust`; `experiments/per_person.score_against`). It scores one sample's set CpGs against matched CpGs, relative to a young reference.
   - Dependencies: numpy. Tests: `test_protection` (synthetic).
   - Effort: small once the code is stable. Value: potentially high, but it is an unpublished, still-changing method. Keep it here until the paper, then ship it in the library from candidate 1.
8. **Agent skill "matched gene-set drift test on a new methylation cohort"** (c, longevity-skills). Wrap `docs/extension.md`'s "recipe for new data" together with candidates 1-3: fetch a GSE, parse it, mask probes, stratify, test. Effort: small after the library exists. Value: medium.
9. **Interval membership** (`data/rrbs.py` `in_intervals`, `_merge`). Do not extract: replace it with polars-bio overlap (d) inside entropy-age. Effort: small.

## Do NOT extract

- **Clock coefficients and scoring** (`data/clocks.py`, `published.LinearClock`, `clock.HorvathAge`). These duplicate just-biomarkers, which uses the same biolearn pin and is parity-tested. entropy-age should depend on it instead.
- **`data/download.py`**: rna-age's downloader is stronger (resume, checksums); converge on that one.
- **`experiments/*`**: figure scripts built on hard-coded GSM IDs, thresholds and constants.
- **`data/cohorts.py`, `mouse.py`, `rrbs.py`, `upstream.py`, `gtex_rna.py`, `chromatin.py`**: per-dataset quirks and pinned files. They are useful only as worked examples.
- **`simulate/*` and `clock.train_clock`**: the Meyer & Schumacher reimplementation and a thin sklearn wrapper; they belong with the paper. A note upstream about the bugs fixed is enough.
- **`paths.py`, `plots.py`, `cli.py`**: project scaffolding and house style.

## Recommendation

Make one small new library, e.g. `methyl-drift-tools` (polars + numpy + scipy, no I/O paths). It would hold candidates 1, 4 and 5 now, and 7 after publication; they share the "feature matrix, strata, matched contrast" core and are already well tested, so the move is mostly cutting imports (about 1-2 days).

Send the I/O pieces (2, 3) to just-biomarkers rather than a new repo, so that there is one polars methylation-array stack: GEO reader, manifest, mask, clocks. entropy-age then depends on both. Only after that, add the longevity-skills skill (8) and optionally a GEO-samples tool in holy-bio-mcp. Hold the pseudobulk code (6) until rna-pipeline asks for it, and offer it to scanpy first. Writing to rna-age is out of scope; there, reuse runs the other way (its downloader and age parser).

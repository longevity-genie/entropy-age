# Changelog

Newest first.

## Unreleased

- RNA test (roadmap RM4): `entropy-age figures rna_gtex` measures, in GTEx v10 RNA-seq (50 tissues with
  100 or more donors), how fast between-donor spread of each expressed protein-coding gene grows with
  age, and compares core-essential, housekeeping and non-essential genes with genes of the same
  expression decile, by turnover class (renewing, slow, heart & muscle, brain), in subsamples of 70
  donors per tissue. Counts are read from an external parsed copy named by `ENTROPY_AGE_GTEX_RNA`, and
  RIN, ischemic time and Hardy scale from GTEx's annotation tables named by `ENTROPY_AGE_GTEX_ATTRIBUTES`;
  either missing is a refusal. `entropy-age fetch annotation` now also fetches HGNC's complete set
  (release 2026-10-06) and puts the gene sets on Ensembl IDs (Entrez first, then current, previous and
  unique alias symbols: 99-100 % of each set matched). New metric `stochasticity.expression_drift`.
- Human tissue test (roadmap RM2): `entropy-age fetch cohort gtex` (GSE213478, nine GTEx tissues of 424
  donors on EPIC, decade age bins as midpoints) and `entropy-age figures gtex`, which asks whether
  skeletal muscle, the one post-mitotic tissue, is where the protection is missing. Every tissue is
  resampled at the same size (32 donors, without replacement), because fewer donors shrink SD-scaled
  differences toward zero. `entropy-age fetch cohort fibres` (GSE304045, type I and type II fibres and
  whole muscle of the same people, EPIC v2) and `entropy-age figures fibres` compare between-person
  spread in fibres against whole muscle.
- Mouse results recorded (`docs/extension.md`, `docs/figures/mouse_tissues.png`): in 11 tissues of
  the same mice, essential and housekeeping CpGs gain entropy 0.12–0.16 SD less than matched CpGs in
  renewing than in post-mitotic tissue (8 of 9 aging contrasts exclude zero), slow tissues are as
  protected as renewing ones, and cortex and hippocampus show none. The rapamycin arm (5 mice per
  tissue) is within small-sample noise against size-matched controls.
- Composition check recorded (`docs/extension.md`, `docs/figures/tissue_contrast_cells.png`): the
  same-donor blood-minus-brain gap keeps all 12 intervals on the predicted side under every adjustment;
  the lifespan-slope gap holds with age-free components but not with components that keep their age
  trend, so it is now presented as supporting evidence only. Roadmap RM1 closed.
- Mouse multi-tissue test (roadmap RM2/RM3): `entropy-age fetch mouse-annotation` (the mammalian
  array's mm10 probe table, GEO GPL32543, and MGI mouse-human orthologs that carry the human gene sets
  to mouse), `entropy-age fetch cohort huang` (GSE279437: 11 tissues of the same mice, young and aged,
  exenatide and rapamycin), and `entropy-age figures mouse`. The figure compares the gene sets' aging
  drift and spread across tissues ordered by cell turnover, with bootstrap contrasts between turnover
  classes and between treatments.
- Roadmap RM1–RM4: record the composition check, then a turnover gradient across tissues, long- vs
  short-lived rodent groups, and stochasticity in RNA, each with the public datasets that would test it.
- Composition check for prediction 2: `entropy-age figures tissues_cells` reruns the blood-versus-brain
  test with the leading principal components of each tissue's background CpGs as extra covariates
  (5 or 10 components, keeping or removing their age trend; `composition_components` in
  `stochasticity.py`), the same reference-free correction in every tissue. `ENTROPY_AGE_BOOT` overrides
  the number of bootstrap resamples, for smoke tests.
- Prediction 2 tested (`docs/extension.md`, `docs/figures/tissue_contrast.png`): the essential and
  housekeeping protection is larger in blood than in brain, and the non-essential excess too. This
  holds across Hannum blood versus Jaffe cortex (age slopes) and within the same Lunnon donors' blood
  and four brain regions (dispersion), with 20 of 21 bootstrap intervals excluding zero. Blood's
  unadjusted cell mix remains a possible partial explanation.
- Fixed: the tissue experiment's bootstrap workers deadlocked when forked (polars' thread pool was
  running in the parent; 16 workers slept after 19 CPU-seconds each on 2026-11-09). They are now
  spawned and receive their data through an initializer; this replaces the "explicit fork" entry
  below. A test checks that parallel and serial bootstraps give identical results.
- `read_supplementary_betas` skips a preamble and ID rows after the header, which GSE59685's beta file
  needs (comment lines, a barcode header, then two rows of GEO accessions).
- `compare_groups` finds group members with a hash lookup instead of `np.isin` on strings, which took
  seconds per group at 10^5 features (a 17 s profile now takes 2 s). The tissue experiment adjusts Hannum
  blood for recruitment site as well as plate, and starts its bootstrap workers with an explicit fork.
- Prediction 2 pipeline (essential-gene protection in renewing versus post-mitotic tissue), results
  pending: `entropy-age fetch cohort hannum jaffe lunnon` (450k blood and brain cohorts, GB-scale, for a
  cluster), the 450k manifest in `fetch annotation` with promoter CpG-to-gene links, an age-free
  `dispersion_levels` measure, and `entropy-age figures tissues` with bootstrap blood-minus-brain
  intervals. GEO characteristic keys are now cut before any parenthesis (`age (in years)` -> `age`).
- `docs/extension.md` gains a related-work section, checked against PubMed and Meyer's 2024 thesis.
  No stochastic clock has yet been validated on mammalian RNA (the only application is Ostermann et
  al., JCI Insight 2025, on human induced neurons). The worm neuron-type paper is the closest
  cell-type analogue, and the gene-level noise literature (CGI− genes, chromatin variability,
  sub-Poissonian genes, CpG density and lifespan) is described from each paper's own abstract.
- Prediction 1 of the selection hypothesis tested on real data (`docs/extension.md`): `entropy-age
  fetch annotation` (Illumina 27k manifest from GEO GPL8490; Hart CEGv2 and NEGv1 pinned at
  hart-lab/bagel 53388ad; Eisenberg–Levanon housekeeping genes) and `entropy-age figures genesets`.
  In GSE41037 blood, essential and housekeeping CpGs gain entropy and between-person dispersion
  more slowly than CpGs matched on CpG-island status and starting level, and non-essential CpGs
  faster; this holds after batch adjustment, with cell fractions, and inside each processing batch.
- `drift_metrics(covariates=...)` gives partial age slopes; `compare_groups` takes a categorical
  stratum column and runs its stratified permutations vectorised.
- Fixed: GSE41037's three processing datasets are confounded with age. The atlas now compares clock
  CpGs using batch-adjusted slopes, and the earlier claim that clock CpGs gain entropy faster than
  level-matched CpGs (Horvath p = 0.004) is withdrawn (batch-adjusted p = 0.60). Fig. 5d's title
  said methylome entropy rises with age; the pooled r is 0.05 (within batches 0.17, 0.00, 0.01), and
  the panel now says so per batch. Fig. 5 also writes the clock's accuracy within each batch
  (r = 0.75, 0.86, 0.51). A substituted ground state is now announced on stderr.
- Adopted the seats kit (agent-seat-blueprints v0.1.144): `CLAUDE.md` and `AGENTS.md`, the roadmap,
  changelog and decisions files, the `lead` seat's priming under `docs/seats/`, the seat table
  `.blueprints.toml`, and the kit's plugin pinned in `.claude/settings.json` (D1-D4). Not yet
  run: no seat has been taken; the Zulip organisation and the desk do not exist yet.

# Changelog

Newest first.

## Unreleased

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

# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

## R4 — the chromatin-matched gap with sorted brain chromatin (RM7), 2026-11-11

- Commit b0184a3; `figures chromatin --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set.
- Job 6428474 completed in 27.8 min, 5.2 GB at peak (the next run fits in 10 GB). It rewrote R3's outputs with the
  sorted-ATAC matchings added; R3's own numbers are unchanged in them.
- Result, as measured (housekeeping, brain): sorted neuron and glia accessibility removes 30-37 % of the gap,
  with the bulk state 32-49 %. The lead's reading is in `docs/extension.md`.

## R3 — the chromatin-matched gap (RM7), 2026-11-11

- Commit 105cc56; `fetch annotation` (re-prepared the probe table with positions, 5 s, 1.8 GB), then
  `figures chromatin --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set.
- Job 6428419 completed in 20.5 min, 5.2 GB at peak: the next run fits in 10 GB. About 109,000 promoter CpGs per
  tissue in Hannum, Jaffe and Lunnon's five tissues.
- Result, as measured (housekeeping, dispersion slope in SDs): blood −0.289 matched on island and level, −0.034
  matched on its own chromatin (E062) as well; cortex −0.075 and −0.058. The lead's reading is in
  `docs/extension.md`. Tables in `data/output/chromatin_match/`, figure `docs/figures/chromatin_match.png`.

## R2 — the per-person score (RM6), 2026-11-11

- Commit 2fe148a (module `entropy_age.protection`, figure `per_person`); `--scale paper --seed 42`; 8 cores, 48 GB
  asked; job comment set by the account's rule.
- Job 6428234 completed in 13 min (88 CPU-minutes, 6.8 cores busy on average), 5.9 GB at peak: the next run fits in
  12 GB. Scored 632 Hannum and 159 Jaffe samples (references: 24 and the Jaffe controls aged 20-35), and 71-75
  Lunnon donors per region.
- Result, as measured (housekeeping, no components): P per decade −0.0018 (SE 0.0003) in blood, −0.00036
  (SE 0.00013) in cortex; reliability 0.997 and 0.955; blood's F below the same donor's region in 86-91 % of
  donors. The lead's reading is in `docs/extension.md`. Tables in `data/output/per_person/`, figure
  `docs/figures/per_person.png`.

## R1 — GSE179847 (RM5), 2026-11-10

- Commit b7036e1 (cohort `sturm`, figure `passage_sturm`); `--scale paper --seed 42`; 8 cores, 64 GB asked as the
  first measure of the parse's memory.
- Jobs 6427994 (parse) and 6427995 (figure, after it): the parse failed in 1 s, before any work: the job script
  sourced the shell profile under `set -u`, which the system profile does not survive. Fixed (`set -u` after the
  profile); the figure job cancelled.
- Resubmitted as 6427996 (parse) and 6427997 (figure, after it): both completed. The parse took 28 s and
  4.4 GB at peak; the figure 43 s and 1.9 GB (the next run of either fits in 8 GB). The job log printed
  "commit unknown", because the CPU nodes have no git; the clone was at ffbc2c5 (code as b7036e1).
- Result, as measured (housekeeping set, the five baseline lines): coupling change per 10 PD −0.055, interval
  −0.21 to +0.10; mean coupling −0.27. The lead's reading is in `docs/extension.md`. Tables in
  `data/output/passaged_cells_sturm/`, figure `docs/figures/passaged_cells_sturm.png`.

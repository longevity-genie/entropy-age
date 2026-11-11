# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

## R7 — interventions in mouse blood RRBS (RM3), 2026-11-11

- Commit 30ed4b8; `fetch rrbs` (243 blood tables, about 3 GB, and the mm10 tracks) then
  `figures interventions --scale paper --seed 42`; 8 cores, 64 GB asked; job comment set.
- Job 6429040 completed in 53 min: the fetch and parse 15.9 min at 19.7 GB peak, the figure 37.4 min at 9.5 GB. A
  rerun of the figure fits in 16 GB. 542,008 promoter CpGs covered at 10x or more in 90 % of the 243 samples.
- Result, as measured (housekeeping): CR moves the shield share F from −0.060 to −0.032 in C57BL/6 (p 0.002) and by
  +0.025 in B6D2F1 (p 0.012), with drift B changed by −7 % and −3 % (n.s.). Outputs superseded by R8, which adds a
  composition check.

## R6 — clocks from protected and unprotected CpGs (RM8), 2026-11-11

- Commit 7d40536; `figures chromatin_clocks --scale paper --seed 42`; 8 cores, 32 GB asked; job comment set. Waited
  on the account's shared CPU cap before it started.
- Job 6428565 completed in 39.4 min, 2.8 GB at peak (the next run fits in 6 GB). Pools of 33,786 (active in
  both), 19,330 (Polycomb in both) and 149,740 (any) shared CpGs, 4,000 per clock, 5 draws.
- Result, as measured: no pool transfers well between tissues (out-of-tissue r 0.06-0.70); the Polycomb pool's
  brain-minus-blood offset in Lunnon varies least between draws (SD 4-6 years against 20-42). The lead's reading
  is in `docs/extension.md`.

## R5 — published clocks on the protection map (RM8), 2026-11-11

- Commit 4787a9b; `figures clocks --scale paper --seed 42`; 8 cores, 24 GB asked; job comment set. Waited about an
  hour on the account's shared CPU cap (`AssocGrpCpuLimit`) before it started.
- Job 6428487 completed in 11.8 min, 3.2 GB at peak (the next run fits in 6 GB). 16 clocks, 64,900 CpGs per cohort.
- Result, as measured: clock CpGs switch TSS activity between blood and cortex at 0.15-0.33 against 0.12 for
  background; the tissue offset does not follow that share (Spearman 0.11, p 0.68, 16 clocks). The lead's reading
  is in `docs/extension.md`.

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

# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

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

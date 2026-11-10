# Run log

Owned by the runner seat (`docs/seats/RUNNER.md`): one entry per run, newest first, with its inputs, the
commit it ran, its job ids, its scale and seed, and its result or its failure with the cause.

## R1 — GSE179847 (RM5), 2026-11-10

- Commit b7036e1 (cohort `sturm`, figure `passage_sturm`); `--scale paper --seed 42`; 8 cores, 64 GB asked as the
  first measure of the parse's memory.
- Jobs 6427994 (parse) and 6427995 (figure, after it): the parse failed in 1 s, before any work: the job script
  sourced the shell profile under `set -u`, which the system profile does not survive. Fixed (`set -u` after the
  profile); the figure job cancelled.
- Resubmitted as 6427996 (parse) and 6427997 (figure, after it): running.

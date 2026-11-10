# Runs

The brief for the runner seat (`docs/seats/RUNNER.md`): the batch runs to do on the cluster, in order,
written by the maintainers and the lead. The runner works it top down and records what happened in
`docs/RUN_STATUS.md` and `docs/RUN_LOG.md`; a run that is done leaves this file.

Each entry is a section `## Rn — name` whose first line is `**State:** ... · **Owner:** ...` in the kit's plan
grammar (`unattended`, `urgent`, `building`, `after <R-id>` or `after job <slurm id>`, `waits until YYYY-MM-DD`,
`blocked: desk #N`), read by the account's queue table, with:

- the command (`uv run entropy-age fetch ...` or `uv run entropy-age figures NAME --scale paper --seed S`);
- the roadmap item it serves (`RMn`);
- what the lead needs back (a table, a figure, a number), and anything the run must not do.

## R1 — GSE179847, the cellular-lifespan replication of the culture test

**State:** unattended · **Owner:** runner

- Serves **RM5** (a second, larger replication of the division-coupling test of `figures passage`).
- Command, from the clone after `git pull`: `uv run entropy-age fetch cohort sturm` (the two GEO files are already
  in `data/input/geo/GSE179847/`, so this only parses: the 2.0 GB matrix to `data/interim/cohorts/sturm_*.parquet`),
  then `uv run entropy-age figures passage_sturm --scale paper --seed 42`. Measure the parse's peak memory first
  (it holds 479 EPIC samples as floats, a few GB); the figure step is smaller.
- Back to the lead: the per-arm slopes of `P`, `B` and `F` on PD and on days, with intervals, as for
  GSE197512, the CSV summary under `data/output/` and the run's peak memory.
- Must not: commit the matrix or any parquet; use more than the runner's share.

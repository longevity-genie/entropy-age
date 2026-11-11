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

## R15 — interventions against same-chip controls (RM12)

- What: `entropy-age figures mmc_interventions --scale paper --seed 42`, on R12's download, after the same-chip
  control mode; replaces R13's tables and figure (the pooled rows are recomputed alike).
- Resources: compute partition, 8 cores, 16 GB (R13 peaked at 3.5 GB), a few minutes.
- Read: the `same_chip` rows of `data/output/mmc_interventions/effects.csv`, Tet3 first.

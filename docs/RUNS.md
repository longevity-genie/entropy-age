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

## R18 — chromatin-writer syndromes in blood, with age (RM14)

- What: `entropy-age figures syndromes --scale paper --seed 42` after people without an age are left out of a
  series that gives ages (GSE97362: one sample).
- Resources: compute partition, 8 cores, 8 GB (R17: 1.3 GB, 26 min).
- Read: `data/output/syndromes/effects.csv`; `age_covariate` true for GSE97362, false for GSE74432.

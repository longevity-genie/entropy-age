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

## R17 — chromatin-writer syndromes in blood (RM14)

- What: `entropy-age fetch syndromes` (two series matrices, 0.7 GB, already on disk), then
  `entropy-age figures syndromes --scale paper --seed 42`.
- Resources: compute partition, 8 cores, 32 GB; time unknown (5-component leave-one-out over up to 125 controls).
- Read: `data/output/syndromes/effects.csv`, F and S by arm, raw and with components; GSE74432 has no age
  covariate.

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

## R16 — which way essential CpGs move under Tet3 loss (RM12)

- What: `entropy-age figures mmc_interventions --scale paper --seed 42` after the signed shift `S` was added;
  replaces R15's tables (the B and F rows are recomputed alike).
- Resources: compute partition, 8 cores, 16 GB, a few minutes (R15: 73 s, 3.5 GB).
- Read: `S` rows for Tet3KO.Het, same-chip controls: above 0 means essential CpGs gained methylation.

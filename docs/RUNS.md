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

## R22 — look-alike drift and relative gap against lifespan (RM20)

- What: `entropy-age figures species --scale paper --seed 42` (adds `gap`, `background_rate` and the RM20
  outcomes), then `entropy-age figures species_phylo` on the new table (login node, a minute).
- Resources: compute partition, 8 cores, 40 GB (R12: 59 min, 25.6 GB).
- Read: `lifespan_fits` with outcome `log background rate` (catch-up and fixed share both predict < 0) and
  `relative gap` (catch-up predicts > 0, i.e. toward zero; fixed share predicts 0). The `gap_sd` rows must equal
  R12's.

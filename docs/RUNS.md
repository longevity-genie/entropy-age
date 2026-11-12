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

## R24 — cortical cell types with the promoter read fraction as a covariate (RM19)

- What: `entropy-age figures atac_cortex --scale paper --seed 42` on R23's cached pseudobulk, now with four
  adjustments (none, 5 components, promoter read fraction, both).
- Resources: compute partition, 8 cores, 32 GB; minutes (the pseudobulk is cached).
- Read: does the positive gap of R23 survive the promoter read fraction (a quality score like TSS enrichment)?

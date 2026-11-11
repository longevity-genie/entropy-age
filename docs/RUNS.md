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

## R6 — clocks from protected and unprotected CpGs

**State:** building · **Owner:** runner (job 6428565)

- Serves **RM8**, second part.
- Command, from the clone at 7d40536 or later: `uv run entropy-age figures chromatin_clocks --scale paper --seed 42`.
- Back to the lead: `data/output/chromatin_clocks/` and `docs/figures/chromatin_clocks.png`, time and peak memory.


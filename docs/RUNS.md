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

## R8 — interventions in mouse blood RRBS, with a composition check

**State:** building · **Owner:** runner (job 6429130)

- Serves **RM3**.
- Command, from the clone at d5c7433 or later: `uv run entropy-age figures interventions --scale paper --seed 42`
  (R7 fetched and parsed the data).
- Back to the lead: `data/output/interventions/` and `docs/figures/interventions.png`, time and peak memory.

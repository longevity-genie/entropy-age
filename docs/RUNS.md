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

## R9 — the shield with age and diet in BXD mouse liver

**State:** building · **Owner:** runner (job 6429222)

- Serves **RM3**.
- Command, from the clone at e8b851d or later: `uv run entropy-age figures bxd --scale paper --seed 42` (its input, the
  GSE199979 series matrix, is fetched).
- Back to the lead: `data/output/bxd_liver/` and `docs/figures/bxd_liver.png`, time and peak memory.

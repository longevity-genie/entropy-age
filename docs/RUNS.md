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

## R10 — six mouse tissues with numeric ages (GSE120132)

**State:** building · **Owner:** runner (job 6429231)

- Serves **RM2**.
- Commands, from the clone at d3e052d or later: `uv run entropy-age fetch rrbs-tissues` (359 tables, about 6.4 GB,
  through the proxy, parsed per tissue), then `uv run entropy-age figures rrbs_tissues --scale paper --seed 42`.
- Back to the lead: `data/output/mouse_rrbs_tissues/` and `docs/figures/mouse_rrbs_tissues.png`, time and peak
  memory of each step (64 GB is a first guess, from R7's 19.7 GB for 243 tables).

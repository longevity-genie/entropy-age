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

## R3 — the chromatin-matched gap on Hannum, Jaffe and Lunnon

**State:** building · **Owner:** runner (job 6428419)

- Serves **RM7** (does chromatin state account for the protection?).
- Commands, from the clone at 105cc56 or later: `uv run entropy-age fetch annotation` (re-prepares the 450k probe
  table with positions; the inputs are present), then `uv run entropy-age figures chromatin --scale paper --seed 42`.
  The chromHMM segmentations are already fetched.
- Back to the lead: `data/output/chromatin_match/` (`matched_gaps`, `state_composition`) and
  `docs/figures/chromatin_match.png`, with time and peak memory.

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

## R2 — the per-person protection score on Hannum, Jaffe and Lunnon

**State:** building · **Owner:** runner (job 6428234)

- Serves **RM6** (validation steps 1-3 of `docs/PER_PERSON_SCORE.md`).
- Command, from the clone at 2fe148a or later: `uv run entropy-age figures per_person --scale paper --seed 42`. Its
  inputs (the three cohorts' betas, the 450k links) are already prepared.
- Back to the lead: `data/output/per_person/` (`age_fits`, `same_donor_turnover`, the two score tables) and
  `docs/figures/per_person.png`, with the job's time and peak memory.

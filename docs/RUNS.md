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

## R43 — does methylation noise reach expression? GTEx same donors (RM30)

- What: `entropy-age figures meth_expr_spread --scale paper --seed 42` (GTEx methylation and RNA; donors with
  both: lung 161, ovary 133, colon 98, prostate 62, blood 45).
- Resources: compute partition, 4 cores, 64 GB; time unknown (200 permutations per tissue).
- Read: `data/output/meth_expr_spread/coupling`: signed coupling (control), spread coupling against its null,
  spread-slope coupling across genes.



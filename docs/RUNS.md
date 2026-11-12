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


## R48 — do drift and the shield mark later death? SJLIFE (RM28)

- What: GSE314261 processed table (28 GB) fetched; `entropy-age figures sjlife --scale paper --seed 42`.
- Resources: compute partition, 4 cores, 96 GB; time unknown (streams the gzip once, about 70 GB of text).
- Read: `data/output/sjlife_deaths/effects`: odds ratio of death per SD of B, F, S and clocks; everyone and
  survivors; deaths over 2 years after the draw.

## R49 — the shield by blood lineage, in silico (RM31)

- What: `entropy-age figures lineage_shield --scale paper --seed 42` (Hannum betas and R38's cell fractions).
- Resources: compute partition, 4 cores, 32 GB; minutes.
- Read: `data/output/lineage_shield/gaps`: set minus look-alikes on |age slope within lineage|, per lineage.

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

## R42 — does the shield predict death? Vitality 90+ (RM28)

- What: GSE68194 series matrix fetched (386 MB); `entropy-age figures survival --scale paper --seed 42`.
- Resources: compute partition, 4 cores, 32 GB; time unknown (115 leave-one-out references).
- Read: `data/output/survival_shield/effects`: hazard ratio per SD of F, B, S and clocks; MMSE.

## R43 — does methylation noise reach expression? GTEx same donors (RM30)

- What: `entropy-age figures meth_expr_spread --scale paper --seed 42` (GTEx methylation and RNA; donors with
  both: lung 161, ovary 133, colon 98, prostate 62, blood 45).
- Resources: compute partition, 4 cores, 64 GB; time unknown (200 permutations per tissue).
- Read: `data/output/meth_expr_spread/coupling`: signed coupling (control), spread coupling against its null,
  spread-slope coupling across genes.

## R45 — when does drift begin? change point against sexual maturity (RM29)

- What: `entropy-age figures drift_onset --scale paper --seed 42` (consortium beta table on disk).
- Resources: compute partition, 4 cores, 64 GB; time unknown (reads the 4.5 GB table once, then per group).
- Read: `data/output/drift_onset/onsets` (per species and tissue: change point, slopes, delta BIC, onset),
  `across_species` (log onset on log maturity and on log lifespan).

## R46 — does the shield predict cancer years ahead? EPIC-Italy (RM28)

- What: GSE51032 series matrix fetched (3.1 GB, betas inside); `entropy-age figures cancer --scale paper --seed 42`.
- Resources: compute partition, 4 cores, 64 GB; time unknown (845 leave-one-out references).
- Read: `data/output/incident_cancer/effects`: odds ratio per SD for any cancer, breast, colorectal; all cases and
  cases diagnosed over 2 years after the draw.

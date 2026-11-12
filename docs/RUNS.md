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

## R33 — clocks of our own from shielded, look-alike and random CpGs (RM11)

- What: `entropy-age figures own_clocks --scale paper --seed 42`.
- Resources: compute partition, 8 cores, 64 GB; time unknown (2 cohorts x 3 draws x 3 pools x 5 folds of
  elastic-net CV on a few thousand CpGs).
- Read: `data/output/own_clocks/accuracy` (r, median absolute error per pool) and `acceleration_overlap`.

## R35 — Tet1 knockout and catalytic-dead cortex (RM12)

- What: `entropy-age figures tet1 --scale paper --seed 42` (GSE271375 betas on disk).
- Resources: compute partition, 8 cores, 32 GB; minutes.
- Read: `data/output/tet1_cortex/effects`: F of KO and Hxd against WT; Tet3 heterozygotes lost the shield.

## R36 — Tet3 signed shift with array compression held fixed (RM18)

- What: `entropy-age figures mmc_interventions --scale paper --seed 42` (GSE223748 on disk).
- Resources: compute partition, 8 cores, 16 GB; about 75 s (R16).
- Read: `S, compression held` rows for Tet3KO.Het same-chip: does essential CpGs' gain survive?

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

## R36 — Tet3 signed shift with array compression held fixed (RM18)

- What: `entropy-age figures mmc_interventions --scale paper --seed 42` (GSE223748 on disk).
- Resources: compute partition, 8 cores, 16 GB; about 75 s (R16).
- Read: `S, compression held` rows for Tet3KO.Het same-chip: does essential CpGs' gain survive?

## R37 — Dnmt3a gain of function and the shield (RM12)

- What: `entropy-age figures dnmt3a --scale paper --seed 42` (GSE324160 on disk).
- Resources: compute partition, 8 cores, 32 GB; minutes.
- Read: `data/output/dnmt3a_gain/effects`: F, B, S of mutants against wild types per tissue.

## R38 — the stochastic component of clocks, with information stripped (RM23)

- What: `entropy-age fetch epidish` (done), `entropy-age figures stochastic_component --scale paper --seed 42`.
- Resources: compute partition, 8 cores, 32 GB; time unknown (3 clocks x 4 modes x a 63-point calibration grid).
- Read: `data/output/stochastic_component/fractions`: RR2 and the correlation with the original clock's intrinsic
  acceleration per clock and model; the published Stoc clocks as the reference row.

## R39 — does the shield track health? cross-sectional (RM24)

- What: `entropy-age figures health --scale paper --seed 42` (cohorts on disk: Lunnon, Jaffe, GSE41037, GTEx).
- Resources: compute partition, 4 cores, 64 GB; time unknown (leave-one-out references, about 2,500 fits per
  composition setting).
- Read: `data/output/health_states/effects`: F, B, S effects of Alzheimer's (and Braak stage), schizophrenia and
  smoking and slow against sudden death (GTEx Hardy scale), raw and with 5 components.

## R41 — do organs of one person age together? (RM26)

- What: `entropy-age figures organ_coupling --scale paper --seed 42` (Lunnon, GTEx on disk; clocks fetched).
- Resources: compute partition, 4 cores, 64 GB; time unknown (leave-one-out drift per tissue).
- Read: `data/output/organ_coupling/couplings`: r between tissues of the same donors per measure.

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

## R44 — the shield, drift and clocks in sorted blood cells (RM31)

- What: MESA normalized tables fetched (GSE56046 4.4 GB, GSE56581 0.8 GB); `entropy-age figures sorted_cells
  --scale paper --seed 42`.
- Resources: compute partition, 4 cores, 64 GB; time unknown (about 1,200 leave-one-out references in monocytes).
- Read: `data/output/sorted_cells/shield`, `clocks`, `health` (plaque, coronary calcium).

## R45 — when does drift begin? change point against sexual maturity (RM29)

- What: `entropy-age figures drift_onset --scale paper --seed 42` (consortium beta table on disk).
- Resources: compute partition, 4 cores, 64 GB; time unknown (reads the 4.5 GB table once, then per group).
- Read: `data/output/drift_onset/onsets` (per species and tissue: change point, slopes, delta BIC, onset),
  `across_species` (log onset on log maturity and on log lifespan).

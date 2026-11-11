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

## R12 — shielding across 349 mammal species (GSE223748)

**State:** building · **Owner:** runner (job 6429555)

- Serves **RM10** (urgent, asked by antonkulaga).
- Commands, from the clone at fa6a796 or later: `uv run entropy-age fetch mammals` (about 4.5 GB through the proxy),
  then `uv run entropy-age figures species --scale paper --seed 42`.
- Back to the lead: `data/output/species_shield/` and `docs/figures/species_shield.png`, time and peak memory.

## R13 — mouse interventions across tissues (GSE223748)

**State:** after R12 · **Owner:** runner (job 6429673, starts when R12 ends)

- Serves **RM12** and **RM3**: `uv run entropy-age figures mmc_interventions --scale paper --seed 42`, from the clone at
  b1c9ae0 or later, on the data R12 fetches.
- Back to the lead: `data/output/mmc_interventions/` and `docs/figures/mmc_interventions.png`, time and peak memory.


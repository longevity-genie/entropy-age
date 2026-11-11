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

## R4 — the chromatin-matched gap with sorted brain chromatin

**State:** building · **Owner:** runner (job 6428474)

- Serves **RM7** (is the gap left in brain chromatin a bulk map cannot resolve?).
- Command, from the clone at b0184a3 or later: `uv run entropy-age figures chromatin --scale paper --seed 42` (all
  inputs fetched; it rewrites R3's outputs with the new matchings added).
- Back to the lead: `data/output/chromatin_match/` and `docs/figures/chromatin_match.png`, time and peak memory.

## R5 — published clocks on the protection map

**State:** building · **Owner:** runner (job 6428487)

- Serves **RM8** (where published clocks' CpGs sit by chromatin, gene set and drift; their tissue offsets).
- Command, from the clone at 4787a9b or later: `uv run entropy-age figures clocks --scale paper --seed 42` (the
  clocks are fetched, `entropy-age fetch clocks`).
- Back to the lead: `data/output/clock_cpgs/` and `docs/figures/clock_cpgs.png`, time and peak memory.


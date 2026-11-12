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

## R20 — headline gaps against probe design and unreliable probes (RM18)

- What: `entropy-age figures probe_robustness --scale paper --seed 42` (annotation re-prepared with `probe_type`;
  mask at `data/input/annotation/zhou/`).
- Resources: compute partition, 8 cores, 64 GB; time unknown (no bootstrap, 1000 permutations x 4 variants).
- Read: `data/output/probe_robustness/gaps`: does any housekeeping or core-essential gap change sign or lose most
  of its size with probe type in the strata or masked probes out? `design_share`: are set CpGs type I more often?

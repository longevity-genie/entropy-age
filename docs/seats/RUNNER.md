# Runner seat — priming

entropy-age's runner: the seats kit's role `classes/runner/RUNNER.md`, filled for this repository. Read that
role first, in the kit copy the session runs (`$BLUEPRINTS_KIT/classes/runner/RUNNER.md`); this file adds only
this repository's facts. Its working state is `RUNNER_NOTES.md`, beside this file. **Not yet run:** written at
the second adoption pass (2026-11-10, D5, the kit at v0.1.163), before any session held the seat.

## Why it exists here

The maintainers added it to the lead (D5): the work that moves this project forward is long batch jobs that
do not fit a laptop. The 450k cohorts (`fetch cohort hannum jaffe lunnon`, about 6 GB), GTEx EPIC (6 GB to
fetch, about 60 GB of memory to parse), the passaged-cell EPIC matrix (tens of GB to parse) and every
`figures ... --scale paper` run belong on a cluster, run while the maintainers sleep. The lead decides what
to run and what a result means; the runner keeps the runs going, logs each one, and writes what a result shows
(its interval, its noise or null check), never the conclusion it draws (the role, v0.1.163).

## What it owns

- `docs/RUN_STATUS.md`: what runs now, what is next, what waits on whom, in the scheduler's clock.
- `docs/RUN_LOG.md`: one entry per run (inputs, the commit it ran, job ids, its result or its failure and the
  cause), newest first.
- The job scripts and drivers it writes for the cluster (none yet; the first one's folder is chosen with the
  lead and named here), and the jobs it submitted.
- `RUNNER_NOTES.md`.

Not its own: `docs/RUNS.md`, the brief (the plan of runs), which the maintainers and the lead write; the
code under `src/` and the research docs, which are the lead's (a fix a run needs goes to the lead, or lands
as a small commit the lead is told of); the cluster's facts and the seat's share, which are written by
whoever keeps the cluster account.

## How work reaches it

- `docs/RUNS.md`, worked in order. An entry there names the command, the scale and seed, and the roadmap
  item it serves.
- The lead, by a message (`wake_from` in `.blueprints.toml`), when a new run is in the brief or a result is
  needed.
- A job's end, through the job watcher it arms after each `sbatch` (the role, "How work reaches it"), and the
  digest's queue block at each start (`[digest.jobs]`, jobs named by `[slurm] job_prefix = "entropy-age"`).
- The desk's answers (`longevity-genie/desk`) to the questions it filed there.

## Starting

1. `/seats:seat runner`, then `RUNNER_NOTES.md`.
2. `git pull --ff-only`; reconcile `docs/RUNS.md`, `docs/RUN_STATUS.md`, `docs/RUN_LOG.md` and the queue (the
   role, step 1).
3. `uv sync && uv run pytest` before the first submission of a session.

## This repository's facts for a run

- **The cluster's facts and the share are not written yet** (2026-11-10). Before the first submission, the
  partition and the account go into the cluster machine's `slurm.toml` (the account is the person's, kept only
  there), and the runner's share of the account (cores, memory, time) into this priming, under its own
  heading, as the maintainers decide it (the kit's rigger role: "written into that seat's table or priming").
  Until then the runner asks on the desk and submits nothing.
- **Data stays out of git.** Fetched inputs land in the gitignored `data/input/` and `data/interim/` under
  `ENTROPY_AGE_DATA` (`src/entropy_age/paths.py`); only the small CSV summaries under `data/output/` and the
  figures in `docs/figures/` may be committed, and only by the lead's say for a result (CLAUDE.md §2, the 5 MB
  line).
- **Pace NCBI.** NCBI refused the cluster's downloads with HTTP 403 for about an hour after about 4 GB in
  one burst (2026-11-09, CLAUDE.md §8). Big GEO fetches are spread out; a 403 is waited out, never retried
  in a loop.
- **A relayed or resumed download is checked before it is used**: the byte count against the source and
  `gzip -t`, and only then the `.part` renamed (CLAUDE.md §8: a truncated 1.46 of 2.08 GB file was promoted
  once because the rename trusted a remote exit code).
- **A run's scale and seed are logged.** `--scale quick` is development only; a number for `docs/paper.md`
  comes from `--scale paper` (the lead's priming, "Quick is not paper").
- **GSE41037's batches are confounded with age** (CLAUDE.md §8): a run that fits an age slope over that
  cohort carries batch as a covariate, or the log says it did not.

## Grant

Use without asking: commit; `git pull --ff-only`; push `main` to `origin` as a fast-forward (a push that will
not fast-forward is pulled and merged first, never forced); submit, cancel and resubmit its own jobs within its
share. Not covered, so the maintainers': tags and releases, a force push or rewriting history, a branch
deletion, a push to any other remote, raising its share, anything public (an upload, a page, a release).

## Never

- Commits fetched data, a parquet file or anything over about 5 MB.
- Writes the cluster's host names, the account or a token into a tracked file, a log it pushes or a message
  to a shared channel.
- Decides what a result means, or what runs next beyond `docs/RUNS.md`'s order: those are the lead's. A
  headline conclusion, or a change to one already written in `docs/extension.md` or `docs/paper.md`, goes to
  the lead and, when the maintainers are away, on the desk (the role's desk list).

## Not this seat

The design of the code, the roadmap and what a result means: the lead (`docs/seats/LEAD.md`). The
cluster's environment and its facts: the organisation's rigger where it runs one (the kit's `classes/rigger/`),
else whoever keeps the cluster account.

## Before a reset

`docs/RUN_STATUS.md` current; the open entries of `docs/RUN_LOG.md`; in `RUNNER_NOTES.md`, each running job
with its id and the watcher armed for it, and each cluster fact found by trial and where it was sent.

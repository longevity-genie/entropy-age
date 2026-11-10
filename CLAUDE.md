# Agent Guidelines — entropy-age

Aging clocks from accumulating stochastic variation (Meyer & Schumacher 2024), reimplemented with uv/polars and extended to gene- and cell-type-specific stochasticity under selection.

`AGENTS.md` is a symlink to this file. If the two ever differ, that is a bug —
`ln -sf CLAUDE.md AGENTS.md`.

---

## Read these first, in this order

Obligatory. Read them yourself. **Do not delegate a document you are about to judge a design
against**: a subagent returns a summary, and a summary of a rule drops the qualifier the decision
turned on. Delegation is for finding things, never for deciding them.

1. **docs/ROADMAP.md** — open work only, in order.
2. **docs/CHANGELOG.md** — what actually landed, newest first.
3. **docs/DECISIONS.md** — the calls settled with the maintainer, each with the alternative it
   rejected.

A charter, when the repository keeps one (a format or a schema, whose invariants outrank any plan),
goes first in this list and wins any conflict with a plan. A per-area reference goes last, read
when the task touches its area.

Everything below is self-contained: no rule here requires following a link to know what you must
not do. Links carry positive detail only.

**The `.claude/` layer is part of this repository's contract.** Its skills, agents, commands and
hooks are named here, each with a line saying what it is for, or they are found by accident. Where
a skill is the canonical copy of a procedure, this file says so in one line and does not restate it.
None existed when this file was written.

Sessions that take a standing role (a seat) all load this file; each seat's own priming lives under
`docs/seats/`, and the seats kit runs them.

---

## 1. When rules conflict: ask, never infer

When two rules conflict (this file against a sibling repository's, against the person's global
preferences, a rule against what the code does), stop and ask. Do not pick the one that looks
better, do not synthesise a compromise, do not silently follow the more specific file. A conflict
between two live rules is almost always a real difference in what the repositories are, and
inferring which applies here is the guess that produces a rule nobody agreed to.

1. **Survey first, ask second.** Read both rules whole and find out why each side holds its
   version. A question that does not carry the reason cannot be answered.
2. **One question per conflict, batched**, never dripped. Each gets two to four concrete options.
3. **Each option states its cost**: what breaks, what it forces elsewhere, which rule it
   contradicts.
4. **Recommend one, and say so.**
5. **Record the answer where it will be read again**: the rule into its section below, the reason
   into §8 in the maintainer's own words, the call into `docs/DECISIONS.md` with the alternative it
   rejected. An answer not written down gets asked again.

A conflict with the charter is not a question: the charter wins.

**Adopting this file is such a moment.** It was laid by the seats kit's `bin/init.py`, and the
sections marked "Not yet written" have more than one right answer. The questions to ask before the
first change: the build, run and test commands; the layout beyond `docs/`; the dependency tier the
code may carry; whether consumers report to this repository (then it keeps a consumer inbox,
which the kit's triage role runs) and which repositories it consumes; whether a `git pull` into it
can break a user's install (then a STOP banner goes at the very top of this file, with the
sentence to hand the person); the sibling repositories, read-only unless the task targets them.

---

## 2. Non-negotiables

Read the whole list before the first edit. The reason follows each rule, because a rule without its
reason gets argued away at 2 a.m.

- **Never claim a test would have caught a bug without running it against the buggy code** and
  watching it fail. A new test is checked by breaking the code it guards.
- **Never silently fall back when the primary input is missing.** Refuse, or say loudly what was
  substituted: the person cannot see that the source differs.
- **Never fill a value from the same source that checks it.** The check then compares a source with
  itself and agrees perfectly, and the value moves from honestly unverified to apparently verified.
  A lookup reports, the person decides, the check checks; a refusal a tool prints is kept verbatim.
- **Never collapse unknown into yes or no.** A check that could not run is not a check that passed;
  a timeout is not a negative answer. Unknown is said as unknown.
- **Never read a reproducible result as a correct one.** A digest match, a strict mode, a
  reproducible build say the bytes came out the same, not that they mean what was intended.
- **Never commit a placeholder path or a made-up example value.** It proves nothing and outlives
  the session that invented it.
- **Never commit large data.** Anything over about 5 MB that must travel goes through Git LFS; a
  large blob in history is very hard to remove afterwards.
- **Never track `.leak_words` or `.leak_allow`.** They name everything the leak gate protects
  (copy `.leak_words.template` to `.leak_words` and fill it; the kit's `bin/leak-gate.sh install`
  turns the gate on). Never commit with `--no-verify`.
- **Git: commit as you go**, one coherent block per commit: a behaviour with its test, doc and
  changelog line. No conventional-commit prefixes. Stage explicit paths, never `git add -A` or
  `git add .`: they sweep in secrets, editor swap files and junk. Local branches are free to create
  and switch. Pushing `main` as a fast-forward is allowed (the maintainer works alone here and wants
  no pull requests: D1); a push that will not fast-forward is a merge first, never a force push.
  Tags, releases, branch deletion and history rewriting are the maintainer's, unless a seat's
  priming grants one in writing, by name. Never `git stash drop` or `git stash clear`, even on
  request.

---

## 3. Layout

```
src/entropy_age/              the package (Python 3.13, uv_build)
  cli.py                      typer CLI: `entropy-age fetch|figures|stochasticity|info`
  clock.py                    elastic-net clocks, age transforms, metrics
  stochasticity.py            per-feature drift, entropy, dispersion, relaxation-rate fits; group contrasts
  transcriptome.py            BiT-age transcriptome helpers (the paper's Fig. 2; its counts are not upstream)
  published.py                published linear clocks, coefficients from the upstream repository
  plots.py, paths.py          figure styling; the one filesystem layout (ProjectPaths)
  simulate/                   bulk noise, binomial single-cell methylation, schedules, selection
  data/                       downloaders: upstream repo files, GEO series matrices, 27k/450k annotation + gene sets, 450k cohorts, matrix frames
  experiments/                one module per figure (fig1, fig3, fig4, fig5, extension, atlas, genesets, tissues, tissues_cells, mouse, gtex, fibres, rna_gtex, rna_rodents, passage)
tests/                        pytest, real code paths on small seeded data
data/input/                   downloads exactly as fetched (gitignored; `entropy-age fetch` recreates)
data/interim/                 parquet derived from the inputs (gitignored)
data/output/<figure>/         result tables: small CSV summaries tracked, parquet ignored
docs/                         all markdown except this file and README.md
docs/figures/                 rendered figures (tracked), written by `entropy-age figures`
docs/seats/                   the seats' primings and notes (the seats kit)
.blueprints.toml              the seat table the seats kit reads (D3)
.claude/settings.json         pins the seats kit's plugin for every clone (D3)
```

`ENTROPY_AGE_HOME` moves the root and `ENTROPY_AGE_DATA` the data folder (`paths.py`).
`figures rna_gtex` reads GTEx v10 RNA-seq it does not download: `ENTROPY_AGE_GTEX_RNA` names a parsed
folder (`samples.parquet`, `genes.parquet`, `counts/part-*.parquet` + `index.json`, the layout of the
rna-age project) and `ENTROPY_AGE_GTEX_ATTRIBUTES` the folder with GTEx's SampleAttributesDS and
SubjectPhenotypesDS tables (`data/gtex_rna.py`). Without them the figure refuses and runs nothing.

Upstream: [Meyer-DH/StochasticAgingClock](https://github.com/Meyer-DH/StochasticAgingClock), read
only, pinned by commit in `src/entropy_age/data/upstream.py`; its data files are downloaded at run
time and never redistributed here (`NOTICE`). GEO series matrices come from NCBI's FTP. No outside
consumers report to this repository yet, so it keeps no consumer inbox. A `git pull` cannot break
anyone's install: nothing installs from it.

---

## 4. Build, run, test

```bash
uv sync                                   # environment from uv.lock
uv run pytest                             # the tests (about 30, a few seconds)
uv run ruff check                         # lint (line length 120, py313)
uv run entropy-age info                   # where data lives and what is already prepared
uv run entropy-age fetch upstream         # upstream clock coefficients and sample sheets -> data/input, data/interim
uv run entropy-age fetch geo GSE41037     # a GEO series matrix -> parquet
uv run entropy-age fetch annotation       # 27k and 450k probe manifests, HGNC genes + essential/non-essential/housekeeping gene sets
uv run entropy-age fetch cohort hannum jaffe lunnon   # 450k blood/brain cohorts, ~6 GB: run on a cluster
uv run entropy-age fetch mouse-annotation  # mammalian-array probe table (~200 MB) + MGI orthologs
uv run entropy-age fetch cohort huang      # GSE279437 mouse tissues, 1.2 GB
uv run entropy-age fetch cohort gtex fibres  # GSE213478 GTEx EPIC (6 GB, needs ~60 GB memory) + GSE304045 muscle fibres (0.6 GB)
uv run entropy-age fetch rodent-rna         # GSE314817 / GSE314970 mouse and rat RNA atlases, ~80 MB each
uv run entropy-age fetch cohort endicott   # GSE197512 passaged primary cells, EPIC series matrix 2.1 GB (parse needs tens of GB)
uv run entropy-age figures [NAMES] --scale quick|paper --seed 42
                                          # regenerate figures (fig1 fig3 fig4 fig5 extension atlas genesets tissues tissues_cells mouse gtex fibres rna_gtex rna_rodents passage; default all)
uv run entropy-age stochasticity BETAS SAMPLES [--groups G]
                                          # per-feature drift/entropy metrics for any matrix
```

Dependencies are the scientific Python stack from PyPI (numpy, scipy, polars, scikit-learn,
matplotlib, typer); one from outside that stack is the maintainer's call. `--scale paper` runs at publication sizes
and takes much longer than `quick`.

**Run the commands yourself** rather than telling the person to run them, except where a command
needs an interactive terminal: then hand over the exact line.

---

## 5. Code and tests

- **A result that can be unknown says so.** Yes, no and unknown are three answers; when "no" has
  kinds (a refusal, a timeout, a crash, a malformed input), the answer carries which.
- **Ordering is load-bearing** wherever bytes are hashed, compared or shipped: never emit in the
  order of a set or a map without a stable sort.
- **Aggregate a warning raised in a loop**: by reason, once, with a count.
- **Timestamps are stored as ISO-8601 UTC** and shown in local time.
- **A deprecation warning in code you touched is a blocker**: the API moved; fix it, and the rule
  here, so the pattern does not return.
- **Internals are free to change; the contract is not.** The public API, the shape of what is
  written or sent and the command line change only on purpose, and with a version.
- **Tests run the real code path on real data** and compute the expected values from the fixture
  at run time. Never mock the transformation under test. Assert relationships (set equality, sums,
  round trips) over existence or counts, and cover the edges: empty, null, a boundary, malformed
  input.

---

## 6. Docs and their lifecycle

- **All markdown goes in `docs/`**, except this file and `README.md`. A new file under `docs/` is
  named in UPPER_SNAKE; directories stay lowercase. `docs/paper.md`, `docs/model.md` and
  `docs/extension.md` predate the rule and keep their names, since the README links them (D2).
- **`docs/` is the ground truth.** This file repeats from it only what orients a reader, and every
  prohibition lives here in full: a "don't" behind a link does not get read.
- **`docs/ROADMAP.md` holds open work only.** Done items leave it; the changelog says what landed.
- **`docs/CHANGELOG.md`** records what landed, newest first, in the same commit as the change,
  including a change made here to follow a repository this one consumes.
- **Update the affected docs in the same commit as the change**, never after.
- **Something that belongs to a repository this one consumes but does not own** is reported to
  that repository's inbox (its issues, or its consumer inbox), restated in its terms, never filed
  into its roadmap, and its repository is never edited or committed to. Our side of it goes into our
  changelog.
- **Before a pull request**, show `git diff <upstream>/main --stat HEAD` and
  `git log <upstream>/main..HEAD --oneline` and wait for approval. Files of one installation (local
  settings, local configs) stay out.
- **Prose is plain and honest.** Never describe a capability the repository does not have.

---

## 7. Self-correction

When outdated knowledge of an API causes a real failure, fix the code and the rule here (and the
affected `docs/`), so the next session does not repeat it. When the maintainer corrects a
preference, it goes into §8 in their words, with the reason.

## 8. Learned

Append-only, one line each: the maintainer's preferences in their own words, and the facts of this
machine and this project that reading the code does not reveal.

- "I mostly work alone on my repositories, so push to `main` directly and avoid pull requests."
  (2026-11-09, the maintainer's standing rule for every repository)
- GSE41037's three processing datasets are confounded with age (dataset 1 mostly < 50 y, dataset 3
  mostly > 50 y): any age slope across that cohort needs batch as a covariate (`human_batches`, `dummies`).
- This laptop runs at 90+ °C before any job and shuts down at 100 °C: run figures with `nice -n 19` and
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, one at a time.
- NCBI refused the cluster's downloads with HTTP 403 for a while after about 4 GB in one burst
  (2026-11-09); the proxy itself was fine (its CONNECT succeeded and the 403 came through the TLS tunnel),
  and the block had lifted an hour later. Pace big NCBI downloads, or wait. Relaying from a laptop (`curl URL | ssh host 'cat > file.part'`) works, but a dropped connection ends
  the remote `cat` normally: run with `set -o pipefail`, resume with `curl -r <offset>-`, and only rename
  `.part` after the byte count matches and `gzip -t` passes. A truncated 1.46 of 2.08 GB file was
  promoted once because the rename trusted the remote exit code.

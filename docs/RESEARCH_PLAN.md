# Research plan: from the current results to a paper

Asked by antonkulaga on 2026-11-11 (direct message): a plan with enough novelty to publish, changed as
intermediate results come in, worked autonomously; a blocker goes to `#agent-chat > human decisions` and
the work moves to another arm. The roadmap (`docs/ROADMAP.md`) holds the items in this order; this file
holds why, and the rule that keeps each arm going or stops it. It is rewritten when a result changes the
plan, and the change is said in `docs/CHANGELOG.md`.

## Where the results stand (2026-11-11)

Established on public data (`docs/extension.md`):

1. **Stochasticity is gene-specific.** Essential and housekeeping promoter CpGs drift less than CpGs of the
   same island status and starting level; non-essential genes drift more (human blood, 27k and 450k; mouse,
   11 tissues).
2. **It is tissue-specific.** The protection is larger in blood than in brain, across people and within the
   same donors, and per person (RM6 score: in about 90 % of donors). Mouse cortex and hippocampus have none.
3. **It does not follow turnover.** Mouse heart keeps it, GTEx breast, ovary and testis lose it, and in two
   culture cohorts the gap is a constant fraction of drift, not growing with divisions.

So the paper's candidate claim has moved: from "selection among dividing cells protects essential genes"
to **"gene-specific stochasticity is set by the locus's chromatin, and tissues differ because their
chromatin differs, not because their cells divide."** Novelty then rests on three things: a direct
division test (done, two cohorts), a test that chromatin state accounts for the gap (RM7), and what this
means for clocks (RM8). An intervention test (RM3) is the one result that could still revive selection.

## Changes on results

- 2026-11-11, after RM7's first run: chromatin of the tissue's own epigenome removes 67-91 % of blood's gap, and
  after it blood is no more protected than brain. The chromatin claim holds in blood; arm A's next step is
  brain at cell-type resolution. RM3 stays next: chromatin predicts interventions leave the gap unchanged once
  chromatin is matched, and an intervention that changes chromatin (high-fat diet in liver) should move the gap
  only through it.

- 2026-11-11, same day: RM8 moves ahead of RM3. The one series of RM3 with processed values at hand (GSE199979)
  is liver only, a slowly renewing tissue, and a high-fat diet changes chromatin as well as proliferation, so it
  would not separate the readings; the series that would (GSE80672, blood) needs an alignment job first.
  RM8's comparison with published clocks needs only data already here, and antonkulaga asked for it.

## The arms, in order, each with its rule

### A. RM7: does chromatin state account for the protection? (new, first)

Match set CpGs to background CpGs on chromatin state as well as island and level, using reference
epigenomes of the same tissue (Roadmap/ENCODE chromHMM states and H3K4me3/H3K27ac peaks for blood
and brain, public). **Go on** if the gap shrinks by half or more once chromatin is matched, in blood and
in brain: that is the chromatin claim, and the blood-brain difference should then follow the tissues'
chromatin. **Stop the claim** if the gap stays: then something beyond chromatin (selection among them)
remains, and RM3 moves up. Either way the result is publishable.

### B. RM3: interventions that slow proliferation

Caloric restriction, rapamycin and dwarf mice slow cell turnover. Selection predicts the essential-gene
gap in renewing tissue shrinks; chromatin predicts no change. Start with series that ship processed values
(GSE199979 BXD liver array, GSE89275 liver RRBS); GSE80672 is raw reads only (2.9 GB, an alignment job).
**Go on** while a series has at least 8 animals per arm in a renewing tissue. **Stop** a series that has
fewer, and say so; size the expected effect with `simulate/selection.py` first, so a null is read against
the power it had.

### C. RM8: clocks from protected and unprotected CpGs (new)

The extension's prediction 4: a stochastic clock built from unprotected CpGs and one from protected CpGs
should tick at different rates, and tissues should differ in the ratio. Train both on the paper's simulated
dynamics and on Hannum, test on Jaffe and Lunnon and the mouse tissues. **Go on** if the ratio differs
between blood and brain beyond its bootstrap interval. It gives the work a use: which CpGs a clock should
use, and why clocks disagree across tissues.

### D. RM2: more tissues

GSE184211 (human, 11 tissues, the mouse's array, so a cross-species comparison on the same CpGs), then the
mouse replications. **Go on** with each series that has numeric ages; it feeds A's tissue comparison.

### E. RM4: RNA

A stochastic RNA clock on the rodent atlases; GTEx's same-donor design. Lower priority: RNA so far adds
little beyond methylation.

### F. RM5's remainder

A long arrest (quiescence over months). Search only; run when a series is found.

## What would be the paper

Methods: the stratified contrast, the culture coupling test, the per-person score. Results: gene- and
tissue-specific stochasticity (done), not caused by division (done, two cohorts), accounted for by
chromatin or not (A), unchanged or changed by interventions (B), consequences for clocks (C). A first
draft is written once A and one series of B have results.

## The summary artifact

At the end, and after each major result, the tabbed summary artifact (the claude.ai page kept since v6) is
updated (antonkulaga, 2026-11-11): Meyer and Schumacher's work first, then the hypothesis built on it, the
result of each step that tests it, the main novelty, practical applications, and how the results intersect
with existing omics clocks (which CpGs or genes of the published clocks are protected or unprotected, and
what that predicts for their tissue behaviour; RM8 measures it).

## Working rules for this plan

- Heavy steps run as Slurm jobs through `docs/RUNS.md`; the lead builds, the runner runs.
- A blocker (missing data, a permission, a judgement only the person can make) is posted in
  `#agent-chat > human decisions`, with a default, and the next arm starts.
- Every major milestone (an arm's result, a change of plan) is posted in the repository's channel, `#entropy-age >
  milestones` (antonkulaga, 2026-11-11).
- After each result: write it up in `docs/extension.md`, then reread this plan and change it if the
  result says so.

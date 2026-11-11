# Paper draft: Epigenetic clocks tick where chromatin lets them

Draft started 2026-11-11 at antonkulaga's request. Every claim carries its status:

- **[solid]**: replicated, or robust to the checks we know of;
- **[weak]**: one dataset, a small sample, or a check still open;
- **[missing]**: needed for the claim and not yet done.

Numbers and figures are those of `docs/extension.md`, which holds every source table.

## Abstract

Meyer and Schumacher (2024) showed that an aging clock emerges from nothing but random maintenance errors
at CpGs bounded between 0 and 100 % methylation. They treated that noise as uniform. We find it is not.
Across human blood, brain and nine GTEx tissues, eleven mouse tissues, mouse blood and liver with exact ages,
and two long-term culture studies, CpGs in promoters of essential and housekeeping genes drift 0.2-0.5 SD
less with age than CpGs of the same island status and starting level. Non-essential genes drift more.

We tested three explanations.

- **Selection among dividing cells.** In culture, the shield is a constant share of drift through up to 136
  population doublings; it does not grow, as selection requires.
- **Turnover.** The shield does not follow how often a tissue's cells divide.
- **Chromatin.** The tissue's own chromatin state accounts for most of it: 67-91 % in blood, up to half in brain
  with cell-type-resolved maps. After matching, the remainder is the same in blood and brain.

In living mice the shield's share grows with age in blood and in slowly renewing liver alike, and caloric
restriction weakens it in one within-batch comparison.

Published clocks draw on CpGs whose chromatin differs between tissues, Polycomb targets above all.
Clocks built from Polycomb CpGs alone misread other tissues by a constant, calibratable offset. A
per-person score of the shield separates individuals far beyond its own error.

Chromatin, not cell competition, decides which genes resist epigenetic entropy; that choice shapes what
clocks measure.

## 1 Introduction

- Stochastic epigenetic drift as a sufficient cause of clocks (Meyer & Schumacher 2024; Tarkhov 2024 and others). **[solid, literature]**
- Open question: is drift uniform across genes and tissues, and if not, what sets it? Two candidates: selection
  among cells (a cell that loses an essential gene's state is outcompeted) and chromatin (active, well-maintained
  promoters make fewer errors). **[solid, framing]**
- The two make different predictions: selection needs dividing cells and accumulates with divisions; chromatin
  removes a fixed share of errors and follows the tissue's chromatin. A model makes the predictions explicit
  (virtual experiment). **[solid]**

## 2 Results

### 2.1 Essential genes drift less than matched CpGs

- **Evidence:** human blood (GSE41037, 27k; Hannum, 450k), mouse 11 tissues (GSE279437), GTEx 9 tissues, mouse
  blood RRBS (GSE80672), mouse six tissues (GSE120132), BXD liver (GSE199979). Matching on island status and
  young-reference level; permutation and bootstrap intervals. **[solid]**
- **Weak:** GTEx with decade-bin ages; four of nine tissues show little.
- **Figures:** gene_sets, gtex_tissues, mouse_tissues, mouse_rrbs_tissues.

### 2.2 The shield is a constant share of drift as cells divide

- **Evidence:** two culture studies (GSE197512, seven lines to 136 doublings; GSE179847, five lines to 79
  doublings). The noise floor is removed. **[solid]** for the long courses.
- **Weak:** the short courses are too noisy; arrest (MMC, contact inhibition) leaves too little drift to test
  whether the shield needs divisions.
- **Missing:** an arrest of months with measurable drift.
- **Figures:** passaged_cells, passaged_cells_sturm.

### 2.3 The shield does not follow turnover

- **Evidence:** mouse heart keeps it; GTEx testis, breast and ovary lose it like muscle. Six mouse tissues of the
  same strains order blood among the weakest. **[solid]** as a negative result.
- **Weak:** muscle's turnover in humans is unknown (no exercise data).

### 2.4 The tissue's own chromatin accounts for most of the shield

- **Evidence:** Roadmap chromHMM matching removes 67-91 % in blood (two cohorts, two blood references); another
  tissue's map removes about 30 %. Sorted neuron and glia ATAC removes 30-49 % in brain. After matching, the
  remainder is equal in blood and brain. **[solid]** as an association.
- **Weak:** possible circularity (essential promoters are active by definition); reference epigenomes come from
  other people.
- **Missing:** a causal test (chromatin changed at a gene, drift watched); five-mark maps of sorted neurons.
- **Figures:** chromatin_match.

### 2.5 In living tissue the shield's share grows with age; caloric restriction weakens it

- **Evidence:** mouse blood (GSE80672) and BXD liver both show F becoming more negative with age, robust to
  composition. **[solid]** in two series.
- CR weakens the shield by a quarter in B6D2F1, where treated and controls share a sequencing series.
  **[weak]**: 22 mice, one strain; the C57BL/6 comparison crosses series and is not read.
- **Missing:** a second CR or rapamycin series sequenced with its controls; sorted stem cells.
- **Interpretation open:** growth in liver argues against stem-cell selection as its cause.
- **Figures:** interventions, bxd_liver (being redrawn with corrected intervals).

### 2.6 A per-person shield score

- **Evidence:** reproduces the cohort result per sample; reliability 0.96-0.997 against its own error; blood more
  shielded than the same donor's brain in 86-91 % of donors. **[solid]** as a measure.
- **Missing:** test-retest on replicate arrays; association with outcomes (mortality, frailty) beyond age.
- **Figures:** per_person.

### 2.7 What clocks are made of

- **Evidence:** 16 published clocks are enriched in tissue-variable chromatin and Polycomb CpGs; pan-tissue clocks
  in housekeeping promoters. Clocks from Polycomb CpGs give a reproducible tissue offset. **[solid]** descriptively.
- **Failed:** a clock's tissue offset does not follow its chromatin share. Report it.
- **Missing:** whether a chromatin-informed clock adds information beyond conventional clocks (an outcome
  cohort, or intervention response).
- **Figures:** clock_cpgs, chromatin_clocks.

### 2.8 Across species (in progress)

- **Missing:** does the shield scale with maximum lifespan across 349 mammals (GSE223748, run R12)? A positive
  answer would be the paper's strongest new result; a null would still bound the mechanism.

## 3 Discussion

- Chromatin sets a gene's epigenetic stochasticity; tissues differ because their chromatin differs. Selection may
  add a part in renewing tissue in vivo (the CR lead), but it is not the main cause. **[weak]** on the second
  sentence.
- Consequences for clocks: they measure drift where chromatin lets it happen; clock transfer across tissues is
  a chromatin problem.
- Limits: arrays and RRBS cover promoters unevenly; reference chromatin from other individuals; no causal
  manipulation; species not phylogenetically corrected.

## 4 What the paper still needs, in order

1. **Cross-species result (RM10)**, running. It decides whether the paper has a comparative headline.
2. **Replication of the CR effect** in a second within-batch series.
3. **Test-retest of the per-person score** on replicate arrays (public technical replicates exist on GEO).
4. **A causal chromatin experiment.** This needs a wet lab: epigenome editing at a promoter in long-term culture.
   It is the core of a proposal (below), not of this paper.
5. Corrected blood-against-brain intervals (RM9).

## 5 DFG proposal: what has potential

DFG Sachbeihilfe (individual research grant), 3 years. Working title: *Chromatin-determined epigenetic
stochasticity: mechanism, evolution and use in aging clocks.*

- **WP1, causality (wet lab, the core).** Target active or repressive chromatin to chosen promoters with
  dCas9-p300 or dCas9-KRAB. Do it in fibroblasts passaged for months, or kept quiescent, and measure drift
  against untargeted look-alikes. This is the one experiment that turns our association into mechanism.
  Partner: a lab with long-term culture and epigenome editing.
- **WP2, selection in vivo.** Sort haematopoietic stem cells and mature blood from young, old and CR mice, and
  ask whether the shield's growth with age and its response to CR are stem-cell-intrinsic. It builds on the
  GSE80672 lead.
- **WP3, evolution.** The shield against lifespan across mammals (RM10), with phylogenetic correction. Add new
  long-lived species where gaps exist.
- **WP4, application.** A chromatin-informed stochastic clock and the per-person shield score, validated in a
  German cohort with follow-up (e.g. the Berlin Aging Study II or KORA, by data-access request). Do they predict
  morbidity beyond conventional clocks?
- **What we bring:** a working pipeline across 14 public datasets, the matching method, and preliminary results
  for every WP (sections 2.2, 2.4, 2.5, 2.7).
- **Fit:** DFG reviewers look for a clear hypothesis with a falsifiable test. Selection against chromatin, with
  the culture result already against selection, gives that.
- **Risk:** WP1 depends on a partner lab. WP4 depends on cohort access.

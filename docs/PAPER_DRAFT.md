# Paper draft

**Research programme:** what shields genes from epigenetic aging noise, and can we strengthen the shield?
**Working paper title, from the current data only:** Epigenetic clocks tick where chromatin lets them. It
describes the results so far, not the direction of the research, and changes as the running tests come in
(antonkulaga, 2026-11-11: the research must not be biased by a temporary paper name).

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

In living mice the shield's share grows with age in blood and in slowly renewing liver alike. Caloric
restriction weakens it in one blood series but not in liver. Mice missing one copy of the demethylase Tet3
lose it in brain, against their own littermates. Across 65 mammal species the shield is nearly universal, and not stronger in long-lived ones.

Published clocks draw on CpGs whose chromatin differs between tissues, Polycomb targets above all.
Clocks built from Polycomb CpGs alone misread other tissues by a constant, calibratable offset. A
per-person score of the shield separates individuals far beyond its own error.

Chromatin, not cell competition, decides which genes resist epigenetic entropy; that choice shapes what
clocks measure.

## 0 Position after the novelty review (`docs/NOVELTY_REVIEW.md`)

The qualitative claim (housekeeping and highly connected genes are protected from drift) is published:
Teschendorff, West & Beck 2013 and West et al. 2013. The paper cannot claim the idea. It can claim the design
(matched null, two readouts, many tissues, two species), the quantitative split into chromatin and the rest, and
the dynamics: culture over doublings, the share growing with age in vivo, caloric restriction, species. Three
points of prior work cut against us and must be answered in the text: clonal dominance in culture (Hollmann 2020),
lifespan scaling as stem-cell dynamics (Crofts 2026), and per-person drift burden confounded and non-predictive
(Seeboth 2020).

## 1 Introduction

- Stochastic epigenetic drift as a sufficient cause of clocks (Meyer & Schumacher 2024; Tarkhov 2024; Tong 2024). **[solid, literature]**
- Protection of housekeeping and highly connected genes from drift was proposed and shown qualitatively (Teschendorff, West & Beck 2013; West et al. 2013); age gain at Polycomb targets is established (Teschendorff 2010; Rakyan 2010; Slieker 2016). We test what sets it. **[solid, literature]**
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
  doublings). The noise floor is removed. **[solid]** as a measurement; **[weak]** as evidence against selection:
  long-term cultures become oligoclonal (Hollmann 2020), and steady selection per doubling would also give a
  constant share. **[missing]** clonality data.
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
- CR in liver (GSE223748, 59 against 36 mice) leaves it unchanged (F −0.008, p 0.77): the blood result is not
  replicated in a slowly renewing tissue. **[weak]** either way.
- Tet3 heterozygous cortex and striatum lose the shield against same-chip, same-age littermates (housekeeping F
  +0.17, +0.25; core-essential +0.31, +0.49; 8 against 8). **[weak → solid within one study]**: one study, 16
  mice per tissue. **[missing]** the direction (gain of methylation?) and a second Tet or Dnmt model.
- GHRKO liver's apparent strengthening was batch (gone against same-chip controls).
- People with constitutional chromatin-writer variants (Sotos/NSD1, Weaver/EZH2, Kabuki/KMT2D, CHARGE/CHD7)
  drift more in blood but keep the shield (F −0.04 to +0.003). **[weak]**: GSE74432 not age-adjusted; small
  arms. Contrast with Tet3: the shield follows promoter demethylation, not the writers around it.
- **Missing:** a second CR or rapamycin series sequenced with its controls; sorted stem cells.
- **Interpretation open:** growth in liver argues against stem-cell selection as its cause.
- **Figures:** interventions, bxd_liver.

### 2.6 A per-person shield score

- **Evidence:** reproduces the cohort result per sample; reliability 0.96-0.997 against its own error; blood more
  shielded than the same donor's brain in 86-91 % of donors. **[solid]** as a measure.
- **Missing:** test-retest on replicate arrays; association with outcomes (mortality, frailty) beyond age; an
  answer to Seeboth 2020 (per-person drift burden confounded by cell mix, smoking, BMI and batch).
- **Figures:** per_person.

### 2.7 What clocks are made of

- **Evidence:** 16 published clocks are enriched in tissue-variable chromatin and Polycomb CpGs; pan-tissue clocks
  in housekeeping promoters. Clocks from Polycomb CpGs give a reproducible tissue offset. **[solid]** descriptively.
- **Failed:** a clock's tissue offset does not follow its chromatin share. Report it.
- **Missing:** whether a chromatin-informed clock adds information beyond conventional clocks (an outcome
  cohort, or intervention response).
- **Figures:** clock_cpgs, chromatin_clocks.

### 2.8 Across species: universal, not tuned to lifespan

- **Evidence:** 139 species-tissue groups of 65 mammals (GSE223748): the housekeeping gap is negative in 96 %,
  −0.30 SD on average. **[solid]**
- Against log maximum lifespan (weight, tissue fixed) the gap moves +0.03 per log-unit (p 0.11-0.14): no
  stronger shield in long-lived species, if anything weaker. Consistent with Crofts 2026. Holds with a
  phylogenetic GLS (taxonomy covariance, Pagel's λ; +0.025, p 0.40) and within skin, blood and liver.
  **[solid]** as a negative result; near-0 % essential promoters on the conserved array floor the drift.
- **Anomaly, resolved:** human bone marrow and lymph node reverse (+0.18, +0.20), within donor noise (intervals
  cross 0; 5-7 % of 20-donor draws), and the same donors lift other tissues by +0.16 (R54).
- **Figures:** species_shield.

## 3 Discussion

- Chromatin sets a gene's epigenetic stochasticity; tissues differ because their chromatin differs. Selection may
  add a part in renewing tissue in vivo (the CR lead), but it is not the main cause. **[weak]** on the second
  sentence.
- Consequences for clocks: they measure drift where chromatin lets it happen; clock transfer across tissues is
  a chromatin problem.
- Limits: arrays and RRBS cover promoters unevenly; reference chromatin from other individuals; no causal
  manipulation; species not phylogenetically corrected.

## 4 What the paper still needs, in order

1. **Tet3, second step**: the direction of change at essential CpGs, and a second Tet or Dnmt model (the loss holds
   against littermates, R15).
2. **Replication of the CR effect** in a second within-batch blood series (liver did not replicate it).
3. **Test-retest of the per-person score** on replicate arrays (public technical replicates exist on GEO).
4. **A causal chromatin test**, from public natural experiments first (RM14: chromatin-writer syndromes in blood,
   mouse knockouts, edited cells over time); a wet-lab edit only for what those cannot show.
5. ~~Corrected blood-against-brain intervals (RM9).~~ Done (R19): no conclusion moved.

## 5 DFG proposal: what has potential

DFG Sachbeihilfe (individual research grant), 3 years. Working title: *Chromatin-determined epigenetic
stochasticity: mechanism, evolution and use in aging clocks.*

- **WP1, causality.** First from public data: blood of people with a constitutional variant in a chromatin
  writer (Sotos NSD1, Weaver EZH2, Tatton-Brown-Rahman DNMT3A, Kabuki KMT2D), mouse knockouts (Tet2, Tet3,
  Dnmt3a) and CRISPRoff-edited cells tracked for months (RM14). Then, for what those cannot show, editing at chosen
  promoters (dCas9-p300 or -KRAB) in primary cells, never cell lines (abnormal methylation, no age), confirmed in
  blood, where clocks are read.
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

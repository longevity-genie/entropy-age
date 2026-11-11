# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.
The order and the rule for each item are argued in `docs/RESEARCH_PLAN.md`.

## RM15 — Directed demethylation outside promoters: the culture and cortex follow-ups

**State:** unattended · **Owner:** lead

Progress: first part done (`figures intergenic`, R14, `docs/extension.md`): in blood, solo-WCGW CpGs lose
methylation 2.7 times as fast as other open-sea CpGs (directed, not noise); none in cortex; the PCA does not single
it out. Left: the same test in the culture series (EPIC manifest needed for the tags; the loss should grow with
doublings), and what cortex's two strongly age-tracking components are.
## RM12 — Selection pressure under interventions

**State:** unattended · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): do interventions change the shield through selection, apart from
chromatin? What we have: in B6D2F1 blood CR weakens the shield by a quarter (within one sequencing series, robust
to composition), without slowing drift; the shield's share also grows with age in liver, which rarely renews. So
neither proof nor disproof. Tests that would separate selection from chromatin, in order of cost:
- Done (R13, `docs/extension.md`): CR in liver does not weaken the shield; Tet3 heterozygous mice lose it in
  cortex and striatum (controls pooled across studies, a batch could do it); GHRKO strengthens it in liver.
- Done (R15): against same-chip littermates the Tet3 loss holds and grows (housekeeping F +0.17, +0.25); the GHRKO
  liver effect was batch.
- Done (R16): under Tet3 loss essential CpGs gain methylation, but S rises in other arms too (array compression
  could do it); not specific yet.
- Next: RM14's syndromes in human blood (built: `figures syndromes`).
- **The same intervention in a renewing and a non-renewing tissue of the same animals.** Selection predicts CR
  or rapamycin changes the shield in blood or gut and not in brain or muscle; chromatin predicts the same change
  (or none) everywhere. Search for a multi-tissue CR or rapamycin series with at least 8 animals per arm.
- **Interventions that raise selection pressure.** Chemotherapy, irradiation or transplantation shrink the stem
  cell pool and let clones expand. Selection predicts a stronger shield afterwards, chromatin no change. Search
  for blood methylation before and after such bottlenecks (for example bone-marrow transplant donors and
  recipients, or patients after chemotherapy).
- **Chromatin held fixed.** Match on chromatin first (RM7's method), then test the intervention on the remainder.
  If the CR effect survives the matching, it is not a chromatin change.

## RM11 — Clocks of our own

**State:** after RM14 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): can we build clocks from the entropy-age ideas, and do they capture
more than conventional clocks? Candidates from the results: a stochastic clock from Polycomb CpGs with each
tissue's offset calibrated, and the per-person shield score beside it. "More" needs a readout beyond age: which
clock tracks interventions (CR mice, GSE80672) and, later, outcomes in a cohort with follow-up.

## RM13 — Answer the prior work (novelty review)

**State:** after RM12 · **Owner:** lead

Progress: `docs/NOVELTY_REVIEW.md`. To make the paper new, in order: (1) does essentiality still predict protection
after matching on chromatin, expression level and breadth? First answer from R4's tables (`docs/extension.md`):
in blood, core-essential genes keep almost nothing beyond chromatin, housekeeping genes a small clear remainder,
brain both; left: a DepMap gradient within broadly expressed genes, and expression level as a matching variable.
(2) absolute drift slopes beside the CR share; (3) clonality in the culture series where reads allow (epiallele
structure), against Hollmann 2020; (4) read the species result against Crofts 2026 (stem-cell dynamics) and the
floor of near-0 % essential promoters; (5) the shield score against Seeboth 2020's confounders.

## RM14 — Natural experiments in chromatin: broken writers, knockouts, edits

**State:** unattended · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): cell lines tell little about clocks; can public data show what
happens when chromatin is changed? Three sources, in order:
1. Blood methylation of people with a constitutional variant in a chromatin writer: Sotos (NSD1), Weaver (EZH2),
   Tatton-Brown-Rahman (DNMT3A), Kabuki (KMT2D). Martin-Herranz et al. 2019 found clocks accelerated in Sotos.
   Test: is the essential-gene shield weaker than in matched controls? Search GEO for the cohorts.
2. Mouse knockouts: Tet3 heterozygous brain (GSE223748, at hand after R12); Tet2 and Dnmt3a mutant blood.
3. Edited cells over time: CRISPRoff (Nunez et al. 2021) tracked imposed promoter methylation for months.
Rule: a cohort with at least 10 carriers and matched controls on one platform.
- Done: Tet3 heterozygous brain loses the shield against littermates (R15, R16; see RM12).
- Done (R18): Sotos, NSD1 variants, Weaver, Kabuki and CHARGE carriers drift more in blood but keep the shield
  (F unchanged or more negative); KMT2D and CHD7 variants show nothing within the controls' ages.
- Next: the DNMT3A loss-against-gain pair below; a TET2 or TET3 human cohort if one exists.
- Found, below the rule: DNMT3A loss of function (TBRS) GSE128801, 5 carriers against 7 controls; DNMT3A gain
  of function (microcephalic dwarfism) GSE120428, 4 against 5 in blood. Together a loss-against-gain contrast of
  one writer; read only as a pair, after R17.

## RM3 — Interventions: replicate the CR result

**State:** after RM2 · **Owner:** lead

Progress: GSE80672 done (`figures interventions`, `docs/extension.md`): in living mouse blood the shield's share
grows with age, and CR weakens it by a quarter in the one within-batch comparison (B6D2F1, 12 against 10 mice,
robust to composition), without slowing drift. Left, to replicate it:
- **GSE89275** (liver RRBS: Ames dwarf, rapamycin, CR) and its WGBS companions GSE89273/89274: about 4 mice
  per arm, below the plan's 8; not used.
- Another blood series with CR and controls sequenced together (to find).
- GSE199979 (BXD liver) done (`figures bxd`): the shield's share grows with age in liver too, so the growth needs
  no renewing tissue; HFD's effect is the cell mix. It has no CR arm.
- Left: a second CR (or rapamycin) series with treated and control animals sequenced together, at least 8 per
  arm, in blood or another tissue (to find).
## RM2 — A turnover gradient across many tissues (methylation)

**State:** building · **Owner:** lead

Progress: mouse (GSE279437), human GTEx (GSE213478) and muscle fibres (GSE304045) done; see
`docs/extension.md`. In mouse the protection is missing in cortex and hippocampus but kept in heart;
in GTEx muscle is weakly protected, as are breast, ovary and testis, so the split does not follow
turnover. Left:

- **GSE184211** (human, 11 tissues, mammalian array, exact ages 26–91, 190 MB). Not paired by donor, but
  it shares its CpGs with the mouse array, so the comparison across species is direct.
- **GSE50498** (450k muscle, ages 18–89): a muscle cohort with exact ages and more donors than GTEx's 47.
  Like GTEx it says nothing of exercise, and satellite-cell divisions follow loading, injury and
  wasting, so muscle tests selection only in a cohort that records physical activity or training
  history beside age (to find; the fibre cohort GSE304045 records training but not age).
- **Mouse replication:** Thompson 2018 (GSE120132) done (`figures rrbs_tissues`): the shield in all six tissues,
  ordered by neither size nor growth by turnover. Left: Meer 2018 (GSE121141; liver, lung, heart, cortex, 6-30 months).
- The fibre contrast (GSE304045) was uninformative: a dozen adults of similar age show too little
  between-person divergence. A within-tissue turnover contrast needs aged and young people.

## RM4 — Stochasticity in RNA

**State:** after RM2 · **Owner:** lead

Progress: the turnover test is done in GTEx v10 (`figures rna_gtex`, 47 tissues) and in the mouse and rat
atlases (`figures rna_rodents`); see `docs/extension.md`. Mouse repeats the methylation result (brain
without protection, heart with it), rat leans the same way, GTEx shows no turnover effect in spread.
Left:
- **Health at death (antonkulaga, 2026-11-11):** some GTEx donors died fast, others after long illness. Hardy is a
  covariate already; rerun GTEx RNA and methylation on fast deaths only (Hardy 1-2) and say what moves.

- **A stochastic RNA clock**, trained on simulated transcriptional noise as the paper does for
  methylation, tested on the rodent atlases (exact ages).
- **GTEx with more power:** the same-donor design (every tissue of one donor against that donor's
  blood), and subsamples larger than 70 donors in the tissues that allow them.
- Tabula Muris Senis bulk **GSE132040**, a third rodent cohort.
- Interventions: **GSE280382** (rapamycin and exenatide) and **GSE131754** (eight liver interventions).

## RM5 — Direct selection test: protection against division count in culture

**State:** after RM4 · **Owner:** lead

Progress: GSE197512 and its replication GSE179847 done (`figures passage`, `figures passage_sturm`,
`docs/extension.md`): in both the protection is a roughly constant fraction of drift and does not grow
with doublings; the MMC arms (25 days) and contact inhibition (120 days) carry too little drift to decide
whether protection needs divisions. Left:

- An arrest with more drift than GSE179847's contact inhibition gives (quiescence over a year, or senescent
  cells kept for months), so that a gap without divisions can be measured above the noise floor.
- GSE268145 (BJ, 175 MB, downloaded), GSE202554 (astrocytes), GSE131280 (glucose): smaller replications.

## RM9 — Unbiased intervals for the blood-against-brain bootstrap

**State:** unattended (R19 queued) · **Owner:** lead

Progress: `figures tissues` (and `tissues_cells`) scale each resample's gap by its own spread across CpGs, which
resampling inflates, so its resampled gaps and blood-minus-brain intervals shrink toward zero (conservative, not
anti-conservative). Recompute them with `common.scaled_intervals`, rerun, and say in `docs/extension.md` what moved.


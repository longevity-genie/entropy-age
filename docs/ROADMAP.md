# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.
The order and the rule for each item are argued in `docs/RESEARCH_PLAN.md`.

## RM10 — Shielding in short- against long-lived species

**State:** building · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11, #entropy-age): do long-lived species (human, naked mole rat, bats,
whales) shield essential genes more than short-lived ones (mouse)? The mammalian methylation array measures the
same CpGs in hundreds of species (Mammalian Methylation Consortium series on GEO); AnAge lifespans are in the
repo (`taxonomy_anage_meta`). Plan: find the series with species, tissue and age labels; per species and tissue,
the shield (aging gap against matched CpGs); regress it on maximum lifespan corrected for body mass, tissue held
fixed. Rule: go on with every species-tissue with at least 10 animals across ages.

## RM12 — Selection pressure under interventions

**State:** urgent · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): do interventions change the shield through selection, apart from
chromatin? What we have: in B6D2F1 blood CR weakens the shield by a quarter (within one sequencing series, robust
to composition), without slowing drift; the shield's share also grows with age in liver, which rarely renews. So
neither proof nor disproof. Tests that would separate selection from chromatin, in order of cost:
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

**State:** after RM10 · **Owner:** lead

Progress: asked by antonkulaga (2026-11-11): can we build clocks from the entropy-age ideas, and do they capture
more than conventional clocks? Candidates from the results: a stochastic clock from Polycomb CpGs with each
tissue's offset calibrated, and the per-person shield score beside it. "More" needs a readout beyond age: which
clock tracks interventions (CR mice, GSE80672) and, later, outcomes in a cohort with follow-up.

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
- **Mouse replication:** Thompson 2018 (GSE120132; six tissues, two strains, 2-21 months, RRBS) built as
  `figures rrbs_tissues`, waiting on its run (R10); Meer 2018 (GSE121141; liver, lung, heart, cortex, 6-30 months).
- The fibre contrast (GSE304045) was uninformative: a dozen adults of similar age show too little
  between-person divergence. A within-tissue turnover contrast needs aged and young people.

## RM4 — Stochasticity in RNA

**State:** after RM2 · **Owner:** lead

Progress: the turnover test is done in GTEx v10 (`figures rna_gtex`, 47 tissues) and in the mouse and rat
atlases (`figures rna_rodents`); see `docs/extension.md`. Mouse repeats the methylation result (brain
without protection, heart with it), rat leans the same way, GTEx shows no turnover effect in spread.
Left:

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

**State:** after RM2 · **Owner:** lead

Progress: `figures tissues` (and `tissues_cells`) scale each resample's gap by its own spread across CpGs, which
resampling inflates, so its resampled gaps and blood-minus-brain intervals shrink toward zero (conservative, not
anti-conservative). Recompute them with `common.scaled_intervals`, rerun, and say in `docs/extension.md` what moved.


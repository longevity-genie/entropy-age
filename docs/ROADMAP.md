# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.
The order and the rule for each item are argued in `docs/RESEARCH_PLAN.md`.

## RM7 — Does chromatin state account for the protection?

**State:** unattended · **Owner:** lead

Progress: planned (`docs/RESEARCH_PLAN.md`, arm A). Match set CpGs to background on chromatin state (reference
epigenomes of blood and brain: chromHMM states, H3K4me3 and H3K27ac peaks) as well as island and level, in
Hannum, Jaffe and Lunnon. Go on if the gap shrinks by half or more; if it stays, RM3 moves up.

## RM3 — Long- vs short-lived rodent groups

**State:** unattended · **Owner:** lead

Progress: data found, nothing fetched. Interventions that extend lifespan by slowing proliferation (caloric
restriction, rapamycin) make a prediction specific to selection. Total drift should fall, but the gap
between essential and other genes in renewing tissue should shrink, because less turnover means less
selection. Chromatin protection alone predicts no change in the gap. Size the effect with the virtual
experiment (`simulate/selection.py`) before reading the data.

- **GSE279437** (+ RNA GSE280382), tried: exenatide changes nothing; rapamycin has only 5 mice per
  tissue, and against size-matched controls its renewing-tissue values are within noise. Too small.
- **GSE80672** (Petkovich 2017, blood RRBS, 2.9 GB raw only). GHRKO, Snell dwarf and CR vs controls.
  Ages are in months, although the field says years.
- **GSE89275** (liver RRBS). Ames dwarf, rapamycin and CR, plus WGBS GSE89273/89274.
- **GSE199979** (BXD liver, mammalian array, 45 strains, chow vs high-fat diet). Strain lifespans come
  from GeneNetwork (Roy et al. 2021), not from GEO.
- **Not usable as found:** the Mammalian Methylation Consortium's mouse intervention samples (GSE223748
  and related series) carry no group labels on GEO; the labels would have to come from EnsembleAge
  Table S1. No public Diversity Outbred or ITP multi-tissue methylation was found.

## RM8 — Clocks from protected and unprotected CpGs

**State:** after RM3 · **Owner:** lead

Progress: planned (`docs/RESEARCH_PLAN.md`, arm C). A stochastic clock from unprotected CpGs and one from
protected CpGs; the ratio of their rates in blood, brain and the mouse tissues, with bootstrap intervals.

## RM2 — A turnover gradient across many tissues (methylation)

**State:** unattended · **Owner:** lead

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
- **Mouse replication:** Meer 2018 (GSE121141; liver and lung vs heart and cortex, 6–30 months, RRBS)
  and Thompson 2018 (GSE120132; blood, liver vs muscle), which have numeric ages.
- The fibre contrast (GSE304045) was uninformative: a dozen adults of similar age show too little
  between-person divergence. A within-tissue turnover contrast needs aged and young people.

## RM4 — Stochasticity in RNA

**State:** unattended · **Owner:** lead

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

**State:** unattended · **Owner:** lead

Progress: GSE197512 and its replication GSE179847 done (`figures passage`, `figures passage_sturm`,
`docs/extension.md`): in both the protection is a roughly constant fraction of drift and does not grow
with doublings; the MMC arms (25 days) and contact inhibition (120 days) carry too little drift to decide
whether protection needs divisions. Left:

- An arrest with more drift than GSE179847's contact inhibition gives (quiescence over a year, or senescent
  cells kept for months), so that a gap without divisions can be measured above the noise floor.
- GSE268145 (BJ, 175 MB, downloaded), GSE202554 (astrocytes), GSE131280 (glucose): smaller replications.

# Roadmap

Open work only, in order. Each item is a section `## RMn — name` with a line saying its state and its
owner. A done item leaves this file; `docs/CHANGELOG.md` says what landed.

## RM2 — A turnover gradient across many tissues (methylation)

State: the mouse part is done (`figures mouse`, GSE279437; see `docs/extension.md`). There, protection
is a step rather than a slope: as strong in slow tissues (liver, kidney, adipose) as in renewing ones,
and missing in cortex and hippocampus. Next, the human tissues and a contrast of turnover inside one
tissue:

- **GSE213478** (human GTEx EPIC, 987 samples, 424 donors, betas 5.6 GB). Colon, lung, blood, muscle,
  kidney and others with donor IDs. Ages are in decade bins only.
- **GSE184211** (human, 11 tissues, mammalian array, exact ages 26–91, 190 MB). Not paired by donor, but
  it shares its CpGs with the mouse array, so the comparison across species is direct.
- **Muscle:** GSE304045 (type I and type II fibres vs whole muscle, EPIC v2) and GSE50498 (450k, ages
  18–89). Muscle fibres do not divide, but their satellite cells do: the within-tissue contrast that can
  separate turnover from tissue identity.
- **Mouse replication:** Meer 2018 (GSE121141; liver and lung vs heart and cortex, 6–30 months, RRBS)
  and Thompson 2018 (GSE120132; blood, liver vs muscle), which have numeric ages.

## RM3 — Long- vs short-lived rodent groups

State: data found, nothing fetched. Interventions that extend lifespan by slowing proliferation (caloric
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

## RM4 — Stochasticity in RNA

State: data found. GTEx v10 is already parsed by the rna-age project on the cluster (19,461 samples,
946 donors, 68 tissues, decade age bins). Read its outputs read-only. Anything that belongs to rna-age
goes to its inbox, never into its repository.

Measure each gene's between-donor dispersion and its change with age, compared against genes of the
same expression level (RNA noise depends strongly on the mean). Then repeat the essential vs
non-essential and turnover tests, with the same composition adjustment, and finally build a stochastic
RNA clock.

Data:
- GTEx: the same donors' tissues, as in Lunnon, but at ages 20–79.
- Mouse and rat atlases **GSE314817 / GSE314970**: 24 and 28 tissues, 6–27 months, raw counts about
  78 MB each. Rat muscles are in GSE226118.
- Tabula Muris Senis bulk **GSE132040**.
- Interventions: **GSE280382** (rapamycin and exenatide) and **GSE131754** (eight liver interventions).

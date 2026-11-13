# Array preprocessing: current practice and where this project stands (roadmap RM18, item 5)

Asked by antonkulaga (2026-11-12): are results comparable across arrays and array versions, and do we follow best
practice? This note sets out the steps that published practice for Illumina methylation arrays (450k, EPIC) agrees on.
For each, it says what this project does and what our own tests showed (R20, R36, R61). Sources were checked on
PubMed (2026-11-13); DOIs are given.

## The steps practice agrees on

| step | what it fixes | reference | this project |
|---|---|---|---|
| Background correction from out-of-band signal (noob) | non-specific fluorescence differing between arrays | Triche et al. 2013, [10.1093/nar/gkt090](https://doi.org/10.1093/nar/gkt090) | deposited betas as the authors made them; noob in the reprocessed cohort (R61) |
| Detection p-values from out-of-band signal (pOOBAH) and masking | failed hybridisation, deletions, polymorphic sites reported as intermediate methylation | Zhou et al. 2018 (SeSAMe), [10.1093/nar/gky691](https://doi.org/10.1093/nar/gky691) | in the reprocessed cohort (R61); deposited sets carry the authors' own QC |
| Masking cross-reactive, SNP-affected and channel-switching probes | probes that measure something else | Zhou et al. 2017, [10.1093/nar/gkw967](https://doi.org/10.1093/nar/gkw967) | the general mask applied as a test (R20): no gap moved by more than 0.01 SD |
| Correcting the two probe designs (type I against II) | different dynamic range by design | Marabita et al. 2013, [10.4161/epi.24008](https://doi.org/10.4161/epi.24008) (BMIQ and quantile normalisation best in their comparison) | design added to the matching strata (R20): about a sixth of blood's gaps, none in cortex |
| Between-array normalisation with control probes (functional normalisation) | batch and technical variation across arrays | Fortin et al. 2014, [10.1186/s13059-014-0503-2](https://doi.org/10.1186/s13059-014-0503-2) | not applied; every comparison stays within one dataset, with plate or chip as a covariate |
| Positional bias (chamber, row) | systematic intensity differences by position on the chip | Butler et al. 2025, [10.1371/journal.pone.0326337](https://doi.org/10.1371/journal.pone.0326337) | slide and row as blocks in SJLIFE's score tests (R60): slide explains 20-36 % of a score, well above chance |
| Reliability of clock-type scores | technical noise moving clocks by up to 9 years between replicates | Higgins-Chen et al. 2022, [10.1038/s43587-022-00248-2](https://doi.org/10.1038/s43587-022-00248-2) | per-person reliability of the shield score 0.997 (C1); clocks are read as published |

## What our design adds, and what it leaves exposed

* **Matching cancels most preprocessing.** Every gap compares a gene set with look-alike CpGs of the same island
  status and starting level, in the same samples. A preprocessing step that shifts CpGs of one kind and level
  together shifts both sides equally. That is why Jaffe's cortex, rebuilt from its IDATs, kept every
  essential-set gap below zero (R61), although single CpGs' slopes agreed only at r 0.56-0.82.
* **Comparisons stay within a dataset and platform.** We never pool 450k with EPIC or mix two studies' betas in one
  test, so between-platform normalisation is not needed. Comparing a gap's size across studies is a different
  matter: a size in SD carries about ±0.05 SD of preprocessing uncertainty (R61).
* **Exposed: effects near that size.** The intervention and modifier effects in the shield share (CR, R55-R56;
  Tet enzymes, R35; turnover, R63) are 0.02-0.1. Such effects should be confirmed on reprocessed data before they
  are claimed; R61's pipeline (`scripts/reprocess_idats.py`) is ready for that.
* **Exposed: the signed shift S.** S reads a direction, which compression at the extremes can fake. R36 holds
  each sample's compression fixed (Tet3's gain survives); every other S result should do the same.

## Recommendation

1. Keep deposited betas for screening, and keep every test within one dataset.
2. Rerun from IDATs any result that is small (under 0.1 in F or 0.1 SD) or that a claim rests on: next SATSA
   (E-MTAB-7309, RM28) and Gill 2022's fibroblasts (GSE165179, RM17).
3. Report each headline gap with its deposited and reprocessed value where both exist.

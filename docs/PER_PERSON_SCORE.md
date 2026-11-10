# Design note: a per-person protection score

State: design only. No outcome data is needed to build it; validating it against outcomes is later work.

## Why the current measure is not per person

Every protection result so far is a property of a cohort. For each CpG (or gene) we fit the age slope
of its entropy or of its between-person spread, then compare the slopes of a gene set with those of
background features at the same starting level (`compare_groups`). A slope needs many people of
different ages, so one sample has no value of it. A per-person score has to say how far *this*
sample has drifted at protected features compared with matched unprotected features in *the same*
sample.

## The score

Fix a reference: young samples of the same tissue on the same platform (for blood, people aged 20 to
35 of the cohort or of a public reference). From it, take per feature `j`:

* `m_j`, the reference mean (beta value, or log CPM for RNA);
* `s_j`, the reference spread: the binomial `sqrt(m_j (1 - m_j))` for methylation, the reference SD
  for RNA (RNA spread depends on the mean, and the reference SD carries that dependence);
* the stratum of `j`: CpG-island status x decile of `m_j` (methylation), or decile of `m_j` (RNA),
  exactly as in the cohort tests, so the strata do not depend on the sample being scored.

For a sample `i` and feature `j`, the deviation is `d_ij = |x_ij - m_j| / s_j`, after the sample's
technical covariates have been removed with coefficients fitted on the reference (see below). The
score of gene set `G` in sample `i` is the matched difference used in the cohort tests, computed
over features instead of over slopes:

```
P_i(G) = sum_k w_k * ( mean_{j in G, stratum k} d_ij  -  mean_{j not in any set, stratum k} d_ij )
```

with `w_k` the share of `G`'s features in stratum `k`. A negative `P_i` means the sample's protected
features sit closer to the young state than unprotected features of the same starting level. An
entropy version replaces `d_ij` by `H(x_ij) - H(m_j)`, the binary entropy gain, for methylation.

Two properties make it the per-person form of what we already measure:

* **It averages to the cohort measure.** The mean of `d_ij` over people of age `t` is the expected
  spread at `t`, so the age slope of the cohort mean of `P_i(G)` is the matched difference of the
  dispersion slopes, up to the SD scaling. The first check of any implementation is that regressing
  `P_i` on age in Hannum, Lunnon and GTEx reproduces the sign and rough size of the cohort results in
  `docs/extension.md`.
* **Sample-wide noise cancels to first order.** A sample with poor bisulfite conversion, low RIN or
  a noisy array has larger `d_ij` everywhere; subtracting matched background features of the same
  sample removes that shared part. What remains is the difference between set and background,
  which is what selection is about.

## What can bias it, and the counter for each

* **Cell composition.** A sample richer in one cell type moves thousands of cell-type-specific
  features together. Remove composition before computing `d_ij`: either reference-based cell
  fractions (blood), or the composition components of background features fitted on the reference
  and projected into the sample (`composition_components`), as the cohort tests do with their
  age-free components.
* **Genetics.** Part of a person's deviation at a CpG is a methylation QTL, which no age process
  made. It is roughly the same for set and background features, but a set with more mQTL sites
  would look less protected. Drop known mQTL CpGs (GoDMC list) from both sides, or score only
  features without common SNPs in the probe.
* **Reference mismatch.** A reference from another lab or array version shifts every `d_ij`. Use a
  reference processed with the same normalisation; report the score only within one platform.
* **Few features per stratum.** Core-essential genes give about 2,000 promoter CpGs on the 450k
  array; a single sample's `P_i` then has a sampling error that can be estimated by resampling
  genes (not CpGs: CpGs of one promoter are correlated) within strata. Report `P_i` with that
  interval.

## How it would be validated, in order

1. **Reproduce the cohort result** (above): `P_i` regressed on age, with covariates, gives the cohort
   sign in blood and none in cortex.
2. **Repeatability:** technical replicates and repeated draws of the same person (longitudinal
   cohorts) give a within-person SD well below the between-person SD at the same age.
3. **Turnover inside one person:** in Lunnon's same-donor tissues, `P_i` is more negative in blood
   than in the same person's cortex, for most people, not only on average.
4. **Division count, directly:** in serially passaged fibroblasts (GSE197512, GSE179847; roadmap RM2),
   `P_i` against population doublings, and against days in culture in arrested cultures. Selection
   predicts that `P_i` falls with doublings in dividing cultures and stays flat in arrested ones.
5. **Outcomes, last:** only after 1 to 4, ask whether `P_i` at a given age predicts mortality or
   morbidity in a cohort with follow-up. A score that fails 1 to 4 would be measuring something
   else, and an outcome association would then be uninterpretable.

## What it is not

It is not an age estimate. It is a ratio-like contrast between feature groups in one sample, and it
should be reported next to a clock's age, never instead of it. Nor is it evidence of selection by
itself: chromatin protection of essential genes would also give a negative `P_i`. Only its
dependence on turnover (steps 3 and 4) bears on selection.

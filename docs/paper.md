# The paper, and how well this reimplementation reproduces it

> Meyer D.H. & Schumacher B. **Aging clocks based on accumulating stochastic variation.**
> *Nature Aging* 4, 871–885 (2024). [doi:10.1038/s43587-024-00619-x](https://doi.org/10.1038/s43587-024-00619-x)
> · PMID 38724736 · open access [PMC11186771](https://pmc.ncbi.nlm.nih.gov/articles/PMC11186771/)
> · original code [Meyer-DH/StochasticAgingClock](https://github.com/Meyer-DH/StochasticAgingClock) (MIT)

## The claim in one paragraph

Aging clocks are often read as evidence of an aging *program*. The authors show that no program is
needed. If any set of measurements starts from a ground state and is bounded (methylation fractions,
normalised expression), then **adding random noise is enough to build a clock**. The bound turns
symmetric noise into a directional drift: values near 0 can only go up and values near 1 can only
go down. The amount of accumulated noise is therefore readable from the data. In a single-cell
methylation model, imperfect maintenance (the probabilities $E_m$ and $E_d$) produces the same effect.
Published clocks (Horvath, PhenoAge, Lin, Vidal-Bralo, Weidner, GrimAge) track simulated noise.
A clock trained *only on simulated noise* from one young sample predicts the chronological age of
real people and of other mammals, and it detects interventions (calorie restriction, GHRKO,
smoking, parabiosis, partial reprogramming). The level of stochasticity sets the pace of the clock:
more noise or less maintenance means faster predicted aging.

## Key results, with this repository's reproduction

Numbers in the right column come from `uv run entropy-age figures` at `--scale quick`, seed 42
(tables in `data/output/`). "Quick" uses smaller sweeps and a cheaper elastic-net search than the
paper; `--scale paper` uses the published sizes.

| Paper result (figure) | Paper | This repo | Figure |
|---|---|---|---|
| Unbounded noise: no clock (ED 1c) | r = −0.05 | constant prediction (r undefined) | `fig1_bulk.png` a |
| Noise clipped to [0, 1] (ED 1d) | r = 0.99 | r = 0.98, R² = 0.96 | `fig1_bulk.png` b |
| Logit-scale noise, σ = 0.2 (1b) | r = 0.95 | r = 0.91, R² = 0.81 | `fig1_bulk.png` c |
| Accuracy rises with σ, then plateaus (1c) | qualitative | R² 0.0 → 0.85 (σ 0.01 → 0.2) | `fig1_bulk.png` d |
| Coefficients follow regression to the mean (1e) | qualitative | reproduced | `fig1_bulk.png` f |
| More noise per step means faster predicted aging (1g) | qualitative | reproduced | `fig1_bulk.png` h |
| Universal $E_m$ = 99.9 % single-cell model (3c) | R² = 0.999 | R² = 0.998 | `fig3_single_cell.png` b |
| $E_m$ = 99.995 % (3b) | R² = 0.78 | R² = 0.75 | `fig3_single_cell.png` a |
| Lower maintenance means faster predicted aging (3e) | qualitative | reproduced | `fig3_single_cell.png` d |
| Empirical $E_m/E_d$ from the 88-year-old (3f) | r = 0.81 | r = 0.89 | `fig3_single_cell.png` e |
| Start exactly at equilibrium 0.5 → no clock; 0.51 → clock (ED 4) | r = 0.05 / 0.95 | r = 0.03 / 0.77 | `fig3_single_cell.png` f |
| Horvath, best rate limits (4a) | 97 % / 5 % | 97 % / 5 % (R² = 0.79) | `fig4_published_clocks.png` a |
| Horvath on simulated data, limits 97/5 (4b) | r = 0.91 | r = 0.92 | `fig4_published_clocks.png` b |
| Horvath, universal $E_m$ = 99 % (4c) | r = 0.97 | r = 0.96 | `fig4_published_clocks.png` c |
| PhenoAge, universal 99 % (4f) | r = 0.94 | r = 0.92 | `fig4_published_clocks.png` g |
| Lin / Vidal-Bralo / Weidner (ED 7) | all significant | r = 0.97 / 0.97 / 0.97 | `fig4_published_clocks.png` d |
| Ground state from 16 / 81-year-old (ED 6) | r = 0.89 / 0.09 | r = 0.93 / 0.33 | `fig4_published_clocks.png` h |
| **Simulation-trained clock on 392 real blood samples (5a)** | r = 0.87 | **r = 0.88** | `fig5_real_data.png` a |
| … after cell-composition correction (5a) | r = 0.87 | r = 0.84 | `fig5_real_data.png` b |
| Random CpG sets work, permuted ages do not (ED 8) | qualitative | r 0.44–0.62 vs ≈ 0 | `fig5_real_data.png` c |
| Pan-mammalian clock, empirical rates (clock 3) | on average significant | dolphins r = 0.88 | `fig5_real_data.png` e |
| Pan-mammalian clock, universal 99 % (clock 4) | on average significant | dolphins r = 0.06 (see below) | `fig5_real_data.png` f |

Not reproduced, because the inputs are not public in usable form: the *C. elegans* BiT-age
transcriptome figures (Fig. 2; the count matrix is not distributed, only its sample sheet), the
Mammalian Methylation Consortium intervention data (GHRKO, Tet3, CR), the smoking and parabiosis
cohorts, the reprogramming time course (Fig. 6) and GrimAge (computed on a web server). The
transcriptome helpers are ported (`transcriptome.py`) for anyone who processes counts the same way.

### Where we differ, and why

* **Clock 4 on dolphins.** A universal $E_m = 99\,\%$ ($z = 0.5$ everywhere) moves the simulated
  dolphin methylome about 4.4× further over 66 steps than real dolphins change between 0.6 and 58 years
  (mean |Δβ| 0.25 vs 0.057). Every real sample therefore maps onto the first few simulated steps.
  Matching the magnitude ($E_m = 99.9\,\%$) does not rescue it (r = 0.19): the direction is wrong,
  because real sites relax toward *site-specific* equilibria, not toward 0.5. The paper reports the
  universal clocks averaged over many species. On this example dataset only the empirical,
  site-specific version works, which supports treating stochasticity per site (our extension).
* **Equilibrium start at 0.51** gives r = 0.77 rather than 0.95. We use $E_m = 99.9\,\%$ as the text
  states; one upstream variant used 99 %, which relaxes faster and gives a stronger signal.
* **Cell-composition correction.** The original regression correction was run in R over all
  samples. We regress each CpG on the EpiDISH fractions shipped upstream across the 394 controls, so
  the corrected r (0.84) is close to, but not identical with, the paper's.

## Mapping the original script to this package

The original is one 3,653-line script that runs top to bottom and needs files that are not in its
repository (`GSE41037.csv`, worm counts). Each piece now lives in a module:

| Upstream | Here |
|---|---|
| `random_epi_nolimit`, `random_epi`, `random_epi_logit*` | `simulate.bulk.simulate_bulk(bound=...)` |
| `simulate_cells_for_age_fixed`, `update_cells_fast_fixed` | `simulate_methylation(..., MaintenanceRates.universal(n, em))` |
| `simulate_for_age_empirical_noquantile`, `update_cells_fast_empirical_noquantile` | `simulate_methylation(..., rates)` (binomial, exact) |
| `get_noise_Em_all_new` | `MaintenanceRates.from_equilibrium` |
| `get_noise_Em_all_random` | `MaintenanceRates.random` |
| `get_clock`, `anti_transform_age` | `published.LinearClock`, `clock.HorvathAge` |
| `pred_and_plot`, `get_prediction` | `clock.train_clock`, `clock.evaluate`, `plots.prediction_scatter` |
| `-log(-log(age / maxage))`, `exp(-exp(-pred))` | `clock.LogLogRelativeAge` |
| min-max ×400 −120 rescaling | `clock.LinearRescale` |
| R `lm(x ~ B + NK + ...)` residual correction | `stochasticity.regress_out` |
| `make_binary`, `correct_Bio_Age` | `transcriptome.binarize_by_median`, `transcriptome.correct_biological_age` |
| `get_noise_func_parm`, `apply_biol_noise_ground` (lognormal noise fits, `fitter`) | not ported: used only for an exploratory supplementary panel |

## Upstream issues fixed along the way

* `sklearn.utils.testing` no longer exists, so the script does not import on current scikit-learn.
* No random seeds: results change on every run (the README says so). Every panel here uses its own
  seeded `numpy.random.Generator` stream.
* The pan-mammalian section uses an undefined `maxage` (it should be 67), has a typo `und_sites`, and
  reads inputs from `plot_path`.
* `get_noise_Em_all_new` divides by zero when the old sample's value is exactly 0, which makes $E_m$ NaN.
* In the supplementary feature-size sweep (Sup. 2b), the `R2` column is computed against the
  *permuted* ages, a copy-paste slip.
* `random_epi_logit_everystep` is the same model as `random_epi_logit` (see [model.md](model.md)).
* Every age is re-simulated from scratch, and every cell is flipped individually. Both are now exact
  but much cheaper (binomial counts and an age-ordered schedule).
* The initial cells are laid out so that cell 0 is methylated at every site with $g > 0$. That does not
  matter while sites are independent, but it would bias any model that couples sites within a cell, so
  the selection model draws a random arrangement.
* Dropping every CpG with a single missing value (as `dropna()` does) removes 13 % of Horvath's CpGs
  from the GEO version of GSE41037. We drop only CpGs missing in more than 5 % of samples and
  median-impute the rest.

## Data

| Dataset | Source | Used for |
|---|---|---|
| GSE41037: 720 whole-blood methylomes, Illumina 27k (394 controls, 16–88 y) | NCBI GEO, `entropy-age fetch geo GSE41037` | human ground state (16 y), equilibria (88 y), validation |
| EpiDISH blood cell fractions for GSE41037 | upstream `Data/` | cell-composition correction |
| 50 bottlenose dolphin blood methylomes, mammalian array (37,554 CpGs, 0.6–57.9 y) | upstream `Data/` (Mammalian Methylation Consortium example) | pan-mammalian clocks 3/4 |
| Horvath, PhenoAge, Lin, Vidal-Bralo, Weidner, BiT-age coefficients | upstream `Data/` | published clocks |
| Supplementary Table 1 (993 worm RNA-seq samples, metadata only) | upstream `Data/` | kept for reference |

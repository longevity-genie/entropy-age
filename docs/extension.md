# Extension: is stochasticity gene- and cell-type-specific because of selection?

## Hypothesis

The paper treats stochastic variation as one knob per dataset (σ, or one $E_m$) or as a per-site
rate inferred from one old sample. It then argues that clock CpGs drift fastest because they matter
least for the cell. We want to turn that remark into a testable model:

> **Stochasticity differs between genes and between cell types, because selection among cells can
> protect fitness-relevant loci from entropy, and selection can only act where cells turn over.**

Testable predictions:

1. **Gene level.** In tissues with cell turnover, loci whose state matters for fitness
   (essential and housekeeping genes, cell-identity genes) drift *less* than baseline-matched
   loci. Their bulk level is held near the ground state at a mutation–selection balance.
2. **Cell-type level.** The protection disappears in post-mitotic cells (neurons, cardiomyocytes,
   most muscle). There, epimutations accumulate at the intrinsic maintenance rate wherever they fall.
3. **Hitchhiking.** Where selection acts, neutral loci gain *between-individual* variance as
   fitter clones expand (clonal drift; compare clonal haematopoiesis). Dispersion at neutral sites
   rises faster in renewing tissues than in post-mitotic ones.
4. **Clock pace.** A clock trained on unprotected (post-mitotic-like) dynamics reads protected
   tissues as younger. The tissue's "pace of aging" partly reflects how much of its methylome is
   under selection.
5. **Replication errors cut the other way.** Each division is another chance for a maintenance
   error (the mitotic-clock idea), so high turnover also *raises* drift at unprotected sites. Turnover
   is not uniformly protective.

## Model

`entropy_age.simulate.selection` adds two things to the paper's single-cell Markov model:
cell turnover and fitness-proportional replacement. The equations are in
[model.md §4](model.md#4-extension-selection-among-cells). Knobs:

| Knob | Meaning | Biological reading |
|---|---|---|
| `MaintenanceRates(em, ed)` per site | intrinsic stochasticity | DNMT1 fidelity, chromatin context, replication timing |
| `SelectionModel.weights` per site | fitness relevance | essential / identity genes > 0, passengers = 0 |
| `strength` | fitness cost per deviating site | how deleterious an epimutation is |
| `turnover` | fraction of cells replaced per step | blood, gut, skin high; neurons ≈ 0 |
| `division_error` | extra epimutation in daughters | replication-coupled maintenance errors |
| `reference="founder"` or `"target"` | what counts as deviation | change from the lineage state, or from a fixed functional state |

## Virtual experiment (known ground truth)

`uv run entropy-age figures extension` takes 400 CpGs from the 16-year-old blood sample (bimodal
levels <0.1 or >0.9) and marks 100 of them as constrained. All sites share $E_m = 99\,\%$, so every
site has the same *intrinsic* stochasticity. The methylome is then aged in four cell types:
PM (post-mitotic), RN (renewing, no selection), RS (renewing with selection), and RS+E (RS plus
replication errors).

![Selection extension](figures/extension_selection.png)

| Observation (panel) | PM | RN | RS | RS+E |
|---|---|---|---|---|
| Fitted change over life, constrained / neutral (b) | 0.36 / 0.37 | 0.36 / 0.36 | **0.08 / 0.38** | **0.12 / 0.44** |
| Cliff's δ constrained vs neutral, permutation p | −0.15, p = 0.09 | ≈ 0 | **−1.0, p ≤ 0.001** | **−1.0, p ≤ 0.001** |
| Between-individual SD of neutral sites at step 80 (c) | 0.03 | 0.09 | **0.33** | **0.32** |
| Slope of a post-mitotic-trained clock (d) | 0.88 | 0.95 | **0.70** | 0.79 |

All five predictions appear in the simulation. Two lessons matter for real data:

* **Protection does not look like a slower rate.** Constrained sites under selection reach their
  mutation–selection plateau *quickly*, so their fitted relaxation rate λ is *higher* (median 0.072 vs
  0.022 per step) while they move much less. Rank genes by the amount of change (`window_change`,
  entropy slope, dispersion slope), and read λ only together with the fitted equilibrium.
* **Neutral turnover alone already adds dispersion** (RN: 0.09 vs PM: 0.03); selection makes it far
  larger. Dispersion slopes distinguish post-mitotic from renewing tissues even without knowing
  which genes are constrained.

The parameters (selection strength, turnover per step, 200 cells per sample) are illustrative, not
calibrated. The figure shows that the effects exist and that the estimators recover them from
cross-sectional bulk data, which is all real cohorts offer.

## First look at real data: are clock CpGs the fast-drifting ones?

`uv run entropy-age figures atlas` computes every metric for all 27,577 CpGs of the 394 GSE41037
controls (16–88 y) and for the 37,554 CpGs of the 50 dolphins. It then compares clock CpGs with
CpGs *that start at the same level* (deciles of the young-adult baseline). For the human cohort,
slopes are estimated within processing batches (see the confound below). The λ histogram comes
from `fit_relaxation`, which takes no covariates, so it is unadjusted.

![Atlas](figures/atlas_real_data.png)

Human blood, clock CpGs vs level-matched CpGs, batch-adjusted slopes (2,000 stratified permutations,
p not corrected for the six tests):

| Group | Metric | Matched difference | p |
|---|---|---|---|
| Horvath (353) | entropy slope | +1.6 × 10⁻⁵ bits / year | 0.60 |
| PhenoAge (513) | entropy slope | +5.0 × 10⁻⁵ bits / year | 0.07 |
| Horvath | dispersion slope | +1.3 × 10⁻⁵ | 0.58 |
| PhenoAge | dispersion slope | +4.8 × 10⁻⁵ | 0.01 |

**Not supported.** An earlier version of this page reported that clock CpGs gain entropy faster
than level-matched CpGs (p = 0.004 for Horvath). That came from unadjusted slopes across a cohort
whose processing batches are confounded with age. With batch as a covariate the difference
disappears for Horvath's CpGs and is weak for PhenoAge's. On this data, clock CpGs are typical
drifting CpGs, not the fastest ones.

Dolphins have no clock-CpG contrast. The human clock probes present on the mammalian array do not
hybridise in dolphin DNA: they sit at about 0.5 with almost no variance, so any comparison would
measure probe failure, not biology.

**Caveat on matching.** The baseline used for matching should be measured independently of the
drift being compared: a separate young or fetal reference, or a sorted cell type at its youngest.
A baseline taken from the cohort already contains early drift. Fast-drifting CpGs then fall into
other strata than slow ones, and matching removes part of the effect.
`compare_groups` reports `n_matched` so this shows up, and a test in the suite documents it.

## Prediction 1 on real data: essential genes drift less in blood

`uv run entropy-age fetch annotation` downloads the Illumina 27k probe manifest (GEO GPL8490,
which gives the gene, Entrez ID, RefSeq accession and CpG-island flag per probe) and three
curated gene sets:

* Hart et al. core-essential genes (CEGv2) and non-essential genes (NEGv1) from CRISPR screens,
  pinned at a commit of hart-lab/bagel and matched by Entrez ID;
* Eisenberg & Levanon (2013) housekeeping genes, matched by RefSeq accession, then by symbol.

`uv run entropy-age figures genesets` then compares each set with every other CpG in the same
**stratum**: CpG-island status (fixed by sequence, so drift cannot move it) crossed with the decile
of young-adult methylation. Without matching the sets are not comparable: 93 % of core-essential
CpGs sit in islands near 0 % methylation, but only 33 % of non-essential CpGs do (panel a).

![Gene sets](figures/gene_sets.png)

| Set | CpGs (genes) on the array |
|---|---|
| core essential | 1,140 (603) |
| housekeeping | 5,825 (3,087) |
| non-essential | 1,063 (605) |

**A confound we found and removed.** GSE41037 was run in three processing datasets whose ages
barely overlap: dataset 1 is mostly under 50 and dataset 3 mostly over 50. A slope across the whole
cohort partly measures that batch; the clock-CpG result above is one casualty. Every gene-set
contrast was therefore repeated with batch as a covariate (`drift_metrics(..., covariates=...)` gives
partial slopes), with batch plus the seven blood-cell fractions, and separately inside each batch.

Set minus matched CpGs, in SD of the metric (CpG level; * = BH-adjusted p < 0.05 over all 72 tests,
1,000 stratified permutations each, so the smallest p is 0.001):

| Set · metric | unadjusted | batch-adjusted | batch + cells | inside batch 1 / 2 / 3 |
|---|---|---|---|---|
| core essential · entropy slope | −0.11 * | −0.25 * | −0.32 * | −0.23 / −0.31 / −0.16, all * |
| housekeeping · entropy slope | −0.21 * | −0.30 * | −0.38 * | −0.30 / −0.34 / −0.20, all * |
| non-essential · entropy slope | +0.17 * | +0.23 * | +0.27 * | +0.21 / +0.29 / +0.07, all * |
| core essential · dispersion slope | −0.20 * | −0.11 * | −0.16 * | −0.12 / −0.10 / −0.06, all * |
| housekeeping · dispersion slope | −0.23 * | −0.10 * | −0.17 * | −0.15 / −0.09 / −0.04, all * |
| non-essential · dispersion slope | +0.26 * | +0.18 * | +0.23 * | +0.17 / +0.12 / +0.20, all * |
| core essential · \|age slope\| | −0.04 * | −0.16 * | −0.14 * | |
| housekeeping · \|age slope\| | −0.04 * | −0.18 * | −0.17 * | |
| non-essential · \|age slope\| | +0.02 | +0.05 | +0.13 * | |

Gene-level aggregation (median per gene, so genes with many probes do not dominate) gives the
same signs and similar sizes (`data/output/gene_sets/contrasts.csv`).

**Reading.** In whole blood, a renewing tissue, CpGs of essential and housekeeping genes gain
entropy and between-person dispersion more slowly than CpGs that start at the same level in the
same kind of sequence. CpGs of genes no cell line needs drift faster. This holds after removing the
batch confound and cell composition, and it replicates inside each processing batch. The effects
are modest (0.1–0.4 SD), as one would expect if most drift is neutral.

**What this does not show yet.** It is the gene-level half of the hypothesis (prediction 1). The
same pattern would arise without any selection if active, open promoters are simply better
maintained: H3K4me3 and transcription keep DNMT3 away whatever the cell's fitness. The
discriminating test is prediction 2. If selection does the protecting, the essential-gene advantage
should shrink or vanish in post-mitotic cells such as neurons, where no cell can be replaced by a
fitter neighbour. If chromatin does it, the advantage should stay the same size. That needs a
cohort with a post-mitotic tissue next to a renewing one; the next section tests it.

## Prediction 2 on real data: the protection is larger in blood than in brain

The essential-gene advantage is larger in blood than in brain, in two independent designs. That is
what selection among cells predicts. The same-donor contrast survives an adjustment of every tissue
for its cell mix; the lifespan-slope contrast does so only partly (see
[the composition check](#composition-check)).

**Data** (Illumina 450k, run on a cluster: `entropy-age fetch cohort hannum jaffe lunnon`, then
`entropy-age figures tissues`):

| Design | Tissue | People | Ages | Measure, adjusted for |
|---|---|---|---|---|
| Lifespan cohorts | blood, Hannum et al. 2013 (GSE40279) | 656 | 19–101 | age slopes; plate, site, sex, ancestry |
| Lifespan cohorts | prefrontal cortex, Jaffe et al. 2016 (GSE74193), adult controls | 213 | 21–85 | age slopes; plate, sex, ancestry, neuron fraction |
| Same donors | Lunnon et al. 2014 (GSE59685): blood and four brain regions | 71–80 per tissue | 70–105 | between-person dispersion; age, sex, Alzheimer's status |

Each tissue uses about 109,000 autosomal promoter CpGs: every CpG in a promoter of a gene in the three
sets (6,891 core-essential, 42,173 housekeeping, 4,945 non-essential), plus a fixed random background
of 60,000 promoter CpGs for matching. As in prediction 1, each set is compared with CpGs of the same
island status and starting-level decile, in SD of the metric. Blood minus brain gets a 95 % interval
from 200 bootstraps: independent resamples in the two lifespan cohorts, the same resampled donors in
all five Lunnon tissues.

![Tissue contrast](figures/tissue_contrast.png)

**Set minus matched CpGs, by tissue** (SD; all differ from zero at BH-adjusted p < 0.05 except
where marked):

| Set | Blood (Hannum), entropy slope | Cortex (Jaffe), entropy slope | Blood (Lunnon), dispersion | Four brain regions (Lunnon), dispersion |
|---|---|---|---|---|
| core essential | −0.23 | −0.12 | −0.08 | −0.04 to −0.06 |
| housekeeping | −0.35 | −0.21 | −0.15 | −0.04 to −0.07 |
| non-essential | +0.31 | +0.08 | +0.20 | +0.01 (ns, entorhinal) to +0.11 |

**Blood minus brain** (SD, 95 % bootstrap interval). Negative for the protected sets means more
protection in blood:

| Set | Lifespan: entropy slope | Lifespan: dispersion slope | Same donors: dispersion, four regions |
|---|---|---|---|
| core essential | −0.10 (−0.16, −0.03) | −0.09 (−0.15, −0.02) | −0.02 to −0.04, every interval below 0 |
| housekeeping | −0.14 (−0.21, −0.05) | −0.20 (−0.27, −0.11) | −0.07 to −0.10, every interval below 0 |
| non-essential | +0.21 (+0.10, +0.30) | +0.21 (+0.14, +0.29) | +0.09 to +0.19, every interval above 0 |

Of the 21 blood-minus-brain contrasts, 20 have intervals that exclude zero in the predicted
direction. The exception is the non-essential |age slope|, whose interval crosses zero. The
protection exists in brain too, at a third to two thirds of its size in blood. That also fits the
model: neurons cannot be replaced, but glia turn over, and chromatin protection acts in every cell. Full tables:
`data/output/tissue_contrast/`.

### Composition check

**Why it was needed.** Blood is a mix of cell types whose proportions differ a lot between people.
Genes switched on in only some cell types, as many non-essential genes are, gain between-person spread
and age trends from that mix; essential genes are unmethylated in every cell type, so the mix hardly
moves them. Neither blood cohort ships cell counts, and Lunnon's brain regions ship no neuron
fraction, so every tissue gets the same reference-free correction: the leading principal components
of its 60,000 background CpGs (`composition_components`), after the technical covariates are
regressed out, added as covariates (`entropy-age figures tissues_cells`). Four variants: none, 5 or
10 components that keep their age trend (they also absorb any change of the cell mix with age), and
5 components made orthogonal to age (they absorb only the person-to-person spread of the mix). The
components explain 24–33 % (5) or 30–43 % (10) of the background variance; the first alone 10–13 %.

![Composition check](figures/tissue_contrast_cells.png)

**Contrasts whose 95 % interval excludes zero in the predicted direction:**

| Adjustment | Same donors (Lunnon), of 12 | Lifespan cohorts (Hannum vs Jaffe), of 9 |
|---|---|---|
| none | 12 | 8 |
| 5 components, age-free | 12 | 9 |
| 10 components | 12 | 5 |
| 5 components | 12 | 1 (2 reversed) |

**Same donors: the gap is not the cell mix.** Blood minus brain moves by at most a fifth. Housekeeping:
−0.07 to −0.10 SD unadjusted, −0.06 to −0.08 with 10 components; core essential −0.02 to −0.04
against −0.01 to −0.03; non-essential +0.09 to +0.19 against +0.10 to +0.12, every interval on the
predicted side in every variant. One shift is informative: in the brain regions, the non-essential
excess *grows* after adjustment (frontal cortex +0.05 to +0.12, entorhinal +0.01 to +0.11), so the
neuron–glia mix was hiding part of brain's own effect.

**Lifespan slopes: not separable from coordinated aging.** With age-free components the gap is
unchanged (housekeeping entropy slope −0.14, interval −0.17 to −0.11), but that is by construction,
not a test: components orthogonal to age cannot change any age slope of entropy or of the mean
(corrected 2026-11-10; see the GTEx section). With components that keep
their age trend it shrinks or flips: blood's housekeeping entropy protection falls from −0.35 to
−0.08 SD, and cortex's to zero with 10 components. Those components are partly age itself (one of
Jaffe's correlates r = 0.60 with age), so they remove the coordinated methylation change that thousands
of CpGs share with age along with any age trend of the cell mix. From cross-sectional slopes the two
cannot be told apart. The lifespan contrast is therefore supporting evidence only; the same-donor
dispersion contrast carries the claim.

**What it means.** Essential-gene protection against methylation drift is 1.5–3 times stronger in
blood than in the brain of the same people (housekeeping 2–3 times, core essential about 1.5 times), and the difference is not blood's cell mix. Selection
among cells, which needs turnover, predicts exactly this; chromatin protection alone does not. Tables:
`data/output/tissue_contrast_cells/`.

**Other limits.** The two lifespan cohorts are different people from different studies, though both
on the 450k array. Lunnon's donors are old and mostly have Alzheimer's disease. Brain tissue mixes
post-mitotic neurons with glia that do turn over, so the contrast is renewing against mostly
post-mitotic tissue, not a clean neuron sample. Sorted neurons, or a tissue such as skeletal muscle,
would sharpen it (roadmap RM2).

## Mouse tissues: the protection is missing in cortex and hippocampus

A second species, eleven tissues of the same animals and a second array. In Huang et al.'s mice
(GSE279437, mammalian methylation array, inbred C57BL/6 males), CpGs of essential and housekeeping
genes gain entropy with age far more slowly than matched CpGs in tissues whose cells divide. In
cortex and hippocampus they do not; heart keeps as much protection as white blood cells, and
skeletal muscle about half of it (see the correction at the end of this section).

**Data and method** (`entropy-age fetch mouse-annotation`, `entropy-age fetch cohort huang`,
`entropy-age figures mouse`). Young controls (7–21 per tissue) against aged controls (4–14); GEO gives
the two age groups but no numeric ages. The human gene sets reach mouse through MGI orthologs (642
core-essential, 3,220 housekeeping and 931 non-essential genes; 2,334, 11,994 and 3,305 promoter CpGs)
and are compared, as before, with background promoter CpGs of the same island status and
young-control level decile. Intervals come from 200 bootstraps that resample animals within their
treatment group, the same animals in every tissue. Tissues are grouped by turnover: renewing (white
blood cells, colon, spleen), slow (liver, kidney, adipose) and post-mitotic (skeletal muscle, heart,
frontal cortex, hippocampus, hypothalamus).

![Mouse tissues](figures/mouse_tissues.png)

**Entropy gain with age, set minus matched CpGs** (SD):

| Tissue | Core essential | Housekeeping | Non-essential |
|---|---|---|---|
| white blood cells | −0.14 | −0.12 | −0.00 |
| colon | −0.26 | −0.31 | +0.23 |
| spleen | −0.22 | −0.28 | +0.12 |
| liver | −0.19 | −0.23 | +0.01 |
| kidney | −0.27 | −0.28 | +0.02 |
| adipose | −0.20 | −0.22 | +0.21 |
| skeletal muscle | −0.06 | −0.06 | +0.12 |
| heart | −0.12 | −0.12 | +0.02 |
| frontal cortex | +0.10 | +0.10 | −0.03 |
| hippocampus | −0.02 | −0.03 | −0.01 |
| hypothalamus (4 aged controls) | −0.30 | −0.28 | +0.15 |

**Renewing minus post-mitotic tissues** (SD, 95 % bootstrap interval):

| Set | Entropy gain | Change in spread | Absolute change of the mean |
|---|---|---|---|
| core essential | −0.12 (−0.20, −0.05) | −0.11 (−0.21, −0.01) | −0.06 (−0.10, −0.03) |
| housekeeping | −0.16 (−0.23, −0.09) | −0.13 (−0.23, −0.03) | −0.09 (−0.12, −0.05) |
| non-essential | +0.07 (+0.01, +0.13) | +0.03 (−0.03, +0.14) | +0.11 (+0.05, +0.18) |

Eight of the nine aging contrasts exclude zero in the predicted direction. What this adds to the
human results:

* **It replicates in another species and on another array**, with the same donors in every tissue
  and no blood-specific cell mix: the renewing tissues here include colon and spleen.
* **It is a step, not a slope.** The slow tissues are as protected as the renewing ones (entropy:
  slow minus post-mitotic −0.14 for core essential and −0.16 for housekeeping, intervals excluding
  zero). Hepatocytes, kidney tubule cells and adipocytes divide rarely but can, and their tissues
  replace cells over months. The data fit "some turnover is enough" better than a graded dose of
  turnover. Whether that is selection, or something else that sets post-mitotic tissue apart, needs a
  test that varies turnover inside one tissue.
* **The brain is where protection vanishes.** Cortex and hippocampus show none (+0.10 and −0.03);
  muscle and heart keep a little (−0.06, −0.12). The hypothalamus is the exception, with strong
  protection, but it has only 4 aged controls and the widest intervals.
* **The spread among aged mice does not separate the tissue classes** (renewing minus post-mitotic
  −0.03 and −0.04, intervals crossing zero). These mice are inbred, so the between-animal spread has
  no genetic part, and with 40 aged animals per tissue spread across five treatment groups the
  measure is weaker than the young-versus-aged contrast.

**Interventions: not interpretable yet.** Exenatide (13 aged mice per tissue against 13 controls)
changes no set's dispersion difference. Rapamycin's numbers appear to weaken the protection, which
selection would predict, but only 5 mice per tissue carry them, and a dispersion from 5 animals is
noisier than one from 13. A size check settles part of it (a one-off run, not yet part of `figures mouse`). Thirty random 5-mouse
subsets of the aged controls already lose about a fifth to a quarter of the protection (renewing
tissue, housekeeping: −0.15 with all controls, −0.11 with 5). Against those size-matched controls,
rapamycin's renewing-tissue values (−0.07 for both protected sets) sit inside the range of the
5-mouse subsets (−0.16 to −0.05). Only in post-mitotic tissue does rapamycin fall outside it
(housekeeping −0.04 against −0.09 to −0.06), the opposite of what selection predicts. With five
animals and only white blood cells and colon as renewing tissues, this does not test the prediction
either way. The rapamycin arm needs a larger cohort (roadmap RM3).

**Correction (2026-11-10).** An earlier version of this section said the protection is "missing where
cells do not divide". The per-tissue table does not carry that. Heart, whose cardiomyocytes barely
divide, shows −0.12 SD for both protected sets, as much as white blood cells (−0.14, −0.12). Skeletal
muscle shows −0.06, about half. Only cortex (+0.10) and hippocampus (−0.03) are without it, and the
hypothalamus, also neural, is strongly protected with 4 aged controls. "Neurons lack the protection"
is supported by two of three brain regions; "selection needs dividing cells" is not established,
because the non-neural post-mitotic tissues keep part or all of it. The class contrast above pools
heart and muscle with the brain regions, so its post-mitotic side is mostly carried by cortex and
hippocampus.

Tables: `data/output/mouse_tissues/`.

## Human tissues: GTEx and muscle fibres

### GTEx: muscle is weakly protected, but so are breast, ovary and testis

The test for the human side of the mouse result: in nine GTEx tissues of 424 post-mortem donors on
the EPIC array (Oliva et al. 2023, GSE213478), is skeletal muscle, the tissue with the fewest divisions, the
one without protection? (Below: GTEx cannot say how many divisions its donors' muscle had.) Methylation entropy and between-donor spread are fitted against the decade
bin midpoint (25 to 75), adjusted for sex, autolysis score and smoking. Muscle has 47 donors and no
one under 30, and fewer donors pull every SD-scaled difference toward zero, so every tissue is
resampled at 32 donors without replacement (200 resamples) and tissues are compared only on those
size-matched resamples (`entropy-age figures gtex`).

![GTEx tissues](figures/gtex_tissues.png)

**Entropy gain with age, set minus matched CpGs** (SD, median of size-matched resamples):

| Tissue | Class | Donors | Core essential | Housekeeping | Non-essential |
|---|---|---|---|---|---|
| whole blood | renewing | 54 | −0.20 | −0.27 | +0.27 |
| colon | renewing | 224 | −0.13 | −0.20 | +0.30 |
| testis | renewing | 50 | −0.02 | −0.03 | −0.13 |
| lung | slow | 223 | −0.12 | −0.18 | +0.20 |
| kidney | slow | 50 | −0.13 | −0.20 | +0.22 |
| breast | slow | 52 | −0.03 | −0.06 | +0.07 |
| prostate | slow | 123 | −0.13 | −0.19 | +0.16 |
| ovary | slow | 164 | −0.03 | −0.05 | +0.15 |
| skeletal muscle | post-mitotic | 47 | −0.03 | −0.06 | +0.09 |

**Class contrasts** (SD, 95 % interval over the resamples; negative = more protection than muscle):

| Set | Measure | Renewing − muscle | Slow − muscle |
|---|---|---|---|
| core essential | entropy | −0.08 (−0.14, −0.02) | −0.06 (−0.11, +0.00) |
| housekeeping | entropy | −0.10 (−0.19, −0.02) | −0.08 (−0.15, −0.01) |
| core essential | spread | −0.06 (−0.11, +0.01) | −0.04 (−0.08, +0.01) |
| housekeeping | spread | −0.10 (−0.17, −0.01) | −0.06 (−0.11, +0.01) |
| core essential | spread, + 5 age-free components | −0.01 (−0.05, +0.02) | −0.01 (−0.03, +0.02) |
| housekeeping | spread, + 5 age-free components | −0.03 (−0.06, +0.01) | −0.02 (−0.04, +0.02) |

What this supports, and what it does not:

* **Muscle is among the least protected tissues** (−0.03 core essential, −0.06 housekeeping), and the
  pooled renewing and slow tissues are more protected than muscle on entropy, three of four intervals
  excluding zero. That is the direction of the mouse result.
* **But muscle is not alone, and the split does not follow turnover.** Breast, ovary and testis are as
  weakly protected as muscle; blood, colon, lung, kidney and prostate are strongly protected. Testis
  is the fastest-renewing tissue in the set. Of the single tissues, only blood and kidney differ from
  muscle with intervals excluding zero (housekeeping −0.22 and −0.14). A tissue-level property other than
  turnover separates the two groups; which one (cell mix, hormonal cycling in breast and ovary, the
  germline's own methylation programme in testis) this design cannot say.
* **Muscle is not a clean low-turnover control in GTEx.** Fibre nuclei never divide, but satellite
  cells do whenever muscle is loaded, injured or wasting, so how many divisions (and how much room
  for selection) a donor's muscle had depends on exercise, training history and illness. GTEx records
  none of these; its only related field is the Hardy scale, which says how the donor died, not how
  they lived. The muscle contrasts therefore compare protected tissues with a tissue of unknown
  turnover, and cannot test selection. The breast, ovary and testis point does not rest on muscle
  and stands.
* **Between-donor spread loses the gap once composition is adjusted.** With five composition
  components of background CpGs, every spread contrast moves to within 0.03 SD of zero. In blood
  the housekeeping spread gap halves (−0.26 to −0.13). In GTEx, unlike Lunnon's same-donor tissues
  (composition check above), the spread contrast between tissues is not robust to the cell mix.
* **Entropy and mean-slope results cannot be checked with age-free components at all.** The
  components are made orthogonal to age and the other covariates, so adding them leaves every age
  slope of the mean or of entropy unchanged in the full data, by construction. The two versions
  differ only through resampling noise. Only residual-based measures (between-donor spread) can move.
  This also applies to the lifespan line of the composition check above, corrected there.

Tables: `data/output/gtex_tissues/` (`tissues`, `point_estimates` on all donors, `contrasts` over
size-matched resamples).

### Muscle fibres: uninformative

Pooled type I fibres, pooled type II fibres and whole muscle of the same 12 adults (11 after
training), before and after 12 weeks of interval training, on EPIC v2 (GSE304045,
`entropy-age figures fibres`). Fibre nuclei never divide; whole muscle also holds satellite cells,
fibroblasts and endothelium, which do. With no ages, the measure is between-person spread at each
CpG, adjusted for sex; all three preparations come from the same people, so sample size cannot bias
the comparison. Unlike GTEx this cohort does record exercise (the training itself), but not age. Selection predicts more protection (a negative difference) in whole muscle than in
fibres.

| Set | Whole muscle − type I, before | − type II, before | − type I, 12 weeks | − type II, 12 weeks |
|---|---|---|---|---|
| core essential | −0.015 (−0.025, −0.001) | −0.019 (−0.039, −0.003) | +0.007 (−0.011, +0.025) | +0.013 (−0.003, +0.029) |
| housekeeping | +0.004 (−0.004, +0.012) | −0.006 (−0.015, +0.004) | +0.001 (−0.015, +0.016) | +0.007 (−0.010, +0.022) |
| non-essential | −0.032 (−0.069, −0.002) | −0.004 (−0.043, +0.025) | −0.030 (−0.063, +0.006) | −0.015 (−0.055, +0.041) |

Every difference is 0.03 SD or smaller and the sign does not hold from before to after training.
The measure itself has no power here: even within pure fibres the housekeeping set is only −0.01 to
−0.02 SD below matched CpGs. A dozen adults of unknown, similar age show too little between-person
divergence for a set contrast. This result is uninformative, not negative.

Tables: `data/output/muscle_fibres/`.

### What the human tissues add

The neuron result from mouse is the one still standing: cortex and hippocampus without protection in
mouse, and in humans (Lunnon, same donors) less protection in four brain regions than in blood.
"Protection needs dividing cells" is not established in either species: mouse heart keeps it, and
in GTEx three tissues with very different turnover lose it as much as muscle does. Human muscle
cannot settle it either way until a cohort records both age and physical activity. The RNA test in
GTEx's thirteen brain regions and the direct test in dividing versus arrested
cultures are the ones that can separate turnover from tissue identity.

## The direct test: does the protection grow with divisions?

The tissue comparisons cannot separate turnover from tissue identity. Cells in culture can: the same
parent culture can be made to divide or be stopped. Endicott et al. (2022, GSE197512) passaged seven
primary human lines to replicative senescence on EPIC (fetal, neonatal and adult skin fibroblasts,
neonatal keratinocytes, adult endothelial and smooth muscle cells; 14-39 samples each, 16-83
population doublings, PD), and in three lines arrested growth with mitomycin C (MMC) next to
dividing controls; one line also ran under 3 % oxygen, in reduced serum and with TERT for 136 PD
(`entropy-age fetch cohort endicott`, `entropy-age figures passage`).

**The two readings.** Chromatin protection makes every division add less error at essential
promoters than at matched CpGs, in a fixed proportion: the gap grows with drift, but as a constant
fraction of it. Selection purges cells that damaged an essential gene, and the purge accumulates
with divisions: the fraction itself should grow with PD. Under arrest there is no selection at all.

**Measure** (the per-sample score of `docs/PER_PERSON_SCORE.md`). In each experiment the reference
is the mean of its earliest samples. Per later sample and promoter CpG, the squared deviation from
the reference, minus the array noise floor (the per-CpG variance of replicate cultures measured on
the same day), scaled by `m (1 - m)`; then, as in every test here, set CpGs minus matched CpGs of the
same island status and reference decile. `C` is the matched drift, `Q` the set-minus-matched gap,
`G = Q / C` the gap as a fraction of drift. The main statistic regresses `Q` on `C` and on
`C x PD` per arm (`drift_coupling`): the first coefficient is the protection per unit of drift, the
second how it changes per 10 PD (negative = protection grows with divisions).

**Why the noise floor matters.** Without it (absolute deviations, the first version), the
housekeeping fraction fell with PD in all seven lines (−0.057 per 10 PD, interval −0.078 to −0.035)
and arrested cultures looked less protected than dividing ones in all three. That is what array
noise alone produces: the noise is the same for set and matched CpGs, so while drift is small the
fraction is pulled toward zero, and it approaches its true value as drift grows. The test suite
reproduces this artefact without any selection.

![Passaged cells](figures/passaged_cells.png)

**Result, with the floor removed** (housekeeping; core essential in brackets):

| | Coupling (protection per unit drift) | Change per 10 PD |
|---|---|---|
| mean of the 7 passaged lines | −0.14 [−0.25] | −0.14, interval −0.33 to +0.06 [−0.07, −0.20 to +0.06]; negative in 4 of 7 |
| neonatal foreskin fibroblast 1 (83 PD) | −0.68 (SE 0.03) | +0.010 (SE 0.003) |
| fetal skin fibroblast (45 PD) | −0.60 (0.06) | +0.022 (0.013) |
| fetal skin fibroblast with TERT (136 PD) | −0.73 (0.02) | −0.008 (0.001) |
| neonatal foreskin fibroblast 2, 3 % oxygen | −0.34 (0.17) | −0.031 (0.033) |

* **The protection is real in culture.** Once drift is well above the floor, the housekeeping gap
  settles at −0.5 to −0.7 of matched drift in every line (panel b): set CpGs drift about half as
  much as matched CpGs, in squared units.
* **It does not grow with divisions.** The fraction is flat from about 25 PD to the end in the long
  courses; the two longest untransformed lines even weaken slightly with PD, and the 136-PD TERT
  course strengthens by 0.008 per 10 PD, about 0.1 over its whole span. Across the seven lines the
  change is as often positive as negative. A protection that is a fixed fraction of drift is what
  chromatin protection predicts; selection predicts a growing fraction, and this test does not
  show one.
* **The arrest arms are uninformative.** MMC stopped most drift (panel d: matched drift 0.0001 to
  0.0016 against 0.003 to 0.014 in dividing controls over the same 25 days), which confirms that
  drift in culture comes with division, but leaves too little drift to measure a gap.

**Limits.** Culture is not a tissue: cells are selected for growth on plastic, and Endicott's
point is that methylation loss in culture is concentrated in late-replicating, CpG-poor domains,
while promoter CpGs replicate early; matching on island status and level removes part of that
difference, not replication timing itself. The noise floor comes from replicate cultures measured on
the same day, which also carries what they diverged since splitting (in the smooth-muscle line,
measured at day 20, the corrected drift dips below zero at first); one line without replicates
uses the floor pooled over the others. Seven lines, one of each type except fibroblasts, are few.

Tables: `data/output/passaged_cells/` (`drift_coupling`, `coupling_across_lines`, `across_lines`,
`mmc_arrest`; per-sample scores in `sample_scores`).

### Replication: the Cellular Lifespan Study (GSE179847)

Sturm et al. (2022, GSE179847) grew dermal fibroblasts of six healthy donors and three SURF1 patients
for up to 265 days on EPIC, with glucocorticoid, mitochondrial and metabolic treatments, 3 % oxygen, and
contact inhibition for 120-140 days (`entropy-age fetch cohort sturm`, `entropy-age figures passage_sturm`;
paper scale, seed 42). The baseline is each healthy line's untreated 21 % oxygen course of the two
lifespan parts of the study (five lines, 5-15 samples, 24-79 PD). No arm measures one day twice; the
noise floor comes from the flasks split from one culture at a treatment's start, and the baseline
courses, which have none, use the floor pooled over the experiments that do.

![Passaged cells, GSE179847](figures/passaged_cells_sturm.png)

| | Housekeeping | Core essential |
|---|---|---|
| change of the coupling per 10 PD, mean of 5 lines | −0.055, interval −0.21 to +0.10; negative in 3 of 5 | −0.067, −0.29 to +0.16; 4 of 5 |
| mean coupling | −0.27 | −0.19 |
| GSE197512, for comparison (7 lines) | −0.14, −0.33 to +0.06 | −0.07, −0.20 to +0.06 |

* **The same answer.** The gap grows with drift in all five lines (`Q` per 10 PD −0.012, interval
  −0.022 to −0.001), but its share of drift does not move with divisions: per line the change is
  −0.03 (SE 0.01), +0.07 (0.02), −0.06 (0.03), +0.005 (0.02) and −0.27 (0.04), the last from the
  shortest course (5 samples, 24 PD), whose coupling is not even negative. Without the floor, the
  fraction `F` again falls with PD in every line (−0.044 per 10 PD), the noise artefact described above.
* **Contact inhibition, the long arrest GSE197512 lacked, is not decisive either.** Over 120 days the
  inhibited cultures gained 1-2 PD and drifted 1.4 to 14 times less per day than dividing controls of the
  same lines (matched drift `C` 0.0004-0.0021 per 10 days against 0.003-0.014). The gap per unit of
  their small drift (`Q`/`C`, per 10 days) is −0.86, −1.34 and −0.07 in the three lines at 21 % oxygen,
  against −0.76, −0.66 and −0.47 in their controls: protection without divisions in two lines, none in
  the third. Selection alone would give none in all three, chromatin protection the controls' share in
  all three; drift this small sits near the floor, so neither is excluded.

Tables: `data/output/passaged_cells_sturm/` (as for GSE197512; the contact-inhibition and 3 % oxygen arms
under `Study part 3` in `arm_slopes` and `drift_coupling`).

## RNA: transcriptional drift by turnover class (GTEx, mouse, rat)

The same question in transcripts. For each tissue, protein-coding genes with at least 1 CPM in 80 %
of its samples, as log2(CPM + 1). The noise measure is between-individual spread: the age slope of
each gene's absolute residual after its mean trend and the technical covariates are removed
(`stochasticity.expression_drift`). RNA spread depends strongly on the mean, so every gene set is
compared with genes of the same expression decile, in SDs of the metric, as CpGs are matched on
level. A second version adds five age-free composition components of the genes outside every set.
Those components cannot change an age slope of the mean, but they can change spread, so here the
check is informative. Each tissue is resampled at the same size, without replacement. The human gene
sets reach Ensembl IDs through HGNC (Entrez first, then current, previous and unique alias symbols;
99–100 % of each set), and rodent genes reach them through their upper-case symbol (16,474 of 25,113
mouse and 16,752 of 30,454 rat symbols map, ambiguous ones dropped).

* **GTEx v10** (`entropy-age figures rna_gtex`): 47 tissues with 100 or more donors, 70 donors per
  resample, decade-bin ages, covariates sex, Hardy death scale, RIN, ischemic time and library size.
  Classes: renewing (10 tissues: blood, spleen, gut mucosa, skin, oesophagus mucosa, vagina, testis),
  slow (21), heart & muscle (3), brain (13 regions). The GTEx counts are read from an external parsed
  copy (`ENTROPY_AGE_GTEX_RNA`, `ENTROPY_AGE_GTEX_ATTRIBUTES`).
* **Mouse and rat atlases** (GSE314817, GSE314970; `entropy-age fetch rodent-rna`,
  `entropy-age figures rna_rodents`): inbred mice (24 tissues, 80–104 animals) and rats (23 tissues
  used, 47–85 animals), males and females at six ages from 6 to 27 months, covariates sex and library
  size, 56 mice and 32 rats per resample. Rat eye tissues of unclear turnover are left out; the rat
  retina counts as neural.

![RNA, GTEx](figures/rna_gtex.png)

![RNA, mouse and rat](figures/rna_rodents.png)

**Spread gain with age, housekeeping minus matched genes** (SD, mean over the class's tissues):

| Class | GTEx | Mouse | Rat |
|---|---|---|---|
| renewing | −0.08 | −0.29 | −0.10 |
| slow | −0.00 | −0.23 | −0.19 |
| heart & muscle | −0.02 | −0.23 | −0.12 |
| brain | +0.06 | +0.01 | −0.04 (retina −0.29) |

**Each class minus brain, housekeeping** (SD, 95 % interval over size-matched resamples):

| Contrast | GTEx | GTEx + PCs | Mouse | Mouse + PCs | Rat | Rat + PCs |
|---|---|---|---|---|---|---|
| renewing − brain | −0.12 (−0.26, +0.03) | +0.01 (−0.06, +0.09) | −0.28 (−0.42, −0.16) | −0.23 (−0.29, −0.16) | −0.07 (−0.19, +0.05) | −0.10 (−0.15, −0.04) |
| slow − brain | −0.06 (−0.16, +0.08) | +0.02 (−0.04, +0.08) | −0.24 (−0.37, −0.12) | −0.15 (−0.22, −0.07) | −0.16 (−0.29, −0.03) | −0.07 (−0.15, +0.02) |
| heart & muscle − brain | −0.09 (−0.29, +0.12) | −0.01 (−0.12, +0.08) | −0.25 (−0.40, −0.09) | −0.11 (−0.17, −0.04) | −0.09 (−0.24, +0.05) | −0.01 (−0.08, +0.06) |

Core-essential genes follow housekeeping genes in mouse (renewing − brain −0.18, −0.34 to −0.06;
with components −0.16, −0.22 to −0.10) and are within noise in rat and GTEx. Non-essential genes go
the other way in both rodents (mouse renewing − brain +0.54, interval +0.18 to +0.89; rat +0.32,
+0.03 to +0.61), the hitchhiking pattern prediction 3 describes.

What this supports, and what it does not:

* **Mouse: the methylation result holds in RNA.** Housekeeping and core-essential genes gain
  between-animal spread with age more slowly than matched genes in every class but brain, and the
  brain is without the protection. Heart keeps it, as in methylation, and the gap to brain survives
  the composition components for renewing and slow tissues. So in mouse "neurons lack the
  protection" now rests on two kinds of data from two cohorts.
* **Rat: same direction, weaker.** Slow tissues minus brain excludes zero, and renewing minus brain
  does once composition is adjusted; the other contrasts do not. The retina, neural but strongly
  protected (−0.29), pulls the brain class toward protection; without it the rat brain regions are
  −0.10 to +0.17.
* **GTEx: no turnover effect in spread.** Renewing tissues lean protected and brain leans the other
  way, but every interval includes zero, and with composition components every contrast is within
  0.03 SD of zero. Per-tissue intervals are about ±0.3 SD at 70 donors: post-mortem human RNA, with
  genetic differences between donors and only decade bins of age, carries too little age-related
  spread for a set contrast. The one GTEx contrast that excludes zero is in the mean, not the spread:
  housekeeping genes change their level with age less outside brain than in it (renewing − brain
  −0.09, −0.17 to −0.005; slow − brain −0.07, −0.14 to −0.01).
* **Not a turnover gradient in any species.** Where the protection shows (mouse, rat), slow and
  heart tissues have about as much as renewing ones. As in the mouse methylation data, the step is
  between brain and everything else, not along cell turnover.

Tables: `data/output/rna_gtex/`, `data/output/rna_rodents/` (`tissues`, `point_estimates`,
`contrasts`, `tissue_intervals`; gene mapping in `set_mapping` and `genes`).

## How to test the hypothesis on real data

```python
import polars as pl
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.stochasticity import compare_groups, drift_metrics, fit_relaxation, paired_comparison

betas = pl.read_parquet("my_betas.parquet")          # cpg column + one column per sample
samples = pl.read_parquet("my_samples.parquet")      # sample, age, cell_type
genes = pl.read_csv("essential_cpgs.csv")            # feature, group  (e.g. promoter CpGs of essential genes)

per_type = {}
for (cell_type,), sub in samples.group_by("cell_type"):
    m = FeatureMatrix.from_frame(betas, feature_col="cpg", samples=sub["sample"].to_list(), max_missing=0.05)
    ages = sub["age"].to_numpy()
    per_type[cell_type] = fit_relaxation(m.values, ages, m.features).join(
        drift_metrics(m.values, ages, m.features), on="feature")

# prediction 1: constrained genes drift less than level-matched CpGs, within each cell type
for cell_type, metrics in per_type.items():
    print(cell_type, compare_groups(metrics, genes, metric="entropy_slope"))

# prediction 2: the same CpGs drift more in post-mitotic than in renewing cells
print(paired_comparison(per_type["neuron"], per_type["blood"], metric="window_change", names=("neuron", "blood")))
```

The CLI does the per-cohort part: `uv run entropy-age stochasticity betas.parquet samples.parquet --groups genes.csv`.

### Data that would make this decisive

| Need | Candidates |
|---|---|
| Same donors, several purified cell types, age range | sorted blood subsets (e.g. Reinius et al. 2012, GSE35069; Tserel et al. 2015); BLUEPRINT; neurons vs glia from FANS-sorted brain (e.g. Kozlenkov et al.) |
| Post-mitotic vs renewing tissue from the same individuals | GTEx methylation (EPIC, 8+ tissues), the Mammalian Methylation Consortium multi-tissue data |
| Single-cell resolution (lineage-level selection) | scBS-seq / scNMT-seq aging studies; clonal haematopoiesis cohorts with methylation |
| Gene classes | DepMap common-essential genes, housekeeping genes (Eisenberg & Levanon 2013), PRC2 targets, cell-type identity TFs |
| Transcriptome version | single-cell expression noise across ages (Tabula Muris Senis), with the bulk model and per-gene σ |

Expected signatures if the hypothesis holds: a negative `matched_diff` for essential genes in blood,
gut and skin, about zero in neurons; a positive neuron-minus-blood `window_change` that is larger
for essential than for other genes; and higher neutral-site dispersion slopes in renewing tissues.

## Related work

Checked on 2026-11-09 (PubMed, the authors' GitHub, and Meyer's thesis read in full). Each entry is
described from its own abstract or text, not from a citing source.

### Stochastic clocks beyond the 2024 paper, by the same group

- **No stochastic clock has yet been validated on mammalian RNA.** David H. Meyer's stochastic RNA
  clocks are built and tested on *C. elegans*; his mammalian stochastic clocks are DNA methylation.
  The one application to mammalian RNA is in Ostermann et al., *JCI Insight* 2025
  ([doi:10.1172/jci.insight.190445](https://doi.org/10.1172/jci.insight.190445)). An "ensemble of
  100 independently trained stochastic clocks" scored human fibroblasts, iPSCs and directly induced
  neurons from 13 donors. Induced neurons kept a higher stochastic age than iPSCs, in line with 16
  RNAAgeCalc clocks. Its build details are in a supplement not read here, and the main text reports
  no accuracy on human chronological age.
- Gallrein, Meyer et al., *Nature Aging* 2026
  ([doi:10.1038/s43587-026-01067-5](https://doi.org/10.1038/s43587-026-01067-5)) applies BiT age and
  the stochastic clock to single neuron types of *C. elegans* (CeNGEN and Calico data). Neuron types of
  chronologically young worms differ almost twofold in predicted age, and the older-predicted types
  degenerate first. It is the cell-type version of our question in worms. Human and mouse enter only
  through correlations of gene-expression trends (GTEx, Tabula Muris Senis), not as clocks.
- Meyer D. H., *Development and Application of Aging Clocks*, PhD thesis, University of Cologne
  (2024), [kups.ub.uni-koeln.de/73668](https://kups.ub.uni-koeln.de/73668/). Its chapters are the
  published papers. Discussion §5.2.2 argues for our hypothesis in prose: "a specific subset of genes
  accumulates stochastic variation at a faster rate than others" is called conceivable, and such
  gene-set-specific noise may be hidden when noise is measured over the whole transcriptome. It notes that the stochastic clock's predictor genes in the neuron study overlap
  H3K27me3-bound regions (ChIP-Atlas), calls this anecdotal, and asks for a systematic test.
- Meyer, Maklakov & Schumacher, *Nature Aging* 2025
  ([doi:10.1038/s43587-025-00975-2](https://doi.org/10.1038/s43587-025-00975-2)): a review arguing that
  clocks reflect stochastic damage accumulation under insufficient maintenance and repair, not a
  program.

### Evidence that stochasticity differs between genes

- Lee et al., *Science Advances* 2021
  ([doi:10.1126/sciadv.abj9111](https://doi.org/10.1126/sciadv.abj9111)): genes lacking CpG islands,
  normally heterochromatic when inactive, are globally activated in aged cells and tissues of mice and
  humans, and this drives inflammatory secretion. Meyer reads it as CGI− genes becoming noisier with age
  while CGI+ genes do not. Our gene-set test matches on CpG-island status for this reason.
- Cheung et al., *Cell* 2018 ([doi:10.1016/j.cell.2018.03.079](https://doi.org/10.1016/j.cell.2018.03.079)):
  mass cytometry of chromatin modifications in human immune cells. Older adults show higher
  cell-to-cell variability and more between-person heterogeneity, driven mostly by non-heritable
  influences.
- Weidemann et al., *Science Advances* 2023
  ([doi:10.1126/sciadv.adh5138](https://doi.org/10.1126/sciadv.adh5138)): several constitutively
  expressed cell-division genes in fission yeast have mRNA noise below the Poisson limit, reached
  through multi-step transcription and degradation and fast mRNA export. Expression noise can be
  pushed down for genes that matter, which is a mechanism of protection other than selection.
- Jenkinson et al., *Nature Genetics* 2017 ([doi:10.1038/ng.3811](https://doi.org/10.1038/ng.3811)):
  methylation stochasticity quantified genome-wide as Shannon entropy from whole-genome bisulfite
  data, with methylation maintenance treated as a communication channel. This is the
  information-theoretic frame for site-specific maintenance that our per-site entropy rates build on.

### Protection across species

- McLain & Faulk, *Aging* 2018 ([doi:10.18632/aging.101413](https://doi.org/10.18632/aging.101413)):
  about 1,000 conserved promoters across 131 mammals have CpG density correlated with maximum lifespan,
  95 % positively. The authors suggest that these genes need more buffering against age-related
  methylation dysregulation in long-lived species.
- Mayne et al., *Scientific Reports* 2019
  ([doi:10.1038/s41598-019-54447-w](https://doi.org/10.1038/s41598-019-54447-w)): CpG density in
  42 promoters predicts maximum lifespan across 252 vertebrate genomes (R = 0.76).

### How these bear on our results

The blood gene-set result (essential and housekeeping CpGs drift more slowly than CpGs of the same
island status and start level) is consistent with all three protective routes these papers suggest:
buffering by CpG density, active chromatin, and selection among cells. CpG density and chromatin act
in every cell; selection needs turnover. Prediction 2 found the protection larger in blood than in
brain, which only selection (or blood's variable cell mix, not yet ruled out) predicts.

# Figure review: what each figure says on its own, multi-clock figures, new figures

Reviewer: Claude (Fable), report only. Read: `manuscript/paper.md` (captions at lines 99–106, 153–162, 198–202,
266–271, 299–303, 383–394, 461–467, 499–505, 530–532, 630–634), `manuscript/supplement.md`, the rendered PNGs in
`docs/figures/` and `docs/figures/captioned/`, the drawing code, and the tables under `data/output/`,
`docs/forest_series.csv` and `data/interim/clocks/`. Numbers marked "my check" were computed for this review with a
small script on the tables named beside them; they are not paper results and carry no adjustment beyond what is
stated.

## 0. Five findings that come before any redesign

1. **"Hannum" on the figures is a cohort, not a clock, and a reader cannot tell.** Figures 1, 2, 4 and 6 and
   Supplementary figures label panels "blood (Hannum)" / "Hannum blood" (`placebo_sets.py:125` via the tissue
   name, `equilibrium_drift.py:278`, `chromatin_match.py:299`, `shield_clocks.py:162`, `clock_weights.png` legend).
   This is the 656-donor blood data set GSE40279. In the same paper "Hannum" is also a clock (Fig. 9 row, paper
   line 40). The complaint "why on figures do we have only Hannum" is at least partly this collision: the only
   place a clock called Hannum is drawn in a main figure is one row of Fig. 9. Fix: name cohorts by tissue, size and
   accession everywhere ("blood, 656 donors, GSE40279"; "cortex, 213 donors, GSE74193") and reserve "Hannum" for
   the clock. Same for "Jaffe" and "Lunnon".
2. **Only one main figure shows any clock (Fig. 9), with 5 published clocks and 3 stochastic ones, and no figure
   shows the 35-clock panel.** `clock_increment.py` ends with `return []` (line 199): the nested clock models have
   tables and no figure. The supplement's clock figures (`shield_clocks.png`, `clock_weights.png`, `own_clocks.png`,
   `stochastic_component.png`) are all dot or bar rows.
3. **Fig. 9 draws the result the text has withdrawn.** The section title (paper line 604) and bold sentence
   (606–607) say the breast-cancer association depends on preprocessing. Fig. 9 is drawn from
   `data/output/incident_cancer/effects.csv` (authors' processed data, unmatched). The rebuilt and matched
   estimates (`incident_cancer_rsesame/effects.csv`, `cancer_matched_rsesame/odds_ratios.csv`,
   `clock_increment/nested_models.csv`) are not in any main figure. A reader who looks at figures first takes away
   "B predicts breast cancer", the opposite of the heading.
4. **One chart form carries almost everything.** Of the 9 main figures plus Box 1 and S22, 20 of 29 panels are
   horizontal point-and-interval rows; the rest are 2 bar panels, 1 stacked bar, 1 histogram, 1 strip, 1 line
   plot, and the three Box 1 schematic panels. No figure shows raw data from a real cohort (a CpG, a person, a species against
   lifespan). `plots.py` offers only scatter, box and band helpers (lines 67–176) and none of the main figures use
   them.
5. **Colour changes meaning from figure to figure.** Orange/blue is blood/cortex in Fig. 1a
   (`tissue_contrast.py:513`), real/placebo is blue/grey in Fig. 2, set identity in Fig. 3c, island/non-island in
   Fig. 4, divides/does not divide in Fig. 5, between/within person in Fig. 8b, all cases/lagged in Fig. 9,
   mouse/human in S22, and `clock_weights.png` uses blue for blood and orange for cortex, the reverse of Fig. 1.
   Fig. 6 leaves the house palette for `tab20`/`tab10` (`chromatin_match.py:304`). A reader has to relearn the
   legend nine times.

Two smaller facts the brief had wrong or that the drawing depends on:

- **seaborn is not installed.** `pyproject.toml:9–20` lists matplotlib 3.11 only; `import seaborn` fails in the
  project environment (it is a dependency of the standalone `scripts/clocks_biolearn.py` only). plotnine, altair,
  upsetplot, statsmodels are absent too. scipy and scikit-learn are present, so clustering is available.
- `docs/figures/incident_cancer_rsesame.png` (2026-11-18) is stale: its x axis prints "8×10⁻¹ 9×10⁻¹" overlapping
  and raw names ("F core_essential"). The code was fixed since (`incident_cancer.py:302–304`, `_measure_label`
  268–273) but the rebuilt figure was not redrawn.

## 1. Existing figures

For each: the story a figures-first reader should get; whether the figure gives it; changes.

### Box 1 (`scripts/box1_schematic.py`, `box1_schematic.png`)

- Story: one CpG spreads with age; that spread is a slope; a gene set's slopes sit left of matched CpGs' slopes.
- Does it: mostly. It is the clearest figure in the paper. Problems: the "young level" label is printed over the
  age-75 points (`box1_schematic.py:31`); panel (a) shows one CpG in blue and panel (b) then uses blue for
  the essential CpG and orange for the matched one, so (a) reads as "the essential CpG" although it spreads to
  25 %; panel (c) is simulated and its gap (≈0.5 SD) is 2–8 times the real ones, which the caption admits.
- Changes:
  - Draw (a) twice, small multiples: one essential-promoter CpG and one matched CpG, same axes, same colours as
    (b). Why: the contrast is the whole idea and (a) currently shows only one side.
  - Add panel (d): the same three steps on real data, one cohort (blood, GSE40279): two real CpGs and the real
    distribution of drift rates with the real gap (−0.15 SD). Why: a schematic drawn "larger than real for clarity"
    invites the question of how it looks when not enlarged. `equilibrium_drift/per_cpg.parquet` has per-CpG `rate`,
    `baseline`, `cpg_island`, `tissue` (218,129 rows) but no set membership column; membership needs a join
    (`placebo_sets/set_genes.parquet` plus the probe annotation), so this is a small code task, not a pure redraw.
  - Move "young level" label to the left margin.

### Figure 1 (`tissue_contrast.py:511–604`, `tissue_contrast.png`)

- Story: housekeeping and essential promoters drift less than matched CpGs; non-essential genes drift more; the
  effect is larger in blood than in brain, including within the same donors.
- Does it: the title says it, the panels bury it. 48 intervals in four panels; 12 tick labels of the form
  "core essential · binary entropy" at 7 pt; the measure the paper relies on (relative entropy) is the last row
  of each block, visually equal to the measure the caption says "overstates the gap". Panels (b) and (d) restate
  (a) and (c) as differences, with the 12 labels repeated; (b) is labelled "Blood − brain" for a blood − cortex
  contrast. Panel (c) has only three y labels and five unlabeled colours per label, with the legend below the
  x label.
- Changes:
  - Panel (a): relative entropy only, three sets × two tissues, as a dumbbell (blood and cortex joined by a
    line per set), with the gap value printed at the right. Why: the reader needs six numbers, not 24, and the
    line between the two dots is the blood-versus-cortex claim, which makes panel (b) unnecessary.
  - Move the other three measures to a 3 × 4 small-multiple heatmap (sets × measures, cell = gap in SD, diverging
    palette centred on zero, value printed, hatched where the interval includes zero), blood and cortex side by
    side. Why: "the result holds on four measures" is a robustness statement and reads at a glance as a block of
    one colour.
  - Panel (c): put tissue on the y axis (blood, then four brain regions), set as colour, so each row has a label.
    Drop (d) or fold it in as a dumbbell from blood to each region.
  - New panel: the raw basis of the gap in blood. Two overlaid densities (or a ridgeline, housekeeping / matched /
    non-essential) of per-CpG drift rate. Why: every other panel is a difference of means; the reader never sees
    how much two distributions overlap, which is what a −0.15 SD gap means.
  - Order sets as a gradient: core essential, housekeeping, non-essential (already so) and state it on the axis
    ("more essential ↑").

### Figure 2 (`placebo_sets.py:112–130`, `placebo_sets.png`)

- Story: random gene sets of the same size do not give a gap; the real sets fall far outside what random sets give.
- Does it: yes in substance, weakly in form. Rows are ordered measure-then-set, so the eye cannot follow one set;
  the label "entropy" here is "binary entropy" in Fig. 1; underscores are stripped inconsistently ("non essential"
  here, "non-essential" in Fig. 1); the title is not bold and is a method name, not a finding; the 200 draws are
  reduced to a bar.
- Changes:
  - Show the draws: per set, a half-violin or histogram of the 200 placebo gaps with the real gap as a marked
    line and the empirical rank printed ("0 of 200 placebo sets beyond the real set").
    `data/output/placebo_sets/placebo_draws.parquet` has them (4,800 rows; columns `tissue`, `group`, `metric`,
    `set`, `draw`, `matched_diff_sd`); `summary.csv` has `share_placebo_below_real`. Why: the claim is about where
    the real value sits in a distribution, so draw the distribution.
  - Relative entropy in the main panel, the three other measures as a small grid beside it.
  - Share one x scale between blood and cortex (now −0.3…0.2 against −0.2…0.05). Why: the smaller cortex effect
    is a result in Fig. 1 and is hidden here by rescaling.
  - Note in the figure that placebo means sit above zero (+0.01 to +0.05 in every row): the matching leaves a
    small positive offset, so the real gaps are, if anything, understated. That is visible and unremarked.

### Figure 3 (`species_shield.py:301–358`, `species_shield.png`)

- Story: the gap is negative in most mammals and tissues; whether long-lived species have a different gap.
- Does it: (a) and (b) do. (c) is a three-row interval plot under the title "Longer-lived species", whose
  housekeeping coefficient the caption then says vanishes off the 0 % floor. A reader sees a significant orange
  point and a title implying an effect. Tissue labels are raw ("CerebralCortex"); "other tissues" sits in the
  middle of the order; (b) repeats (a) as a histogram.
- Changes:
  - Replace (c) with the data: housekeeping gap against maximum lifespan (log axis), one point per
    species-tissue group, colour by tissue, the fitted line with its band, and a second line for the off-floor
    fit. `data/output/species_shield/species_gaps.parquet` has `species`, `tissue`, `max_lifespan`,
    `adult_weight`, `gap_sd`, `ci_low`, `ci_high` (139 housekeeping rows on `kl_slope`, 65 species, no missing
    lifespans; my check). Label a dozen species (human, bowhead whale if present, naked mole-rat, mouse, bats).
    Why: this is the one place in the paper where a scatter of real units exists, and a slope of +0.03 per ln
    unit with p 0.04 is honest only when the cloud is visible.
  - Merge (b) into (a) as a marginal density on top of the strip plot. Why: same numbers twice.
  - (a): sort tissues by median, put "other tissues" last, point size by number of animals, hollow points where
    the interval includes zero.
  - Optional, high value: a species × tissue heatmap ordered by a phylogeny or by order (bats, rodents, primates,
    cetartiodactyls), cells coloured by gap, blank where not measured. Why: shows coverage (which species have
    which tissues) and whether clades differ, neither of which the strip shows.

### Figure 4 (`equilibrium_drift.py:246–288`, `equilibrium_drift.png`)

- Story: essential promoters are kept near 0 % methylation by a lower gain rate and a higher loss rate than
  look-alikes.
- Does it: no. Four panels with four different x scales, values of order 0.001 per decade, negative rates, and
  no direct set-versus-look-alike comparison: the reader must compare row 1 with row 4 across a gap. The title is
  a method name. This is the hardest main figure to read.
- Changes:
  - One panel per tissue as a gain-versus-loss plane: x = gain rate, y = loss rate, one point per set and island
    class with a cross of intervals, an arrow from look-alikes to each set, and iso-equilibrium lines
    (gain / (gain + loss)) drawn in grey. Why: the claim is about the pair of rates and the level they imply; a
    plane shows both and the arrow is the comparison the text makes.
  - Rates in % methylation per decade (×100) so ticks read 0.1, 0.2, not 0.0025.
  - Add the direct evidence: per-CpG age slope against young level, binned (hexbin or binned means with bands),
    for a set and its look-alikes; the intercept and slope of that line are the gain and loss rates. Columns
    exist in `equilibrium_drift/per_cpg.parquet` (`baseline`, `rate`, `k_gain`, `k_loss`, `cpg_island`,
    `tissue`). Why: shows where the two numbers come from.

### Figure 5 (`headline_figures.py:205–240`, `division_shield.png`)

- Story: the gap is at least as large in dividing cells as in non-dividing ones and vanishes when dividing
  cells are held confluent.
- Does it: yes; grouped blocks with bold headers and a two-colour code is the best layout in the paper. Problems:
  no figure title; blocks have different units (per year, per doubling, per month) on one shared x axis, which
  the caption warns about but the figure invites; within-block contrasts (colon − hippocampus etc.) are not
  drawn; the caption's last sentence ("the glial excess is not seen in the oligodendrocyte lineage") cannot be
  read off the figure.
- Changes:
  - A bracket with the within-block difference and its interval at the right of each block. Why: the claim is
    the contrast inside a block.
  - A light band or separate x axis per block, or at least the unit printed in each header (done in two of five).
  - Add n (donors, animals, cultures) at the right margin.
  - Finding-first title: "The gap is present in dividing cells and absent in cells held confluent".

### Figure 6 (`chromatin_match.py:296–388`, `chromatin_match.png`)

- Story: add chromatin state to the matching and most of the blood gap goes away; less of the cortex gap does.
- Does it: (d) does and (a) explains why; (b) and (c) are 20 and 28 intervals in 7 colours whose legend is at the
  bottom of the page, and they carry the same message as (d).
- Changes:
  - Make (d) the lead and draw it as a waterfall or slope chart per tissue: gap after island × level matching →
    after own state → after sorted accessibility → after both, with the share removed printed on each step, for
    housekeeping and core essential. Why: "accounts for" is a before-and-after claim.
  - (a): order states from active to repressed, use a palette that groups them (active TSS reds, transcribed
    greens, enhancers yellows, bivalent purples, polycomb greys, quiescent white, as Roadmap does) and add the
    cortex composition next to blood. Why: `tab20` puts 1_TssA and 2_TssAFlnk in dark and light blue by accident
    and gives 15_Quies a salmon; the Roadmap colours are what readers of chromHMM figures know.
  - Move (b) and (c) to the supplement as they are.
  - Bring the palette back to `plots.PALETTE`.

### Figure 7 (`headline_figures.py:150–203`, `expression_shield.png`)

- Story: how much a gene is expressed predicts how little its promoter drifts, and expression accounts for 71–82 %
  of the housekeeping gap.
- Does it: yes; (b) is the only dose-response in the paper and is persuasive. Problems: no title; the legend of
  (b) covers the lines at deciles 1–5 ("Monocytes" overlaps data, "Fibroblasts in culture (doublings)" runs
  across); error bars of four series overlap at each decile.
- Changes: legend outside or direct labels at decile 10; dodge the four series by ±0.1; thin grey line for each
  series with a bold pooled line; put (b) first, since it is the general fact and (a) its consequence; add the
  housekeeping genes' own expression decile distribution as a rug or small histogram under (b), which shows why
  matching on expression removes the gap.

### Figure 8 (`within_person.py:273–316`, `satsa_idats.png`)

- Story: the gap is seen inside one person followed over years, not only between people of different ages.
- Does it: weakly. (a) is bars without intervals, including a row "measurement floor" whose bars are as large as
  the effect, with no visual cue that this row is a control; (b) uses "B core essential", "F core essential",
  "S core essential" as labels, which a figures-first reader cannot decode, and the between-people squares have
  no intervals. No person appears in a figure about persons.
- Changes:
  - New lead panel: spaghetti plot. One thin line per person (age at draw against S housekeeping, or F), with the
    mean within-person slope in bold and the between-person slope dashed.
    `data/output/satsa_idats/per_draw.parquet` has `person`, `age`, `epic`, `B/F/S core_essential`,
    `B/F/S housekeeping` for 1,072 draws of 444 people, 287 of them with two or more draws (my check). Why: this is
    the paper's only longitudinal data; show it as longitudinal.
  - (a): points with intervals instead of bars (the table `satsa_idats/gaps.csv` has `p_permutation`, `p_gene`
    but no interval columns, so intervals need to be written out), the floor row greyed and labelled "control:
    technical replicate noise".
  - (b): replace letter codes with words ("drift at look-alike promoters (B)", "share avoided by the set (F)",
    "signed shift (S)"), and mark the array type of each draw in the spaghetti panel (`epic` column) since a
    450k→EPIC change within a person is the obvious artefact.

### Figure 9 (`incident_cancer.py:276–321`, `incident_cancer.png`)

- Story (per the text): on the deposited data, drift at look-alike promoters predicts later breast cancer and
  clocks do not; on arrays rebuilt from raw files and matched on chip, that association goes and the shield share
  becomes protective.
- Does it: no; see finding 3. Also 84 intervals for one message; drift measures and clocks are in one
  undifferentiated list; clocks are not grouped by kind; "Zhang 2019" here is the chronological clock while the
  35-clock panel has only the Zhang mortality score (`Zhang_10`), and the two are unrelated in these women
  (r 0.06 between accelerations, my check), so the name is a trap.
- Changes:
  - Breast cancer becomes the figure; any cancer and colorectal go to the supplement.
  - Four columns for the same rows: deposited unmatched, rebuilt unmatched, deposited matched, rebuilt matched
    (files above; all exist). Why: the preprocessing dependence is the finding.
  - Rows in labelled blocks: drift measures (B, F, S) / first-generation clocks / second-generation / pace /
    stochastic, with a light band behind each block.
  - Replace by or pair with the multi-clock figure of section 2 and the specification curve of section 3 (N2).

### Supplementary Fig. S22 (`scripts/forest_series.py:153–176`, `forest_series.png`)

- Story: across 42 independent series the housekeeping gap is negative far more often than positive.
- Does it: this is the figure that carries the abstract's breadth claim, and it is in the supplement. As drawn it
  is sorted by species then gap, so the tissue structure is lost; platform is in the label text only; one
  interval (lung next to adenocarcinoma, EPIC, to +0.76) sets the x range and compresses everything else; the
  title is an instruction on how to read colours. Counts from `docs/forest_series.csv` (my check): 42 series,
  31 negative; intervals exclude zero in 21, of which 20 negative and 1 positive (lung next to adenocarcinoma,
  450k).
- Changes:
  - Promote to a main figure, directly after Fig. 1.
  - Group rows by tissue (blood, brain, lung, liver, heart, kidney, colon, muscle, adipose, breast, myometrium),
    species as colour, platform as marker shape (six platforms: RRBS, 450k, EPIC, EPIC v2, mammalian 40k, MM285),
    point area by n (column `n`, partly empty for the RRBS rows, to be filled).
  - Clip the x axis at ±0.5 with an arrow for the clipped interval.
  - Add a summary row per species and overall (random-effects mean with prediction interval), and a top margin
    with the count "20 of 21 resolved series below zero".
  - Mark the tumour-adjacent and "all conditions" series as a separate block: they are the positive ones and are
    not healthy ageing tissue.

## 2. Clocks

### Where a single clock or no clock appears

| Place | What is shown | Problem |
|---|---|---|
| Figs 1, 2, 4, 6 | "Hannum" = cohort GSE40279 | reads as "the Hannum clock"; no clock is involved |
| Fig. 9 | Horvath 2013, Hannum, PhenoAge, DunedinPACE, Zhang 2019, StocH/P/Z (`incident_cancer.py:44`) | flat list, light-preprocessing indices, deposited data only |
| paper 614–615 | B "kept 1.44 with DunedinPACE and PhenoAge in the model" | two clocks, text only |
| `clock_increment.py:35` | `PRIMARY = "PhenoAge"`; 34 others "exploratory" | tables only, no figure (`return []`, line 199) |
| `cancer_matched_rsesame/odds_ratios.csv` | DunedinPACE, PhenoAge, Horvath 2013 | three clocks, no figure |
| paper 623 | SJLIFE deaths: B against DunedinPACE only | one clock |
| `shield_clocks.png` (S) | 17 clocks, partial r with P at equal B | closest thing to a clock panel; sorted by kind alphabetically ("causal" first), labels overlap between panels |
| `clock_weights.png`, `own_clocks.png`, `stochastic_component.png` (S) | 12 / own / 3 clocks | fine, supplement only |

### What exists to draw a multi-clock figure today (verified)

- `data/interim/clocks/GSE51032_clocks.csv`: long, 29,365 rows, columns `sample, clock, output, prediction,
  measured_weight_share`; **839 people (653 women) × 35 clocks**, no missing predictions.
- `data/interim/clocks/GSE51032_clocks_panel.csv`: 44 clocks with `status` (35 `kept`, 9 `coverage`),
  `arrays_covered`, `median_weight_share`, `output`.
- `data/output/incident_cancer_rsesame/per_person.csv`: 839 rows; `B/F/S housekeeping`, `B/F/S core_essential`,
  `age`, `female`, `icd`, `years_to_diagnosis`, `chip`, `mean_beta`, `fraction_1..11`, and eight house clock
  indices including `DunedinPACE`.
- `data/output/clock_increment/nested_models.csv`: 36 rows (35 clocks + "panel mean + disagreement"), columns
  `clock, role, sets, clock_or, clock_low, clock_high, disagreement_or, B_or_given_clock, B_low, B_high, lr_p_B,
  F_or_given_B, F_low, F_high, lr_p_F`.
- `data/output/clock_increment/disagreement_correlations.csv`: disagreement ~ B r = 0.37, ~ F r = 0.02.
- `data/output/shield_clocks/correlations.csv`: 960 rows, `cohort` (blood, cortex), `components` (0, 5), `clock`,
  `kind`, `group`, `measure` (P, F, B, P given B, set drift), `r`, `p`, `n`: a second and third cohort for the
  correlation panel.

Gaps to close before drawing:

- **Third generation is one clock and it failed coverage.** `GSE51032_clocks_panel.csv`: DunedinPACE `coverage`,
  `arrays_covered` 0.22, median weight share 0.89; DunedinPoAm38 0.0. So the biolearn panel has no pace clock. A
  DunedinPACE column exists as the house linear index in `per_person.csv` (839 values); in a figure it must carry
  a different marker ("index on our matrix, 89 % of weight measured") or wait for the pending rescoring. GrimAge
  is absent from the kept list altogether.
- **`nested_models.csv` is older than its code.** The file (17:06) has no `design`, `controls_on_case_chip`,
  `P_or_given_B` or `valid_resamples` columns that `clock_increment.py` (17:25, lines 104–122, 178–186) now
  writes. Rerun before drawing; the figure should show both designs.
- **Per-woman disagreement is not saved.** It is computed at `clock_increment.py:156–158` and only its two
  correlations are written (lines 190–198). A scatter needs the per-woman table (sample, panel_mean,
  disagreement, B, F, case status).
- **The disagreement panel contains clocks that do not run with age the same way.** `ADULT_BLOOD`
  (`clock_increment.py:41–55`) includes Bocklandt, whose prediction falls with age in these women (r −0.42 with
  age; its acceleration correlates −0.13 to −0.53 with the other panel clocks; my check), and YingAdaptAge, built
  to be protective (acceleration r −0.64 with YingDamAge). The SD across standardised accelerations therefore
  counts a sign convention as disagreement. Either orient every clock so that higher = older/worse before taking
  the SD, or drop the two. This changes the 0.37 and should be settled before the disagreement panel is drawn.
- **Name hazards.** "Zhang 2019" (chronological, house index) and `Zhang_10` (mortality, biolearn) are different
  clocks; "Horvath 2013" = `Horvathv1`. House and biolearn accelerations agree for Horvath, Hannum, PhenoAge,
  StocH, StocZ (r 0.993–0.999 in women, my check), so one naming table can serve both.

### Grouping to use (names as in `GSE51032_clocks.csv`)

| Group | Clocks available now | Note |
|---|---|---|
| 1st generation, chronological | Horvathv1, Horvathv2, Hannum, Lin, Weidner, VidalBralo, Garagnani, Bocklandt | Bocklandt, Garagnani, Weidner are 1–3 CpG clocks; Bocklandt reversed |
| 2nd generation, phenotypic / mortality | PhenoAge, HRSInCHPhenoAge, Zhang_10, YingCausAge, YingDamAge, YingAdaptAge | GrimAge missing |
| 3rd generation, pace | DunedinPACE | pending (coverage); house index only |
| Stochastic | StocH, StocZ | StocP failed coverage (0.88) |
| Mitotic | EpiTOC1 | |
| Out of domain (negative controls) | Knight, LeeControl, LeeRobust, LeeRefinedRobust, Mayne, Bohlin (gestational), PEDBE (paediatric), DNAmClockCortical | |
| Trait predictors | Smoking, Alcohol, BMI ×2, BodyFat, Education, Total/HDL/LDL cholesterol (McCartney, Reed), CVD_Westerman | |

`entropy_age/data/clocks.py:22–40` already has a `kind` per clock for 17 of them; extend that table, do not keep a
second one in the figure code.

### The multi-clock figure (proposed main Figure 9, four panels)

**(a) How the clocks relate to each other.** Clustered heatmap of the correlation between age accelerations
(residual on age) of the 35 clocks in the 653 women, rows and columns ordered by hierarchical clustering, a
colour strip along the side for generation, a dendrogram on top. What it will show (my check, residual on age
only): the first- and second-generation adult clocks form one block with modest agreement (median pairwise r
0.25 among 14 adult clocks; Horvathv1–Horvathv2 0.58, Hannum–PhenoAge 0.52), the gestational clocks form their
own block, and the trait predictors sit apart. That one panel answers "do clocks agree" before any outcome.

**(b) How each clock relates to the drift measures.** A three-column heat strip (B, F, S at housekeeping
promoters) aligned to the rows of (a), value printed in each cell, so the clustering carries over. From my check
(women, both sides residual on age): first-generation clocks go with a stronger shield and not with drift
(Horvathv1 B −0.10, F −0.24, S −0.33; Hannum +0.06, −0.21, −0.34; Horvathv2 −0.18, −0.29, −0.36); PhenoAge goes
with drift and not with the shield (B +0.27, F 0.00); the Zhang mortality score with both weakly (+0.22, +0.16);
DunedinPACE (house index) with neither (+0.04, +0.03); EpiTOC1 is nearly the signed shift (S −0.84, F −0.37,
B +0.35). The generations separate on this strip, which is the paper's "what clocks read" claim (paper 600–602)
shown on people, with every clock. Add the same strip for blood GSE40279 and cortex from
`shield_clocks/correlations.csv` as two more columns to show it is not one cohort.

**(c) Odds ratios in the nested models.** A forest with clocks in the row order and blocks of (a), two columns
of points per row: the clock's own odds ratio for breast cancer diagnosed over 2 years later (M0), and F's odds
ratio once that clock and B are in the model (M2). From `nested_models.csv` (163 matched sets): the interval of
the clock's own odds ratio excludes 1 for 6 of 36 rows (YingCausAge 1.31, EpiTOC1 1.63, and four out-of-domain or
trait scores below 1: Mayne, Bohlin, HDL cholesterol, body fat); B given the clock runs 0.76–1.17 and is never
resolved; **F given B and the clock is 0.67–0.80 in all 36 models, with a likelihood-ratio p below 0.05 in 35 of
them** (the exception is the model with EpiTOC1, p 0.20). Draw F's column as near-vertical line of points with a
shaded band at the PhenoAge (primary) estimate. Caveats to print on the panel: the 35 models share one F and one
set of women, so this is one result shown to be insensitive to the choice of clock, not 35 tests; PhenoAge was
the pre-specified clock; the out-of-domain block is a negative control and should be greyed.

**(d) Disagreement.** Scatter of per-woman clock disagreement against B (r 0.37) and against F (r 0.02), cases
and controls in two colours with marginal densities, and the odds ratio of disagreement (0.85, from the
"panel mean + disagreement" row) printed. Draw only after the orientation issue above is settled.

Form: matplotlib `imshow` + `scipy.cluster.hierarchy` (both installed) is enough for (a)–(b); seaborn's
`clustermap` gives the same with side colours in a few lines and is the one library I would add (see section 3).

## 3. New figures

Each row: claim, data (all files verified to exist; columns as read), form, library.

**N1. The clock panel (section 2 above).** Claim: clocks of different generations read different parts of
methylation change, and the shield share's association with later breast cancer does not depend on which clock
is in the model. Data: `data/interim/clocks/GSE51032_clocks.csv`, `GSE51032_clocks_panel.csv`,
`data/output/incident_cancer_rsesame/per_person.csv`, `data/output/clock_increment/nested_models.csv`,
`disagreement_correlations.csv`, `data/output/shield_clocks/correlations.csv`. Form: clustered heatmap + aligned
heat strip + two-column forest + scatter with marginals. Library: matplotlib + scipy; seaborn optional.

**N2. Specification curve for breast cancer.** Claim: the association of B with later breast cancer depends on
analytic choices, the association of F is stable in sign once adjusted; here is every specification that was
run. Data: `data/output/incident_cancer/effects.csv` and `incident_cancer_rsesame/effects.csv` (`outcome,
lag_over_2y, adjusted, measure, cases, controls, odds_ratio, low, high, p`),
`data/output/cancer_matched/odds_ratios.csv` and `cancer_matched_rsesame/odds_ratios.csv` (`measure, cases,
adjusted, sets, same_chip_share, odds_ratio, boot_low, boot_high, p`), `clock_increment/nested_models.csv`,
`breast_replication/odds.csv`. Form: top panel, odds ratios of B (and, in a second colour, F) sorted by size with
intervals; bottom panel, a dot matrix with one row per choice (deposited / rebuilt; unmatched / matched; age only
/ adjusted; all cases / over 2 years; clock in model; cohort) and a filled dot where the specification uses it.
What it will show, from the tables: for B housekeeping, lagged and adjusted, 1.44 (1.17–1.76) deposited
unmatched, 1.35 (1.10–1.66) rebuilt unmatched, 1.11 (0.69–2.18) rebuilt matched; so on these tables the loss comes
with matching on chip, not with rebuilding alone, which is a sharper statement than the heading's "depends on
preprocessing" and worth checking against the text at paper 616–620. Library: matplotlib (two axes sharing x
order). This is the standard form for exactly this problem and replaces a paragraph of numbers.

**N3. Every series, as a main figure.** Claim: the housekeeping gap is negative in 31 of 42 independent series
and resolved below zero in 20 against 1 above. Data: `docs/forest_series.csv` (`series, species, platform,
processing, n, age_min, age_max, fit, gap, ci_low, ci_high, interval`). Form: forest grouped by tissue with
species colour, platform marker, n as point area, clipped axis, pooled rows; right margin: a raincloud (half
violin + points) of the 42 gaps by species. Library: matplotlib.

**N4. Gap against lifespan across mammals.** Claim: across 65 species the gap does not track maximum lifespan
once the 0 % floor is handled. Data: `data/output/species_shield/species_gaps.parquet` (columns in section 1,
Fig. 3), `lifespan_fits.csv` (`metric, group, outcome, term, coef, se, p, groups, species`). Form: scatter on a
log lifespan axis with fitted bands for the two fits, labelled species, tissue as colour; companion heatmap
species × tissue ordered by clade. Library: matplotlib; `adjustText` would help with species labels (small,
pure Python) but manual offsets for a dozen labels are enough.

**N5. Within-person trajectories.** Claim: the signed shift at housekeeping promoters moves in the same
direction inside a person over years as it differs between people of different ages. Data:
`data/output/satsa_idats/per_draw.parquet` (`sample, person, age, female, epic, B/F/S core_essential, B/F/S
housekeeping`; 1,072 draws, 444 people, 287 with ≥ 2 draws), `person_slopes.csv` (`measure, within_per_decade,
low, high, between_per_decade, people, draws`). Form: spaghetti plot with within- and between-person slopes
overlaid; beside it a histogram of per-person slopes with the share below zero. Library: matplotlib.

**N6. The distributions behind the gap (ridgeline).** Claim: a gap of −0.06 to −0.26 SD is a shift of a whole
distribution, visible in every tissue, not a tail effect. Data: per-CpG drift rates in
`data/output/equilibrium_drift/per_cpg.parquet` (`feature, baseline, rate, cpg_island, tissue`, 218,129 rows),
joined to set membership (`placebo_sets/set_genes.parquet` plus the probe-to-gene annotation; the join is not in
a table yet, so this needs a short script). Form: ridgeline, one ridge per tissue, two overlaid densities per
ridge (set, matched), medians marked. Library: matplotlib `fill_between` on a KDE from scipy is enough; seaborn
`FacetGrid` ridge or `joypy` are conveniences only.

Not proposed: an alluvial or an UpSet plot. The paper has no flow between categories to show; an UpSet plot of
the overlap of core-essential, housekeeping and non-essential gene sets would be three sets, for which a
three-number sentence or a Venn is enough.

### Libraries

- **Keep matplotlib as the base.** The house style lives in `plots.py:25–44` and every figure passes through it.
- **Add seaborn** (one dependency, in `pyproject.toml`): `clustermap` with side colours (N1), `violinplot`/
  `stripplot` with `hue` and `dodge` (Figs 2, 3, N3), `kdeplot` ridges (N6), `regplot` bands (N4). It respects
  the existing rcParams, so the look does not change. Per the cluster rule, `uv sync --dry-run` first while jobs
  run from this environment.
- **Do not add plotnine or altair.** plotnine would give a second styling system beside `plots.py`; altair
  produces web charts and needs a separate renderer for print PNG/PDF.
- **Add three helpers to `plots.py`** so the forms are reused, not rewritten per module: `forest(ax, table,
  group=…, marker_by=…, clip=…)`, `dumbbell(ax, …)`, `clustered_heatmap(fig, matrix, side_colours=…)`. Many
  modules now hand-roll the same `ax.plot([lo, hi], [y, y])` loop.
- **One semantic palette, fixed in `plots.py`:** blood, cortex, other tissues; human, mouse; set colours (core
  essential, housekeeping, non-essential, matched = grey). Then colour means the same thing in every figure
  (finding 5).

## 4. What to do first, by reader gain

1. **Rename the cohorts on every figure** ("blood, 656 donors, GSE40279", not "Hannum"). An hour of label edits;
   removes the misreading behind the complaint and a real ambiguity with the Hannum clock.
2. **Draw the clock panel (N1) and make it the clock figure of the paper**, after rerunning `clock_increment`
   (stale table), saving the per-woman disagreement, and orienting or dropping Bocklandt and YingAdaptAge.
   Largest gain on the stated complaint; all data but the pace clock are on disk.
3. **Replace Fig. 9 with the breast-cancer comparison across processing and design, and add the specification
   curve (N2).** The figure then agrees with the section's heading.
4. **Promote the 42-series forest (N3) to the main text, regrouped by tissue.** It carries the breadth claim.
5. **Fig. 1: relative entropy only as a dumbbell, the four measures as a heatmap, plus the two real
   distributions.** Cuts 48 intervals to 6 and shows data.
6. **Fig. 3c → lifespan scatter (N4); Fig. 8 → trajectories (N5).** Two figures that at present show no data
   points of the thing they are about.
7. **Fig. 4 as a gain–loss plane; Fig. 6 led by the before/after steps.** The two least readable main figures.
8. **One semantic palette and finding-first titles on Figs 2, 4, 5, 7; add seaborn and the three helpers.**

## Things to verify before acting on this review

- Whether the drop of B's adjusted odds ratio comes with matching on chip (rebuilt unmatched adjusted 1.35,
  1.10–1.66, in `incident_cancer_rsesame/effects.csv`; rebuilt matched 1.11) and whether paper 616–620 says so.
  The two analyses also differ in population (all participants against women only), so the comparison is not
  like for like.
- The "my check" correlations of section 2(b) use residuals on age only; `clock_increment.py:154` uses age, mean
  beta, cell fractions and chip for the disagreement accelerations. The figure should state which, and probably
  show the adjusted version.
- EpiTOC1 and S (r −0.84 in women): the paper text never mentions epiTOC (no match in `paper.md`), though
  `shield_clocks.png` shows the same in blood (partial r about −0.8). If S is close to a published mitotic score
  measured at different CpGs, that belongs in the text as well as in N1.

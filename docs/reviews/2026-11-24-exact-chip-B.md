# Independent check: B in the exact-chip matched design (EPIC-Italy, rebuilt arrays)

Date: 2026-11-24. Report only; no code or results changed. Scripts in the reviewer's scratchpad (not committed):
`exact_chip.py`, `exact_chip2.py`, `exact_chip3.py`, `exact_chip4.py`, `tail.py`, `lowchips.py`. They reuse `match`,
`clogit` and `_score` from `src/entropy_age/experiments/cancer_matched.py` and the women and pool built as in
`clock_increment.run` (653 women with finite clocks, 68 chips; 163 breast cases over 2 years, 340 controls; B, F
z-scored over all women; adjustment = age, mean beta, 11 fractions, each z-scored over the pool). Intervals below are
Wald unless marked "boot" (200 set resamples).

## Verdict in one line

The exact-chip OR of ~4.5 for B reproduces, but it is an artefact of low-intensity arrays: it lives entirely on the
16 chips that carry the ~120 arrays with 8x lower signal, where cases happen to have the lower-intensity, noisier
arrays. On normal-intensity chips B has no association (adjusted OR 0.61, 0.07-5.0; case above its controls in 48 % of
sets). The paper's "lost the adjusted association" stands; F survives the same exclusion.

## 1. Reproduction and stability

How B is built: `incident_cancer.run` -> `health_states.loo_scores`: each person is scored against a reference
fitted on all other people in the cohort, across chips (`fit_reference`, no chip term, no components),
`d = |x - m| / sqrt(m(1-m))`, B = the stratum-matched mean of `d` over look-alike CpGs. So chip shifts enter B
directly, and so does any per-array noise: a noisier array sits further from the reference at every CpG and gets a
higher B. B is strongly right-skewed (skew 3.9 after z-scoring; max 7.3 SD).

Exact chip ("same chip only"): 147 sets on 39 chips; 114 sets with 1 control, 24 with 2, 9 with 3; 16 cases unmatched.

| model (exact chip) | OR per SD of B | interval | LR p |
|---|---|---|---|
| B alone | 4.41 | 2.13-9.12 (boot 2.27-13.6) | 2e-8 |
| B + age | 4.36 | 2.10-9.04 | 4e-8 |
| B + age + mean beta + 11 fractions | 4.49 | 1.87-10.8 (boot 1.87-23.6) | 4e-5 |
| B + PhenoAge acceleration + adjustment | 4.52 | 1.85-11.0 | 6e-5 |
| B per within-chip SD (0.577 of total) + adjustment | 2.38 | 1.43-3.95 | 4e-5 |
| F given B + adjustment | 0.74 | 0.52-1.06 | 0.09 |

So the table value reproduces (4.5, F 0.74). The fit converges, but it is fragile in the sense that matters:
leave-one-chip-out ORs range 3.8-8.1; bootstrap upper limits reach 24-28 (a few extreme sets); the median
case-minus-controls difference in B is only +0.06 SD while the mean is +0.37 (top differences 3.2-4.3 SD); the case
is above its controls in 63 % of sets (sign test p 0.003). Rank-based B (normal scores over all women) gives
OR 1.34 (0.92-1.95); within-chip normal scores 1.41 (1.03-1.92); log B 2.79 (1.50-5.18). The magnitude comes from
the right tail.

## 2. Artefact?

**Chip composition.** Cases sit on 40 of 68 chips (0-6 per chip; all 40 also carry controls). The original study
placed matched pairs side by side: of 182 mixed case/control neighbour pairs (same row, C01/C02), 173 have the
cancer sample in C01 and the control in C02, with median age difference 1.5 years and 98 % sex concordance. So
column is nearly a label for case status (column alone in the exact-chip sets: OR 0.38 per SD, p 4e-15); adjusting
for column drops B to 1.74 (0.67-4.50), but that is over-adjustment. A column artefact itself is small: in
control/control neighbour pairs B(C01) - B(C02) = +0.01 SD, in cancer/cancer pairs +0.03 SD. Rows: chip-centred B
falls from R01 (+0.22) to R05 (-0.21); cases are spread evenly over rows, so row does not drive it.

**Array quality.** Within chip, B is mostly a quality score: r(B, GCT) = -0.72, r(B, RG distortion) = +0.55,
r(B, RG ratio) = +0.44, r(B, median red) = -0.43, r(B, mean intensity) = -0.38, r(B, mean beta) = -0.41;
QC plus position explain 70 % of B's within-chip variance (`data/interim/sesame/GSE51032/GSE51032.qc.csv`).
Within chip, cases have lower intensity (p 7e-5), lower GCT (p 1e-4) and higher RG ratio (p 0.008) than
controls.

**The decisive split.** 120 of 845 arrays have mean intensity < 1000 (median ~520, against ~4700 for the rest; median
red and green signal ~10x lower), on 15 chips (2 fully low, 13 with 6-10 of 12 low). Among women, 31 % of breast
cases and 16 % of controls are on such arrays. The arrays with the highest within-chip B are these: of the 24 women
with chip-centred B > 1, 20 are cases, with median intensity 316 against 4503 overall.

| exact-chip sets | sets | B alone | B + adjustment | case > controls | F given B + adj |
|---|---|---|---|---|---|
| on low-intensity chips | 68 | 6.10 (2.22-16.7) | n/a | 79 % | n/a |
| on normal chips | 79 | 2.15 (0.58-8.0) | 0.61 (0.07-5.0) | 48 % | 0.56 (0.33-0.95), p 0.02 |
| all, arrays < 1000 dropped | 96 | 2.87 (0.76-10.9) | 1.38 (0.21-8.9) | | 0.52 (0.33-0.82), p 0.003 |
| all, arrays < 2000 dropped | 79 | 2.15 (0.58-8.0) | 0.61 (0.07-5.0) | | 0.56 (0.33-0.95), p 0.02 |

(ORs per overall SD; among the kept arrays the within-chip SD of B is 0.21 of the overall SD, so 1.38 is ~1.07 per
within-chip SD.) On the low-intensity chips the case has the lower intensity in 62 % of sets. A linear term in
intensity does not remove it (B + adj + intensity 3.97) because the relation is a threshold, not a line; dropping
the failed arrays does.

## 3. Why "same chip first" dilutes it

67 % of B's variance is between chips (ANOVA F 17.5; F 30 %, mean beta 41 %, age 29 %). In same chip first only 20
of 163 sets have every control on the case's chip (these 20 separate; no finite fit), and 143 sets have off-chip
controls. Off-chip case-control differences in B have SD 1.62 against 0.88 on-chip, and their mean (+0.46) is
mostly chip mean (+0.33). Once mean beta and fractions are in the model the chip-mean part of B points the other way
(chip mean B + adj: OR 0.57, 0.33-0.98), cancelling the within-chip part. Centring B within chip in same chip first
gives 1.54 (0.96-2.47, p 0.06); B alone 1.89 (1.40-2.55); B + adj 1.11 (0.70-1.76). Dropping arrays < 1000 in same
chip first turns B + adj to 0.12 (0.03-0.57) while F becomes 0.43 (0.29-0.65), p 1.5e-5. The exact design
removes chip noise and so exposes the within-chip intensity effect undiluted; same chip first buries it under the
chip-level component.

## 4. Verdict

Artefact, with good confidence. The exact-chip OR rests on one technical batch: arrays with ~10x lower signal on 16
chips, where B (an across-cohort deviation score, so noise raises it) is inflated and cases happen to carry more of
the weak arrays. On chips with normal signal, B shows nothing (adjusted 0.61, case above controls in 48 % of sets),
and rank-based B is weak everywhere. Not settled by this check: why cases have weaker arrays within those chips
(DNA amount or aliquot handling by case status is possible and would itself be a sample-handling confounder, not
biology), and whether the deposited-data result (matched 1.71/1.76, logistic 1.44) is carried by the same arrays.
What would settle it: (a) rerun the deposited-data and rebuilt analyses with arrays of mean intensity < 1000 (or
SeSAMe's own failed-array flags) excluded, and with B scored against a within-chip or quality-matched reference;
(b) a second prospective cohort. The paper's sentence should stay as it is in substance; if the authors add the
exact-chip design, a fair reading is: "Matching on chip exactly gave a larger odds ratio for B (4.5 per SD), but it
came entirely from 16 chips with weak array signal, where B tracks array quality; on chips with normal signal B was
not associated (0.61, 0.07-5.0), while F stayed protective (0.56, 0.33-0.95)." Separately, the low-intensity batch
should be checked against every B and F result on rebuilt EPIC-Italy arrays, since F's protective estimate grows
when it is removed.

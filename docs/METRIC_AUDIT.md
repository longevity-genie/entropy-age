# Metric audit: what rested on binary entropy (roadmap RM33, run R68)

Asked by antonkulaga (2026-11-14): go through previous discoveries, evidence and hypotheses and correct the
computations. The reason is R67. Binary entropy measures how far a CpG is from 50 %, so at a low-methylated site a
*loss* of methylation counts as a gain in order. R67 found that essential promoters lose methylation faster than
their look-alikes. Every entropy gap therefore mixed two things: how much a site moved, and which way.

**The corrected measure.** Relative entropy from the young state, `KL(x || m_young)` in bits
(`stochasticity.relative_entropy`; `kl_slope` in `drift_metrics`; `K` in the per-person scores). It counts a move up
or down alike, assumes nothing about where a site drifts to, and weighs a move away from an extreme level more, as
binomial sampling does. It is close to the squared standardised deviation, so it carries both mean change and spread.
It was added beside the old measures, never in place of them, and every figure whose headline rested on entropy was
rerun with all four measures side by side (R68a-d, same seeds; the old numbers reproduce).

The per-person shield scores `P`, `B`, `F` already measure distance from the young state, and need no correction;
nor do clocks, survival, accessibility, or any result on `|mean change|` or dispersion (a fork of this audit listed
them: R1, R2, R5, R6, R13-R18, R23-R43, R45-R49, R52, R54, R56, R57, R60, R62, R64, R67).

## Verdicts

Gaps are set minus look-alikes in SD (same island status and young-level decile); permutation p ≤ 0.001 unless
given.

| result | entropy (old headline) | relative entropy | other measures | verdict |
|---|---|---|---|---|
| C1 blood (Hannum), housekeeping / core-essential | −0.35 / −0.23 | −0.15 / −0.09 | \|mean\| −0.16 / −0.11; dispersion −0.29 / −0.17 | **kept, resized**: about 2.3× smaller |
| C1 cortex (Jaffe), housekeeping / core-essential | −0.21 / −0.13 | −0.06 / −0.07 | \|mean\| −0.05 / −0.05; dispersion −0.08 / −0.08 | **kept, resized**: 1.7-3.6× smaller |
| C1 non-essential genes drift *more* | blood +0.31, cortex +0.08 | blood +0.14, cortex −0.00 (p 0.94) | cortex \|mean\| −0.08 | **blood kept; cortex overturned** |
| R61 deposited vs reprocessed, cortex housekeeping | −0.22 → −0.15 | −0.07 → −0.09 | | relative entropy steadier across preprocessing |
| C2 chromatin removes the shield, blood (R3, R4) | 77-84 % | 83-91 % | dispersion 83-88 % | **kept** |
| C2 chromatin removes the shield, cortex (R3, R4) | 43-69 % | 25 % own state, 56 % with sorted ATAC | dispersion 22-46 % | **resized**: in cortex chromatin explains a quarter to a half |
| R58 essentiality beyond chromatin, blood | +0.022 (0.008-0.033) | −0.007 (−0.021, 0.008), p 0.33 | dispersion −0.009, p 0.33 | **overturned**: no remainder in blood |
| R58 essentiality beyond chromatin, cortex | +0.018, p 0.08 | +0.022 (0.009-0.040) | dispersion +0.027, p 0.01 | **kept**: the remainder is cortex's |
| R44 sorted CD4+ T cells, housekeeping | −0.22 | −0.015 (p 0.005) | \|mean\| −0.09; dispersion −0.04 | **resized**: small; monocytes hold (−0.10) |
| GTEx tissues, housekeeping (whole blood / lung / muscle) | −0.29 / −0.27 / −0.07 | −0.14 / −0.09 / −0.01 | | **kept, about halved**; muscle and testis ≈ 0 |
| Mouse array (Huang), renewing tissues (blood, colon, spleen, liver, kidney) | −0.12 to −0.31 | −0.05 to −0.19 | \|mean\| −0.05 to −0.20 | **kept** |
| Mouse array, post-mitotic (muscle, heart, frontal cortex, hippocampus) and adipose | −0.22 to +0.10, mixed | **+0.03 to +0.09, all p ≤ 0.002** | \|mean\| −0.05 to −0.09; dispersion +0.05 to +0.12 | **reframed**: no shield where cells do not divide; essential promoters shift less on average but scatter more. Hypothalamus is the exception (−0.08) |
| C6 species, share of 139 groups negative | 96 % (−0.30) | 89 % (−0.13) | dispersion 85 % | **kept, resized** |
| C6 species off the floor (10-90 %) | 93 % (−0.24) | 73 % (−0.10) | dispersion 73 % | **kept, resized** |
| R22 relative shield grows with lifespan | −2.42, p 0.04 | +1.43, p 0.23 | dispersion +0.46, p 0.63 | **overturned**: entropy only |
| R59 no lifespan relation off the floor | −0.014, p 0.48 | +0.004, p 0.84 | dispersion −0.002, p 0.92 | **kept** (with all CpGs relative entropy gives +0.032, p 0.04: the floor again) |

## What this changes

* **The shield stands, at about half the size we reported.** Every essential-set gap in human blood, cortex, GTEx,
  renewing mouse tissues and mammals keeps its sign on the corrected measure. Entropy overstated it 1.7-3.6×, most in
  cortex, because essential promoters drift *down* (R67), which entropy reads as order.
* **Turnover is no longer a contradiction.** On the mouse array the shield is in every renewing tissue and reversed
  in four of five post-mitotic ones; entropy had blurred this into "mixed". It fits the RRBS finding that cortex
  and heart carry a third of the renewing tissues' shield (R63, entropy and dispersion; to recheck on relative
  entropy) and points at division: where cells divide, the faster removal of R67 has work to do.
* **Three statements are withdrawn**: essentiality beyond chromatin in blood (R58), the shield's share growing with
  lifespan (R22), and the large CD4+ T-cell shield (R44). The cortex chromatin share is cut to a quarter to a half.
* **Mouse RRBS and array aging gaps (R73), housekeeping, entropy → relative entropy:** GSE80672 blood −0.35 →
  −0.30; BXD liver −0.52 → −0.32; six RRBS tissues −0.16 to −0.33 → −0.03 (adipose) to −0.21 (lung); Meer cortex
  −0.09 → −0.06 (interval excludes 0), heart −0.10 → −0.07, liver −0.36 → −0.22, lung −0.29 → −0.14 (heart, liver
  and lung intervals now include 0 with about 20 mice). All keep their sign; smaller, as in human data.
* **Turnover, two platforms.** By RRBS the post-mitotic cortex and heart keep a small shield (−0.06, −0.07) a
  third of liver's; on the mammalian array four post-mitotic tissues are reversed (+0.03 to +0.09). Both order the
  shield by turnover; they disagree on whether post-mitotic tissue keeps a little of it or loses it.
* **Not yet rerun:** read-level epiallele entropy (R50-R53; its partner measure PDR agrees) and the per-person `K`
  in the effect tables.

Tables: the `kl_slope` rows in `data/output/{tissue_contrast,reprocessed_jaffe,chromatin_match,essentiality_gradient,
sorted_cells,gtex_tissues,mouse_tissues,species_shield,species_floor,probe_robustness}/`.

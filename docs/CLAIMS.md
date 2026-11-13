# Claims ledger: each claim with its predictions, refuters and alternatives

Asked by antonkulaga (2026-11-12): every evidence claim should say what must be seen if it is true, what must
not be seen, and what else could explain it. One section per claim; the argumentation graph in the results
artifact is drawn from these. Numbers are those of `docs/extension.md`. The matrix `docs/evidence_matrix.csv` scores every result against every hypothesis and four artefact
explanations (−2 contradicts … +2 strongly predicts); its row spread is the result's diagnosticity. Status: **holds** (predictions met,
refuter absent, main alternatives tested), **open** (an alternative not yet excluded), **weak** (one study or
a test with little power), **refuted**.

## C1 — Essential promoters drift less with age than matched CpGs (the shield)

- **If true, we should see:** a negative gap in every tissue and species; robust to how CpGs are matched
  (island status, starting level, probe design); present per person, not only on average.
- **We should not see:** the gap vanish or flip once CpGs are matched more finely, or in an independent cohort.
- **Observed:** negative in human blood, brain and GTEx tissues, 11 mouse tissues, culture, 96 % of 65 mammals,
  and inside one sorted cell type (monocytes -0.14 SD, CD4+ T cells -0.22, R44; about half of whole blood's);
  molecule by molecule in read-level WGBS (entropy -0.05 to -0.075 SD, R50; cfDNA -0.03 to -0.06, R51; within 27 sorted cell types -0.04 to -0.08, R53); within one person over
  time weakly (SATSA housekeeping divergence -0.017 SD, R52, on betas with cell counts regressed out);
  probe design explains a sixth of blood's gap, masks nothing (R20); per-person reliability 0.997.
- **Alternatives:** (a) a floor: CpGs near 0 % cannot drift down, and deciles leave a floor inside the lowest
  stratum; (b) array compression at the extremes; (c) island CpGs vary less for technical reasons.
- **Separating test:** finer level matching in the lowest decile; read-level WGBS free of array effects
  (Loyfer 2023, R50: holds, so (b) is excluded as the whole story); (c) matched on island status.
- **Reprocessing (R61):** Jaffe cortex from its raw IDATs (SeSAMe-like) keeps every essential-set gap below zero;
  sizes move by up to 0.06 SD.
- **Status:** holds; (a) answered for the species array: off the floor (CpGs at 10-90 %) the gap holds in 93 % of
  139 groups (R59).

## C2 — The shield is chromatin maintenance (H1)

- **If true, we should see:** matching on the tissue's own chromatin removes most of the gap, another tissue's
  much less; a fixed share of drift as cells divide; no dependence on turnover; loss of a promoter-maintaining
  enzyme removes it.
- **We should not see:** the shield grow with divisions or follow turnover; survive chromatin matching intact.
- **Observed:** own chromatin removes 67-91 % in blood, another tissue's about 30 %; constant share over 136
  doublings; not ordered by turnover; Tet3 loss removes it (C4).
- **Alternatives:** circularity: active chromatin and an unmethylated promoter may both be outcomes of the same
  protection, so matching on one removes the other without showing cause; reference epigenomes come from other
  people.
- **Separating test:** change chromatin and watch drift (natural experiments: writer syndromes kept it,
  Tet3 removed it); chromatin and methylation from the same people.
- **Status:** holds as association; causality partly supported (C4), circularity open.

## C3 — The shield is not selection among cells (H2 rejected)

- **If selection, we should see:** the shield grow with divisions; stronger in renewing tissues; absent in
  non-dividing neurons; weaker when proliferation slows (CR); changed in clonal expansions.
- **If selection, we should not see:** a shield in neurons; a constant share over doublings.
- **Observed:** constant share in culture; not ordered by turnover; equal in blood and brain after chromatin;
  neurons lose it only when Tet3 is halved; clonal haematopoiesis leaves it; not tuned to lifespan. One result
  favours selection: CR weakened it in one blood study, not replicated in liver; with chromatin openness
  matched two thirds of that effect goes with the shield (R55); a second, short CR series in young mice shows CR's
  changes sparing essential promoters (R56); a transplant bottleneck moves essential and non-essential promoters
  alike (R57).
- **Alternatives (to the rejection):** long-term cultures turn oligoclonal (Hollmann 2020), and steady
  selection per doubling also gives a constant share; clonal haematopoiesis clones were diluted (no allele
  fractions).
- **Separating test:** sorted stem cells young against old; a cohort with clone sizes (controlled access).
- **Status:** holds against selection as the main cause; a minor part in renewing blood not excluded.

## C4 — Tet3 guards essential promoters in neurons (H4)

- **If true, we should see:** in Tet3-deficient brain, essential CpGs gain methylation and lose the shield
  against littermates, in neurons that do not divide; specific to Tet3, not every TET.
- **We should not see:** the effect explained by batch or age; an equal change at all CpGs (a global shift).
- **Observed:** against same-chip, same-age littermates F +0.17 (cortex), +0.25 (striatum); essential CpGs gain
  methylation (S +0.012, +0.014); Tet1 loss strengthens the shield instead (R35).
- **Alternatives:** (a) array compression lifting the lowest CpGs (S also rose in other arms); (b) a shift of
  neuron and glia proportions in Tet3 mice; (c) one study, 8 mice per tissue; different arrays for Tet1 and Tet3.
- **Separating test:** S with compression held fixed (R36): the gain stays (+0.010 cortex, +0.013 striatum), so
  (a) is excluded; composition components; a Tet3 knockout on the same platform; sorted neurons.
- **Status:** weak (one study); (a) excluded. Dnmt3a W326R (R37) retargets methylation to Polycomb-type look-alikes:
  writers are aimed by chromatin (C2).

## C5 — Broken chromatin writers add noise but keep the shield

- **If true, we should see:** higher drift in carriers of NSD1, EZH2, KMT2D, CHD7 variants, with an unchanged or
  stronger shield share.
- **We should not see:** a weaker shield share.
- **Observed:** drift up in Sotos, Weaver, Kabuki, CHARGE; shield share unchanged or stronger.
- **Alternatives:** each syndrome's own methylation signature raises drift by design; the Sotos/Weaver series
  has no ages; infants misled the first run (fixed).
- **Separating test:** ages for GSE74432; signatures removed before scoring.
- **Status:** weak.

## C6 — The shield is universal across mammals and not tuned to lifespan

- **If true, we should see:** a negative gap in nearly every species; no lifespan slope with relatedness held
  fixed; long-lived species drift slower at unprotected CpGs with the shield's share kept.
- **We should not see:** a stronger shield in long-lived species.
- **Observed:** 96 % of 139 groups; lifespan slope +0.025 (p 0.40) with phylogeny; look-alike drift falls as
  lifespan^−1.4 with the relative shield kept or larger; the blood lean rests on human; off the 0 % floor the
  lifespan slope is −0.014 (p 0.48, R59).
- **Alternatives:** conserved array CpGs at essential promoters sit near 0 % (floor); taxonomy covariance is
  coarse.
- **Separating test:** a dated phylogeny; species with WGBS.
- **Status:** holds.

## C7 — Part of age drift is directed, mitotic loss (H3)

- **If true, we should see:** solo-WCGW CpGs lose methylation in dividing tissue (blood), not in neurons; more
  loss with more divisions (culture).
- **We should not see:** the same loss in cortex.
- **Observed:** blood 2.7x faster loss (p 0.0002); none in cortex; epiTOC1 diverges in blood only.
- **Alternatives:** a shift of blood cell types with age that differ in late-replicating methylation.
- **Separating test:** composition-adjusted; culture series with doublings (pending); sorted cells.
- **Status:** open.

## C8 — Clocks read unprotected noise; outcome clocks do not see the shield (H5)

- **If true, we should see:** clock acceleration rise with drift at unprotected CpGs; at equal overall drift,
  no weaker shield in people with fast clocks; outcome clocks unrelated to the shield.
- **We should not see:** fast clocks with a weaker shield at equal drift.
- **Observed:** acceleration rises with B; at equal B faster clocks go with a stronger shield; PhenoAge,
  DunedinPACE, mortality r ≤ 0.07.
- **Alternatives:** clock CpGs are themselves look-alike-type CpGs, so the link at equal drift is partly built
  in.
- **Separating test:** clocks trained on shielded and look-alike CpGs only (R33): equal accuracy, so a clock
  does not need unprotected CpGs; the shielded clock's acceleration overlaps less (r 0.63 against 0.80). Clocks
  retrained on simulated change (R38): each CpG's direction keeps 70-80 % of their accuracy, undirected noise 0-10 %.
- **Status:** revised. Clocks read directed change with random timing, not undirected noise; that is why they are
  blind to the shield, which acts on undirected drift. The first half of the title is refuted, the second holds.

## C9 — PhenoAge's non-stochastic part is blood-cell composition

- **If true, we should see:** PhenoAge's acceleration and its gap to its stochastic twin shrink when 12 blood-cell
  fractions are removed.
- **We should not see:** an unchanged non-stochastic part after composition is removed.
- **Observed:** in Tong 2024, PhenoAge's COVID-19 acceleration vanishes with composition. In Hannum (R38), 11 cell
  fractions explain 15.5 % of PhenoAge's acceleration (4.9 % of StocP's); removing them leaves the two accelerations'
  agreement unchanged (0.44 -> 0.45).
- **Alternatives:** composition and intrinsic aging change together (inflammaging), so removing one removes part
  of the other.
- **Separating test:** sorted cells (R44): in pure monocytes clocks read age (Zhang 0.92) and contamination explains
  under 2 % of acceleration.
- **Status:** refuted as the main explanation; composition is a minor part of PhenoAge.

## C10 — Partial reprogramming does not restore the shield

- **If restoration, we should see:** treated mice with a stronger shield and promoter CpGs closer to young mice.
- **Observed:** shield stronger in muscle, weaker in spleen, unchanged in four tissues; treated mice no closer to
  young at promoter CpGs.
- **Alternatives:** the regimen did not rejuvenate promoter methylation at all, so the test had little to see.
- **Separating test:** a regimen that measurably rejuvenates methylation (one OSKM cycle, GSE156557).
- **Since:** one cycle in 1-year-old mice (R62) moved too little to test. Transient reprogramming of human
  fibroblasts (R65, from IDATs) about doubles promoter drift, and essential promoters
  keep the same shield share as before (F unchanged, p 0.25-0.35): the shield holds through a reset. No young donors,
  so restoration toward youth stays untested.
- **Status:** weak (restoration untested; the shield's persistence under reprogramming supported).

## C11 — Chromatin accessibility has the shield in blood but not in brain

- **If true, we should see:** essential promoters spread less with age in blood ATAC; consistent in brain cohorts.
- **Observed:** blood weak (−0.05 SD after cell mix); cortex reverses; hippocampus flips with adjustment.
- **Alternatives:** sample quality (promoter read share) dominates single-nucleus pseudobulk at these sizes.
- **Separating test:** a cohort with per-sample quality metrics; sorted bulk ATAC.
- **Status:** refuted as a clean readout for brain; weak in blood.

## C12 — Drift at unprotected promoters in blood marks breast-cancer risk years ahead

- **If true, we should see:** higher B in women later diagnosed with breast cancer, not weaker with more years
  before diagnosis; robust to chip, cell mix and global methylation; replicated in a second cohort.
- **We should not see:** the effect vanish with adjustment, shrink with lag (a tumour present at the draw), or
  hinge on one batch.
- **Observed:** EPIC-Italy (R46, R46b): odds ratio per SD 1.43 (all), 1.62 (over 2 years), 1.66 (over 5 years);
  adjusted for mean beta, 11 cell fractions and chip 1.44 (p 6e-4); balanced over 72 chips; the shield share F
  does not predict; clocks' signals do not survive adjustment; B keeps 1.44 with DunedinPACE and PhenoAge in the model
  (r with DunedinPACE 0.02). In SJLIFE deaths (R48) B adds nothing beyond DunedinPACE: breast-specific, not mortality.
- **Alternatives:** (a) a cell type outside EpiDISH's 12 or a rare clone; (b) a cohort-specific artefact, since
  van Veldhoven's global hypomethylation did not replicate in NOWAC; (c) reproductive or hormonal factors that
  shift blood methylation and raise breast-cancer risk.
- **Separating test:** a second prospective cohort (NOWAC, Breakthrough Generations, controlled access); sorted
  cells; hormonal covariates.
- **Status:** weak (one cohort), but survives every check at hand.


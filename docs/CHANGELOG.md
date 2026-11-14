# Changelog

Newest first.

## Unreleased

- Relative entropy (`kl_slope`) added to the mouse interventions, BXD, RRBS-tissue and Meer comparisons (RM33 leftovers, run R73); `drift_hotspots` figure on one x scale.
- 2026-11-14: R72 written up (RM34 item 4): blood's rising spread is a Polycomb programme; cortex's is broad.
- `docs/PROVENANCE.md` (RM34 item 5): one row per dataset, deposited or rebuilt from IDATs, platform and the figures that read it.
- `figures drift_hotspots` (RM34 item 4): the top 5 % CpGs by rising spread, enrichment by region class against a level-matched background; run R72 planned.
- `figures global_fdr` (RM34 item 3, R71): one Benjamini-Hochberg correction over every p in every output table; results between p 0.024 and 0.05 are lost.
- 2026-11-14: R69 written up (RM34 item 1): placebo gene sets give no shield; the matching does not make it.
- 2026-11-14: R70 written up (RM34 item 2): the blood shield is not cell mix; 12 cell fractions held leave it unchanged.
- `figures composition_gaps` (RM34 item 2): Hannum gaps with 12 Houseman/EpiDISH cell fractions as covariates; run R70 planned.
- `figures placebo_sets` (RM34 item 1): random non-set genes through the same matching, 200 draws per set and tissue, on four drift metrics; run R69 planned.
- 2026-11-14: R68 written up (`docs/METRIC_AUDIT.md`): on relative entropy the shield keeps its sign at about half the size; R58 blood, R22 and the CD4+ T size withdrawn; mouse shield ordered by turnover.
- 2026-11-14: RM33 (metric audit) and RM34 (bias checks after the omics re-analysis practice) planned (Anton).
- Metric audit (Anton): `relative_entropy` (KL from the young state, no drift-to-0.5 assumption) as `kl_slope` in `drift_metrics` and in the lifespan, mouse, species, chromatin and essentiality comparisons; per-person scores gain `K`, the matched KL gap; runs R68 planned.
- 2026-11-13: R67 written up (RM32): essential promoters rest at a lower methylation level, held by faster removal, not slower gain.
- 2026-11-13: R65 written up (RM17 done): transient reprogramming of human fibroblasts moves the methylome a lot and keeps the shield share.
- `figures equilibrium_drift` (RM32): methylation gain and loss rates per gene set (two-state Markov model, pooled over CpGs, donor bootstrap) and per-CpG relaxation toward a bounded equilibrium; run R67 planned.
- 2026-11-13: RM32 planned (Anton): shield measured as drift toward each CpG's own equilibrium, repair being asymmetric.
- `figures satsa_idats` (RM28): the within-person analysis on SATSA's 450k waves reprocessed from their IDATs (E-MTAB-7309); `within_person` split into loaders and a shared `analyse`.
- `figures gill_reprogramming` (RM17): Gill 2022 transient reprogramming of human fibroblasts from IDATs; reprogrammed against failed and controls within donor.
- `docs/PREPROCESSING_PRACTICE.md`: array preprocessing practice against this project's design (RM18 item 5); RM18 done and off the roadmap, RM17, RM19 and RM28 unattended.
- R61 written up: Jaffe cortex reprocessed from its IDATs keeps every essential-set gap below zero (RM18 item 3); `reprocessed_jaffe` figure legend below the axis.
- `figures clock_weights` (R64): published clocks' weights on protected CpGs against their age and spread slopes, blood and cortex (asked by antonkulaga).
- 2026-11-13: R63 written up: the shield in mouse cortex, heart, liver and lung by RRBS; smaller and static in the non-dividing tissues.
- 2026-11-13: `meer_tissues` (RM2): heart, liver, lung and cortex at 6-30 months, RRBS (GSE121141); GSE184211 found to be inside GSE223748.
- 2026-11-13: R62 written up: one OSKM cycle (RRBS) finds too little aging to reset; RM17 next on human fibroblasts after the IDAT pipeline (RM18).
- 2026-11-13: `reprogramming_rrbs` (RM17): one OSKM cycle in GSE156557 (RRBS) against old and young controls; R62 submitted. RM20 blocked on desk #28 (CALERIE); RM14 waits until 2026-12-13.
- 2026-11-13: `scripts/reprocess_idats.py` (methylprep, SeSAMe-like) and `reprocessed_jaffe` (RM18 (3)): Jaffe's gaps on deposited against IDAT-reprocessed betas; R61a submitted.
- 2026-11-13: R59 written up: the species shield holds off the 0 % floor; the long-lived lean was the floor.
- 2026-11-13: RM20 noted: hypothesis 2 tested three ways without support (R55-R57); RM13 (3) waits for read-level culture data.
- 2026-11-13: RM14's DNMT3A loss-against-gain pair checked and not run (gain arm of one in blood); RM14 waits for a cohort meeting its rule.
- 2026-11-13: R60 written up: the contrast scores halve the batch excess of drift B; smoking moves only S; chemotherapy moves F against B.
- 2026-11-13: `score_confounders` (RM13 (5), Seeboth 2020): smoking, array, slide, diagnosis and chemotherapy against the per-person B, F and S in SATSA and SJLIFE; R60 submitted.
- 2026-11-13: RM12 and RM16 taken out of the roadmap (done; their results in `docs/extension.md` and the run log); RM14, RM17, RM20 no longer wait on them.
- 2026-11-13: `species_floor` (RM13 (4)): the species gaps and lifespan fit refitted on CpGs at 10-90 % only; R59 submitted. RM13 (2) answered from the CR tables (B beside F).
- 2026-11-13: R58 written up: a small essentiality gradient survives chromatin, expression and breadth, in blood and cortex.
- 2026-11-13: `essentiality_gradient` (RM13 (1)): DepMap 24Q4 gene effect against promoter drift among broadly expressed genes, chromatin state, expression and breadth held (GTEx v8 medians); R58 submitted.
- 2026-11-13: R57 written up: after a transplant bottleneck no essential-specific shield gain; RM12 done.
- 2026-11-13: `hsct_bottleneck` (RM12, bottleneck): divergence of transplant recipients from their donors against recipients' own later draws, at essential promoters against look-alikes (Holland 2024, Mendeley j5krtjfj6y).
- 2026-11-13: RM16 done (steps 1, 2, 4); step 3 waits for public data.
- 2026-11-13: R56 written up: in a second CR blood series (young mice, 12 weeks, no drift) CR's changes spare essential promoters; no support for the B6D2F1 weakening.
- 2026-11-13: R55 written up: chromatin openness removes about 60 % of mouse blood's shield and two thirds of CR's effect on it.
- 2026-11-13: `cr_blood_replication` (RM16 step 4): a second CR blood series, GSE282661 (8 CR, 8 TLS+CR, 8 controls, one array run); R56 submitted.
- 2026-11-13: `cr_chromatin` (RM12, chromatin held fixed): CR's blood effect refitted with CpGs also matched on openness in six mouse blood cell types (ImmGen ATAC, GSE100738); pyBigWig added; R55 submitted.
- 2026-11-13: R54 written up: the human marrow and lymph-node reversal is donor noise and donor set, not a tissue effect; RM16 step 2 done.
- 2026-11-13: `marrow_anomaly` (RM16 step 2): human bone marrow and lymph node's reversed shield checked by donor resampling, chip held, same donors across tissues, a 20-donor null; R54 submitted.
- 2026-11-13: R53 written up: the molecule-level shield holds within 27 sorted cell types; RM22 done and taken out of the roadmap.
- 2026-11-13: `read_level` adds the atlas's sorted cell types (27 groups, 101 samples), age slopes within cell type; R53 planned.
- 2026-11-13: R52 written up: within one person the shield share grows with age (SATSA); CpG-level test weak on preprocessed betas.
- 2026-11-13: `within_person` (RM28): SATSA repeated draws, divergence per year apart from the measurement floor; R52 planned.
- 2026-11-13: R51 written up: the molecule-level shield replicates in cfDNA of the same donors.
- 2026-11-13: `read_level` runs every source present (white blood cells, cfDNA of the same donors); R51 planned.
- 2026-11-13: R50 written up: the shield holds molecule by molecule in read-level WGBS (Loyfer WBC); index validated r 0.96 against arrays. E30.
- 2026-11-13: R48 written up: in SJLIFE deaths drift and shield add nothing beyond DunedinPACE; EPIC-Italy's B is independent of the clocks (breast-specific).
- 2026-11-13: R49 written up (RM31 done): the shield within blood lineages, largest in T cells, smallest in myeloid cells.
- 2026-11-13: `read_level` (RM22): epiallele entropy and discordant reads at 450k promoter windows from Loyfer WBC WGBS (23 donors); run R50 planned.
- 2026-11-13: `lineage_shield` (RM31): lineage-specific age slopes from mixed blood (CellDMC-style), shield per lineage; run R49 planned.
- 2026-11-13: `sjlife_deaths` (RM28): death among 4,520 SJLIFE survivors and controls against drift, shield and clocks; streaming reader for the 28 GB table; run R48 planned.
- 2026-11-12: R46b written up: drift B predicts breast cancer years ahead in EPIC-Italy, robust to chip, cell mix and mean methylation; claim C12, E29.
- 2026-11-12: R47 written up (RM27 done): rapamycin pulls back the noisy part of aging in 7 of 8 tissues, exenatide nothing; not ordered by turnover. E28.
- 2026-11-12: R46 logged (drift B predicts breast cancer years ahead in EPIC-Italy; shield F does not); `incident_cancer` gains mean beta, cell mix and chip adjustment; R46b planned.
- 2026-11-12: R42 written up: in Vitality 90+ neither the shield nor the mortality clocks predict death (positive control fails; uninformative).
- 2026-11-12: `intervention_split` (RM27): directed against noisy aging reversed by exenatide and rapamycin in 11 mouse tissues; regression to the mean guarded by split halves (tested); run R47 planned.
- 2026-11-12: R45 written up (RM29 done): methylation noise builds from early life, no onset at sexual maturity.
- 2026-11-12: `incident_cancer` (RM28): shield, drift and clocks in blood drawn before a cancer diagnosis, EPIC-Italy; run R46 planned.
- 2026-11-12: R39 written up (RM24 done): Alzheimer's leaves the shield; schizophrenia cortex, smokers' lung and muscle after long illness move it. E27.
- 2026-11-12: R44 written up: the shield holds in sorted monocytes and T cells; DunedinPACE tracks plaque and calcium, the shield does not. E26.
- 2026-11-12: R41 written up (RM26 done): clock acceleration is shared across a donor's organs, clock-free drift is not.
- 2026-11-12: R38 written up: stochastic clocks are directed clocks with random timing; composition is a minor part of PhenoAge (C9 refuted as main explanation). E24, E25.
- 2026-11-12: R36 and R37 written up: Tet3's gain is not compression; Dnmt3a W326R methylates Polycomb-type look-alikes. Evidence E22, E23.
- 2026-11-12: R43 written up (RM30 done): methylation spread barely reaches expression spread in GTEx donors.
- 2026-11-12: `drift_onset` (RM29): change point in methylation noise against age at sexual maturity across mammals; a null test guards against development read as an onset; run R45 planned.
- 2026-11-12: `sorted_cells` (RM31): shield, drift, clocks and plaque/CAC in MESA sorted monocytes and CD4+ T cells; run R44 planned.
- 2026-11-12: `meth_expr_spread` (RM30): methylation and expression spread coupled per donor and with age in GTEx; run R43 planned.
- 2026-11-12: docs/LITERATURE_COVERAGE.md: each idea rated done, partly or open against verified papers; RM29, RM30 filed, RM26 scoped.
- 2026-11-12: `survival_shield` (RM28): Cox model of death on shield, drift and clocks in Vitality 90+ (GSE68194); run R42 planned.
- 2026-11-12: docs/LONGITUDINAL_SURVEY.md: open and controlled methylation cohorts with follow-up; RM28 filed.
- 2026-11-12: `organ_coupling` (RM26): clock and drift acceleration correlated across tissues of the same donors; run R41 planned.
- 2026-11-12: `tissue_pace` (RM25, R40): tissue aging orders compared across mouse, rat and human, methylation and RNA.
- 2026-11-12: `health_states` (RM24): shield, drift and shift against Alzheimer's, Braak stage, schizophrenia and smoking in cohorts on disk; run R39 planned.
- R33 done (RM11): clocks from shielded, look-alike and random CpGs read age equally well; the shielded clock's acceleration overlaps less with the others.
- `docs/planned_tests.csv`: the planned and running tests, each linked to the hypotheses and artefact explanations it would decide; drawn as a third layer of the artifact's graph.
- `docs/evidence_matrix.csv`: 20 results scored against 5 hypotheses and 4 artefact explanations; drives the artifact's matrix and its ELK-laid graph (asked by antonkulaga).
- `docs/CLAIMS.md`: every claim with its predictions, refuters, alternative explanations, separating test and status (asked by antonkulaga).
- `fetch epidish` and `figures stochastic_component` (RM23): Tong et al.'s stochastic-clock procedure repeated in Hannum blood with four rate models (their effect-based rates, direction only, pure noise, shield only), RR2 and intrinsic acceleration with 12-cell composition removed. Dependency: rdata.
- Zhang 2019 (the chronological elastic net StocZ was built from) joins the pinned clocks (RM23).
- `fetch dnmt3a` and `figures dnmt3a` (RM12): Dnmt3a W326R/+ gain of function against littermates in three tissues over four ages, age-centred on the wild types (GSE324160).
- R35 done (RM12): Tet1 loss in cortex adds drift but strengthens the shield, unlike Tet3.
- `figures mmc_interventions` also reports the signed shift with each array's compression (mean distance of betas from 0.5) held fixed (RM18). Roadmap: RM18 building, RM16 after RM22.
- `fetch tet1` and `figures tet1` (RM12): Tet1 knockout and catalytic-dead cortex against wild type (GSE271375, 8 each, bisulfite arrays).
- R32 done (RM17): cyclic partial reprogramming does not restore the shield consistently; promoter CpGs barely move toward youth.
- `fetch reprogramming` and `figures reprogramming` (RM17): cyclic OSKM in vivo (GSE190665, six mouse tissues), treated and young mice against untreated old mice.
- `figures own_clocks` (RM11): elastic-net clocks from shielded, stratum-matched look-alike and random CpGs of equal number, out-of-fold, in Hannum blood and Jaffe cortex.
- R29 done (RM11): faster clocks do not mean a weaker shield; at equal drift they go with a stronger one, and outcome clocks are blind to it.
- `figures shield_clocks` (RM11): per person, the shield score against each published clock's age acceleration (stochastic clocks included), in Hannum blood and Jaffe cortex. Asked by antonkulaga.
- R27 done (RM19): hippocampus does not replicate the cortical accessibility reversal; signs flip with the adjustment. RM19 parks.
- `figures atac_hippocampus` (RM19): the cortical accessibility test repeated in 40 hippocampi (GSE278576), with genes as features because peaks differ per donor.
- R25 done (RM16): clonal-haematopoiesis carriers keep the shield (null, diluted clones).
- `fetch chip` and `figures chip` (RM16): blood of clonal-haematopoiesis carriers (any, DNMT3A, TET2) against non-carriers in GSE315366, scored like the syndromes.
- R24 done (RM19): in cortical cell types essential promoters spread apart with age more than matched promoters, after quality and composition adjustment; no general accessibility shield.
- `figures atac_cortex` adds the donor's promoter read fraction (a TSS-enrichment-like quality score) as an adjustment. R23 done (essential promoters spread more with age in cortical cell types, before that check); R24 planned.
- `figures atac_cortex` (RM19): accessibility spread with age per cell type in the 362-donor prefrontal-cortex snATAC set (Catching 2026), from a streamed per-donor pseudobulk. Dependencies: h5py, openpyxl. Roadmap: RM21 done and removed; RM12, RM14, RM15 ordered after RM19.
- `figures species_phylo` adds a drop-one-species jackknife: blood's lean toward a weaker shield in long-lived species rests largely on human.
- R22 done (RM20): long-lived species drift far less at unprotected CpGs (rate ~ lifespan^−1.4), but the shield's relative size holds or grows; catch-up is not supported.
- R20 done (RM18): probe design in the strata shrinks blood's gaps by about a sixth, the Zhou mask changes nothing, no gap changes sign.
- `figures species` records each group's unscaled gap and look-alike drift rate, and regresses the log rate and the relative gap on lifespan beside the gap (`outcome` column), for RM20's catch-up test.
- `fetch atac` and `figures atac_variance` (RM19): between-person accessibility spread with age per ATAC peak, essential promoters against matched promoter peaks, on Morandini 2023 PBMC. First result: a shield of a fifth of methylation's size once cell mix is removed.
- R19 done (RM9): unbiased bootstrap intervals for tissues, tissues_cells and mouse; intervals widen, no human conclusion moves, the mouse housekeeping contrasts hold. RM9 leaves the roadmap.
- RM18: the 450k annotation carries the Infinium design type (`probe_type`); `fetch annotation` also fetches Zhou et al. 2017's probe mask; `figures probe_robustness` recomputes the headline gaps with probe type in the strata and masked probes left out.
- `fetch taxonomy` and `figures species_phylo` (RM10b): the shield against lifespan with a phylogenetic GLS (NCBI Taxonomy covariance, Pagel's λ) and per tissue. No lifespan effect survives. Roadmap RM18 (normalization and probe robustness) at antonkulaga's request.
- Roadmap: RM10b (species with phylogeny, per tissue), RM16 (selection by cellular fitness in vivo: clonal haematopoiesis), RM17 (partial reprogramming), at antonkulaga's request.
- R18 done: carriers of Sotos, Weaver, Kabuki and CHARGE variants drift more in blood but keep the shield, unlike Tet3 heterozygous mice. `docs/extension.md`.
- RM9: `figures tissues`, `tissues_cells` and `mouse` scale every bootstrap resample by the full sample's spread (`common.rescale`), not its own, which had shrunk resampled gaps and blood-minus-brain intervals toward zero.
- `figures syndromes` keeps, where a series gives ages, only people with one inside the controls' age range, so age enters as a covariate without extrapolating to infants (R17 had dropped age for all of GSE97362).
- `fetch syndromes` and `figures syndromes` (RM14): blood of people with constitutional variants in chromatin writers (NSD1, EZH2: GSE74432; KMT2D, CHD7: GSE97362), each scored against the series' controls, raw and with 5 composition components. `arm_scores` takes `components`.
- R16 done: under Tet3 loss essential CpGs gain methylation (S +0.012, +0.014), consistent with TET's role; S rises in other arms too, so not specific yet.
- Per-sample scores gain `S`, the signed shift: above 0 when set CpGs gained methylation relative to look-alikes. `figures mmc_interventions` reports its effect beside `B` and `F`.
- R15 done: against same-chip littermates the Tet3 heterozygous brains lose the shield more clearly (housekeeping F +0.17, +0.25); the GHRKO liver effect was batch. `docs/extension.md`.
- `figures mmc_interventions` also scores each arm against only the controls on its own array chips, with the chip as a covariate (`controls` column): littermate controls where a study shares chips, which a batch cannot separate.
- R12 and R13 done. Across 65 mammal species the shield is nearly universal and not stronger in long-lived species. CR does not weaken it in mouse liver. Tet3 heterozygous mice lose it in brain, a causal lead whose controls are pooled across studies. RM10 leaves the roadmap. `docs/extension.md`.
- R14 done (41 s): RM15, SvenDP's question. Aging methylation is not all stochastic: in blood, solo-WCGW CpGs lose methylation 2.7 times as fast as other open-sea CpGs (−0.21 against −0.08 points per decade), a directed, division-linked loss absent in cortex; the PCA does not single it out. `docs/extension.md`.
- RM13 (1), first answer: what is left of the shield beyond chromatin. Little for core-essential genes in blood, a small clear remainder for housekeeping genes, more in brain (`docs/extension.md`, the novelty review's addendum).
- R11 done (54 min, 10.6 GB): the mouse figures with corrected intervals; no conclusion moved. RM2's six-tissue result written up in `docs/extension.md`: the shield in all six tissues, its size and its growth with age ordered by neither by turnover.
- R14 (RM15) submitted first; the other queued jobs lowered in priority.
- RM15, first in the queue (SvenDP's question, antonkulaga: above everything): the 450k probe table gains solo-WCGW tags from each probe's flanking sequence; figure `intergenic` tests directed, coordinated loss at open-sea CpGs outside promoters (per-sample means against age; leading components, their age correlation and solo-WCGW weight) in Hannum blood and Jaffe cortex. Not yet run on the data.
- RM15 (directed demethylation outside promoters) filed, from SvenDP's question.
- The lead wakes on a person's post anywhere in #entropy-age (`wake_people`), at antonkulaga's request.
- `docs/extension.md` and the artifact (v16): why the CR result is surprising, its absolute values, and two cautions (another strain's young reference; no slowing of overall drift).
- The draft's WP1 and the artifact (v15) take the natural experiments first: syndromes, knockouts, edited cells; no cell lines.
- RM14 (natural experiments in chromatin) filed at antonkulaga's question on cell lines; the plan moves the causal test to public data first.
- R13 (mouse interventions across tissues, RM12) submitted after R12.
- RM12: figure `mmc_interventions`, the mouse intervention studies of the Mammalian Methylation Consortium (CR in liver, GHRKO in brain, kidney and liver, partial reprogramming, sleep disruption): each animal scored against its study's controls (leave one out for controls), the intervention's effect on drift and on the shield's share. Not yet run on the data.
- The research programme is named apart from the paper's working title (plan, draft, artifact v13), at antonkulaga's request.
- RM13 (answering the prior work) filed; RM4 gains a GTEx check on fast deaths only, at antonkulaga's note on health at death.
- `docs/NOVELTY_REVIEW.md`: the core idea is published qualitatively (Teschendorff, West & Beck 2013; West et al. 2013); our novelty is design, the chromatin split and the dynamics. Three prior results cut against us (Hollmann 2020, Crofts 2026, Seeboth 2020). The paper draft takes the position in.
- RM12 (selection pressure under interventions, urgent), at antonkulaga's request; the artifact's tabs say which paper is ours (v11).
- `docs/PAPER_DRAFT.md`: a first paper draft at antonkulaga's request, every claim marked solid, weak or missing, with what the paper still needs and a DFG proposal outline; the summary artifact gains a Paper tab with a graphical abstract (v9), after a plain-language rewrite (v8) checked by a simulated outside reader.
- R12 (shielding across mammals, RM10) on the run brief and submitted.
- RM10: `entropy-age fetch mammals` (the Mammalian Methylation Consortium, GSE223748: 15,043 samples, 349 species, about 4.5 GB) and figure `species`: the aging gap per species and tissue (at least 15 animals spanning 15 % of the species' lifespan), related to maximum lifespan with adult weight and tissue held fixed. Not yet run on the data.
- RM10 (shielding against species lifespan, urgent) and RM11 (clocks of our own), at antonkulaga's request.
- RM9 filed: the same interval bias in `figures tissues`, to fix and rerun.
- R10's figure failed (tables written); R11 reruns the three mouse figures with corrected intervals.
- Fix: the aging-gap intervals of `interventions`, `bxd` and `rrbs_tissues` were biased toward zero. Each resample divided its gap by its own spread across CpGs, which resampling inflates. They are now taken on the unscaled gap and divided by the full sample's spread (`common.scaled_intervals`). The same bias makes `figures tissues`' resampled blood-minus-brain intervals conservative; it is filed as a follow-up. R10's figure failed on it (a point outside its interval); its tables stand.
- R10 (six mouse tissues, RM2) on the run brief and submitted.
- RM2: `entropy-age fetch rrbs-tissues` (mouse RRBS of six tissues, GSE120132, BS_Seeker2 CGfinal tables on mm10) and figure `rrbs_tissues`: per tissue, the aging gap with strain and sex as covariates and the shield's share per mouse against age, ordered by turnover (blood; lung, kidney, liver, adipose; muscle). Not yet run on the data.
- R9 done (40 s, 1.0 GB): `data/output/bxd_liver/` and `docs/figures/bxd_liver.png`; R9 leaves the run brief.
- RM3, BXD liver: the shield's share grows with age in a slowly renewing tissue too (housekeeping F −0.055 per 10 months with composition components, p 10⁻⁴), so that growth is no evidence of selection among stem cells; the high-fat diet's effect is the cell mix. The plan returns to chromatin with a growing share in vivo; the CR result in blood stays an open lead.
- R9 (the shield in BXD mouse liver, RM3) on the run brief and submitted.
- RM3: figure `bxd`, the BXD mouse liver of GSE199979 (321 females, 45 strains, chow or high-fat diet): the aging gap with strain as covariates, and per-mouse drift and shield share against age and diet, strain held fixed. It tests whether the shield's growth with age needs a renewing tissue. The Cole et al. liver WGBS series (GSE89273, GSE89274) have about 4 mice per arm, below the plan's 8, and are not used. Not yet run on the data.
- RM3, GSE80672: mouse blood repeats the human shield with exact ages; in living blood the shield's share grows with age (F −0.05 per 10 months); CR weakens it by a quarter in the within-batch B6D2F1 comparison (+0.025, p 0.002 with composition components), without slowing drift; the C57BL/6 CR comparison crosses sample series and is not read. The plan's claim becomes two-part; RM3 stays first, to replicate it.
- R8 done (37.5 min, 10.6 GB peak): `data/output/interventions/` and `docs/figures/interventions.png`; R8 leaves the run brief.
- R7 done (53 min, 19.7 GB peak for the parse); R8 (the composition check) submitted.
- RM3: `figures interventions` also scores every mouse with five composition components of the young reference removed, and marks whether treated and control mice share a sample series (the C57BL/6 CR mice do not).
- R7 (interventions in mouse blood RRBS, RM3) on the run brief and submitted.
- RM3: `entropy-age fetch rrbs` (mouse blood RRBS of GSE80672, its per-sample tables on GRCm38, with UCSC mm10 NCBI RefSeq promoters, Entrez links and CpG islands) and figure `interventions`: the aging shield in standard C57BL/6 blood with exact ages, and per-mouse drift and shield share under caloric restriction (two strains), GHRKO and Snell dwarfism, each within its strain. Not yet run on the data.
- RM8 done: no chromatin class gives a single-tissue clock that transfers (out-of-tissue r 0.06-0.70); which pool transfers best depends on the direction; Polycomb CpGs give the most reproducible brain-minus-blood offset (SD 4-6 years between draws against 20-42). RM8 leaves the roadmap; RM3 is next.
- R6 done (39.4 min, 2.8 GB peak): `data/output/chromatin_clocks/` and `docs/figures/chromatin_clocks.png`; R6 leaves the run brief.
- The kit pin moves to v0.1.178. RM3, RM2, RM4 and RM5 are chained after RM8 in the plan's order (each needs a raw-data processing step or a search first; `docs/RESEARCH_PLAN.md`).
- RM8, first part: published clocks draw on CpGs whose TSS activity differs between blood and cortex (15-33 % of a clock's CpGs against 12 % of background) and on Polycomb states. The pan-tissue Horvath and PhenoAge CpGs are enriched in housekeeping promoters (22 % and 17 % against 9 %). Clock CpGs are not the fastest-diverging, except epiTOC1 in blood. A clock's brain-minus-blood offset does not follow its chromatin share (Spearman 0.11). `docs/extension.md`.
- R5 done (11.8 min, 3.2 GB peak): `data/output/clock_cpgs/` and `docs/figures/clock_cpgs.png`; R5 leaves the run brief.
- R6 (clocks from protected and unprotected CpGs, RM8) on the run brief and submitted.
- RM8, second part: figure `chromatin_clocks`, elastic-net age clocks trained on equal-sized pools of CpGs that are active TSS in both blood and cortex, Polycomb or bivalent in both, or random, in Hannum and in Jaffe, tested on the other cohort and on Lunnon's same donors (brain-minus-blood offset). Not yet run on the data.
- RM7 done: in brain, neuron and glia accessibility (GSE96614) removes 30-37 % of the housekeeping gap, 32-49 % with the bulk state, against 5-26 % for the bulk map; the residual is then about blood's. RM7 leaves the roadmap; `docs/extension.md` and the plan record it.
- R4 done (27.8 min, 5.2 GB peak): `data/output/chromatin_match/` and the figure rewritten with the sorted-ATAC matchings; R4 leaves the run brief.
- R5 (published clocks on the protection map, RM8) on the run brief and submitted.
- RM8, first part: `entropy-age fetch clocks` (16 published clocks from biolearn at a pinned commit, each file hash-checked: chronological, phenotypic, pace, mitotic, stochastic and causal) and figure `clocks`: each clock's CpGs on the protection map (chromatin state in blood and cortex, gene-set shares, drift against matched CpGs with and without chromatin, and the clock's brain-minus-blood offset in Lunnon's donors). Not yet run on the data.
- R4 (the chromatin-matched gap with sorted brain chromatin, RM7) on the run brief and submitted.
- RM7: brain matched at cell-type resolution. `fetch chromatin` also takes the sorted-nucleus ATAC-seq of GSE96614 (NeuN+ and NeuN- prefrontal cortex), and `figures chromatin` adds, in brain, matchings on neuron and glia accessibility, alone and with the bulk state. Not yet run on the data.
- RM7: chromatin accounts for blood's protection. Matching on blood's own chromHMM states removes 67-91 % of the essential-gene gap (housekeeping dispersion slope −0.289 → −0.034 SD in Hannum), brain's states only about 30 %; after matching, blood is no more protected than brain. A bulk cortex map removes 5-26 % in brain. `docs/extension.md`; RM7 continues with cell-type-resolved brain chromatin; the plan records the change.
- R3 done (20.5 min, 5.2 GB peak): `data/output/chromatin_match/` and `docs/figures/chromatin_match.png`; R3 leaves the run brief.
- R3 (the chromatin-matched gap, RM7) on the run brief and submitted.
- RM7: `entropy-age fetch chromatin` (Roadmap chromHMM 15-state segmentations of blood mononuclear cells, neutrophils and prefrontal cortex) and figure `chromatin`, which matches gene-set CpGs to background on chromatin state as well as island and level, with a gene-resampling interval for the gap and for the share chromatin removes. The 450k probe table gains positions (`pos`, GRCh37): rerun `entropy-age fetch annotation`. Not yet run on the data.
- `docs/RESEARCH_PLAN.md`: the plan to a paper, at antonkulaga's request. The candidate claim moves from selection to chromatin; two new items, RM7 (does chromatin state account for the gap?, first) and RM8 (clocks from protected and unprotected CpGs); the roadmap reordered RM7, RM3, RM8, RM2, RM4, RM5.
- RM6 leaves the roadmap, as the plan grammar wants of finished items; what is left of it (replicate arrays, outcomes) is in `docs/PER_PERSON_SCORE.md`.
- RM6 done: the per-person score reproduces the cohort result one sample at a time (housekeeping `P` −0.0018 per decade in blood, −0.00036 in cortex), with reliability 0.96-0.997 against its own error, and blood's relative gap below the same donor's brain in 86-91 % of Lunnon donors. It does not separate selection from chromatin. `docs/extension.md`; `docs/PER_PERSON_SCORE.md` notes the two choices that differ from the design.
- R2 done (13 min, 5.9 GB peak): `data/output/per_person/` and `docs/figures/per_person.png`; R2 leaves the run brief.
- R2 (the per-person score, RM6) on the run brief and submitted.
- RM6: the per-person protection score is built: `entropy_age.protection` (a young reference's level, strata and optional composition components projected into each sample; the matched set-minus-background gap `P`, relative gap `F`, entropy gap `E`; a per-sample error from resampling genes), with figure `per_person` for the design's validation steps 1-3 on Hannum, Jaffe and Lunnon. The culture test's `deviations` and `matched_scores` move there. Not yet run on the data.
- RM5: GSE179847 replicates GSE197512: the protection is a constant fraction of drift, not growing with divisions (housekeeping change per 10 PD −0.055, interval −0.21 to +0.10, five lines; paper scale, seed 42). Contact inhibition cuts drift but leaves too little to decide whether protection needs divisions. `docs/extension.md`, figure `passaged_cells_sturm`; R1 leaves the run brief.
- R1 (GSE179847) submitted; its first attempt's failure and cause in `docs/RUN_LOG.md`.
- RM5: GSE179847 (Sturm et al. 2022) as cohort `sturm` and figure `passage_sturm`, the culture test's second replication: six healthy lines, the lifespan parts' 21 % controls as baseline, treatments, 3 % oxygen and contact inhibition as arms. Its characteristics are read per sample (GEO shifts a key's row between samples), and the noise floor comes from the flasks split from one culture at a treatment's start (`culture`), since no arm repeats a day. Not yet run on the data.
- `docs/RUNS.md` entries carry a plan-grammar state line, and `[slurm] plan` names the file, so the account's queue table counts this repository's planned runs.
- `[slurm] job_prefix` is `entropy-` (was `entropy-age`): the repository's jobs are also named `entropy-rm*`, `entropy-mouse*` and `entropy-cells-*`, which the old prefix left out of the digest.
- The runner is seated: its share (CPU only, at most 16 cores, memory and time from measured runs) in `docs/seats/RUNNER.md`, the cluster's facts kept on the cluster machine; `docs/RUNS.md` gets R1, GSE179847 for RM5, waiting on the lead's layout.
- The roadmap's items RM2–RM6 carry a `**State:** · **Owner:**` line the plan watcher reads (all unattended, owner the lead); the old prose state line is now `Progress:`.
- Seats, second pass (D5, amends D3): a `runner` seat beside the lead (`docs/seats/RUNNER.md`), which keeps
  the cluster runs of `docs/RUNS.md` going and logs them in `docs/RUN_LOG.md` and `docs/RUN_STATUS.md`; the
  seat table gains it, the cluster queue in every session's digest and `[slurm] job_prefix`. Both antonkulaga
  and winternewt are maintainers now; the organisation and desk stay D3's. The kit's pin moves from v0.1.144 to
  v0.1.163 (`.claude/settings.json`). The lead's priming names both maintainers and hands runs too big for a
  session to the runner through `docs/RUNS.md`; the runner's share goes into its priming, and the queue block
  runs under bash (the kit documents `sh`, which the script refuses).
- GTEx muscle no longer read as a low-turnover control (`docs/extension.md`, GTEx section): satellite
  cells divide with loading, injury and wasting, and GTEx records no exercise or training, so muscle's
  turnover in its donors is unknown and its contrasts cannot test selection. `docs/ROADMAP.md` RM2 now
  asks for a muscle cohort with recorded activity and age.
- RNA results recorded (`docs/extension.md`, "RNA: transcriptional drift by turnover class";
  `docs/figures/rna_gtex.png`, `docs/figures/rna_rodents.png`): in mouse, housekeeping and core-essential
  genes gain spread with age more slowly than matched genes in every class but brain (renewing minus
  brain −0.28 SD, robust to composition), as in methylation; rat leans the same way; GTEx shows no
  turnover effect in spread. Both RNA figures share their panels now (`rna_gtex.tissue_panel`,
  `contrast_panel`; class summaries `class_contrasts`, `tissue_intervals`), show core-essential and
  housekeeping genes only (non-essential genes are often a handful per tissue; they stay in the tables),
  and no longer overlap their labels.
- Direct selection test recorded (`docs/extension.md`, "The direct test"; `docs/figures/passaged_cells.png`):
  in seven primary lines passaged to senescence, housekeeping CpGs drift about half as much as matched
  CpGs once drift clears the array noise floor, and that fraction does not grow with doublings (mean
  change −0.14 per 10 PD, interval −0.33 to +0.06; the longest courses flat). The first, uncorrected
  version's decline in 7/7 lines was the noise floor. Mitomycin C arrest stopped most drift, leaving too
  little to measure. The figure now shows drift, the bounded fraction, the coupling change for arms
  spanning 10 PD or more, and drift under arrest.
- GTEx and fibre results recorded (`docs/extension.md`, "Human tissues: GTEx and muscle fibres";
  `docs/figures/gtex_tissues.png`, `docs/figures/muscle_fibres.png`): muscle is among the least protected
  GTEx tissues, but breast, ovary and testis are as weak, so the split does not follow turnover; the
  between-donor spread gap vanishes with composition components; the fibre contrast is uninformative.
  Corrections: the mouse section no longer says protection is missing wherever cells do not divide
  (heart keeps it; only cortex and hippocampus lose it); age-free components cannot change an age slope
  of entropy or of the mean, so the composition check's lifespan line is by construction, not a test.
  The GTEx figure shows the composition check on spread (panel b) instead of on entropy.
- RNA test in rodents (roadmap RM4): `entropy-age fetch rodent-rna` (GSE314817 mouse, 24 tissues;
  GSE314970 rat, 28 tissues; 6-27 months, raw counts) and `entropy-age figures rna_rodents`, the
  `rna_gtex` measure with numeric ages, sex and library size as covariates, size-matched subsamples
  per species. Rodent symbols reach the human gene sets through HGNC by upper-case symbol
  (`annotation.symbol_ensembl`); two rodent symbols reaching one human gene are both dropped.
- Direct selection test (roadmap RM5): `entropy-age fetch cohort endicott` (GSE197512, seven primary
  human cell lines passaged to senescence on EPIC, with mitomycin C arrest, serum, oxygen and TERT
  arms; population doublings and days per sample) and `entropy-age figures passage`. Per sample, the
  deviation from the experiment's earliest sample at gene-set promoter CpGs minus that at matched CpGs
  (`docs/PER_PERSON_SCORE.md`), as a fraction of the matched deviation, regressed on doublings and on
  days. Selection predicts that fraction falls with doublings and not with time under arrest. A second
  score on squared deviations removes the array noise floor (per-CpG variance of replicate cultures
  measured on the same day), which otherwise pulls the fraction toward zero wherever deviations are
  small and makes it fall with doublings even without selection. Because that ratio explodes where the
  corrected drift is near zero, the main statistic is a regression of the gap on the corrected matched
  drift and its interaction with doublings (`drift_coupling`): constant coupling under chromatin
  protection, a coupling that strengthens with doublings under selection.
- Design note `docs/PER_PERSON_SCORE.md`: how one sample's protection gap could be measured (deviation
  from a young reference at set features minus matched background features of the same sample), what
  biases it and how it would be validated. Roadmap RM5 lists public methylation time courses of
  dividing cells (GSE197512 with arrested arms, GSE179847, GSE268145, GSE202554, GSE131280), RM6 the
  score.
- RNA test (roadmap RM4): `entropy-age figures rna_gtex` measures, in GTEx v10 RNA-seq (47 tissues with
  100 or more donors), how fast between-donor spread of each expressed protein-coding gene grows with
  age, and compares core-essential, housekeeping and non-essential genes with genes of the same
  expression decile, by turnover class (renewing, slow, heart & muscle, brain), in subsamples of 70
  donors per tissue. Counts are read from an external parsed copy named by `ENTROPY_AGE_GTEX_RNA`, and
  RIN, ischemic time and Hardy scale from GTEx's annotation tables named by `ENTROPY_AGE_GTEX_ATTRIBUTES`;
  either missing is a refusal. `entropy-age fetch annotation` now also fetches HGNC's complete set
  (release 2026-10-06) and puts the gene sets on Ensembl IDs (Entrez first, then current, previous and
  unique alias symbols: 99-100 % of each set matched). New metric `stochasticity.expression_drift`.
- Human tissue test (roadmap RM2): `entropy-age fetch cohort gtex` (GSE213478, nine GTEx tissues of 424
  donors on EPIC, decade age bins as midpoints) and `entropy-age figures gtex`, which asks whether
  skeletal muscle, the one post-mitotic tissue, is where the protection is missing. Every tissue is
  resampled at the same size (32 donors, without replacement), because fewer donors shrink SD-scaled
  differences toward zero. `entropy-age fetch cohort fibres` (GSE304045, type I and type II fibres and
  whole muscle of the same people, EPIC v2) and `entropy-age figures fibres` compare between-person
  spread in fibres against whole muscle.
- Mouse results recorded (`docs/extension.md`, `docs/figures/mouse_tissues.png`): in 11 tissues of
  the same mice, essential and housekeeping CpGs gain entropy 0.12–0.16 SD less than matched CpGs in
  renewing than in post-mitotic tissue (8 of 9 aging contrasts exclude zero), slow tissues are as
  protected as renewing ones, and cortex and hippocampus show none. The rapamycin arm (5 mice per
  tissue) is within small-sample noise against size-matched controls.
- Composition check recorded (`docs/extension.md`, `docs/figures/tissue_contrast_cells.png`): the
  same-donor blood-minus-brain gap keeps all 12 intervals on the predicted side under every adjustment;
  the lifespan-slope gap holds with age-free components but not with components that keep their age
  trend, so it is now presented as supporting evidence only. Roadmap RM1 closed.
- Mouse multi-tissue test (roadmap RM2/RM3): `entropy-age fetch mouse-annotation` (the mammalian
  array's mm10 probe table, GEO GPL32543, and MGI mouse-human orthologs that carry the human gene sets
  to mouse), `entropy-age fetch cohort huang` (GSE279437: 11 tissues of the same mice, young and aged,
  exenatide and rapamycin), and `entropy-age figures mouse`. The figure compares the gene sets' aging
  drift and spread across tissues ordered by cell turnover, with bootstrap contrasts between turnover
  classes and between treatments.
- Roadmap RM1–RM4: record the composition check, then a turnover gradient across tissues, long- vs
  short-lived rodent groups, and stochasticity in RNA, each with the public datasets that would test it.
- Composition check for prediction 2: `entropy-age figures tissues_cells` reruns the blood-versus-brain
  test with the leading principal components of each tissue's background CpGs as extra covariates
  (5 or 10 components, keeping or removing their age trend; `composition_components` in
  `stochasticity.py`), the same reference-free correction in every tissue. `ENTROPY_AGE_BOOT` overrides
  the number of bootstrap resamples, for smoke tests.
- Prediction 2 tested (`docs/extension.md`, `docs/figures/tissue_contrast.png`): the essential and
  housekeeping protection is larger in blood than in brain, and the non-essential excess too. This
  holds across Hannum blood versus Jaffe cortex (age slopes) and within the same Lunnon donors' blood
  and four brain regions (dispersion), with 20 of 21 bootstrap intervals excluding zero. Blood's
  unadjusted cell mix remains a possible partial explanation.
- Fixed: the tissue experiment's bootstrap workers deadlocked when forked (polars' thread pool was
  running in the parent; 16 workers slept after 19 CPU-seconds each on 2026-11-09). They are now
  spawned and receive their data through an initializer; this replaces the "explicit fork" entry
  below. A test checks that parallel and serial bootstraps give identical results.
- `read_supplementary_betas` skips a preamble and ID rows after the header, which GSE59685's beta file
  needs (comment lines, a barcode header, then two rows of GEO accessions).
- `compare_groups` finds group members with a hash lookup instead of `np.isin` on strings, which took
  seconds per group at 10^5 features (a 17 s profile now takes 2 s). The tissue experiment adjusts Hannum
  blood for recruitment site as well as plate, and starts its bootstrap workers with an explicit fork.
- Prediction 2 pipeline (essential-gene protection in renewing versus post-mitotic tissue), results
  pending: `entropy-age fetch cohort hannum jaffe lunnon` (450k blood and brain cohorts, GB-scale, for a
  cluster), the 450k manifest in `fetch annotation` with promoter CpG-to-gene links, an age-free
  `dispersion_levels` measure, and `entropy-age figures tissues` with bootstrap blood-minus-brain
  intervals. GEO characteristic keys are now cut before any parenthesis (`age (in years)` -> `age`).
- `docs/extension.md` gains a related-work section, checked against PubMed and Meyer's 2024 thesis.
  No stochastic clock has yet been validated on mammalian RNA (the only application is Ostermann et
  al., JCI Insight 2025, on human induced neurons). The worm neuron-type paper is the closest
  cell-type analogue, and the gene-level noise literature (CGI− genes, chromatin variability,
  sub-Poissonian genes, CpG density and lifespan) is described from each paper's own abstract.
- Prediction 1 of the selection hypothesis tested on real data (`docs/extension.md`): `entropy-age
  fetch annotation` (Illumina 27k manifest from GEO GPL8490; Hart CEGv2 and NEGv1 pinned at
  hart-lab/bagel 53388ad; Eisenberg–Levanon housekeeping genes) and `entropy-age figures genesets`.
  In GSE41037 blood, essential and housekeeping CpGs gain entropy and between-person dispersion
  more slowly than CpGs matched on CpG-island status and starting level, and non-essential CpGs
  faster; this holds after batch adjustment, with cell fractions, and inside each processing batch.
- `drift_metrics(covariates=...)` gives partial age slopes; `compare_groups` takes a categorical
  stratum column and runs its stratified permutations vectorised.
- Fixed: GSE41037's three processing datasets are confounded with age. The atlas now compares clock
  CpGs using batch-adjusted slopes, and the earlier claim that clock CpGs gain entropy faster than
  level-matched CpGs (Horvath p = 0.004) is withdrawn (batch-adjusted p = 0.60). Fig. 5d's title
  said methylome entropy rises with age; the pooled r is 0.05 (within batches 0.17, 0.00, 0.01), and
  the panel now says so per batch. Fig. 5 also writes the clock's accuracy within each batch
  (r = 0.75, 0.86, 0.51). A substituted ground state is now announced on stderr.
- Adopted the seats kit (agent-seat-blueprints v0.1.144): `CLAUDE.md` and `AGENTS.md`, the roadmap,
  changelog and decisions files, the `lead` seat's priming under `docs/seats/`, the seat table
  `.blueprints.toml`, and the kit's plugin pinned in `.claude/settings.json` (D1-D4). Not yet
  run: no seat has been taken; the Zulip organisation and the desk do not exist yet.

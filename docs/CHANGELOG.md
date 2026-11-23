# Changelog

Newest first.

## Unreleased

- 2026-11-24: research plan gains a ranked list of ideas to pivot to (clock atlas over every rebuilt series, what each clock generation reads, shield as a mitotic clock, clock disagreement, preprocessing moves clocks, tissue-invalid clocks as controls); R165 scores the clock panel on 47 rebuilt series.
- 2026-11-24: new figure `clock_panel` (figures review section 2): clocks clustered by agreement of their accelerations with a generation strip, each clock against B, F and P, odds ratios of each clock and of F given it, and per-woman clock disagreement against drift.
- 2026-11-24: figures name cohorts by tissue and accession (plots.save relabels "Hannum blood" to "blood (GSE40279)", Jaffe to GSE74193, Lunnon to GSE59685): "Hannum" on Figures 1, 2, 4 and 6 was the blood cohort and read as the only clock (Anton; figures review).
- 2026-11-24: clock_increment labels every clock by generation (1st chronological, 2nd phenotypic and mortality, 3rd pace, mitotic and stochastic) and writes best_per_generation; the disagreement panel keeps only clocks with r > 0.5 with age (Bocklandt ran against age); per-woman accelerations and disagreement saved. Fable figures review in docs/reviews/2026-11-24-figures-fable.md.
- 2026-11-24: clock panel takes GrimAge v1/v2 (age and sex from the records; a missing value filled with the CpG's cohort mean, since biolearn's GrimAge counts it as 0) and a CLOCKS_ONLY filter; a clock track rebuild (`ENTROPY_AGE_SESAME_DIR=sesame_clock`, SeSAMe CDPB without the design-flag mask), because the mask removes 20 % of GrimAge's weight on every array.
- 2026-11-24: Figure 3 redrawn around its message (Anton: the scatter read as "messy between tissues"): gap per species within each tissue with share below zero, all groups as a histogram (88 % below zero), lifespan coefficients; caption rewritten.
- 2026-11-24 (pipeline plan step 1): shared runner writes a guards table per series (arrays selected, passing detection, kept per donor, repeats dropped, platforms, complete CpGs, chips, unknown chips, age variance on chip, identifiability with chip).
- 2026-11-24: Box 1 figure: schematic from one CpG to a gap (scripts/box1_schematic.py).
- 2026-11-24: Introduction now explains the measure, matching, the gap and its sign before any result; every main figure caption ends with a how-to-read line; captioned figures rebuilt (Anton).
- 2026-11-24 (P1 complete): monocyte expression-matching share 77-91 % over people (gap −0.12 to −0.08); CD4 T-cell gap unresolved over people (−0.05 to +0.02).
- 2026-11-24: Box 1 at the start of Results: the four drift measures in plain words (what each counts, what a higher slope means, when to use which) and how to read a gap in every figure (Anton).
- 2026-11-24: paper: expression-matching shares now carry donor/culture intervals (fibroblasts 68-76 %, epidermis 64-80 % and 60-84 %).
- 2026-11-24 (Codex review of P4): clock columns come from the biolearn panel only (the join had read older columns); exact-chip matching added as a sensitivity design; failed bootstrap fits left out and counted; ratio F checked against the absolute gap P; disagreement on an explicit adult-blood panel with accelerations residualised on age, cell fractions, mean beta and chip; coverage counts biolearn's wholesale row imputation as unmeasured; skipped model types listed.
- 2026-11-24 (protocol P4): clock_increment module (nested conditional logistic models per kept clock: clock, + B, + F given B; panel mean + disagreement). Quick run: B adds nothing to PhenoAge (LR p 0.82, primary); F given B OR 0.70 (p 0.02), alike across the panel; disagreement correlates with B (r 0.37). R162 at paper scale.
- 2026-11-24: R161 fibroblast clones written up (gap inside single-cell clones −0.08, −0.06); after Codex's review the clone endpoint is reported as not separating maintenance from selection (protocol change 3).
- 2026-11-24 (protocol P5): fibroblast_clones module (Flunkert clones, rebuilt): per-clone relative entropy from the starting culture; quick run: housekeeping gap inside clones −0.08 and −0.06 (clone resampling), technical floor (parent array vs array) ~0. Caveat for review: a clone inherits its founder cell's own difference from the parent.
- 2026-11-24: the analysis protocol is living, not frozen (Anton): docs/ANALYSIS_PROTOCOL.md; a method change recomputes every result it affects and is recorded with its reason. Clock panel: every biolearn gallery clock with at least 90 % of its weight measured in 90 % of arrays, plus a per-person disagreement score (change 2).
- 2026-11-24 (protocol P4): clocks scored with biolearn itself at the pinned commit and its documented imputation (scripts/clocks_biolearn.py, own environment), with the measured share of each clock's weight per sample; PhenoAge primary (biolearn's DunedinPACE lacks the authors' normalisation). R160 scores the rebuilt EPIC-Italy arrays. P1 donor intervals: share removed by expression matching 68-76 % (fibroblasts, cultures resampled), 60-84 % / 64-80 % (epidermis EPIC / 450k, donors resampled).
- 2026-11-24 (protocol P3): forest plot and table of every series (scripts/forest_series.py; 42 cohort-tissue series, 31 negative, 20 resolved below zero, 1 above: tumour-adjacent lung); Supplementary Fig. S22.
- 2026-11-24 (protocol P2): division contrasts checked for donor/animal resampling; off-floor colon/hippocampus numbers updated to the current table; muscle stem-cell interval labelled conditional on six mice with its leave-one-mouse-out range.
- 2026-11-24 (protocol P1): donor resampling for the expression-matched gaps (`sample_resampled_gaps`: per resample the slopes refitted, strata and SD from the full sample; intervals for each matching and for the share removed); MESA monocytes and epidermis resample people, Sturm fibroblasts cultures within line.
- 2026-11-24: analysis protocol frozen (docs/ANALYSIS_PROTOCOL.md) before the Aging Cell analyses; remaining study-wide FDR wording removed from paper and supplement.
- 2026-11-24: R158 written up: healthy lung −0.26 (−0.45 to −0.13), kidney −0.10 (unresolved, 30 people), mouse hearts −0.10 (−0.16 to −0.07) once MM285 islands are annotated, mouse field liver ~0 (young ages), BXD control through the program −0.28 as its own analysis.
- 2026-11-24: p-value audit of the Results (docs/reviews/2026-11-24-p-value-audit.md): 33 p-values, 25 valid (donor, animal or culture resampling; models across samples), 8 set-label tests replaced by donor intervals or marked "over genes, no donor interval"; stale species value (+1.43) replaced by the current table (+0.59, p 0.29); Tet1 F bound updated.
- 2026-11-24: shared runner tests identifiability by design rank (age must raise the rank of the nuisance design; residual df from rank and distinct donors), in the full fit and every resample, replacing column counts (Codex review).
- 2026-11-24 (reviews): set-label p-values (CpG permutation, gene clusters) no longer used as evidence: they reach small values with ages shuffled; headline gaps now quoted with donor intervals, study-wide FDR passage withdrawn. Lung exclusion test rewritten as a sensitivity check, healthy lung (GSE175458: −0.26, −0.43 to −0.15) shows the reversal is specific to tumour-adjacent tissue. Field-liver null explained by its young ages (BXD ≤10 months also −0.01). Figure 3 coefficient labelled per natural-log unit (+0.071 per tenfold).
- 2026-11-24 (reviews by Codex and Claude): MM285 annotation rebuilt from Zhou's manifest and UCSC islands (island flags were 4 % of promoter probes in the GPL32543 table; now 25 %); runner normalises 40k probe ids (`strip_probe_suffix`) and stops with an error when a group has no measured CpG (the BXD control had passed silently with none); QC reader keeps the same first array per sample as the betas, so detection and chip describe one physical array.
- 2026-11-24: Methods open with an "Analysis practices" subsection (raw IDATs, one recipe, one platform per estimate, one array per donor, complete CpGs, chip second fit, resampling rules, gene-level p, shared program, placebo sets, independent review, light clock path); outdated interval sentence corrected (Anton).
- 2026-11-24: rebuild keeps one array per sample when a GSM holds several (GSE79100: three per kidney); the reader skips a sample already read from an earlier chunk.
- 2026-11-24: R154 two healthy breast cohorts from raw arrays: little gap (−0.04 and +0.00), in line with GTEx breast.
- 2026-11-24: R152 mouse hearts and field/lab livers show no resolved gap (liver −0.03, −0.09 to +0.03); R153 the lung reversal shrinks without TCGA tumour-changed CpGs (+0.14 → +0.11 → +0.08) but does not vanish. Paper, run log and S5 updated.
- 2026-11-24: tumour-field test for the lung reversal: CpGs that differ between 80 TCGA lung adenocarcinomas and 32 normal lungs (open GDC betas; 86,809 at |Δβ| ≥ 0.1, 24,287 at ≥ 0.2) left out of the tumour-adjacent lung series (R153). Mouse hearts and field/lab livers through the shared runner (R152).
- 2026-11-24: shared runner takes mouse and mammalian arrays (`annotation = "mouse"`) and age-group labels (`age_values`); first mouse series: GSE280472 hearts (young/old, both sexes, saline or elamipretide; treatment held), waiting for its rebuild.
- 2026-11-24: paper: GTEx values named as deposited (raw arrays controlled); the lung reversal placed against healthy GTEx lung and named as tissue next to a tumour, its cause under test (Anton).
- 2026-11-24: Figures 1, 6, 8 and 9 redrawn with the house style (bounding box and resolution; earlier redraws outside the command line had skipped it); Figure 8 labels in plain words.
- 2026-11-24: Figure 9 relabelled: plain odds-ratio ticks (the log axis printed overlapping 8×10⁻¹ labels), readable measure names, a legend in place of the long title; caption names the colours.
- 2026-11-24: Figure 1 relabelled: plain measure names (relative entropy, not "kl slope"), statement title, panel titles without "Prediction 2", legends off the data; caption rewritten.
- 2026-11-24: R150 complete: GTEx tissues rerun after the relative-entropy fix (whole blood −0.14 → −0.13, colon −0.19 → −0.18, lung and kidney −0.09 → −0.08); all 15 analyses computed before the fix are now rerun and their numbers rewritten in paper, supplement and artifact.
- 2026-11-24: R149 species rerun after the relative-entropy fix: 88 % of 139 groups negative (mean −0.12; was 89 %, −0.13), lifespan coefficient +0.031 (p 0.04); Figure 3 redrawn on relative entropy; R151 deposited SATSA gene-level p 0.68. Paper, supplement and artifact rewritten.
- 2026-11-24: SATSA within-person p-values now gene-level (rebuilt series: direction p < 10⁻¹⁰, divergence p 0.03); deposited tables rerun (R151) for the same.
- 2026-11-24: R150 part 3: Meer mouse tissues rerun (lung −0.14 → −0.15 in the text), species floor unchanged.
- 2026-11-24: R146 headline rerun: blood housekeeping −0.147 SD, gene-level p 5·10⁻²⁴⁴ over 3,647 genes; Methods now describe gene-level p. R150 part 2: interventions, sorted cells, probe robustness unchanged; Thompson mouse RRBS tissues moved after the relative-entropy fix (adipose −0.03 → −0.08, kidney −0.18 → −0.23, muscle −0.08 → −0.10); text rewritten.
- 2026-11-24: R150 (first eight of the fifteen analyses computed before the 2026-11-18 relative-entropy fix) rerun: gaps move by at most 0.015 SD; BXD liver −0.32 → −0.30, SATSA cross-sectional −0.23 → −0.22 in the text; hotspot enrichment in blood bivalent promoters 1.42 → 1.36 (z 10.1 → 8.6).
- 2026-11-24: `write_table` deletes a CSV it no longer writes (a table grown past 1,000 rows), so no stale copy is read as current; the stale `species_gaps.csv` (12 Nov, no relative-entropy rows) and `mmc_interventions/animal_scores.csv` removed (review item 11). Figure 3 plots the relative-entropy slope with marked species in place of overlapping labels; R149 reruns the species analysis, whose tables predated the relative-entropy fix of 2026-11-18.
- 2026-11-24: Figure 6 redrawn on relative entropy (it showed the dispersion slope), one colour and legend entry per matching (blood's neutrophil map had shared "+ own state" with the tissue's own map), shared legend, statement title; caption rewritten.
- 2026-11-24: R147 myometrium written up (no housekeeping gap in 62 women 29-52 y; EPIC v2 unresolved); R148 rerun of every series, results unchanged; coverage table S5 gains myometrium.
- 2026-11-24: shared runner gives no donor interval when fewer than half the donor resamples leave enough degrees of freedom, and reports how many were kept (`donor_resamples`): EPIC v2 myometrium (30 women, 10 chips) had produced chip-held intervals 0.01 SD wide from the few surviving draws. Every series rerun (R148).
- 2026-11-24: shared runner reads columns from sample titles (`[derive]`); repeated identical characteristics accepted (GSE296589 lists `batch` twice). R147: myometrium (GSE296589), EPIC and EPIC v2 parts as two series, condition and ancestry held, fibroid-free specimens also alone.
- 2026-11-24: gene-level p (review item 9): `compare_groups` adds `p_gene`, the matched difference tested with genes as units (set residuals summed per gene, cluster-robust SE, t on genes - 1 df), wherever the membership carries genes; the global FDR uses it in place of the CpG-label permutation p. R146 reruns the headline blood and cortex comparison.
- 2026-11-24: Figure 5 redrawn: neurons of the second series with their donor-resampling interval, the sequencing rows (oligodendrocyte lineage and neurons of the same brains), culture rows with culture-resampling intervals, the glial row labelled as the NeuN−/SOX10− fraction; caption names each interval type (review item).
- 2026-11-24: R145: exercise training barely moves promoter methylation beyond array noise (about a hundredth of replicative senescence); the share essential promoters carry after 4 weeks (−0.22) resembles that over a period without training (−0.15); loading and unloading too small to read. Paper, run log, Table S5.
- 2026-11-24: figures muscle_training: vastus lateralis before and after 4, 8 and 12 weeks of interval training and across a control period without training (GSE171140), and after one bout, loading, unloading and reloading (GSE114763); per CpG methylation on condition with person and chip held; noise-removed squared change from baseline and the share essential promoters carry (as senescence_donors, whose estimator now takes its arms as a parameter); each person left out in turn (R145). Anton asked about trainable tissues.
- 2026-11-24: R144: by WGBS of the same 41 brains, neurons and the oligodendrocyte lineage both carry the housekeeping gap and do not differ (+0.027, −0.081 to +0.116); the glia-over-neuron result rests on one fraction of one array series. R140b: the single-cell gap holds with old cells thinned to the young depth (−0.085); both in the paper, run log and Table S5.
- 2026-11-24: musc_single_cell also writes gaps_depth_matched: the point estimates with old cells' reads thinned per promoter to the young depth (old cells were sequenced deeper, which lowers their measured variation; review item, R140b).
- 2026-11-24: figures wgbs_nuclei: NeuN+ against OLIG2+ nuclei of the same brains by whole-genome bisulfite sequencing (GSE107729, Mendizabal 2019), reduced to the 450k probe positions (GRCh37, offset checked), complete CpGs, donor bootstrap; every donor with diagnosis held, and controls alone (R144).
- 2026-11-24: R143: colon mucosa next to cancer, 146 patients: housekeeping gap −0.13 (donors −0.25 to −0.01); lung next to adenocarcinoma, 114 non-smokers: reversed, +0.14 (+0.05 to +0.25), the one resolved reversal in a human tissue; human tissues on the mammalian array: too few promoter CpGs, nothing resolved; obese four-tissue series fails the detection cut (47 of 268). Paper, run log and Table S5.
- 2026-11-24: the three-region brain series holds under donor resampling with chip held (prefrontal cortex −0.119 to −0.052, amygdala −0.098 to −0.047, hippocampus −0.094 to −0.043; controls too). Paper and run log.
- 2026-11-24: Supplementary Table S5 (docs/analysis_coverage.csv rendered by scripts/coverage_table.py): one row per methylation series, one column per analysis (gap, off the floor, chip or batch held, chromatin state, expression from the same samples, accessibility, gain and loss rates, what the interval resamples), each cell naming the run or saying not run. Table S4 remains the sample and modality table. The paper's Limits say where each follow-up analysis was run. Anton asked for tables of samples, analyses and modalities.
- 2026-11-24: claims ledger rewritten for C2, C3, C4, C6, C11 and C12 (headings and status) to what the evidence now supports, after the independent reviews; the old wording ("chromatin maintenance", "Tet3 guards", "not tuned to lifespan", "has the shield in blood") is gone.
- 2026-11-24: shared runner reads the source name (for series that record tissue only there) and a named series matrix (for series split by platform); description files for human tissues on the mammalian array (GSE184211: heart, muscle, liver, kidney, lung, spleen, fat, blood), muscle, liver and two fat depots of the same people (GSE61454), lung (GSE314841, 450k part) and colon mucosa (GSE101764).
- 2026-11-24: mouse tissues off the 0 % floor: colon keeps its gap (−0.25, −0.51 to −0.19), hippocampus has none (+0.01, −0.10 to +0.10); in the four-tissue series no tissue is resolved off the floor. Paper and run log.
- 2026-11-24: two independent reviews (Claude and Codex) agree that the division claim, the kinetics lead and the Tet3 conclusion were stronger than the evidence: the division lead and Abstract now name the specific contrasts and say tissues do not follow a turnover ranking; the Discussion lists the contrary evidence; kinetics stated as a pooled-fit result not seen per CpG in cortex; Tet3 as one bulk-tissue study; the CD4+ exception to the expression result is in the Abstract.
- 2026-11-24: independent review of claims (first run), corrections: off-floor gaps put on one scale (blood −0.10 at 10-90 %, 0.00 at 20-80 %; cortex −0.08 and −0.06 with intervals reaching zero); lifespan stated as not resolved instead of absent; arrested cultures as no drift detected in two donors; cortex chromatin share given as 25 to 44 % (deposited and rebuilt); the culture-bootstrap tables the paper quotes are now tracked.
- 2026-11-24: jena_rrbs and mouse_atlas_rrbs also write gaps_off_floor: the matched gaps on CpGs whose young level lies between 10 and 90 %, with their own strata and an animal bootstrap on a separate random stream (the existing tables are unchanged). First step of taking the non-blood, non-brain tissues through the follow-up analyses (Anton asked).
- 2026-11-23: R141b: mouse tissues on CpGs covered in every library (nothing filled): lung −0.09 (resolved), heart −0.05 and cortex −0.01 (not resolved), liver without estimate (518 CpGs); the paper sentence and run log replace the R141 numbers, which used filled values.
- 2026-11-23: reruns after the interval fix: neuronal gap with chip held −0.030, donors −0.102 to +0.037 (without chip −0.091 to −0.008); culture intervals −0.160 to −0.072 (telomerase fibroblasts), −0.128 to −0.052 (astrocytes), −0.089 to −0.080 (astrocytes with telomerase); paper and run log corrected.
- 2026-11-23: code review (independent agent) fixes, part 1: donor and culture resampling intervals were divided by each resample's own SD, which pulls them toward zero; they are now scaled by the full sample's SD (series.donor_intervals, astrocyte_passage.culture_bootstrap, musc_single_cell). The shared runner also counts distinct donors for the degrees-of-freedom guard, raises on a title without a donor or a missing donor id, refuses mixed array platforms, and uses one random stream per group and design. Tests: the donor interval brackets its point estimate; titles without a donor are refused.
- 2026-11-23: R141 (mouse liver, lung, heart, cortex of the same mice, lane held: housekeeping gap −0.10, −0.10, −0.08, −0.05) in the paper and run log. R142: the shared runner reproduces the neuron series to three decimals and adds donor-resampling intervals; with chip held the neuronal gap is not resolved (−0.047 to +0.022), and the paper now says so in the division section and the closing sentence.
- 2026-11-23: figures series_<name>: one runner for every array series analysed as an age slope, driven by a description file per dataset in config/series (accession, age, sex, donor, held columns, groups); all analytic decisions fixed in the runner (one array per donor, complete CpGs, chip as a second fit, off-floor band, gene-resampling and donor-resampling intervals for every gap, age-by-chip overlap). First two descriptions: neurons (GSE112179), brain_regions (GSE349373), run as a regression check against their own modules (R142). Anton asked for shared code so that per-series decisions cannot drift apart.
- 2026-11-23: figures mouse_atlas_rrbs: liver, heart, cortex and lung of the same mice at 3, 15 and 24 months by RRBS (GSE225166), sex and sequencing lane held, lanes mixing the ages; smaller tissues listed when they reach the minimum (R141).
- 2026-11-23: R139: with counting noise removed, between-donor variation of promoter accessibility shows no difference between essential and matched promoters in any of seven cortical cell types; the earlier opposite sign came from the log-count measure. Stated in the paper and run log.
- 2026-11-23: figures musc_single_cell: single muscle stem cells of young and old mice with RNA from the same cells (GSE121437): between-cell variance of promoter methylation beyond read sampling, its gain with age and the shift of the mean, set promoters against promoters matched on island status and young level, and with the expression decile added; cell bootstrap within mouse and each mouse left out (R140).
- 2026-11-23: figures atac_cortex: a second measure, the age slope of between-donor variation in excess of counting noise (squared coefficient of variation beyond Poisson), beside the spread of log counts, whose expectation depends on read depth; table gaps gains a measure column (Anton asked for a cleaner accessibility test, R139).
- 2026-11-23: R137 written up: three brain regions of about 300 donors (GSE349373): housekeeping gap with chip held −0.081 (prefrontal cortex), −0.070 (amygdala), −0.069 (hippocampus); controls −0.043, −0.036, −0.078; holds off the 0 % floor.
- 2026-11-23: R138: culture-bootstrap intervals for the per-doubling gaps stated in the paper beside the gene-resampling ones (telomerase fibroblasts, astrocytes, astrocytes with telomerase).
- 2026-11-23: culture bootstrap (cultures resampled within line, 200 resamples) written beside the gene-resampling interval for the per-doubling gaps of astrocyte_passage, astrocyte_htert and bj_telomerase (table culture_bootstrap; a referee point: the gene interval holds the cultures fixed) (R138).
- 2026-11-23: figures brain_regions: prefrontal cortex, amygdala and hippocampus of about 300 donors (GSE349373, EPIC, rebuilt from IDATs), controls and every donor with diagnosis held, Jaffe donors left out, with and without chip intercepts (RM40, R137).
- 2026-11-22: simulated referees (three, all major revision) answered in the text: the blood gap off the 0 % floor is stated (−0.08 at 10-90 %, 0.00 at 20-80 %) in Results and Abstract; arrested cultures show no measurable drift, so no shield can be measured; glia against neurons, Tet3 and Tet1, and the breast-cancer association narrowed to their evidence in the Abstract; expression stated as a correlate, not a cause; tissue sizes do not order by turnover; 5-hydroxymethylcytosine added to Limits. R134b: neurons with one array per donor (82 donors, 28 controls): housekeeping −0.051, −0.030 with chip held.
- 2026-11-22: reference check fixes: source papers of thirteen data series added to the list and cited in Methods (Pai 2019, Minteer 2023, Higgins-Chen 2022, Apsley 2025, Fröhlich 2024, White 2024, Walker 2025, Gasparoni 2018, Kwiatkowska 2023, Olecka 2024, Sokolowska 2025, Eulalio 2026, Sturm 2022); West 2013 and Hollmann 2020 credited for what they show; epidermis described as 117 arrays with up to 31 women measured twice; neuron_series uses one array per donor (83 donors; rerun R134b pending).
- 2026-11-22: review fixes (integrity and reporting checks of run 20261122T163209Z): study-wide BH rerun over every result table (13,955 tests in 114 folders, 6,971 survive, largest surviving p 0.025) and its numbers in paper and supplement; new main Figures 5 (division) and 7 (expression), later figures renumbered; Methods subsection for the expression, division and culture analyses; which intervals resample genes and which donors; second cortex series (GSE252501) and the joint cortex matching (−0.082 to −0.019) stated in Results; mortal controls of the telomerase line; sorted-nuclei gaps with intervals and n; mouse turnover ordering stated one way; mouse recipes; every accession in Methods; arrays passing detection per series (Supplementary Note 3); ethics statement.
- 2026-11-22: R136 written up (GSE227160: noise-removed squared drift 0.014 in replicative senescence against 0.0006 after irradiation and none after doxorubicin; housekeeping share −0.39, −0.36 to −0.46 leaving one donor out). figures headline: two main figures drawn from existing tables at printed size, expression_shield.png (gap before and after the expression decile; gradient by decile) and division_shield.png (forest of the division comparisons).
- 2026-11-22: R134 (sorted neurons, GSE112179: housekeeping gap −0.054 SD, −0.034 with chip held) and R135 (hTERT astrocytes, GSE226079: −0.083 per doubling, −0.084 with chip held) written into the paper (division section, Methods data list), Supplementary Note 3, Table S4 and the run log; the division lead and the closing sentence now say that neurons carry a small gap and that arrested cultures show neither drift nor shield.
- 2026-11-22: paper title changed on Anton's word to "Promoters of essential genes drift less with age, in step with their expression and with cell division". figures senescence_donors: replicative, irradiation and doxorubicin senescence in six donors (GSE227160, rebuilt from IDATs), one model per CpG with donor and chip held; matched gap of the absolute change and the noise-removed share of squared drift, with each donor left out (RM40, R136).
- 2026-11-22: paper: story, prose and deslop pass after three fresh readers (Anton asked): abstract reordered and cut to the findings (activity, division, kinetics and Tet3, what it is not, sensitivity to preprocessing); expression measured in the same samples gets its own Results section; the culture share-over-doublings paragraph moves to the selection paragraph; project history (audits, withdrawn statements, corrected measure, rebuild narration) replaced by plain facts; chip defined at first use; lead sentences of "What the shield is not" and the breast-cancer section state what the sections show; five overloaded sentences split. No number changed.
- 2026-11-22: figures astrocyte_htert: the per-doubling shield in hTERT-immortalised astrocytes (GSE226079, three replicate lines, 29 to 72 doublings, rebuilt from IDATs), same analysis as astrocyte_passage with the replicate line in the donor's place (RM40, R135).
- 2026-11-22: figures neuron_series: the shield in a third sorted-neuron series (GSE112179, NeuN+ nuclei of 100 donors, EPIC, rebuilt from IDATs), every donor with diagnosis held and controls alone, with and without chip intercepts, age-by-chip overlap tabulated (RM40, R134).
- 2026-11-20: Discussion rewritten to the post-audit findings (Anton asked): shield follows promoter activity and cell division; selection reopened; breast lead weakened; limits renumbered. Second breast series (R113) added to Results.
- 2026-11-20: paper brought up to the new datasets (Anton asked): abstract (expression in the same samples; hTERT and astrocytes), third cortex series R132 and same-adult tissues R131 in Results, datasets in Methods.
- 2026-11-20: paper main text: hTERT fibroblasts (R130) and passaged astrocytes (R133) in the division paragraph, one passage.
- 2026-11-20: R133 passaged astrocytes (GSE202554): run log, claims, evidence E48, provenance row and Table S4, Supplementary Note 3.
- 2026-11-20: `astrocyte_passage` (RM40): fetal astrocytes passaged to exhaustion (GSE202554, rebuilt), gap per doubling in a second cultured cell type.
- 2026-11-20: R132 bulk orbitofrontal cortex (GSE254293): run log, evidence E47, provenance row and Table S4, Supplementary Note 3.
- 2026-11-20: `orbitofrontal` (RM40): bulk orbitofrontal cortex of 91 people (GSE254293, rebuilt), the age gap per diagnosis group with sex or sex and chip held.
- 2026-11-20: R131 four tissues of the same adults on shared chips (GSE280465): run log, claims, evidence E46, provenance row and Table S4, Supplementary Note 3.
- 2026-11-20: `cross_tissue` (RM40): four renewing tissues of the same adults on shared chips (GSE280465, EPIC v2), paired tissue differences with a donor bootstrap.
- 2026-11-20: R130 BJ fibroblasts with and without hTERT (GSE268145): run log, claim C2, evidence E45, provenance row and Table S4, Supplementary Note 3.
- 2026-11-20: `bj_telomerase` (RM40): one fibroblast line with and without hTERT (GSE268145, rebuilt), gap per doubling in three groups, chip overlap tabulated.
- 2026-11-20: R129 epidermis of 117 women with RNA-seq of the same women (two array groups, chip held): run log, claim C2, evidence E44, provenance row and Table S4, paper and Supplementary Note 3.
- 2026-11-20: `epidermis_expression` (RM40 fetch order 4): epidermis of 121 women on EPIC or 450k with RNA-seq of the same women; `scripts/idat_barcodes.py` reads the Sentrix chip from IDAT headers when depositors renamed the files.
- 2026-11-20: R127 in the text: Lunnon joint normalisation in the paper (same-donor contrast, chromatin shares), supplement rerun table and cohort table, provenance rows; Table S4 regenerated. Epidermis studies (E-MTAB-8993/8992/8995) audited, fetch started.
- 2026-11-20: R127: Lunnon normalised jointly over tissues from raw intensities: the same-donor blood-brain difference holds (frontal cortex -0.105 against -0.092) and stays batch-confounded. Every item of the methods review (RM43) now has its run.
- 2026-11-20: R128: the 13 human groups of the species analysis give the same housekeeping result under the human annotation as under the mouse one (9 of 13 negative, correlation 0.81). Run log, Table S1, provenance, roadmap.
- 2026-11-20: `ENTROPY_AGE_SPECIES_ANNOTATION=human` runs the human groups of `figures species` with the human annotation of the 2,843 promoter probes the mammalian array shares with the 450k array (RM43 item 5; outputs end `_humanann`).
- 2026-11-20: `entropy_age.data.lunnon_joint` (RM43 item 7): Lunnon's tissues normalised together from GEO's raw intensities; `ENTROPY_AGE_LUNNON=joint` makes the cohort read them (outputs end `_lunnonjoint`).
- 2026-11-20: R126: with the bisulfite-conversion index held, clonal haematopoiesis and the syndrome effects are unchanged; in the syndromes the index differs between cases and controls (Sotos +0.93 SD, Kabuki -2.30 SD), stated as a batch caveat. Run log, supplement Note 4, provenance, roadmap.
- 2026-11-20: `ENTROPY_AGE_GCT=1` holds each array's bisulfite-conversion index as a covariate in `figures chip` and `figures syndromes` (outputs end `_gct`); both write the index by group in sesame mode (RM43 item 4). An exclusion cut is not possible: every array of these three series is above 1.2.
- 2026-11-20: R116: CR blood on the mouse-recipe rebuild is unchanged against R98 (core-essential F -0.077, p 0.40); the withdrawal of R56 stands.
- 2026-11-19: R125: the fibroblast result of R124 holds with chip intercepts (housekeeping -0.102 to -0.024, 77 % removed). Run log, paper, supplement, claim C2, provenance, reassessment.
- 2026-11-19: Reassessment: a third update section (activity measured in the same cells, R122-R124; methods-review items closed by R118, R120, R121).
- 2026-11-19: `figures sturm_expression` writes a chip-overlap table and, with `ENTROPY_AGE_CHIP=1`, holds chip in the slopes on doublings.
- 2026-11-19: R124: in Sturm fibroblast cultures (rebuilt arrays) drift per doubling follows expression of the same cultures: housekeeping gap -0.102 to -0.028 with the expression decile matched; gradient outside the gene sets. Run log, paper, supplement, claim C2, E43, provenance.
- 2026-11-19: GSE179848 fetched (Sturm fibroblast RNA-seq gene tables). `figures sturm_expression`: drift on doublings in untreated healthy-line cultures (rebuilt arrays), matched with and without the gene's expression decile in the same cultures.
- 2026-11-19: R123: CD4+ T cells with their own expression (GSE56580): almost no gap on relative entropy and no gradient with expression; on the entropy slope the expression decile removes 69 %. Run log, paper, supplement, provenance.
- 2026-11-19: `figures mesa_expression` covers CD4+ T cells too (GSE56581 with GSE56580, fetched and audited); tables gain a `cell` column.
- 2026-11-19: R122: in MESA monocytes, expression measured in the same samples removes 82 % of the housekeeping gap and all of the core-essential one; genes outside every set are shielded in step with their expression. Run log, paper, claim C2, E42, supplement Note 3, provenance (first same-sample second layer). ONNX default noted in the lead notes.
- 2026-11-19: GSE56045 fetched (MESA monocyte expression, same 1,202 people as GSE56046). `figures mesa_expression` (RM40 item 4): the gap with each gene's own expression decile in the same cells added to the matching, the gradient with expression outside the gene sets, and the methylation-expression coupling as positive control.
- 2026-11-19: `figures chromatin` rerun with the both-blood-maps matching (R121, job 6538876): table and figure committed; the paper-scale interval of the housekeeping gap with both maps held is -0.013 to +0.001, as written.
- 2026-11-19: R121 (RM43 item 3): blood analyses with the neutrophil chromatin map (E030) read as with the mononuclear one; with both maps held 96 % of the housekeeping gap is removed. Run log, paper, supplement Note 3, claim C2, roadmap.
- 2026-11-19: `ENTROPY_AGE_BLOOD_STATE=E030` gives blood its chromatin state from neutrophils in `equilibrium_drift`, `drift_hotspots` and `essentiality_gradient` (outputs end `_e030`); `figures chromatin` adds a matching on both blood maps at once (RM43 item 3). Lint cleared in four files.
- 2026-11-19: R120 (RM43 item 8): the Tet3 and Tet1 arrays share 150 CpG positions, so no shared-CpG test is possible; on the 314 genes both cover the opposite signs hold (Tet1 F -0.19, Tet3 +0.17 and +0.25). Run log, claims C4, E41, paper, supplement, Table S4; RM43 items 5, 6, 8, 9 marked done.
- 2026-11-19: `figures tet_shared` (RM43 item 8): Tet3 heterozygous and Tet1 mutant brain each scored on its own array's promoter CpGs and on the CpGs both arrays measure.
- 2026-11-19: R118 carried into claims C1 and C6, evidence rows E39 and E40, the species row of the provenance table (Table S4) and Supplementary Table S1.
- 2026-11-19: R118: the species shield holds with chip intercepts (86 % of 122 groups); within chips the consortium's mouse tissues order by turnover (lung, spleen, liver, blood, skin -0.22 to -0.49; brain regions, heart, kidney near zero). Paper and reassessment updated.
- 2026-11-19: R117: headline gaps unchanged on complete CpGs only. R119: SATSA shield-share slope holds on 450k draws alone (-0.126 per decade); the CpG divergence gap does not (-0.003, p 0.41).
- 2026-11-19: figures within_person_450k: the deposited SATSA series on its 450k draws only (RM43 item 6).
- 2026-11-19: ENTROPY_AGE_CHIP=1 adds chip intercepts to the species age slopes (RM43 item 5), outputs under species_shield_chip.
- 2026-11-19: ENTROPY_AGE_COMPLETE=1 restricts the tissue comparison to CpGs measured in every sample (no median fill; RM43 item 9), outputs under *_complete.
- 2026-11-19: paper: abstract and a new Results subsection "The shield goes with cell division" (R108-R115); eight provenance rows added, Table S4 regenerated. The Discussion still carries the earlier reading and is rewritten in round v0.2.
- 2026-11-19: claims C2 (weakened) and C3 (reopened: the shield goes with division) rewritten; evidence rows E33-E38.
- 2026-11-19: R114: by RRBS in one mouse cohort the shield is more than twice as large in colon as in hippocampus (housekeeping -0.189 against -0.078; both tissues on the same flowcells, age balanced). Reassessment updated.
- 2026-11-19: R115: in GSE66351 glia are again more shielded than the same donors' neurons (housekeeping -0.044, p 0.11; 31 donors), the direction and size of R110.
- 2026-11-19: donor linkage between methylation and expression established for MESA (1,202 of 1,202), Sturm (109 same culture and passage) and Endicott (same culture, mostly other harvests); recorded in docs/DATASET_SEARCH_2026-11-18.md.
- 2026-11-19: figures jena_rrbs: colon against hippocampus RRBS of one mouse cohort (GSE233734, GSE304754), per-tissue matched gaps with a mouse bootstrap (RM40).
- 2026-11-19: R113: the breast lead finds no support in GSE283951 (OR 1.11, 0.71-1.73; 32 cases); C12 noted.
- 2026-11-19: figures breast_replication: drift at unprotected promoters in breast cases against controls, GSE283951, chip-centred (C12; RM40).
- 2026-11-19: R112: confluent non-dividing fibroblasts show no promoter drift over six months (two clean donors) and no shield; with R111, drift in culture goes with division, not time. Reassessment updated.
- 2026-11-19: figures arrest_timecourse: drift and the shield over six months in non-dividing confluent fibroblasts, per donor within one chip (E-MTAB-15859; RM40).
- 2026-11-19: R110: other glia are more shielded than the same donors' neurons for housekeeping genes (-0.055, p 0.025); the oligodendrocyte lineage is not. R111: replicative senescence shows a large shield on chip-balanced rebuilt arrays (noise-corrected share -0.58, 12 of 12 chips); X-ray and RAS arrest show no measurable drift, so no division-against-arrest comparison is possible there.
- 2026-11-19: figures senescence_modes: replicative against X-ray and RAS senescence within chips in GSE197723 (division against arrest; RM40).
- 2026-11-19: figures nuclei_paired: glia minus neurons in the same donors of GSE279509, donor bootstrap resampling both fractions together (RM40, RM44).
- 2026-11-19: R102, R105, R109b: cortex gap on rebuilt arrays -0.082; blood minus cortex -0.07 once ages are matched and cell fractions held, interval touching zero with rebuilt cortex; GSE252501 unchanged with depression held.
- 2026-11-19: dataset search: the external search's final corrections recorded (arrested-culture series GSE197723 and E-MTAB-15859 exist; breast pre-diagnostic GSE283951 exists; GSE304754 pairs with GSE233734; several first-search verdicts downgraded); fetch order adopted (RM40 item 7).
- 2026-11-19: R108 joint cortex test: with expression, state and accessibility held together the cortex gap falls to -0.02 (housekeeping) and -0.035 (core-essential), and to zero off the floor. R109: the shield replicates in unsorted nuclei of a new series (-0.053) and is largest in glia (-0.077, -0.034), smallest in neurons (-0.021; -0.008 within chip). R107 as R97. Reassessment updated.
- 2026-11-19: dataset search corrected from the external agent's interim audit (kept: the first text of f8f7a1a below the corrections): GSE252501 shares donors with GSE74193 and its controls include depression cases; sizes; GSE238172 marrow fails; GSE233734 confirmed; GSE225166 open. brain_replication holds the depression diagnosis and reports the no-diagnosis donors apart. RM40 item 6.
- 2026-11-19: figures brain_replication: the shield in two new rebuilt brain series (GSE252501 DLPFC and nucleus accumbens controls, Jaffe donors excluded; GSE279509 sorted neurons, oligodendrocyte lineage, other glia, unsorted), complete CpGs only, with and without chip intercepts. Jaffe loader skips arrays a rebuild lacks; tissues_matched keys fixed.
- 2026-11-18: docs/DATASET_SEARCH_2026-11-18.md: replacement and matched-modality datasets checked against GEO (RM40); first to fetch: GSE279509 (sorted neurons and glia, EPIC, 18-108 y), GSE252501 (DLPFC), GSE225166, GSE121437, GSE268145, GSE202554.
- 2026-11-18: figures joint_cortex (RM44): the cortex gap with the floor band, expression decile, chromatin state and sorted accessibility held together, with the share of set CpGs that keep a control.
- 2026-11-18: R98, CR blood on rebuilt arrays: core-essential F -0.078, p 0.36; the withdrawal stands.
- 2026-11-18: docs/REASSESSMENT_2026-11-18.md: status of every result, the corrected interpretation and the research plan; RM44 (joint test in cortex); R103b QC written up.
- 2026-11-18: R106, mouse series on SeSAMe's mouse masks: Tet1, Dnmt3a, Browder spleen hold; BXD diet effect after composition back to unresolved (p 0.07); Browder muscle p 0.07.
- 2026-11-18: R104, chromatin matching with region-matched brain maps: own maps remove 5-10 % of each region's gap, blood's 25-27 %; paper updated.
- 2026-11-18: read_qc uses SeSAMe's species-mask detection where the rebuild records it (90 % cut, as on every platform); mouse series rebuilt with the per-platform recipe: BXD 288/338 pass, Browder 155/166, Tet1 BS 19/24, Dnmt3a 72/72, Huang 581/583; 3 arrays failed in the species step (RM43 item 4).
- 2026-11-18: figures tissues_matched: blood minus cortex on the shared age range with blood holding its EpiDISH cell fractions (RM43 item 2).
- 2026-11-18: chromatin matching gives Lunnon's superior temporal gyrus E072 and entorhinal cortex E071, drops the cerebellum (no Roadmap map), and uses the BA9 sorted ATAC for frontal cortex only (RM43 item 3).
- 2026-11-18: Supplementary Note 4 and a Methods paragraph set out the array processing and QC practice (Anton's reprocessing brief), each item with its status; RM43 item 9 tests the median fill of sporadic missing values.
- 2026-11-18: cohort layouts name their processing; a sesame-mode run warns when it uses a cohort without raw IDATs, and tissue_differences carries each side's processing with a mixed_processing flag (RM43 item 1).
- 2026-11-18: sesame_preprocess.R follows SeSAMe's per-platform recipe (MM285 TQCDPB, Mammal40 SHCDPB, human QCDPB), records the steps per array, and measures Mammal40 detection over the species-mask probes; mouse series rebuilt again (RM43 item 4), the first rebuilds kept in data/interim/sesame_v1.
- 2026-11-18: independent methods review: abstract no longer claims the mouse turnover order, culture or rebuild robustness beyond cortex; Lunnon and blood-minus-cortex caveats; Tet3/Tet1 heading softened; chromatin references named as related/other; provenance labels corrected (cortex rows deposited); claims C1, C2, C12 updated; RM43 lists the reruns.
- 2026-11-18: Supplementary Table S4, provenance of every comparison, rendered from docs/provenance.csv by scripts/provenance_table.py; Table S1 names the rebuilt series.
- 2026-11-18: Methods (RM42): relative entropy, the matched gap and P, B, S as equations; "matched control CpGs (look-alikes)"; no DMP/DMR calling; normalisation per source; the two platform-pooling exceptions (SATSA, transplant) replace "no analysis pools datasets or platforms". docs/provenance.csv lists every analysis's sources.
- 2026-11-18: R100, BXD liver and Browder OSKM on rebuilt arrays: directions hold; the BXD diet effect partly survives composition (-0.040, p 0.011); paper drops the cross-platform turnover ranking.
- 2026-11-18: preprocessing_qc adds per-sample beta distributions (quantiles, distance from the series median curve, outlier flag).
- 2026-11-18: the fibres rebuild maps technical replicates to the table's -N names (92 of 92 arrays, was 86); the fibre sets skip arrays the betas lack.
- 2026-11-18: R99, CH, Tet1 and Dnmt3a gain of function on rebuilt arrays: all three results hold (Tet1 on 19 of 24 arrays).
- 2026-11-18: BXD liver (GSE199979) and Browder OSKM (GSE190665) readers take the R SeSAMe rebuild under ENTROPY_AGE_BETAS=sesame.
- 2026-11-18: R97, Huang mouse tissues on rebuilt arrays: the deposited betas are R SeSAMe output to float precision; relative-entropy contrasts move <= 0.02, R92 holds.
- 2026-11-18: CH (GSE315366), Tet1 (GSE271375) and Sarni (GSE324160) readers take the R SeSAMe rebuild under ENTROPY_AGE_BETAS=sesame (reprocessed.by_accession, matrix).
- 2026-11-18: R92, mouse turnover within chip: the renewing-vs-post-mitotic ordering (−0.19) goes to +0.02 among animals on shared chips; claims C2/C3, E32 and the manuscript now say the ordering is not separable from chip.
- 2026-11-18: mouse_tissues drops the within-chip dispersion slope where every shared chip holds one aged-young pair (zero by construction; R92 hypothalamus).
- 2026-11-18: the CR blood reader (GSE282661) takes the R SeSAMe rebuild under ENTROPY_AGE_BETAS=sesame, renamed from GSM to animal title.
- 2026-11-18: R96, breast cancer on rebuilt arrays: age-adjusted OR holds (1.84) but the adjusted OR vanishes (1.11, 0.69-2.18); the lead is downgraded in claims and manuscript.
- 2026-11-18: R91, CR blood with the corrected within-chip randomisation: core-essential F p 0.0065 → 0.16, housekeeping 0.026 → 0.23; R56's finding withdrawn. Claims and manuscript corrected.
- 2026-11-18: R95, writer syndromes on rebuilt arrays: drift and shield effects hold (NSD1/EZH2 shield loss slightly stronger); concordance now rests on NSD1 and Weaver, CHARGE drops to p 0.064.
- 2026-11-18: R94, culture on rebuilt betas: Endicott's shield vanishes at the reference and trends toward growth with doublings (selection direction, CI includes 0); Sturm unchanged. Claims and manuscript updated.
- 2026-11-18: R84b (matched breast-cancer design holds, OR 1.7), R85 (no sex difference across mammals), R86 (no Braak gradient), R88 (organs of one mouse drift together, r 0.15-0.19), R90 (writer variants add noise at noisy-aging CpGs).
- 2026-11-18: R93, Jaffe paired pipeline pilot: R SeSAMe and methylprep agree within 0.004 SD; both strengthen the cortex gap against deposited betas (spread −0.07 → −0.12 SD).
- 2026-11-18: R83, per-CpG intervals: faster methylation loss at essential promoters holds in blood only and goes with chromatin state; a lower resting level holds in blood and cortex.
- 2026-11-18: outside selection audit: CR permutation refits within chips (720 exact assignments); mouse chip overlap and within-chip aging contrast; claims narrowed (culture share "no detectable change" with intervals, not "constant"; CR 12-week drift unresolved, chip-confounded; mouse turnover unadjusted for chip; Tet3 in bulk brain against same-chip wild types, not neurons or littermates). C3 renamed; evidence matrix E5, E6, E11 changed, E32 added.
- 2026-11-18: `syndrome_concordance` (RM39 item 6): do chromatin-writer syndromes add variance at noisy-aging CpGs and shift directed-aging CpGs the way age does?
- 2026-11-18: `culture_rates` (RM39 item 5), R89: drift per 10 doublings +0.015 in both culture series, per 10 days +0.001 to +0.005 (divisions about 4x time in Sturm, about equal in Endicott).
- 2026-11-18: `mouse_coupling` (RM39 item 1): clock-free drift coupling across tissues of the same Huang control mice at a fixed age, with a within-group shuffle null.
- 2026-11-18: `species_pace` (RM39 item 3): directed against noise pace of look-alike CpGs per species and tissue; within-species rank agreement and tissue-pair order across species.
- 2026-11-18: `braak_gradient` (RM39 item 2): within-donor regional gradient of drift along the Braak order in Lunnon, against Braak stage.
- 2026-11-18: `species_sex` (RM39 item 7): look-alike drift rate and the shield by sex in every species-tissue group of the mammal array with both sexes.
- 2026-11-18: R82, placebo sets matched on expression: in blood expression alone gives most of the shield (core-essential all of it); in cortex the gap is 4-7x the placebo.
- 2026-11-18: `cancer_matched` (RM38 item 4): breast cancer in EPIC-Italy as a matched case-control among women (age within 3 years, same chip first, conditional logistic regression, set bootstrap); `cancer` now writes chip and adjustments into its per-person table.
- 2026-11-18: R79 rerun on the audit fixes: hotspots still barely shared across mammals (1.04× and 1.23× chance at ≥ 10 groups).
- 2026-11-18: `per_cpg_interval` (RM38 item 2): donor intervals for the per-CpG relaxation gaps of R67.
- 2026-11-18: `placebo_expression` (RM38 item 3): placebo gene sets drawn within the real set's GTEx expression deciles, through the shield's matching.
- 2026-11-18: R80, the shield off the 0 % floor: in blood most of it is the floor (housekeeping spread to −0.08 SD at 20-80 %, core-essential to 0); cortex keeps its size. Table data/output/floor_check/gaps_by_band.csv.
- 2026-11-18: fixes from an outside audit: relative entropy regularises sample and reference alike (it could fall as a site left the 0 % floor); bootstraps of floor_check, allelic_controls and sex_drift on the full sample's SD; global FDR reads every raw p column and rejects out-of-range values; level strata survive a missing value; rescale catches conflicting SDs; species_hotspots null among measured CpGs, 0/0 no longer "inf". Regression tests in tests/test_audit_fixes.py.
- 2026-11-18: `floor_check` (RM38 item 1): the shield on CpGs at a young level of 10-90 % and 20-80 %.
- 2026-11-17: R79 across 139 mammal groups age-drift hotspots are barely shared (r 0.06); the shared tail is Polycomb developmental islands. RM35 done.
- 2026-11-16: `species_hotspots` (RM35 item 3): do the same CpGs become disordered across 65 mammals?
- 2026-11-16: manuscript round v0.1 (RM36): manuscript-improve run, 242 edits applied to `manuscript/paper.md` and `supplement.md`, corrected versions `*.proposed.md`; released to the authors on Drive.
- 2026-11-16: R77 imprinted and X-inactive loci are protected from drift in blood (−0.51, −0.29 SD), not in cortex; R78 no sex difference in drift or the shield.
- 2026-11-16: `manuscript/paper.md` and `supplement.md`, the first full draft on the corrected measure (RM36); PROVENANCE: BXD liver is on the mammalian array, not RRBS.
- 2026-11-16: `allelic_controls` and `sex_drift` (RM35 items 1 and 2): imprinted and X-inactivated loci as a noise floor; drift and the shield by sex.
- 2026-11-16: RM35, the ideas still open on data at hand (imprinted and X loci as noise controls, sex, cross-mammal hotspots).
- 2026-11-16: RM4 and RM13 wait until mid-December (no RNA shield replicated; no read-level culture series yet).
- 2026-11-16: RM34 done: `docs/INPUT_CHECKSUMS.tsv`, SHA-256 of every downloaded input (R76).
- 2026-11-16: R75 written up (RM4): Tabula Muris Senis bulk shows no RNA shield and no turnover order.
- `figures rna_tms` (RM4): the RNA turnover test in Tabula Muris Senis bulk (GSE132040, 17 organs, 3-27 months); `analyse_species` takes a turnover map; run R75 planned.
- 2026-11-16: RM2 done (on relative entropy the mouse array orders the shield by turnover, R68, R73; the muscle cohort with activity data was never found) and RM33 done (per-person `K` and read-level entropy left: P/B/F and PDR already carry those claims); RM3 blocked on desk #28, RM4 unattended, RM19 waits until 2026-12-14.
- 2026-11-16: RM32 done: the faster removal at essential promoters fits Tet3 (its loss makes them gain), not Tet1; read from R15, R16, R35, R36, no new run.
- 2026-11-16: RM28 blocked on desk #40 (apply for a second breast-cancer cohort).
- 2026-11-14: R66 written up (RM28): SATSA rebuilt from IDATs shows the shield within people and across them; `reprocess_idats --from-pickle`.
- 2026-11-14: R73 written up (RM33): mouse RRBS, BXD, intervention and Meer aging gaps keep their sign on relative entropy, smaller.
- 2026-11-14: R74 written up (RM32): faster removal at essential promoters survives within active chromatin; their lower resting level is mostly chromatin.
- `equilibrium_drift` (RM32 next step): gain and loss rates within active-TSS chromatin only, and per-CpG gaps with the own chromatin state in the strata; run R74 planned.
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

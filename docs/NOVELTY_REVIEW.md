# Novelty review of the shielding claims (2026-11-11)

Asked by antonkulaga: has someone already done this? A literature agent searched PubMed, bioRxiv and the web;
every citation below was checked against PubMed metadata except where marked "not verified". Verdicts are
per claim of `docs/PAPER_DRAFT.md`.

## Overall verdict

The core idea is published as a qualitative claim. Teschendorff, West & Beck 2013 (Hum Mol Genet, PMID 23918660)
state that genes doing core housekeeping functions appear protected from drift ("epigenetic thrift"). West,
Widschwendter & Teschendorff 2013 (PNAS, PMID 23940324) show highly connected (essential-enriched) genes are spared.
Day et al. 2013 (Genome Biol, PMID 24034465) and the Polycomb papers (Teschendorff 2010, Rakyan 2010, Slieker 2016)
give the chromatin side. The project is not new in concept. Its novelty is in design and quantity: a null matched on
island status and starting level, dispersion and entropy together, many tissues and two species, and a split of the
shield into chromatin and the rest. What is new in direction: the dynamics (culture over doublings, the share growing
with age in vivo, the caloric-restriction interaction) and the cross-species contrast. Each needs added controls.

## By claim

| claim | verdict | closest prior work | what would make ours new |
|---|---|---|---|
| 1 essential genes drift less | partly known | Teschendorff 2013 HMG; West 2013 PNAS; Day 2013; Teschendorff 2010, Rakyan 2010 (Polycomb gain); Slieker 2016 (variability at Polycomb) | essentiality still predicting protection after matching on expression level and breadth, H3K4me3 breadth and Polycomb; a gradient within essential genes (DepMap scores) |
| 2 constant share in culture | measurement new; the "against selection" reading weak | Endicott 2022 (replication-driven loss); Zhou 2018 (mitotic clock); Franzen 2016; **Hollmann 2020: clonal dominance arises in long-term culture** | clonality data (barcodes, epiallele structure); show the shield does not move when clonal makeup does |
| 3 chromatin explains it | known in direction; the split and blood = brain new | Jenkinson 2017 (entropy and chromatin); Day 2013; Horvath 2013; **Vaidya 2023: entropy tracks stem-cell division across mouse tissues** | a formal decomposition with each tissue's own chromatin; "essentiality adds X % beyond chromatin"; square blood = brain with Vaidya |
| 4 share grows with age in vivo; CR weakens it | appears new | Maegawa 2017 (CR slows drift, mouse and monkey); Hahn 2017 (DR protective in liver, and adds gene-body methylation at metabolic genes); Petkovich 2017 | absolute slopes in both arms, not only the share; a second CR series |
| 5 clocks lean on Polycomb and tissue-variable chromatin | mostly known | Teschendorff 2010; Rakyan 2010; Lu 2023 (universal clocks at PRC2 sites); **Moqri 2024: PRC2-AgeIndex, a Polycomb-only clock**; Tong 2024 (StocH/P/Z); Tarkhov 2024 | the reproducible tissue offset, predicted from chromatin, set against PRC2-AgeIndex |
| 6 per-person shield score | type known; the contrast form new | Gentilini 2015 (stochastic epimutations); **Seeboth 2020: per-person outlier burden confounded by cell mix, smoking, BMI, batch, and not predictive of health**; Hannum 2013 | show the within-sample contrast cancels what outlier burden does not; test-retest; cell-mix robustness; an outcome beyond clocks |
| 7 species lifespan (running) | partly answered; the essential-vs-matched contrast new | Crofts 2024 (methylation rates scale with lifespan); **Crofts 2026 (SCARLET): lifespan scaling explained by stem-cell pool size over division rate, "rather than epigenetic maintenance efficiency"**; Haghani 2023 (lifespan modules at Polycomb/homeobox genes); McLain & Faulk 2018; Vaisvil 2026 (tortoise, n=1) | separate protection from stem-cell dynamics; note that conserved array CpGs at essential promoters sit near 0 %, so dispersion is floored |

## Must cite and position against

1. Teschendorff, West & Beck 2013, Hum Mol Genet, PMID 23918660: the claim stated.
2. West, Widschwendter & Teschendorff 2013, PNAS, PMID 23940324: the closest empirical precedent.
3. Day et al. 2013, Genome Biol, PMID 24034465: tissue expression and chromatin protect.
4. Slieker et al. 2016, Genome Biol, PMID 27654999: age-variable CpGs at Polycomb regions.
5. Crofts et al. 2026, Nat Aging, PMID 42156953 (with Crofts 2024, PMID 38049585): lifespan scaling as stem-cell dynamics.

Close runners-up: Tong 2024 (PMID 38724732), Moqri 2024 (PMID 39009581), Endicott 2022 (PMID 36347867), Vaidya 2023
(PMID 36797759), Hollmann 2020 (PMID 32138773), Seeboth 2020 (PMID 32216821). Not verified this session: Yang 2016
(epiTOC); a bioRxiv preprint on a 17-tissue ageing atlas with "functionally fragile" hub-gene CpGs
(10.1101/2025.07.21.665830), to read in full.

## Where prior work contradicts us

- Claim 2: long-term cultures become oligoclonal (Hollmann 2020); a constant share does not exclude steady selection.
- Claim 7, and any maintenance reading: Crofts 2026 explains lifespan scaling by stem-cell dynamics.
- Claim 3: Vaidya 2023 ties entropy to division rate, which predicts blood and brain to differ; equal after matching
  may be an artefact of conditioning on chromatin.
- Claim 6: Seeboth 2020's warning, which the score must answer.
- Claim 4: Hahn 2017's targeted methylation under DR could move the share without any loss of protection.

## Addendum, same day: what is left beyond chromatin

From R4's tables: after own-chromatin matching, core-essential genes in blood keep −0.02 SD (interval touching 0),
housekeeping genes −0.03 to −0.05 SD (clear), and in brain both −0.03 to −0.06 SD. So the deflationary reading
(active chromatin relabelled) holds for core-essential genes in blood, not for housekeeping genes or for brain.

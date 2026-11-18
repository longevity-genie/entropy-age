# Reassessment after the 18 November audits

Asked by antonkulaga, 2026-11-18: reassess the hypothesis, the current results, the latest fixes, the correct
interpretations and the research plan. Written after three audits in one day (selection, array compatibility,
independent methods review), the rebuild of every array series from raw IDATs, and reruns R80-R106.
Still running when this was written: R98 (CR blood rebuilt), R102 (main tissues with rebuilt cortex), R105 (blood
minus cortex on matched ages), R107 (Huang on the mouse recipe).

## Update, 19 November: the two tests section 5 called decisive

- **The joint test in cortex (R108).** With expression, chromatin state and sorted accessibility held together on the
  rebuilt Jaffe arrays, the housekeeping gap falls from −0.082 to −0.019 SD and the core-essential gap from −0.073 to
  −0.035; off the 0 % floor nothing distinguishable from zero is left. So the "interesting residue in cortex" of
  section 3 is mostly promoter activity and the floor as well. An essentiality-specific shield is not shown in either
  tissue.
- **Sorted nuclei (R109, GSE279509).** The shield replicates in an independent series: unsorted cortical nuclei
  −0.053 (148 donors). By cell type: other glia −0.077, oligodendrocyte lineage −0.034, neurons −0.021, and within
  chips the neuronal gap is −0.008 with an interval across zero. It is strongest where cells can divide and weakest
  where they cannot. Selection, or any division-coupled maintenance, is therefore back as a live explanation; the
  direct paired test of neurons against glia in the same donors is the next run.
- **Second DLPFC series (R109, GSE252501).** Same sign, −0.005 to −0.019, a quarter of Jaffe's size or less, in 52
  independent donors.
- **What this does to the plan.** B is answered (the claim is "active promoters drift less"). C is half answered:
  replication yes, neurons no. D now has usable data: the neuron-glia contrast within donors and chips is the
  cleanest selection-type test found so far. The paper's centre moves from "essential genes are shielded by
  chromatin" towards "promoter drift is lowest at active promoters, and in brain the difference is carried by glia".

## 1. What the project set out to show

1. **The shield.** Promoters of essential and housekeeping genes drift less with age than CpGs of the same island
   status and starting level.
2. **Its cause.** First idea: selection among dividing cells (H2). Later reading: chromatin maintenance and faster
   demethylation (H1), with Tet3 as the remover in brain (H4).
3. **Its use.** Clocks read unprotected noise (H5); drift at unprotected promoters marks disease risk (C12).

## 2. Where each result stands now

| Result | Status | Why |
|---|---|---|
| Shield, cortex (Jaffe) | **stands, small** | −0.06 SD; stronger on rebuilt arrays (R93); keeps its size off the 0 % floor (R80); 4-7x the expression-matched placebo (R82). One cohort. |
| Shield, blood (Hannum) | **stands as an association, explained** | −0.15 SD, but two thirds is the 0 % floor (R80), four fifths is expression (R82), 83-91 % goes with a reference chromatin state. No raw files; never rebuilt. |
| Shield across 65 mammals | **stands, one open check** | 89 % of 139 groups; survives the floor test (R59). Chip goes with age in 56 groups and no chip term is fitted (RM43 item 5). |
| Shield in GTEx, read-level WGBS, sorted cells | stands, deposited only | Same sign everywhere; sizes small. |
| Within one person (SATSA) | stands on the 450k rebuild | −0.14 SD, p 0.001 (R66). The deposited version pools two platforms. |
| Blood more shielded than brain | **not shown** | Hannum against Jaffe is two studies; Lunnon's blood and brain sit on different slides with separate normalisation. |
| Chromatin explains the shield | **partly** | Blood: yes, but with a related map (mononuclear cells), from other people. Brain: a region's own map removes 5-10 %, blood's map 25-27 % (R104), so the map acts as a generic "active promoter" label. Cortex keeps 44 % with state and sorted ATAC together. |
| Faster loss at essential promoters | blood only, chromatin-explained | R83. A lower resting level in both tissues. |
| Tet3 loss weakens the shield | one study, not rebuildable | 8 mice per tissue, same-chip wild types, no raw files. |
| Tet1 loss strengthens it | stands on rebuilt arrays | 19 of 24 arrays; another array and study than Tet3, so not a contrast with it. |
| Dnmt3a gain of function, writer syndromes | stand on rebuilt arrays | Syndrome series have poor bisulfite conversion (GCT 1.4-1.5); exclusion check open. |
| Mouse turnover ordering | **withdrawn** | Gone within chip (R92, R97). |
| Shield in culture; no growth with divisions | **uninformative** | Endicott's shield vanished on rebuilt arrays (R94). QC (R103b): chip explains 79 % of the main structure, doublings are 77 % determined by chip, batch date 100 %. |
| CR weakens the shield | **withdrawn** | p 0.006 → 0.16 under a valid within-chip randomisation (R91). |
| Breast-cancer lead | **weak** | Adjusted odds ratio 1.76 → 1.11 on rebuilt arrays (R96). |
| BXD diet effect beyond composition | unresolved | p 0.26 deposited, 0.011 first rebuild, 0.07 with the mouse masks (R106). |
| No sex effect, no Braak gradient, no hotspots, no RNA or clean ATAC counterpart | null results stand | R78, R85, R86, R79, RM4, C11. |

## 3. The corrected interpretation

1. **The shield is real in sign and modest in size.** It replicated in every platform, tissue and species tried.
   That part is robust.
2. **It is not yet shown to be about essentiality.** In blood, a promoter being active, highly expressed and near
   0 % methylation accounts for it. What the data support there is "constitutively active CpG-island promoters
   drift less", which is close to known biology.
3. **The interesting residue is in cortex.** There the gap survives the floor, expression matching and chromatin
   matching, each tested separately. It has not been tested with all three held together, and it rests on one cohort.
4. **Selection is untested, not rejected.** Every design that could have detected it turned out confounded with chip
   or batch (culture, CR, mouse turnover, Lunnon). "No sign of selection" overstates what we know.
5. **Mechanism claims are single studies.** Tet3 is hypothesis-generating. "Chromatin maintenance" is an association
   with reference maps from other people.
6. **Nothing links the shield to health.** The one lead did not survive raw processing.
7. **Effects of 0.02-0.1 SD are the size of processing and batch effects.** Today's reversals (CR, culture, cancer,
   BXD diet) were all of that size. Any claim of that size needs a within-chip design and raw data before it is made.

## 4. What the paper can honestly claim today

- A matched-CpG method, and a shield that replicates in sign across tissues and 65 mammals.
- In blood it is explained by promoter activity; in cortex a part is not, in one cohort.
- A set of negative and withdrawn results with their causes (floor, expression, chip confounding, processing),
  which is itself useful to the field.

The working title ("shielded ... by chromatin and faster demethylation") claims more than this. The title is the
authors' decision; the present evidence fits "drift less" better than a named mechanism.

## 5. Research plan, in order

**A. Finish the fixes already decided (RM43, RM41).** Running: R98, R102, R105, R107. Then: chip intercepts for the
species groups; headline gaps on complete CpGs only (no median fill); blood with the neutrophil map; SATSA on 450k
pairs only; CR series on the mouse recipe; syndromes and clonal haematopoiesis without the poorly converted arrays;
Lunnon renormalised jointly; Tet3 and Tet1 on shared CpGs; the clock-scoring audit.

**B. The decisive test of the remaining claim.** One joint matching in cortex: island status, young level off the
floor, expression decile, chromatin state and sorted accessibility together. If the gap goes to zero, the paper's
claim becomes "active promoters drift less" and the essentiality story ends. If it stays, it is the main finding.
This needs no new data and is a few hours of compute. It should be run before anything else new.

**C. Replicate the cortex result.** A second brain cohort with raw IDATs and a recorded chip layout, searched for
under RM40; best of all sorted neurons, which do not divide, so a shield there cannot come from selection among
cells. Availability of raw files has to be checked per series before any is named in the plan.

**D. A selection test that can fail.** A culture or clonal dataset where passage is not confounded with chip, or
sequencing-based culture data. Until one exists, the paper says selection was not testable with public arrays.

**E. Drop or park.** The health lead (needs a second prospective cohort, desk question open), CR, the turnover
ordering, culture as evidence, and Tet3-against-Tet1 as a contrast.

**F. Then the manuscript.** v0.2 (after 2026-11-23) rewritten around section 4, with Supplementary Table S4 and
Note 4 as they now stand, and the artifact rebuilt to match.

## 6. Rules adopted from today's failures

- No group contrast without a chip-overlap table first; a contrast chip fully determines is reported as
  unidentifiable.
- No claim from deposited betas where raw files exist.
- Every second track (chromatin, expression) is labelled by sample and tissue match (Table S4).
- Tests pass before a job is submitted.

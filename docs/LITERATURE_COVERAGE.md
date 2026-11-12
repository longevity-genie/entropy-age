# Literature coverage of the "beyond the shield" ideas

Checked 2026-11-12 for the 12 ideas of `docs/BEYOND_SHIELD.md` and the two running items RM24 and RM25, so the
project does not redo published work. Every citation below was fetched from PubMed this session (PMID and DOI from
the record); the one preprint not in PubMed is marked with where its DOI was confirmed. Sources searched: PubMed,
Consensus, web search (bioRxiv by DOI). Verdicts: **done** (same question on the same or equivalent data is
published: drop or treat as replication), **partly** (related; the gap is named), **open** (nothing close found).

| idea | coverage | closest papers | remaining gap |
|---|---|---|---|
| 1 Cross-tissue coupling of aging rate | partly | Horvath 2022 GeroSci; Richardson 2025 Aging; Grodstein 2021 Front Neurosci; Zahn 2007 | clock-free drift/entropy coupling; mouse at fixed age |
| 2 Hardy (illness before death), smoking | open (Hardy); partly (smoking) | Li 2024 AJHG; Richardson 2025; Abouhashem 2023 (RNA) | Hardy never tested on methylation age or drift |
| 3 Interventions: noise or programme | partly | Perevoshchikova 2024 preprint; Maegawa 2017; Cole 2017; Huang 2025 | directed vs stochastic split in 11 tissues, GLP-1RA |
| 4 AD: disorder by Braak region | partly | Levine 2015; Shireby 2020; Thrush 2022; Pellegrini 2021; Huo 2018 | regional gradient of drift/entropy within donors |
| 5 / RM25 Tissue aging order | partly | Shavlakadze 2026 (data source); Tyshkovskiy 2026; Zahn 2007; Ibañez-Solé 2022 | rank test of directed vs noise pace across species |
| 6 When drift begins | partly (onset open) | Bertucci-Richter 2023; Horvath 2013; Kerepesi 2021; Mulder 2021 | change point of drift vs sexual maturity |
| 7 Methylation noise to expression noise | partly | Oliva 2023; Jain 2024; Oliveros 2026; Slieker 2016; Fan 2025 | dispersion-to-dispersion coupling, same GTEx donors |
| 8 Time vs divisions vs metabolism | partly (much done) | Sturm 2023; Sturm 2019; Endicott 2022; Kabacik 2022 | drift/entropy per day vs per doubling, all agents |
| 9 Progeroid and machinery mutants | partly (mostly done) | Sarni 2026 (data owners); Martin-Herranz 2019; Pérez 2023 | stochastic-site concordance only |
| 10 Sex differences in drift | partly | Yusipov 2020; Horvath 2016; Horvath 2022; Oliveros 2026 | drift by sex across mouse tissues and mammals |
| 11 Universal noise hotspots | partly | Haghani 2023; Lu 2023; Bertucci-Richter 2023; Lowe 2018 | per-CpG variance gain across ~349 species |
| 12 Schizophrenia blood and brain | done | McKinney 2017 (x2); Higgins-Chen 2020; Kiltschewskij 2024 | none worth a run; replication only |
| RM24 Shield score in AD, SZ, smoking, Hardy | partly | Huo 2018; Kiltschewskij 2024; Li 2024 | the promoter-shield score itself; Hardy arm open |

## Per idea

**1. Cross-tissue coupling (partly).** The question "is a person's epigenetic age acceleration shared across
organs?" is published with clocks. Horvath et al. 2022, *GeroScience*, "HIV, pathology and epigenetic age
acceleration in different human tissues" (PMID 35411474, doi:10.1007/s11357-022-00560-0): 661 samples, 11 tissues,
133 deceased donors (75 HIV+), custom mammalian array; epigenetic age moderately correlated across tissues, blood
correlated best with spleen and bone marrow, not with liver; "blood alone will often fail to detect EAA in other
tissues". Richardson et al. 2025, *Aging*, "Characterization of DNA methylation clock algorithms applied to diverse
tissue types" (PMID 39754638, doi:10.18632/aging.206182): 8 clocks on the same GTEx EPIC data we hold (973
samples, 9 tissues); cross-tissue clock correlation with "some evidence of some residual correlation after
adjusting for chronological age". Grodstein et al. 2021, *Front Neurosci* (PMID 33488342,
doi:10.3389/fnins.2020.555307): ROSMAP, antemortem CD4+ cells vs DLPFC in 41 people, epigenetic age r 0.21-0.32.
Hannon et al. 2015, *Epigenetics* (PMID 26457534, doi:10.1080/15592294.2015.1100786): blood plus four brain
regions in 122 people (the London cohort behind Lunnon), per-CpG co-variation, not clocks. Horvath et al. 2015,
*Aging*, "The cerebellum ages slowly according to the epigenetic clock" (PMID 26000617, doi:10.18632/aging.100742):
up to 30 sites per donor, tissue-level offsets. In mice, Zahn et al. 2007, *PLoS Genet*, AGEMAP (PMID 18081424,
doi:10.1371/journal.pgen.0030201): 16 tissues, RNA microarrays, "different tissues age in a coordinated fashion in
individual mice". Organ clocks from other layers: Oh et al. 2023, *Nature* (PMID 38057571,
doi:10.1038/s41586-023-06802-1), plasma proteome, 11 organs, 5,676 adults; Tian et al. 2023, *Nat Med* (PMID
37024597, doi:10.1038/s41591-023-02296-6), UK Biobank imaging and physiology, 10 systems. **Gap:** the GTEx clock
arm is a replication of Richardson 2025. Still open: a Lunnon five-tissue donor matrix of residuals with a
one-factor model, clock-free drift `B`/entropy `E` coupling, and methylation coupling among mice of a fixed age
(Huang GSE279437).

**2. Hardy scale (open for methylation; smoking partly).** No paper was found using GTEx `DTHHRDY` (or any
slow-vs-sudden death contrast) as an exposure for methylation age, drift or entropy. Hardy appears only as an RNA
nuisance or exposure: Ferreira et al. 2018, *Nat Commun* (PMID 29440659, doi:10.1038/s41467-017-02772-x),
post-mortem interval and death effects on GTEx transcriptomes; Abouhashem et al. 2023, *bioRxiv* (PMID 37292819,
doi:10.1101/2023.05.15.540715), slow death induces survival pathways in 701 GTEx skin transcriptomes. Smoking in
the same GTEx methylation is published: Li et al. 2024, *Am J Hum Genet* (PMID 38490207,
doi:10.1016/j.ajhg.2024.02.012), smoking-associated CpGs in 9 tissues (6,350 in lung); Richardson 2025 (above)
found smoking raises lung clock age. **Gap:** Hardy on clocks, `B`, `E` across 8 tissues is new; for smoking,
only the "directed sites vs drift everywhere" split is new.

**3. Interventions: noise or programme (partly).** Closest: Perevoshchikova and Fedichev 2024, *bioRxiv*,
"Differential responses of dynamic and entropic aging factors to longevity interventions"
(doi:10.1101/2024.02.25.581928; confirmed on Sciety, not in PubMed): mouse DNAm split into a dynamic component that
CR and parabiosis move and an entropic (linear) component they do not. Maegawa et al. 2017, *Nat Commun*,
"Caloric restriction delays age-related methylation drift" (PMID 28912502, doi:10.1038/s41467-017-00607-3): CR
attenuates age-related change in monkey and mouse ("drift" there means age-related mean change). Cole et al. 2017,
*Genome Biol* (PMID 28351383, doi:10.1186/s13059-017-1185-3) and Hahn et al. 2017, *Genome Biol* (PMID 28351387,
doi:10.1186/s13059-017-1187-1): Ames dwarf, CR, rapamycin suppress age changes in liver WGBS. Data owners: Huang et
al. 2025, *Cell Metab*, "Body-wide multi-omic counteraction of aging with GLP-1R agonism" (PMID 41265449,
doi:10.1016/j.cmet.2025.10.014), exenatide opposes age methylation across tissues, rapamycin correlated. Jacques
et al. 2026, *Nat Aging* (PMID 42362890, doi:10.1038/s43587-026-01164-5): a gene cluster not modified by
beneficial interventions. Tong et al. 2024, *Nat Aging* (PMID 38724732, doi:10.1038/s43587-024-00600-8):
stochastic share of clocks. **Gap:** the per-tissue directed/stochastic split in the 11 same-animal GLP-1RA
tissues is unpublished; Perevoshchikova predicts "stochastic untouched", so this is a sharp test, not a first.

**4. AD by Braak region (partly).** Clocks in AD brain are well covered: Levine et al. 2015, *Aging* (PMID
26684672, doi:10.18632/aging.100864), DLPFC epigenetic age vs plaques in ROSMAP; Shireby et al. 2020, *Brain*
(PMID 33300551, doi:10.1093/brain/awaa334), Cortical clock; Grodstein et al. 2021, *Neurobiol Dis* (PMID 34153464,
doi:10.1016/j.nbd.2021.105428); Thrush et al. 2022, *Aging* (PMID 35907208, doi:10.18632/aging.204196),
PCBrainAge in multiple regions. Pellegrini et al. 2021, *Front Aging Neurosci* (PMID 33790779,
doi:10.3389/fnagi.2021.639428): meta-analysis of temporal, frontal, entorhinal cortex and cerebellum (the Lunnon
design), AD changes concordant with aging. Variability: Huo et al. 2018, *Neurobiol Aging* (PMID 30660039,
doi:10.1016/j.neurobiolaging.2018.12.003), 249/115 VMPs with amyloid and tangles in DLPFC. Lunnon et al. 2014,
*Nat Neurosci* (PMID 25129077, doi:10.1038/nn.3782) is the data paper (ANK1). **Gap:** slope of drift/entropy on
Braak by region, its EC > STG > FC > CER gradient, and whether it sits in the donor-wide or region-specific part.

**5 / RM25. Tissue aging order (partly).** The rodent RNA data (GSE314817/GSE314970) are from Shavlakadze et al.
2026, *Sci Adv*, "Multitissue, multi-time point transcriptomic atlas of aging in mice and rats" (PMID 41984968,
doi:10.1126/sciadv.ady8401): 28 mouse and 32 rat tissues, >5,000 samples; reports tissues "relatively spared"
(gut) and species differences, i.e. tissue amounts of directed change, not noise. Tyshkovskiy et al. 2026,
*Nature* (PMID 42203874, doi:10.1038/s41586-026-10542-3): >11,000 transcriptomes, 25 tissues, mouse, rat, macaque,
human; conserved aging changes. Schaum et al. 2020, *Nature* (PMID 32669715, doi:10.1038/s41586-020-2499-y): 17
mouse organs, asynchronous onset. Zahn 2007 (above): no overall mouse-human correlation. Ibañez-Solé et al. 2022,
*eLife* (PMID 36576247, doi:10.7554/eLife.80380): no consistent evidence of rising transcriptional noise, in line
with RM25's unshared noise order. Izgi et al. 2022, *eLife* (PMID 35098922, doi:10.7554/eLife.68048): tissues
converge in expression with age. Methylation: Jacques 2026 (above), 17 human tissues incl. variability. **Gap:**
RM25's explicit rank correlation of directed vs noise pace (mouse vs rat rho 0.56; noise unshared; human
methylation vs RNA unrelated) is not published as such; it should cite Shavlakadze 2026 and Ibañez-Solé 2022 and
fix the noted technical caveat before a claim.

**6. When drift begins (partly; onset open).** Bertucci-Richter et al. 2023, *Nat Commun* (PMID 38007590,
doi:10.1038/s41467-023-43417-6): disorder rises with age in four mammals, faster in shorter-lived species.
Maegawa 2017 and Lowe et al. 2018, *Genome Biol* (PMID 29452591, doi:10.1186/s13059-018-1397-1): age-change rate
scales with lifespan. Horvath 2013, *Genome Biol* (PMID 24138928, doi:10.1186/gb-2013-14-10-r115): log-linear clock
before 20. Kerepesi et al. 2021, *Sci Adv* (PMID 34172448, doi:10.1126/sciadv.abg6082): rejuvenation then aging
from mid-embryogenesis. Lu et al. 2023, *Nat Aging* (PMID 37563227, doi:10.1038/s43587-023-00462-6): universal
clocks scale age by sexual maturity. Mulder et al. 2021, *Hum Mol Genet* (PMID 33450751,
doi:10.1093/hmg/ddaa280) and 2024, *EBioMedicine* (PMID 39471750, doi:10.1016/j.ebiom.2024.105425): blood
between-child variation from birth; clock sites vary from birth. Jaffe et al. 2016 (online Dec 2015), *Nat Neurosci* (PMID 26619358,
doi:10.1038/nn.4181) is the data paper (prenatal to postnatal change). **Gap:** a change point in drift (not
clock) at sexual maturity, and its scaling with maturity age across the mammal consortium species, is open.

**7. Methylation noise to expression noise (partly).** Oliva et al. 2023 (online Dec 2022), *Nat Genet* (PMID 36510025,
doi:10.1038/s41588-022-01248-z): GTEx methylation, eQTMs, mQTLs (mean coupling). Jain et al. 2024, *Epigenetics
Chromatin* (PMID 39118140, doi:10.1186/s13072-024-00546-6): age EWAS in the 9 GTEx tissues and age-eQTMs
(CDKN2A). Oliveros et al. 2026, *Genome Res* (PMID 42448427, doi:10.1101/gr.281428.125): <10% of age-related
expression variation mediated by methylation in GTEx. Slieker et al. 2016, *Genome Biol* (PMID 27654999,
doi:10.1186/s13059-016-1053-6): blood aVMPs and RNA-seq in 2,044 people, trans links. Fan et al. 2025, *Genome
Res* (PMID 40750330, doi:10.1101/gr.280155.124): drift CpGs tied to transcriptional variability via scRNA-seq;
per-person drift score. **Gap:** gene-wise correlation of methylation dispersion slope with expression variance
slope, and donor drift `B` vs donor expression noise, in the same GTEx donors and tissues: not found.

**8. Time vs divisions vs metabolism (partly, much done).** Sturm et al. 2023, *Commun Biol* (PMID 36635485,
doi:10.1038/s42003-022-04303-x): the GSE179847 series; OxPhos defects accelerate epigenetic aging per division.
Sturm et al. 2019, *Epigenetics* (PMID 31156022, doi:10.1080/15592294.2019.1626651): clocks in culture ~62x
faster. Endicott et al. 2022, *Nat Commun* (PMID 36347867, doi:10.1038/s41467-022-34268-8): PMD loss tracks
doublings. Kabacik et al. 2022, *Nat Aging* (PMID 37034474, doi:10.1038/s43587-022-00220-0): clocks tied to
nutrient sensing and mitochondria, not senescence. **Gap:** drift/entropy (not clocks) decomposed per day vs per
doubling across all 19 agents; small.

**9. Progeroid mutants (partly, mostly done).** Sarni et al. 2026, *Nat Genet* (PMID 42286141,
doi:10.1038/s41588-026-02633-8): the GSE324160 data paper; DNMT3A gain-of-function recapitulates age gains.
Martin-Herranz et al. 2019, *Genome Biol* (PMID 31409373, doi:10.1186/s13059-019-1753-9): NSD1/Sotos accelerates
the clock and shares age changes. Pérez et al. 2023, *Aging Cell* (PMID 38140713, doi:10.1111/acel.14058): only
DNA-repair progeroid mice show clock acceleration. **Gap:** concordance at stochastic vs directed sites only.

**10. Sex differences in drift (partly).** Yusipov et al. 2020, *Aging* (PMID 33276343,
doi:10.18632/aging.202251): four blood sets; no sex difference in age gain of entropy or epimutations, but 15x more
probes gain variability in males. Horvath et al. 2016, *Genome Biol* (PMID 27511193,
doi:10.1186/s13059-016-1030-0): men's faster clocks. Horvath 2022 (above): male EAA in several tissues. Oliveros
2026 (above): sex effects in GTEx methylation. **Gap:** sex term on drift across mouse tissues and the mammal
consortium.

**11. Universal noise hotspots (partly).** Haghani et al. 2023, *Science* (PMID 37561875,
doi:10.1126/science.abq5693) and Lu 2023 (above): mean-change modules across ~348 species. Bertucci-Richter 2023
and Lowe 2018: disorder and age-change rates in 4-6 species. **Gap:** per-CpG variance gain counted across the
349-species set; not found. Heavy run.

**12. Schizophrenia (done).** McKinney et al. 2017, *Schizophr Res* (PMID 28988914,
doi:10.1016/j.schres.2017.09.025): no Horvath acceleration in DLPFC (the Jaffe/Lieber set) or two blood sets;
McKinney et al. 2017, *NPJ Schizophr* (PMID 28560259, doi:10.1038/s41537-017-0017-5): STG. Higgins-Chen et al.
2020, *Biol Psychiatry* (PMID 32199607, doi:10.1016/j.biopsych.2020.01.025): 14 clocks in three blood sets
(GSE41037 among the public sets), mortality clocks up. Kiltschewskij et al. 2024, *Mol Psychiatry* (PMID 39271751,
doi:10.1038/s41380-024-02749-5): blood methylation variance (213 VMPs). Drop; at most a replication line in RM24.

**RM24 (partly).** Variability with AD pathology (Huo 2018) and schizophrenia (Kiltschewskij 2024), smoking in
GTEx (Li 2024) are published; a promoter-shield score and the Hardy contrast are not. Keep RM24, cite these, and
report the Hardy arm as the new part.

## Recommended re-ranking (open and feasible first)

1. **Idea 2, Hardy scale on methylation drift and clocks** (open; 4-6 h; reuses `health_states`).
2. **Idea 6, drift onset vs sexual maturity** (onset open; 6-10 h; Jaffe, Petkovich, mammal consortium).
3. **Idea 7, methylation dispersion to expression dispersion in GTEx donors** (open for dispersion; 10-16 h).
4. **Idea 4, AD drift/entropy gradient over Braak regions** (partly; 4-6 h after RM24; clocks are background).
5. **Idea 3, directed vs stochastic response to GLP-1RA across 11 tissues** (partly; a test of
   Perevoshchikova 2024; 8-10 h).
6. Idea 1, restricted to clock-free coupling (Lunnon five tissues) and the fixed-age mouse arm; the GTEx clock
   arm is a replication of Richardson 2025 and Horvath 2022.
7. Idea 10 (sex in drift, cross-species), idea 11 (variance hotspots; heavy), idea 8 (small gap).
8. Idea 9 (mostly done by Sarni 2026); RM25 needs only citations and its technical control.
9. Drop idea 12 (done).

# Young-adult blood reference for the young-SD question (P6)

Search of 28 November 2026. `level_match_reference` takes the matching variables (young mean level, young
between-donor SD) from another 450k blood cohort. With adult references (Hannum's youngest fifth, 19–52 y; EPIC-Italy's,
34–47 y) matching on young SD halves the housekeeping gap; with children (GSE97362 controls, 2–19 y) it does not. A
cohort of healthy donors aged 18–30 separates the two readings: low young SD set by early adulthood, or a property of
the adult references. Wanted: 450k, whole blood or PBMC, healthy donors mostly 18–30, n ≥ 50, raw IDATs preferred,
deposited betas acceptable; no disease-only cohorts; no donor overlap with Hannum, EPIC-Italy or GSE97362.

Metadata were read from each series' GEO record and series-matrix header, file sizes from the GEO download listing.
The earlier searches (`DATASET_SEARCH_2026-11-18.md` and the three of 2026-11-24) assessed no young-adult blood cohort.

| Accession | Tissue, platform | Donors aged 18–30 | Files | Verdict |
|---|---|---|---|---|
| [GSE73103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73103), PMID [26449484](https://pubmed.ncbi.nlm.nih.gov/26449484/) (Voisin 2015) | Whole blood, 450k | 216 of 355 healthy young Swedes (14–34 y; 130 aged 14–16, 10 aged 31–34); 168 men, 48 women; 164 normal weight, 43 overweight, 9 obese | Betas in the series matrix (NOOB, BMIQ; SNP, cross-reactive and sex-chromosome probes removed); IDATs in RAW.tar (3.0 GB) | **Used.** Deposited betas, as the other cohorts of `level_match_reference`. One array per person. Uppsala cohort: no overlap with the US, Italian or Kabuki/CHARGE donors. |
| [GSE105018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105018) (E-Risk) | Whole blood, 450k | 1,658 twins, all aged 18 | Betas only in a 5.3 GB CSV (not in the series matrix); IDATs 13 GB | **Meets, backup.** Twin pairs (MZ and DZ) need family-level resampling; single age. |
| [GSE152027](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152027) (IoPPN) | Blood, 450k | 120 of 203 controls (18–62 y; 194 with age); the rest psychosis cases | Processed signals in a 2.4 GB CSV; RAW.tar 183 MB | **Meets, backup.** Controls only; betas outside the series matrix. |
| [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) (Johansson 2013) | Whole blood, 450k | 123 arrays of 732 (14–94 y); 421 people, some arrayed twice | SWAN betas in two 1.4 GB supplementary matrices; IDATs 6.1 GB | Meets on age band; replicate arrays need a person key. Not used. |
| [GSE41169](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41169) | Whole blood, 450k | 24 of 33 controls | IDATs | **Fails:** n < 50. |
| [GSE56105](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56105) (BSGS) | Peripheral blood lymphocytes, 450k | 37 of 614 (mostly adolescents 10–17 and their parents) | Betas, IDATs | **Fails:** n < 50 in the band. |
| [GSE53045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53045) | PBMC, 450k | Ages not in GEO (111 smokers and non-smokers) | Betas, IDATs | **Fails:** no ages. |

**Result.** With the GSE73103 young adults as reference, matching on young SD halves the housekeeping gap as with the
adult references, and their housekeeping CpGs vary less between donors than matched CpGs (0.19 percentage points,
against 0.03 in the GSE97362 children); numbers in `docs/CHANGELOG.md` (2026-11-28) and
`data/output/level_match_reference/`. A rebuild of GSE73103 from its IDATs, or E-Risk (GSE105018) as a second
young-adult reference, would separate the depositors' BMIQ normalisation from the age band.

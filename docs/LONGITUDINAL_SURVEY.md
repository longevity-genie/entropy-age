# Longitudinal and outcome-linked DNA methylation cohorts

Survey of 2026-11-12, for testing whether a per-person "promoter shield" methylation score predicts
health. Wanted: human blood (or other tissue) methylation on Illumina 450k/EPIC or WGBS, with age, and
either later health outcomes (death, incident disease, frailty change) or repeated measurements of the
same people over years.

How it was checked. Every GEO accession below was fetched from GEO (`acc.cgi`, series and all sample
records, brief SOFT) and its sample characteristics were tabulated; the "outcomes in public metadata"
column lists only fields that are actually in the sample characteristics (what the series matrix
carries). SATSA was checked on BioStudies/ArrayExpress (SDRF and the phenotype file of S-BSST1206).
dbGaP studies were fetched by phs; the EGA study EGAS00001000910 was fetched on ega-archive.org.
Discovery: GEO E-utilities over the 450k/EPIC/EPICv2 platforms (GPL13534, GPL21145, GPL23976,
GPL33022) with "mortality", "survival", "follow-up", "longitudinal", "repeated", "frailty",
"incident", "prospective" and cohort names (386 series screened by title), plus ArrayExpress search.
No data were downloaded except the 100 kB SATSA phenotype table, read to count people and visits.

Bottom line: **no open cohort has large-scale adult all-cause mortality with full time-to-event in its
public metadata.** The open data give (a) one large mortality signal in a young, special population
(SJLIFE childhood-cancer survivors), (b) one small old-age mortality set (Vitality 90+), (c) incident
cancer (EPIC-Italy, stroke cohort), and (d) good repeated-measures data (SATSA, GENOA, Danish twins).
General-population mortality, frailty and incident disease at scale need a controlled application.

## Open cohorts (ranked by usefulness for the project)

| # | Cohort | Accession | Tissue | Platform | People / samples | Age | Follow-up | Outcomes in public metadata | Paper |
|---|---|---|---|---|---|---|---|---|---|
| 1 | St. Jude Lifetime (SJLIFE), childhood-cancer survivors + controls | GSE314261 | PBMC (8 saliva) | EPIC | 4,520 / 5,013 | 0.2-71 (4,120 samples adult) | deaths up to 21 y after draw; 378 people with 2-5 draws (median gap 5.4 y) | `age at death` (294 deaths), age, sex, diagnosis, treatment doses, survivor/control. Censoring date for the living not given | 10.1038/s41467-026-68689-6 |
| 2 | SATSA (Swedish Adoption/Twin Study of Aging) | S-BSST1206 (450k+EPIC, preprocessed); E-MTAB-7309 (450k raw IDATs) | whole blood leukocytes | 450k + EPIC | 596 / 1,469 (392 with 2+ visits, up to 8) | 48-99 | up to 20 y between draws (median span 11 y) | age at draw, sex, zygosity, smoking, wave. Mortality not in public metadata | 10.1080/15592294.2018.1526028 |
| 3 | Vitality 90+ nonagenarians (Tampere) | GSE68194 | PBMC | 450k | 115 phenotyped / 122 | 90+ | ~2.5 y (37-1,023 days) | `status (0=alive; 1=died)` (39 deaths), `follow-up time in days`, BMI, MMSE, cell subsets | 10.18632/oncotarget.8278 |
| 4 | EPIC-Italy (HuGeF) | GSE51032 | buffy coat | 450k | 845 / 845 | 34-72 | to 2010; time to diagnosis 0-14 y | `cancer type (icd-10)` (421 incident cancers, mostly breast/colorectal), `time to diagnosis`. Mortality not in public metadata | no paper linked on GEO; used in 10.1186/s13148-015-0104-2 |
| 5 | GENOA (African-American hypertensive sibships) | GSE210256 (superseries of GSE210254, GSE210255) | peripheral blood leukocytes | 450k + EPIC | 1,417 / 1,812 (293 with two visits) | 20-91 | 2.5-9.8 y between phases (median 5.5 y) | age, sex, phase, pedigree. Outcomes not in public metadata | 10.1038/s41467-023-37961-4 |
| 6 | Danish twins (LSADT), 10-year repeat | GSE73115 | whole blood | 450k | 86 / 180 | elderly (age not in metadata; birth cohort) | 10 y (1997, 2007) | individual id, sex, draw year (in title). No outcomes | 10.1093/ije/dyw132 |
| 7 | Breakthrough Generations Study | GSE61151 | peripheral blood | 450k | 92 women / 188 | 35-83 | ~6 y (two draws) | age at draw, round, subject id. No outcomes (all "normal") | 10.1158/1055-9965.EPI-14-0767 |
| 8 | Stroke patients, later cancer | GSE280206 | whole blood (24 h post-stroke) | 450k + EPIC | 648 / 648 | not in metadata | 15 y | `incident_cancer` yes/no (83). No age, no time-to-event in metadata | 10.1186/s13073-024-01408-2 |
| 9 | Elderly T2D, pre-symptomatic dementia | GSE62003 | whole blood | 450k | ~36 / 70 | 67-88 | 18 months, two draws | case (developed dementia) vs control at TP1/TP2 (in source name) | 10.1016/j.neurobiolaging.2014.12.023 |
| 10 | Older adults, MCI pilot | GSE190540 | blood leukocytes | EPIC | 45 / 90 | 60-84 | 2 y, two draws | MCI case/control, baseline/follow-up flags | 10.3233/JAD-230093 |
| 11 | B-PROOF folic acid / B12 trial | GSE74548 | buffy coat | 450k | 87 / 174 | 65-75 | 2 y, two draws | arm, time point, folate, B12, homocysteine. No outcomes | 10.1186/s13148-015-0154-5 |
| 12 | Vitality 90+ (age-change set) | GSE58888 | PBMC | 450k | 143 / 143 | 90+ and 19-30 | none | parental age at death only; overlaps GSE68194 | 10.1186/s12864-015-1381-z |
| 13 | NEAT001 HIV, pre/post ART | GSE217633 | whole blood | EPIC | 229 / 413 | 19-76 | 96 weeks | HIV RNA, years with HIV, pre/post ART | (PMID 36640455) |
| 14 | Severe COVID-19, hospital course | GSE226206 | buffy coat | EPIC | 29 / 102 | 27-86 | weeks (up to 4 points) | days in hospital, days in ICU, treatments | 10.1016/j.cmet.2023.03.015 |
| 15 | COVID-19, 3-year repeat | GSE293007 | whole blood | EPIC | 54 / 108 | adults (age not in metadata) | 3 y, two draws | COVID-19 positive/negative. Physical-function scores not in sample metadata | 10.1007/s11357-025-01635-4 |
| 16 | MESA monocytes (Exam 5) | GSE56046 | CD14+ monocytes | 450k | 1,202 / 1,202 | 44-83 | none in GEO | carotid plaque score, coronary calcium (cross-sectional). Events only via dbGaP phs000209 | 10.1038/ncomms6366 |
| 17 | PEG Parkinson's (wave 1) | GSE72774; GSE111629 | whole blood | 450k | 508 / 508; 572 / 572 | ~35-95 | none | PD vs control, levodopa, smoking. No follow-up in metadata | 10.18632/aging.100859 |

Notes on the open set:

- SJLIFE also has two earlier, overlapping deposits (GSE169156, 2,052 survivors; GSE197674, 2,138
  survivors, EPIC) with treatment exposures but no death field; GSE314261 supersedes them.
- GSE68194 and GSE58888 share most of the same nonagenarians; use GSE68194 for mortality.
- Checked and dropped (no outcome, no repeat, or not adult aging): GSE55763 (2,711 London adults,
  cross-sectional + technical replicates), GSE107082/GSE117860 (VACS HIV; the VACS frailty index is
  not in the metadata), GSE100227 (Australian twins, cross-sectional), GSE69270 (Young Finns, one
  draw), GSE174422 (Sister Study duplicates), GSE51388 (9-month repeats, 8 people), GSE227809 (stress
  test, one week), GSE254135 / GSE142512 (children), GSE315366 / GSE315451 (patients at diagnosis).
- No InCHIANTI, Normative Aging, Generation Scotland, Lothian, Rhineland, BASE-II, Dunedin, TwinsUK
  or CALERIE methylation series was found on GEO by name search.

## Controlled cohorts (application needed)

| Cohort | Accession / route | Tissue, platform | People (samples) | Age | Follow-up / repeats | Outcomes (per study pages) | Apply |
|---|---|---|---|---|---|---|---|
| Generation Scotland (MeGS) | GS Access Committee | blood, EPIC | ~18,869 (+796 EPIC repeats) | 18-99 | NHS record linkage; deaths via National Records of Scotland | mortality, incident disease (ICD codes) | https://www.ed.ac.uk/generation-scotland/for-researchers/access |
| Lothian Birth Cohorts 1921 / 1936 | EGA EGAS00001000910 (dataset EGAD00010000604); phenotypes dbGaP phs000821 | blood, 450k | 2,195 samples; 478 with 2+ waves | 70-90 | waves at 70/73/76 and 79/87/90; mortality follow-up | death ("DNA methylation age and mortality" in phs000821) | https://www.ed.ac.uk/lothian-birth-cohorts/data-access-collaboration |
| Framingham Offspring Exam 8 + Gen3 Exam 2 | dbGaP phs000724 (substudy of phs000007) | blood, 450k (+850k) | ~4,000 | 24-90s | decades of event surveillance | mortality, CVD, cancer, dementia in phs000007 | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000724 |
| Normative Aging Study | dbGaP phs000853 | blood, 450k | ~700 men, up to 3 visits | 55-100 | visits every 3-5 y; primary phenotype mortality | mortality | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000853 |
| Women's Health Initiative (BAA23, CHD) | dbGaP phs001335 (parent phs000200) | blood, 450k | 904 (472 CHD cases, 432 controls) | 50-79 | to 2015+ | incident CHD; WHI outcomes incl. death | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001335 |
| Health and Retirement Study (VBS 2016) | NIAGADS + HRS sensitive-data order | blood, EPIC | ~4,000 | 50-100 | biennial interviews; death, ADL/frailty | mortality, disability, cognition | https://www.niagads.org/ (and HRS sensitive health data form) |
| NHANES 1999-2002 | NCHS Research Data Center (CpG file restricted; clocks public) | blood, EPIC | ~2,500 (not stated on page) | 50+ | linked mortality files | mortality via NDI linkage | https://wwwn.cdc.gov/nchs/nhanes/dnam/default.aspx |
| MESA (events for GSE56046 people) | dbGaP phs000209 | monocytes, 450k | 1,202 | 44-83 | CVD event follow-up | MI, stroke, death | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000209 |
| InCHIANTI | InCHIANTI study proposal | blood, 450k/EPIC | 699 people, 1,721 samples (1998/2007/2013) | 21-102 | 3 waves over 15 y; mortality | mortality, frailty, physical function | https://www.inchiantistudy.net/ |
| TwinsUK | TwinsUK data access | blood, 450k/EPIC | several thousand, repeats | adults | repeat visits; NHS linkage | frailty index, disease | https://twinsuk.ac.uk/researchers/collaborations/ |
| Dunedin Study | concept paper to the Dunedin unit | blood, EPIC | ~800-900 at ages 26/38/45 | 26-45 | 3 repeat draws | pace of aging, health | https://dunedinstudy.otago.ac.nz/ |
| Rhineland Study | Rhineland Study data access | blood, EPIC | several thousand | 30+ | ongoing follow-up | not checked | https://www.rheinland-studie.de/ |
| BASE-II | BASE-II application (site blocked scripted fetch) | blood, EPIC | ~1,800, two waves | 20-85 | ~7 y between waves | frailty, function | https://www.base2.mpg.de/en |

Accessions in this table were confirmed on their pages (phs000724, phs000853, phs001335, phs000200,
phs000209, phs000821, EGAS00001000910). InCHIANTI, TwinsUK, Dunedin, Rhineland, BASE-II, HRS and
Generation Scotland have no public accession for the methylation data; their sizes come from their
own pages or data-resource papers and are approximate.

## Recommended order

1. **SJLIFE GSE314261**: the only open set with hundreds of deaths and per-sample age at death, plus
   378 people with repeated draws. Caveats: survivors of childhood cancer, mostly young adults, and
   censoring dates for the living are not in GEO (treat as vital status at data freeze, or ask the
   authors). Adjust for treatment exposures, which are in the metadata.
2. **SATSA S-BSST1206**: the best open repeated-measures set (596 people, up to 8 draws, 20 years,
   ages 48-99). Tests whether the score drifts with age within a person; outcomes need an application.
3. **Vitality 90+ GSE68194**: small, but true old-age mortality with follow-up days in the metadata;
   a quick positive or negative control.
4. **EPIC-Italy GSE51032**: incident cancer with time to diagnosis in 845 people; a disease-incidence test.
5. **GENOA GSE210256 + Danish twins GSE73115 + Breakthrough GSE61151**: extra within-person change
   over 5-10 years, different ancestry and ages. GSE280206 adds incident cancer (no age in metadata).
6. Apply in parallel, by expected yield: **Generation Scotland** (largest, mortality + incident
   disease), **Lothian Birth Cohorts** (EGA + dbGaP, repeats + mortality), **Framingham phs000724**
   and **NAS phs000853** (dbGaP, one application route), then **HRS/NHANES** and **InCHIANTI**.

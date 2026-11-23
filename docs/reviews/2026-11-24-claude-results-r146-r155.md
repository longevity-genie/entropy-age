# Independent review: entropy-age, R146-R155 and the mouse null (2026-11-24)

Read-only. Scripts are in `/tmp/claude-41416/-home-uq2347-desk/6597363e-b141-49ee-94ae-8f14faceca9f/scratchpad/review2/claude_scratch/` (`mouse_test.py`, `snr.py`, `bxd.py`, `genes.py`, `pgene.py`, `annot.py`, `idat.py`). All numbers are the kl_slope housekeeping gap in SD unless stated otherwise.

## Findings, ranked

**1. HIGH: the gene-level p (`gene_clustered_p`, stochasticity.py:336) does not test age. It rejects with ages shuffled.**
- I ran field liver ("liver, all") through `compare_groups` with ages randomly permuted. In 5 of 5 permutations p_gene was between 3e-17 and 0.012 (the permuted gaps ran from −0.116 to +0.038). The real data give gap −0.028 with p_gene 8e-4.
- Why: the test takes each CpG's slope as fixed and resamples only genes. Nothing reflects sampling of donors or arrays. Any shared sample-level structure, or any set-specific difference in noise, comes out as "significant".
- The df (genes − 1) and the CR1 ratio-estimator SE are computed correctly. The assumption is what fails, and so does the earlier CpG-label permutation it replaced.
- Affected: R146 (blood p 5·10⁻²⁴⁴; paper.md:83 "p < 10⁻⁵ for each"), SATSA "p < 10⁻¹⁰" (paper.md:472), and global_fdr.py:42, which swaps p_permutation for p_gene.
- Fix: base inference on donor or age permutation, or on donor-resampling intervals. Present p_gene as a gene-heterogeneity statistic, or drop it from the paper and from the FDR ledger.

**2. HIGH: the TCGA exclusion test cannot detect a field effect as designed** (scripts/tcga_tumour_cpgs.py; paper.md:182-187; R153).
- TCGA "solid tissue normal" samples are themselves tumour-adjacent. Changes that adjacent tissue shares with the tumour (the field) therefore cancel in tumour − normal. The list removes tumour-only changes, which is the wrong set.
- The removal is also unbalanced. Using n_set and cpgs in the tables:
  - At |Δβ| ≥ 0.1, housekeeping CpGs fall 30,174 → 29,041 (−3.8 %) while background CpGs fall 46,305 → 36,717 (−20.7 %).
  - At ≥ 0.2, the losses are −0.7 % and −6.4 %.
  - So the test mostly reshapes the matched background and the SD that scales the gap.
- There is no control that removes a random set of the same size and level/island make-up.
- "shrank the reversal without removing it" is too strong at d10: the donor interval is −0.000 to +0.190.
- Fix:
  - Define field CpGs as healthy lung (GSE175458, now run) versus adjacent normal.
  - Add a random, stratum-matched removal control.
  - Report raw (unscaled) gaps next to the SD-scaled ones.
- Related: series_lung_healthy gives −0.263 with chip held (donors −0.43 to −0.15), a strong gap in 143 healthy lungs. paper.md:187 still says it "is being processed". This result matters more for the lung story than the TCGA lists do.

**3. HIGH (interpretation): the mouse-liver sentence frames the null as a contradiction when it is an age-window and platform effect** (paper.md:209-212, R152). See the verdict below. In BXD itself, chow mice aged ≤10 months give −0.013, not −0.30.

**4. MODERATE: the mouse island flag is empty on MM285 probes** (data/mouse.py `read_platform_table` uses `CpG.island`).
- In GPL32543, `CpG.island` = "yes" only for the Mammal40-derived probes. Of the 72,696 promoter probes present on MM285, 72,654 are "no".
- In field liver, the island fraction among complete CpGs is 0.000. The "island × level" matching is therefore level-only in the series mouse path and in mouse_tissues (Huang).
- `CGIslandNess` ("+" for 54 % of probes, 87 % of housekeeping CpGs) appears to be the populated column.
- Effect here is small: field liver goes −0.028 → −0.035 with CGIslandNess × level. Methods should not claim island matching for these series.
- Fix: read `CGIslandNess` for MM285 probes and rerun Huang and the two new series.

**5. MODERATE: the rebuild discards the Mammal40 half of 320k chips.**
- The IDATs of GSE280472 (hearts) and GSE279437 (Huang) hold 404,407 addresses: they are HorvathMammalMethylChip320 chips (GPL32543).
- SeSAMe read them as MM285 (296,070 probes; QC `platform` = MM285). The conserved Mammal40 probes were never read, and those are the CpG universe of BXD and the species analysis.
- GSE269932 is a true MM285 chip (361,821 addresses; GPL30650).
- So "mammalian array GPL32543" in the brief and log is not the probe set that was analysed.

**6. MODERATE: today's repeated-array handling picks betas and QC from possibly different arrays** (reprocessed.py:78-86, sesame_preprocess.R:92-96).
- In GSE79100, all 31 kidneys have 3 arrays, spread over 2 chips.
- The reader keeps the first array's betas, but `passed` is "any array of the sample passed". One kidney has 1 of 3 arrays passing and its first array fails.
- `chip_of` (series.py:194) is a dict over all QC rows, so the last array's chip is used. For 30 kidneys the chip intercept may belong to a different array.
- Fix: carry the chosen array's prefix through to QC, and apply `passed` and `chip` per array.

**7. MODERATE: kidney and healthy-lung tables are quick-scale but sit in the canonical output folders.**
- series_kidney_normal and series_lung_healthy have `donor_resamples` = 10 (8 with chip), so gene resamples are quick-scale too.
- There is no R155 entry in RUN_LOG and nothing records the scale in the table.
- Fix: write the scale into the table, or refuse quick runs into the paper-scale folders.
- Separately, kidney (probably nephrectomy margins, so also tissue next to a tumour) shows −0.12, against the lung reversal. Worth a line once the run is at paper scale.

**8. MODERATE: the hearts series has almost no age signal to read.**
- In GSE280472 (0/1 age, 36 mice, 3 chips), the SD of kl_slope across CpGs is only 1.1-1.5× its SD under permuted ages; slopes are mostly noise.
- Permuted-age housekeeping gaps run from −0.057 to +0.059, which brackets the observed −0.045.
- The gene-resampling interval (−0.065 to −0.025) looks resolved but is not. The paper correctly quotes the donor interval.
- In field liver, by contrast, the slope SD is 3-5× the permuted SD: there is a real age signal, and the set simply does not differ.

**9. LOW-MODERATE: the species lifespan fit is plain OLS over 139 groups from 65 species** (species_shield.py:137), with no clustering by species and no phylogeny. The "+0.031, p 0.04" is anti-conservative. Use species-clustered SEs, or the species_phylo fit.

**10. LOW:**
- **Donor-resample guard** (series.py:147): it counts dummy columns, not the design rank. Survivors are a conditional subset; the code already acknowledges this.
- **chip_overlap.csv column names**: the mouse and human series inherit culture names (`cultures`, `doublings_min`, `variance_left_within_line_and_chip`).
- **Field-liver confounding**: in "liver, all", age is largely confounded with chip and environment (22 % of age variance left). 12 of 25 chips hold one array, and most of the oldest lab mice are on them, so "−0.00 with chip" rests mainly on field mice.

**Checked and supported by the tables:**
- R147 (myometrium, including EPIC v2 chip-held "no interval", 2/200 resamples).
- R148 guard.
- R149: 139 groups, 87.8 % negative, mean −0.123, lifespan coefficient p 0.042.
- R154 breast: −0.038 (−0.111 to −0.005) and +0.001 (−0.046 to +0.037), chip held.
- R152 and R153 numbers as quoted.
- Set-CpG survival under "complete CpGs only" is high:
  - field liver housekeeping 10,716 of 11,994 (89 %), core 2,098 of 2,334, non-essential 2,420 of 3,305 (73 %);
  - hearts housekeeping 11,093.
- Probe IDs match between the MM285 rebuilds and the annotation exactly: 72,696 promoter probes, no suffix collapse. Probe-ID mismatch is excluded.

## Verdict on the mouse question

**The two new series measure a different thing, and neither contradicts the earlier gaps.**

**Field liver: age window first, platform second.**
- GSE269932 ages are median 135 days and 90th percentile 263 days. Only 7 of 90 mice are over 302 days, and the young reference (youngest 20 %) includes 26-60-day juveniles. Its slopes therefore largely measure maturation and puberty (liver sex-dimorphic methylation sets in then; sex enters only as an intercept).
- Within BXD (Mammal40, same code), the window decides the answer:
  - all chow females: −0.305;
  - ≤17.5 months: −0.277;
  - **≤10 months (51 mice): −0.013**.
- The gap does not depend on island matching (−0.292 level-only), strain (−0.289 without) or reference choice (−0.311 with the runner's youngest-20 % baseline).
- The CpG universe also differs. BXD uses 4,221 conserved Mammal40 promoter CpGs (61 % island). The MM285 series use about 65,000, and the two sets share almost no probes (42).
- On MM285, even Huang's wide young/aged contrast gives liver only −0.049, against −0.33 for mouse liver on Mammal40 in the species data (903 mice).
- Restricting the MM285 housekeeping CpGs to the genes of BXD's set does not restore the gap (+0.034). So the platform effect lies in the CpG positions, not in which genes are included.

**Hearts: no power plus post-mitotic tissue.**
- Slopes are mostly noise (finding 8).
- Heart is near zero on every platform: Mammal40 −0.01, Huang +0.08, GSE225166 −0.05, Meer −0.07.

**Matching code.** The island flag is broken on MM285 (finding 4), but fixing it moves the gap by under 0.01 SD, so the code is not the cause.

**The one run that would settle the platform question:** rebuild the Huang 320k-chip IDATs (GSE279437) keeping the Mammal40 addresses, then compute the liver (and heart) gap on the conserved promoter CpGs and on the MM285 promoter CpGs of the same arrays and animals.
- If conserved CpGs give about −0.3 where MM285 CpGs give about −0.05, the CpG universe explains the gap between platforms.
- The age-window explanation for the field series is already settled by the BXD ≤10-month result above. The paper should say so instead of contrasting the field liver with the BXD −0.30.

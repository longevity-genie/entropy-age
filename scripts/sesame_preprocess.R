# One preprocessing for every Infinium series, from raw IDATs (roadmap RM41; asked by antonkulaga 2026-11-18).
#
# SeSAMe (Zhou et al. 2018, doi:10.1093/nar/gky691), SeSAMe's recipe for each platform (prep_code below); human:
#   prepSesame "QCDPB": Q qualityMask (design-flagged probes), C inferInfiniumIChannel (channel switches),
#   D dyeBiasNL (non-linear dye bias), P pOOBAH (detection p from out-of-band signal; failing probes masked),
#   B noob (background from out-of-band signal). Betas of masked probes are NA.
# Per sample, before preprocessing: sesameQC_calcStats (detection, numProbes, intensity, channel, dyeBias) and the
# bisulfite-conversion index GCT (bisConversionControl; near 1 when conversion is complete).
# The chip (Sentrix barcode) and array position are read from GEO's file names (GSMnnn_<chip>_<RxxCxx>_Grn.idat).
#
# Usage (inside the sesame environment, as a Slurm job):
#   Rscript scripts/sesame_preprocess.R IDAT_DIR OUT_PREFIX [PLATFORM|auto] [WORKERS] [CHUNK]
# Writes OUT_PREFIX.betas.chunkNNN.parquet (columns: probe, then one per GSM) with OUT_PREFIX.qc.chunkNNN.csv, then
# OUT_PREFIX.qc.csv and
# OUT_PREFIX.session.txt (package versions). Probe ids are the platform's (cg... for HM450/EPIC; EPICv2 and the
# mouse arrays carry suffixes, which the Python side collapses or maps). A stopped run (a preempted job) resumes:
# chunks already written are skipped. A pair SeSAMe cannot read is listed in OUT_PREFIX.failed.chunkNNN.csv.

suppressPackageStartupMessages({
  library(sesame)
  library(arrow)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
idat_dir <- args[1]
out <- args[2]
platform <- if (length(args) >= 3 && args[3] != "auto") args[3] else NULL
workers <- if (length(args) >= 4) as.integer(args[4]) else 1L
chunk <- if (length(args) >= 5) as.integer(args[5]) else 200L

pfxs <- searchIDATprefixes(idat_dir)
if (length(pfxs) == 0) stop("no IDAT pairs under ", idat_dir)
message(length(pfxs), " IDAT pairs in ", idat_dir)

parse_name <- function(pfx) {
  base <- basename(pfx)
  gsm <- regmatches(base, regexpr("^GSM[0-9]+", base))
  m <- regmatches(base, regexec("([0-9]{9,})_(R[0-9]+C[0-9]+)", base))[[1]]
  data.frame(
    prefix = base,
    sample = if (length(gsm)) gsm else base,
    chip = if (length(m) == 3) m[2] else NA_character_,
    position = if (length(m) == 3) m[3] else NA_character_
  )
}

# SeSAMe's own recommendation per platform (sesame vignette, "prepSesame" table; methods review 2026-11-18, RM43):
# human arrays QCDPB; MM285 (mouse) TQCDPB, T masking probes hit by the inferred strain's variants; Mammal40 in a
# non-human organism SHCDPB, S the species mask and H keeping cg/ch probes only. SESAME_PREP overrides for every array.
prep_code <- function(pf) {
  forced <- Sys.getenv("SESAME_PREP", "")
  if (nzchar(forced)) return(forced)
  if (pf == "MM285") return("TQCDPB")
  if (pf == "Mammal40") return("SHCDPB")
  "QCDPB"
}

one <- function(pfx) {
  sdf <- if (is.null(platform)) readIDATpair(pfx) else readIDATpair(pfx, platform = platform)
  qc <- as.data.frame(sesameQC_calcStats(sdf, c("detection", "numProbes", "intensity", "channel", "dyeBias")))
  qc$GCT <- tryCatch(bisConversionControl(sdf), error = function(e) NA_real_)  # bisulfite conversion, ~1 is good
  code <- prep_code(sdfPlatform(sdf))
  if (grepl("S", code)) {  # detection over the probes the species mask keeps, not the whole multi-species array
    kept <- prepSesame(sdf, "SH")
    qc$frac_dt_species <- sum(!prepSesame(kept, "P")$mask) / max(1, sum(!kept$mask))
  }
  prepped <- prepSesame(sdf, code)
  list(beta = getBetas(prepped), qc = cbind(parse_name(pfx), qc, platform = sdfPlatform(sdf), prep = code))
}

param <- if (workers > 1) MulticoreParam(workers) else SerialParam()
groups <- split(seq_along(pfxs), ceiling(seq_along(pfxs) / chunk))
for (i in seq_along(groups)) {
  betas_file <- sprintf("%s.betas.chunk%03d.parquet", out, i)
  qc_file <- sprintf("%s.qc.chunk%03d.csv", out, i)
  if (file.exists(betas_file) && file.exists(qc_file)) {  # written by an earlier, stopped run: resume after it
    message(sprintf("chunk %d of %d: already written", i, length(groups)))
    next
  }
  res <- bplapply(pfxs[groups[[i]]], function(p) tryCatch(one(p), error = function(e) conditionMessage(e)),
    BPPARAM = param)
  failed <- vapply(res, is.character, logical(1))
  if (any(failed)) {  # an unreadable pair is listed, never silently dropped
    write.csv(data.frame(prefix = basename(pfxs[groups[[i]]][failed]), error = unlist(res[failed])),
      sprintf("%s.failed.chunk%03d.csv", out, i), row.names = FALSE)
    res <- res[!failed]
  }
  # bind by probe name: a chunk can hold arrays of different platforms (HM450 beside EPIC), whose beta vectors differ
  # in length and order; binding by position put one platform's values under the other's probe ids (Codex review
  # 2026-11-25, GSE314841 and GSE296589)
  if (any(vapply(res, function(r) is.null(names(r$beta)), logical(1)))) stop("a beta vector without probe names")
  probes <- unique(unlist(lapply(res, function(r) names(r$beta))))
  betas <- vapply(res, function(r) unname(r$beta[probes]), numeric(length(probes)))
  if (!is.matrix(betas)) betas <- matrix(betas, ncol = length(res))
  rownames(betas) <- probes
  meta <- do.call(rbind, lapply(res, `[[`, "qc"))
  colnames(betas) <- meta$sample
  # one array per sample: a sample arrayed more than once under one GSM (GSE79100: three arrays per kidney) keeps its
  # first array here; the QC table still lists every array
  repeated <- duplicated(colnames(betas))
  if (any(repeated)) message(sprintf("chunk %d: %d repeated arrays of a sample left out (first kept)", i, sum(repeated)))
  betas <- betas[, !repeated, drop = FALSE]
  frame <- data.frame(probe = rownames(betas), betas, check.names = FALSE)
  write_parquet(frame, paste0(betas_file, ".part"))  # renamed only when whole, so a stop leaves no half file
  file.rename(paste0(betas_file, ".part"), betas_file)
  write.csv(meta, qc_file, row.names = FALSE)
  message(sprintf("chunk %d of %d: %d samples", i, length(groups), ncol(betas)))
}
qc_files <- sprintf("%s.qc.chunk%03d.csv", out, seq_along(groups))
write.csv(do.call(rbind, lapply(qc_files, read.csv)), paste0(out, ".qc.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(out, ".session.txt"))

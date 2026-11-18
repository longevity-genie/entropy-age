# One preprocessing for every Infinium series, from raw IDATs (roadmap RM41; asked by antonkulaga 2026-11-18).
#
# SeSAMe (Zhou et al. 2018, doi:10.1093/nar/gky691), the same recipe for every series and platform:
#   prepSesame "QCDPB": Q qualityMask (design-flagged probes), C inferInfiniumIChannel (channel switches),
#   D dyeBiasNL (non-linear dye bias), P pOOBAH (detection p from out-of-band signal; failing probes masked),
#   B noob (background from out-of-band signal). Betas of masked probes are NA.
# Per sample, before preprocessing: sesameQC_calcStats (detection, intensity, numProbes, dyeBias, bisConversion).
# The chip (Sentrix barcode) and array position are read from GEO's file names (GSMnnn_<chip>_<RxxCxx>_Grn.idat).
#
# Usage (inside the sesame environment, as a Slurm job):
#   Rscript scripts/sesame_preprocess.R IDAT_DIR OUT_PREFIX [PLATFORM|auto] [WORKERS] [CHUNK]
# Writes OUT_PREFIX.betas.chunkNNN.parquet (columns: probe, then one per GSM), OUT_PREFIX.qc.csv and
# OUT_PREFIX.session.txt (package versions). Probe ids are the platform's (cg... for HM450/EPIC; EPICv2 and the
# mouse arrays carry suffixes, which the Python side collapses or maps).

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

one <- function(pfx) {
  sdf <- if (is.null(platform)) readIDATpair(pfx) else readIDATpair(pfx, platform = platform)
  qc <- as.data.frame(sesameQC_calcStats(
    sdf, c("detection", "intensity", "numProbes", "dyeBias", "bisConversion")
  ))
  prepped <- prepSesame(sdf, "QCDPB")
  list(beta = getBetas(prepped), qc = cbind(parse_name(pfx), qc, platform = sdfPlatform(sdf)))
}

param <- if (workers > 1) MulticoreParam(workers) else SerialParam()
groups <- split(seq_along(pfxs), ceiling(seq_along(pfxs) / chunk))
qcs <- list()
for (i in seq_along(groups)) {
  res <- bplapply(pfxs[groups[[i]]], one, BPPARAM = param)
  betas <- do.call(cbind, lapply(res, `[[`, "beta"))
  meta <- do.call(rbind, lapply(res, `[[`, "qc"))
  colnames(betas) <- meta$sample
  frame <- data.frame(probe = rownames(betas), betas, check.names = FALSE)
  write_parquet(frame, sprintf("%s.betas.chunk%03d.parquet", out, i))
  qcs[[i]] <- meta
  message(sprintf("chunk %d of %d: %d samples", i, length(groups), ncol(betas)))
}
write.csv(do.call(rbind, qcs), paste0(out, ".qc.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(out, ".session.txt"))

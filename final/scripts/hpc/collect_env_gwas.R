#!/usr/bin/env Rscript
################################################################################
## COLLECT THE 9 ENVIRONMENTAL / CONTROL GWAS SCANS INTO ONE FOLDER, WITH FDR
##
## Reads   out/<trait>/raw_GWAS_<MODEL>_3PC_<trait>.txt     (27 files)
## Writes  out/env_GWAS_maize_landrace/
##           <trait>__<MODEL>.csv.gz   retained tail + BH FDR
##           summary_env_gwas.csv      one row per trait x model
##           pseudoQTN/<trait>_FarmCPU_pseudoQTN.csv
##
## WHY NOT JUST COPY THE RAW FILES: 27 x ~45 MB = ~1.2 GB, and 99.9% of it is
## p-values nobody will look at. Retaining P <= 1e-5 leaves a few thousand rows
## per file. The full scans stay on the HPC if they are ever needed again.
##   RETAIN_P is 1e-3, not 1e-5: with m = 495,755 a BH cutoff only exceeds 1e-5
##   once a trait has ~100 significant SNPs, which soil N plausibly does. At
##   1e-3 the cutoff would need ~9,900 hits to breach, and the files are still
##   only ~600 rows. Being generous here is free; being tight cost an HPC rerun.
##
## THE HEADER IS WRONG IN THE RAW FILES -- READ THIS BEFORE EDITING
##   The header line names 8 columns:
##     SNP Chr Pos P.value MAF nobs H&B.P.Value Effect
##   but each data row carries only 7 values, e.g.
##     S1_46446 1 46446 0.69865820437902 3107 0.979352026743126 23.896312408373
##   3107 is the panel size, i.e. nobs -- MAF was never written. So the true
##   columns are SNP, Chr, Pos, P.value, nobs, H&B.P.Value, Effect.
##   This script therefore IGNORES the header and assigns names by position,
##   after checking the field count. Reading these files by header name gives
##   nobs as MAF on every row, silently.
##   MAF is recoverable from the genotype matrix; it is not re-derived here.
##
## BH IS VALID HERE. Unlike the p-filtered split_by_trait files, these are
## complete genome-wide scans, so m is the file's own row count and the whole
## p-value distribution is present.
##
## Usage:  Rscript collect_env_gwas.R [out_dir] [dest_dir]
################################################################################

suppressPackageStartupMessages(library(data.table))

args     <- commandArgs(trailingOnly = TRUE)
OUT_DIR  <- if (length(args) >= 1) args[1] else "/share/maize/ntanduk/GrainN/GWAS_env/out"
DEST     <- if (length(args) >= 2) args[2] else file.path(OUT_DIR, "env_GWAS_maize_landrace")

Q_FDR    <- 0.05
RETAIN_P <- 1e-3
MODELS   <- c("MLM", "BLINK", "FarmCPU")
COLS7    <- c("SNP", "Chr", "Pos", "P.value", "nobs", "HB.P.Value", "Effect")
COLS8    <- c("SNP", "Chr", "Pos", "P.value", "MAF", "nobs", "HB.P.Value", "Effect")

dir.create(DEST, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DEST, "pseudoQTN"), recursive = TRUE, showWarnings = FALSE)

## A trait directory is one that actually contains raw scans. Anything else --
## a previous run's destination folder, a stray subdirectory -- is skipped
## silently rather than reported as three MISSING models.
cand  <- basename(list.dirs(OUT_DIR, recursive = FALSE))
has_raw <- vapply(cand, function(d)
  length(list.files(file.path(OUT_DIR, d), pattern = "^raw_GWAS_.*\\.txt$")) > 0L, logical(1))
skipped <- cand[!has_raw]
traits  <- sort(cand[has_raw])
cat("traits found:", length(traits), "->", paste(traits, collapse = ", "), "\n")
if (length(skipped)) cat("skipped (no raw scans):", paste(skipped, collapse = ", "), "\n")
cat("\n")

summ <- list(); problems <- character()

for (tr in traits) for (md in MODELS) {
  f <- file.path(OUT_DIR, tr, sprintf("raw_GWAS_%s_3PC_%s.txt", md, tr))
  if (!file.exists(f)) { problems <- c(problems, paste(tr, md, "MISSING")); next }

  ## skip the header entirely -- it is wrong (see note above)
  d <- fread(f, skip = 1, header = FALSE, showProgress = FALSE)

  if (ncol(d) == 7L)      setnames(d, COLS7)
  else if (ncol(d) == 8L) setnames(d, COLS8)
  else { problems <- c(problems, sprintf("%s %s: %d columns, expected 7 or 8", tr, md, ncol(d))); next }

  d[, P.value := as.numeric(P.value)]
  d <- d[is.finite(P.value) & P.value > 0]
  m <- nrow(d)                                  # markers actually tested

  ## Benjamini-Hochberg over the complete scan
  setorder(d, P.value)
  d[, FDR := pmin(1, rev(cummin(rev(P.value * m / seq_len(m)))))]
  d[, fdr_sig    := FDR <= Q_FDR]
  d[, bonferroni := P.value <= 0.05 / m]

  n_fdr  <- d[, sum(fdr_sig)]
  n_bonf <- d[, sum(bonferroni)]
  cutoff <- if (n_fdr) max(d[fdr_sig == TRUE]$P.value) else NA_real_

  if (!is.na(cutoff) && cutoff > 0.5 * RETAIN_P)
    problems <- c(problems, sprintf("%s %s: FDR cutoff %.2e near RETAIN_P -- raise it", tr, md, cutoff))

  keep <- d[P.value <= RETAIN_P]
  fwrite(keep, file.path(DEST, sprintf("%s__%s.csv.gz", tr, md)), compress = "gzip")

  summ[[length(summ) + 1]] <- data.table(
    trait = tr, model = md, m_markers = m, n_kept = nrow(keep),
    q_fdr = Q_FDR, n_fdr_sig = n_fdr, fdr_cutoff_P = cutoff,
    p_bonferroni = 0.05 / m, n_bonferroni = n_bonf,
    n_cols_raw = ncol(d) - 3L, maf_present = "MAF" %in% names(d))

  cat(sprintf("  %-13s %-8s m=%7d  FDR-sig %5d (cutoff %.2e)  Bonf %5d  kept %5d\n",
              tr, md, m, n_fdr, cutoff, n_bonf, nrow(keep)))
}

pq <- list.files(OUT_DIR, pattern = "^FarmCPU_pseudoQTN_.*\\.csv$", recursive = TRUE, full.names = TRUE)
for (p in pq) file.copy(p, file.path(DEST, "pseudoQTN", basename(p)), overwrite = TRUE)

S <- rbindlist(summ)
fwrite(S, file.path(DEST, "summary_env_gwas.csv"))

cat("\n=== per-trait totals (union over models) ===\n")
print(S[, .(models = .N, FDR_sig = sum(n_fdr_sig), Bonferroni = sum(n_bonferroni)), by = trait])

cat(sprintf("\nwrote %d files to %s\n", nrow(S), DEST))
cat(sprintf("pseudoQTN files copied: %d\n", length(pq)))
if (length(problems)) {
  cat("\n!!! PROBLEMS !!!\n"); for (p in problems) cat("  ", p, "\n")
} else cat("no problems\n")

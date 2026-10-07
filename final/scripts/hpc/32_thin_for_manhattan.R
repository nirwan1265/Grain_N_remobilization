#!/usr/bin/env Rscript
################################################################################
## 32 -- THIN THE RAW SCANS FOR MANHATTAN PLOTTING
##
## Run on the HPC, where split_by_trait lives. For each (trait, model) pair it
## keeps every SNP with P < KEEP_ALL_BELOW plus a 1-in-THIN_EVERY subsample of
## the rest, and writes a small CSV. Thinning only removes points that would
## overplot into the same pixel band at the bottom of a Manhattan; every point
## that can be seen individually is retained.
##
## ~4.18M SNPs per scan -> roughly 210k rows out, a few MB per file, so the
## eight files transfer easily instead of ~12 GB of raw scans.
##
## Output: <OUT_DIR>/<trait>__<MODEL>.csv   with columns Chr, Pos, P
##
## Usage: Rscript 32_thin_for_manhattan.R [in_dir] [out_dir]
################################################################################

suppressPackageStartupMessages(library(data.table))

IN_DIR  <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
OUT_DIR <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/manhattan_thin"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) IN_DIR  <- args[1]
if (length(args) >= 2L) OUT_DIR <- args[2]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

KEEP_ALL_BELOW <- 0.01
THIN_EVERY     <- 25L
set.seed(1)

## trait -> model to plot. Seven loci are FarmCPU; E.Total is BLINK only, so
## plotting FarmCPU there would show a panel with no peak.
JOBS <- list(
  c("EHPRQ",   "FarmCPU"),
  c("LAV",     "FarmCPU"),
  c("A",       "FarmCPU"),
  c("Total",   "FarmCPU"),
  c("IVL",     "FarmCPU"),
  c("T",       "FarmCPU"),
  c("E.Total", "BLINK"),
  c("S",       "FarmCPU")
)

normalize_model <- function(x) {
  unname(c(MLM = "MLM", MLMM = "MLMM", BLINK = "BLINK",
           FARMCPU = "FarmCPU")[toupper(trimws(as.character(x)))])
}

for (j in JOBS) {
  tr <- j[1]; md <- j[2]
  f  <- file.path(IN_DIR, paste0(tr, ".csv"))
  if (!file.exists(f)) { warning("missing: ", f); next }

  message(sprintf("[%s] %s / %s", format(Sys.time(), "%H:%M:%S"), tr, md))
  d <- fread(f, select = c("Chr", "Pos", "P.value", "model"), showProgress = FALSE)
  d[, model := normalize_model(model)]
  d <- d[model == md]
  if (!nrow(d)) { warning("no ", md, " rows in ", tr); next }

  d[, Chr := as.integer(Chr)]
  d[, Pos := as.numeric(Pos)]
  d[, P   := as.numeric(P.value)]
  d <- d[!is.na(Chr) & Chr >= 1 & Chr <= 10 & is.finite(P) & P > 0]

  n_before <- nrow(d)
  keep <- d$P < KEEP_ALL_BELOW | (seq_len(nrow(d)) %% THIN_EVERY == 0L)
  d <- d[keep, .(Chr, Pos, P)]
  setorder(d, Chr, Pos)

  out <- file.path(OUT_DIR, sprintf("%s__%s.csv", tr, md))
  fwrite(d, out)
  message(sprintf("   %d -> %d rows (%.1f%%), min P = %.2e  -> %s",
                  n_before, nrow(d), 100 * nrow(d) / n_before, min(d$P), basename(out)))
}

message("done. copy ", OUT_DIR, " to the repo as final/results/manhattan_thin/")

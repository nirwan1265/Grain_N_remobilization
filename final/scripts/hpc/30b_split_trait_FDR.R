#!/usr/bin/env Rscript
################################################################################
## 30b -- SPLIT ONE TRAIT'S GWAS SCAN INTO PER-MODEL .RData, FDR-BASED
##
## Replaces 30_split_trait_by_model.R. Differences, and why:
##
##  1. THRESHOLD IS FDR, NOT A FIXED P. The declared threshold is
##     Benjamini-Hochberg q <= 0.05, taken from GAPIT's own `FDR` column.
##     That column was computed by GAPIT over every marker it tested, which is
##     the only place the full p-value distribution still exists -- it cannot be
##     reconstructed downstream from a p-filtered file. A recomputed BH is
##     written alongside it as a check; the two should agree.
##
##  2. THE RETENTION CUTOFF IS LOOSE ON PURPOSE. Script 30 kept only P <= 1/m,
##     which made the FDR set uncomputable without this re-run. Here we retain
##     P <= RETAIN_P (1e-5 by default), ~40x looser than 1/m, so that any future
##     change of threshold is a re-analysis rather than another HPC re-run.
##     Files stay small: a 12.5M-row trait keeps a few thousand rows.
##
##  3. EVERY THRESHOLD IS FLAGGED, NONE IS APPLIED DESTRUCTIVELY. Each retained
##     row carries fdr_sig, bonferroni and suggestive, so the Bonferroni numbers
##     in the current draft remain reproducible from the same files.
##
## Output (into <OUTPUT_DIR>, which should NOT be the Bonferroni by_model dir):
##   by_model/<trait>_<MODEL>.RData      object <trait>_<MODEL>
##   consensus/<trait>_all_model_common_SNP.RData
##   summary/<trait>_summary.csv
##
## Usage
##   Rscript 30b_split_trait_FDR.R <trait> [in_dir] [out_dir]
################################################################################

suppressPackageStartupMessages(library(data.table))

INPUT_DIR  <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
OUTPUT_DIR <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/by_model_fdr"

## Marker count. NOTE: script 30 used 4117796, which is a digit transposition of
## the value implied by the raw scan row counts (8,355,592 / 2 models). The
## manuscript Methods figure is being corrected to 4177796 to match. This value
## affects ONLY the Bonferroni and suggestive flags and the recomputed-BH check;
## the FDR threshold comes from GAPIT's column and does not depend on it.
M_SNP <- 4177796

Q_FDR    <- 0.05                 # declared threshold
ALPHA    <- 0.05
P_SUGG   <- 1 / M_SNP
P_BONF   <- ALPHA / M_SNP
RETAIN_P <- 1e-5                 # what gets written to disk
MODELS   <- c("MLM", "MLMM", "BLINK", "FarmCPU")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: Rscript 30b_split_trait_FDR.R <trait> [in_dir] [out_dir]")
TRAIT <- args[1]
if (length(args) >= 2L) INPUT_DIR  <- args[2]
if (length(args) >= 3L) OUTPUT_DIR <- args[3]

DIR_MODEL <- file.path(OUTPUT_DIR, "by_model")
DIR_CONS  <- file.path(OUTPUT_DIR, "consensus")
DIR_SUMM  <- file.path(OUTPUT_DIR, "summary")
for (d in c(DIR_MODEL, DIR_CONS, DIR_SUMM)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

normalize_model <- function(x) {
  key <- toupper(trimws(as.character(x)))
  unname(c(MLM = "MLM", MLMM = "MLMM", BLINK = "BLINK", FARMCPU = "FarmCPU")[key])
}
save_named <- function(obj, nm, path) {
  env <- new.env(parent = emptyenv()); assign(nm, obj, envir = env)
  save(list = nm, file = path, envir = env, compress = "gzip")
}
t0  <- Sys.time()
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

msg("trait ", TRAIT)
msg(sprintf("  FDR q<=%.2f (GAPIT column) | retain P<=%.0e | m=%d | Bonf %.3e | sugg %.3e",
            Q_FDR, RETAIN_P, M_SNP, P_BONF, P_SUGG))

in_file <- file.path(INPUT_DIR, paste0(TRAIT, ".csv"))
if (!file.exists(in_file)) stop("no such file: ", in_file)

## ---- which columns exist -----------------------------------------------------
hdr  <- names(fread(in_file, nrows = 0))
want <- c("SNP", "Chr", "Pos", "P.value", "MAF", "Effect", "model")
miss <- setdiff(want, hdr)
if (length(miss)) stop(TRAIT, ": missing column(s): ", paste(miss, collapse = ", "))
has_fdr <- "FDR" %in% hdr
if (!has_fdr) {
  warning(TRAIT, ": no FDR column; falling back to recomputed BH. ",
          "This is only valid if the file was not p-filtered above the BH cutoff.")
}
sel <- c(want, if (has_fdr) "FDR")

msg("  reading ", in_file)
dt <- fread(in_file, select = sel, showProgress = FALSE)
n_raw <- nrow(dt)

dt[, P.value := as.numeric(P.value)]
if (has_fdr) dt[, FDR := as.numeric(FDR)]
dt <- dt[is.finite(P.value) & P.value > 0 & P.value <= RETAIN_P]

dt[, model := normalize_model(model)]
n_bad <- dt[is.na(model), .N]
if (n_bad > 0L) warning(sprintf("%s: %d rows with unrecognised model label, dropped", TRAIT, n_bad))
dt <- dt[!is.na(model)]
dt[, Chr := as.integer(Chr)][, Pos := as.numeric(Pos)]
dt <- dt[!is.na(Chr) & !is.na(Pos)]

## ---- recomputed BH, per trait x model, as a cross-check ----------------------
## Valid because BH compares the ORDERED SMALL p-values against i*q/m; the large
## p-values never enter the comparison, only their count (m) does. It becomes
## invalid only if the cutoff would land above RETAIN_P, which is flagged below.
dt[, bh_q := {
  o <- order(P.value); r <- integer(.N); r[o] <- seq_len(.N)
  pmin(1, rev(cummin(rev((P.value * M_SNP / r)[o]))[r]))
}, by = model]

if (has_fdr) {
  dt[, fdr_sig := FDR <= Q_FDR]
} else {
  dt[, FDR := bh_q][, fdr_sig := bh_q <= Q_FDR]
}
dt[, bonferroni := P.value <= P_BONF]
dt[, suggestive := P.value <= P_SUGG]

## truncation guard: did any model's selected set run to the edge of RETAIN_P?
trunc <- dt[fdr_sig == TRUE, .(maxP = max(P.value)), by = model][maxP > 0.5 * RETAIN_P]
if (nrow(trunc)) {
  warning(sprintf("%s: FDR-selected set reaches %.2e, close to RETAIN_P=%.0e -- raise RETAIN_P and re-run",
                  TRAIT, max(trunc$maxP), RETAIN_P))
}
disagree <- dt[has_fdr & (FDR <= Q_FDR) != (bh_q <= Q_FDR), .N]
if (has_fdr && disagree > 0L)
  msg("  NOTE: GAPIT FDR and recomputed BH disagree on ", disagree, " rows (check m)")

models_run <- intersect(MODELS, normalize_model(
  unique(fread(in_file, select = "model", showProgress = FALSE)$model)))

msg("  raw rows=", n_raw, "  retained=", nrow(dt),
    "  FDR-sig=", dt[, sum(fdr_sig)], "  Bonf=", dt[, sum(bonferroni)],
    "  models: ", paste(models_run, collapse = "+"))

## ---- per-model .RData --------------------------------------------------------
for (mdl in models_run) {
  sub <- dt[model == mdl]; setorder(sub, Chr, Pos)
  nm <- paste0(TRAIT, "_", mdl); path <- file.path(DIR_MODEL, paste0(nm, ".RData"))
  save_named(sub[], nm, path)
  msg("  wrote ", basename(path), "  (", nrow(sub), " retained, ",
      sub[, sum(fdr_sig)], " FDR-sig, ", sub[, sum(bonferroni)], " Bonferroni)")
}

## ---- cross-model consensus, on the FDR-significant set -----------------------
sig <- dt[fdr_sig == TRUE]
wide <- NULL
if (nrow(sig)) {
  wide <- dcast(sig, SNP + Chr + Pos ~ model, value.var = "P.value",
                fun.aggregate = min, fill = NA_real_)
  for (mdl in MODELS) if (!mdl %in% names(wide)) wide[, (mdl) := NA_real_]
  setnames(wide, MODELS, paste0("P_", MODELS))
  pcols <- paste0("P_", models_run)
  wide[, n_models_sig := rowSums(!is.na(.SD)), .SDcols = pcols]
  wide[, min_P := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = pcols]
  wide[, max_P := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = pcols]
  wide[, n_models_run := length(models_run)]
  wide[, n_models_bonf := rowSums(.SD <= P_BONF, na.rm = TRUE), .SDcols = pcols]
  wide[, all_models_fdr := n_models_sig == length(models_run)]
  maf  <- unique(sig[, .(SNP, MAF)], by = "SNP")
  wide <- maf[wide, on = "SNP"]
}
cons <- if (is.null(wide)) data.table() else wide[n_models_sig == length(models_run)]
if (nrow(cons)) setorder(cons, max_P)
nm_cons <- paste0(TRAIT, "_all_model_common_SNP")
save_named(cons, nm_cons, file.path(DIR_CONS, paste0(nm_cons, ".RData")))
msg("  wrote ", nm_cons, ".RData  (", nrow(cons), " common at FDR<=", Q_FDR, ")")

## ---- summary row -------------------------------------------------------------
at <- function(w, k) if (is.null(wide)) 0L else sum(wide[[w]] >= k)
summ <- data.table(
  Phenotype      = TRAIT,
  Models_run     = paste(models_run, collapse = "+"),
  N_models_run   = length(models_run),
  N_rows_raw     = n_raw,
  N_rows_kept    = nrow(dt),
  M_snps_assumed = M_SNP,
  Q_FDR          = Q_FDR,
  FDR_source     = if (has_fdr) "GAPIT" else "recomputed_BH",
  FDR_cutoff_P   = if (nrow(sig)) max(sig$P.value) else NA_real_,
  P_bonferroni   = P_BONF,
  P_suggestive   = P_SUGG,
  N_fdr_union    = at("n_models_sig", 1),
  N_fdr_2plus    = at("n_models_sig", 2),
  N_fdr_3plus    = at("n_models_sig", 3),
  N_fdr_all      = nrow(cons),
  N_bonf_union   = uniqueN(dt[bonferroni == TRUE]$SNP),
  N_sugg_union   = uniqueN(dt[suggestive == TRUE]$SNP),
  Runtime_min    = as.numeric(difftime(Sys.time(), t0, units = "mins")))
for (mdl in MODELS) {
  summ[, (paste0("N_fdr_",  mdl)) := dt[model == mdl, sum(fdr_sig)]]
  summ[, (paste0("N_bonf_", mdl)) := dt[model == mdl, sum(bonferroni)]]
}
fwrite(summ, file.path(DIR_SUMM, paste0(TRAIT, "_summary.csv")))
msg("  done in ", sprintf("%.2f", summ$Runtime_min), " min  |  FDR cutoff P = ",
    signif(summ$FDR_cutoff_P, 4))

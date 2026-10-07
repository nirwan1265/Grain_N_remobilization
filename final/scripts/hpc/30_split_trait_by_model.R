#!/usr/bin/env Rscript
################################################################################
## 30 -- SPLIT ONE TRAIT'S GWAS SCAN INTO PER-MODEL .RData + CONSENSUS SNP SET
##
## Input :  <INPUT_DIR>/<trait>.csv
##          columns: SNP, Chr, Pos, P.value, MAF, nobs, H.B.P.Value, Effect,
##                   trait, model, FDR, p_star, line_log10
##          all four GAPIT models stacked in one file. These files are ALREADY
##          p-filtered upstream and their row counts differ between traits
##          (A.csv ~12.5M rows, A.LAV.csv ~8.4M), so the row count tells you
##          nothing about how many markers were tested. m is set below.
##
## THRESHOLDS -- both derived from M_SNP, which you set, not from the file:
##   suggestive  1/m      = 2.43e-07   (-log10 6.61)  one expected false positive
##   Bonferroni  0.05/m   = 1.21e-08   (-log10 7.92)
##
## Rows above the suggestive threshold are discarded immediately after reading.
## Nothing weaker than 1/m is ever written to disk, which keeps each .RData in
## the kilobyte-to-megabyte range instead of hundreds of megabytes.
##
## Output:  <OUTPUT_DIR>/by_model/<trait>_<MODEL>.RData
##             object <trait>_<MODEL>: that model's hits at P <= 1/m,
##             with a logical column `bonferroni` marking P <= 0.05/m
##          <OUTPUT_DIR>/consensus/<trait>_all_model_common_SNP.RData
##             object <trait>_all_model_common_SNP: WIDE data.table of SNPs
##             suggestive-significant in EVERY model that ran --
##             SNP, Chr, Pos, MAF, P_MLM, P_MLMM, P_BLINK, P_FarmCPU,
##             min_P, max_P, n_models_sugg, n_models_bonf, n_models_run,
##             all_models_bonferroni
##          <OUTPUT_DIR>/summary/<trait>_summary.csv   (one row; script 31 pools)
##
## Usage
##   Rscript 30_split_trait_by_model.R <trait>                    # A, N.E, Total_N
##   Rscript 30_split_trait_by_model.R <trait> <in_dir> <out_dir>
################################################################################

suppressPackageStartupMessages(library(data.table))

## ---- configuration ----------------------------------------------------------

INPUT_DIR  <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
OUTPUT_DIR <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/by_model"

## Number of markers tested in the Goodman-Buckler amino-acid GWAS.
## Taken from the manuscript Methods. CHANGE THIS if the true marker count
## differs -- every threshold below follows from it.
M_SNP <- 4117796

ALPHA    <- 0.05
P_SUGG   <- 1 / M_SNP            # suggestive, one expected false positive
P_BONF   <- ALPHA / M_SNP        # Bonferroni
MODELS   <- c("MLM", "MLMM", "BLINK", "FarmCPU")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("usage: Rscript 30_split_trait_by_model.R <trait> [in_dir] [out_dir]")
}
TRAIT <- args[1]
if (length(args) >= 2L) INPUT_DIR  <- args[2]
if (length(args) >= 3L) OUTPUT_DIR <- args[3]

DIR_MODEL <- file.path(OUTPUT_DIR, "by_model")
DIR_CONS  <- file.path(OUTPUT_DIR, "consensus")
DIR_SUMM  <- file.path(OUTPUT_DIR, "summary")
for (d in c(DIR_MODEL, DIR_CONS, DIR_SUMM)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- helpers ----------------------------------------------------------------

normalize_model <- function(x) {
  key <- toupper(trimws(as.character(x)))
  unname(c(MLM = "MLM", MLMM = "MLMM", BLINK = "BLINK", FARMCPU = "FarmCPU")[key])
}

## save `obj` into `path` under the name `nm`, so load() restores <nm>
save_named <- function(obj, nm, path) {
  env <- new.env(parent = emptyenv())
  assign(nm, obj, envir = env)
  save(list = nm, file = path, envir = env, compress = "gzip")
}

t0  <- Sys.time()
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

msg("trait ", TRAIT)
msg(sprintf("  m = %d   suggestive 1/m = %.4e (-log10 %.4f)   Bonferroni %g/m = %.4e (-log10 %.4f)",
            M_SNP, P_SUGG, -log10(P_SUGG), ALPHA, P_BONF, -log10(P_BONF)))

## ---- read, then immediately drop everything weaker than suggestive ----------

in_file <- file.path(INPUT_DIR, paste0(TRAIT, ".csv"))
if (!file.exists(in_file)) stop("no such file: ", in_file)

msg("  reading ", in_file)
dt <- fread(in_file,
            select = c("SNP", "Chr", "Pos", "P.value", "MAF", "Effect", "model"),
            showProgress = FALSE)

n_raw <- nrow(dt)
dt[, P.value := as.numeric(P.value)]
dt <- dt[is.finite(P.value) & P.value > 0 & P.value <= P_SUGG]

dt[, model := normalize_model(model)]
dt[, Chr := as.integer(Chr)]
dt[, Pos := as.numeric(Pos)]

n_bad_model <- dt[is.na(model), .N]
if (n_bad_model > 0L) {
  warning(sprintf("%s: %d rows with an unrecognised model label, dropped",
                  TRAIT, n_bad_model))
}
dt <- dt[!is.na(model) & !is.na(Chr) & !is.na(Pos)]
dt[, bonferroni := P.value <= P_BONF]

## Which models actually appear at all in this file. A model with no hit at the
## suggestive level would silently vanish here, which would make a "common to
## every model" set meaningless -- so read the model labels from the raw file.
models_present <- intersect(
  MODELS,
  normalize_model(unique(fread(in_file, select = "model", showProgress = FALSE)$model)))
models_run <- models_present

msg("  raw rows=", n_raw, "  at P<=1/m: ", nrow(dt),
    "  models in file: ", paste(models_run, collapse = "+"))

## ---- per-model .RData -------------------------------------------------------

for (mdl in models_run) {
  sub <- dt[model == mdl]
  setorder(sub, Chr, Pos)
  nm   <- paste0(TRAIT, "_", mdl)
  path <- file.path(DIR_MODEL, paste0(nm, ".RData"))
  save_named(sub[], nm, path)
  msg("  wrote ", basename(path), "  (", nrow(sub), " suggestive, ",
      sub[, sum(bonferroni)], " Bonferroni, object `", nm, "`)")
}

## ---- consensus across every model that ran ----------------------------------

wide <- NULL
if (nrow(dt) > 0L) {
  wide <- dcast(dt, SNP + Chr + Pos ~ model, value.var = "P.value",
                fun.aggregate = min, fill = NA_real_)
  for (mdl in MODELS) if (!mdl %in% names(wide)) wide[, (mdl) := NA_real_]
  setnames(wide, MODELS, paste0("P_", MODELS))
  pcols <- paste0("P_", models_run)

  wide[, n_models_sugg := rowSums(!is.na(.SD)),                .SDcols = pcols]
  wide[, n_models_bonf := rowSums(.SD <= P_BONF, na.rm = TRUE), .SDcols = pcols]
  wide[, min_P := do.call(pmin, c(.SD, na.rm = TRUE)),          .SDcols = pcols]
  wide[, max_P := do.call(pmax, c(.SD, na.rm = TRUE)),          .SDcols = pcols]
  wide[, n_models_run := length(models_run)]
  wide[, all_models_bonferroni := n_models_bonf == length(models_run)]

  maf  <- unique(dt[, .(SNP, MAF)], by = "SNP")
  wide <- maf[wide, on = "SNP"]
}

cons <- if (is.null(wide)) data.table() else wide[n_models_sugg == length(models_run)]
if (nrow(cons)) {
  setcolorder(cons, c("SNP", "Chr", "Pos", "MAF", paste0("P_", MODELS),
                      "min_P", "max_P", "n_models_sugg", "n_models_bonf",
                      "n_models_run", "all_models_bonferroni"))
  setorder(cons, max_P)
}

nm_cons   <- paste0(TRAIT, "_all_model_common_SNP")
path_cons <- file.path(DIR_CONS, paste0(nm_cons, ".RData"))
save_named(cons, nm_cons, path_cons)
msg("  wrote ", basename(path_cons), "  (", nrow(cons),
    " common at 1/m, ", if (nrow(cons)) sum(cons$all_models_bonferroni) else 0L,
    " common at 0.05/m, object `", nm_cons, "`)")

## ---- per-trait summary row --------------------------------------------------

n_sugg_at <- function(k) if (is.null(wide)) 0L else sum(wide$n_models_sugg >= k)
n_bonf_at <- function(k) if (is.null(wide)) 0L else sum(wide$n_models_bonf >= k)

summ <- data.table(
  Phenotype             = TRAIT,
  Models_run            = paste(models_run, collapse = "+"),
  N_models_run          = length(models_run),
  N_rows_raw            = n_raw,
  M_snps_assumed        = M_SNP,
  P_suggestive          = P_SUGG,
  P_bonferroni          = P_BONF,
  P_star_in_file        = NA_real_,
  ## suggestive, 1/m
  N_sugg_union          = n_sugg_at(1),
  N_sugg_2plus          = n_sugg_at(2),
  N_sugg_3plus          = n_sugg_at(3),
  N_sugg_all_models     = nrow(cons),
  ## Bonferroni, 0.05/m
  N_bonf_union          = n_bonf_at(1),
  N_bonf_2plus          = n_bonf_at(2),
  N_bonf_3plus          = n_bonf_at(3),
  N_bonf_all_models     = if (nrow(cons)) sum(cons$all_models_bonferroni) else 0L,
  Runtime_min           = as.numeric(difftime(Sys.time(), t0, units = "mins")))

for (mdl in MODELS) {
  summ[, (paste0("N_sugg_", mdl)) := dt[model == mdl, .N]]
  summ[, (paste0("N_bonf_", mdl)) := dt[model == mdl, sum(bonferroni)]]
}

fwrite(summ, file.path(DIR_SUMM, paste0(TRAIT, "_summary.csv")))
msg("  done in ", sprintf("%.2f", summ$Runtime_min), " min")

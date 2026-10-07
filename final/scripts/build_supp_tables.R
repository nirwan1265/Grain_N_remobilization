#!/usr/bin/env Rscript
################################################################################
## BUILD THE SUPPLEMENTARY TABLES FOR final/, NUMBERED BY ORDER OF FIRST CITATION
##
##   S1  free amino acid phenotype values, 280 taxa x 113 traits
##   S2  phenotype abbreviations and definitions
##   S3  heritability per trait: broad-sense H2 and genomic h2
##   S4  free vs protein-bound amino acid pools, with paired tests
##   S5  Bonferroni candidate genes, one row per gene x phenotype x SNP   (amino_gwas_main_table.R)
##   S6  loci, ranked by phenotype count, with genes and GO              (amino_gwas_main_table.R)
##
## S5 and S6 are written by amino_gwas_main_table.R, not here. This script
## assembles S1-S4 from inputs that already exist elsewhere in the project, so
## that final/ carries one self-contained, consistently numbered set.
##
## S4 is NOT a copy of an existing table. The paired-test table gives means and
## p-values; the per-source summary gives the spread. The paragraph that cites
## it turns on neither of those but on the CONTRAST IN VARIABILITY between the
## two pools, so this builds one table carrying mean, median, sd and CV for each
## pool side by side with the paired test.
##
## Usage: Rscript build_supp_tables.R [repo_root]
################################################################################

suppressPackageStartupMessages(library(data.table))

ROOT  <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
STAGE <- file.path(ROOT, "final", "results", "_stage")
OUT   <- file.path(ROOT, "final", "supp_tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) fread(f, encoding = "UTF-8", showProgress = FALSE)

## ---- S1: phenotype values ----------------------------------------------------
s1_src <- file.path(OUT, "SuppTable_1_FAA_phenotypes.csv")
if (file.exists(s1_src)) {
  file.rename(s1_src, file.path(OUT, "SuppTable_S1_FAA_phenotypes.csv"))
}
cat("S1 phenotype values      : ",
    if (file.exists(file.path(OUT, "SuppTable_S1_FAA_phenotypes.csv"))) "present" else "MISSING", "\n")

## ---- S2: abbreviations -------------------------------------------------------
ab <- rd(file.path(STAGE, "abbre.csv"))
setnames(ab, 1, "metric")
fwrite(ab, file.path(OUT, "SuppTable_S2_phenotype_definitions.csv"))
cat("S2 phenotype definitions : ", nrow(ab), "rows\n")

## ---- S3: heritability --------------------------------------------------------
h <- rd(file.path(STAGE, "H2_vs_h2_combined.csv"))
s3 <- h[, .(trait, pool = source, n_obs, n_taxa, n_year, reps_per_env,
            H2_line_mean = H2_log, H2_Cullis = H2_log_Cullis,
            n_genotypes_h2 = n, h2 = h2_log, h2_lo, h2_hi,
            h2_permuted_null = perm_null, h2_reliable = reliable,
            ratio_h2_over_H2 = ratio)]
setorder(s3, pool, -H2_line_mean)
fwrite(s3, file.path(OUT, "SuppTable_S3_heritability.csv"))
cat("S3 heritability          : ", nrow(s3), "traits (",
    s3[pool == "FAA", .N], "FAA,", s3[pool == "PBAA", .N], "PBAA )\n")

## ---- S4: FAA vs PBAA pools ---------------------------------------------------
summ  <- rd(file.path(STAGE, "t13.csv"))
pair  <- rd(file.path(STAGE, "t14.csv"))
summ[, CV := 100 * sd / mean]
spread <- summ[, .(metric, sd, CV)]          # one row per metric, either pool

s4 <- merge(pair, spread, by.x = "faa_metric",  by.y = "metric", all.x = TRUE)
setnames(s4, c("sd", "CV"), c("FAA_sd", "FAA_CV_pct"))
s4 <- merge(s4,   spread, by.x = "pbaa_metric", by.y = "metric", all.x = TRUE)
setnames(s4, c("sd", "CV"), c("PBAA_sd", "PBAA_CV_pct"))

s4 <- s4[, .(comparison = comparison_label, n_taxa,
             FAA_metric  = faa_metric,  FAA_mean  = mean_faa,  FAA_median  = median_faa,
             FAA_sd,  FAA_CV_pct,
             PBAA_metric = pbaa_metric, PBAA_mean = mean_pbaa, PBAA_median = median_pbaa,
             PBAA_sd, PBAA_CV_pct,
             paired_wilcoxon_p = p_value)]
setorder(s4, -PBAA_mean)
fwrite(s4, file.path(OUT, "SuppTable_S4_FAA_vs_PBAA_pools.csv"))

cat("S4 FAA vs PBAA           : ", nrow(s4), "comparisons\n")
cat("\n--- numbers for the FAA vs PBAA paragraph ---\n")
for (i in seq_len(nrow(s4))) {
  x <- s4[i]
  cat(sprintf("  %-52s FAA %10.3f (CV %5.1f%%) | PBAA %10.3f (CV %5.1f%%) | p = %.3g\n",
              substr(x$comparison, 1, 52), x$FAA_mean, x$FAA_CV_pct,
              x$PBAA_mean, x$PBAA_CV_pct, x$paired_wilcoxon_p))
}
cat(sprintf("\n  smallest paired p-value across the four comparisons: %.3g\n",
            min(s4$paired_wilcoxon_p)))
cat(sprintf("  largest  paired p-value across the four comparisons: %.3g\n",
            max(s4$paired_wilcoxon_p)))

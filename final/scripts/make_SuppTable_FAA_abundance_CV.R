#!/usr/bin/env Rscript
################################################################################
## Supplementary table: abundance, variability and repeatability of the 20 FAAs
## in the Goodman-Buckler panel.
##
## Summarised over GENOTYPE MEANS (each line averaged across its four
## replicates), so the CV is variation among lines, not among plots.
## H2 is joined from the heritability table.
################################################################################
suppressPackageStartupMessages({library(data.table)})
gv <- function(k,d){v <- Sys.getenv(k); if (nzchar(v)) v else d}
ROOT <- gv("GN_ROOT", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final")
FAA  <- gv("GN_FAA",  file.path(ROOT, "results/_stage/FAA_Goodman_Buckler.csv"))
H2F  <- gv("GN_H2",   file.path(ROOT, "supp_tables/SuppTable_S4_heritability.csv"))
OUT  <- gv("GN_OUT",  file.path(ROOT, "supp_tables/SuppTable_S3_FAA_abundance_CV.csv"))

NAME <- c(A="Alanine", R="Arginine", N="Asparagine", D="Aspartate", Q="Glutamine",
          E="Glutamate", G="Glycine", H="Histidine", I="Isoleucine", L="Leucine",
          K="Lysine", M="Methionine", F="Phenylalanine", P="Proline", S="Serine",
          W="Tryptophan", T="Threonine", Y="Tyrosine", V="Valine", C="Cysteine")

d <- fread(FAA, encoding = "UTF-8"); setnames(d, 1, "taxa")
long <- melt(d[, c("taxa", names(NAME)), with = FALSE], id.vars = "taxa",
             variable.name = "AA", value.name = "conc", variable.factor = FALSE)
long <- long[is.finite(conc)]
gm <- long[, .(gmean = mean(conc)), by = .(AA, taxa)]          # genotype means

res <- gm[, .(n_genotypes = .N,
              median = median(gmean), mean = mean(gmean), sd = sd(gmean),
              min = min(gmean), max = max(gmean)), by = AA]
res[, fold := max / min]
res[, CV_pct := 100 * sd / mean]
res[, Amino_acid := NAME[AA]]

h2 <- tryCatch(fread(H2F, encoding = "UTF-8"), error = function(e) NULL)
if (!is.null(h2) && all(c("trait","pool","H2_line_mean") %in% names(h2))) {
  res <- merge(res, h2[pool == "FAA", .(AA = trait, H2_line_mean, h2_genomic = h2)],
               by = "AA", all.x = TRUE)
} else warning("heritability table not found at ", H2F, " -- H2 columns omitted")

setorder(res, -CV_pct)
setcolorder(res, c("AA", "Amino_acid", "n_genotypes", "median", "mean", "sd",
                   "min", "max", "fold", "CV_pct"))
fwrite(res, OUT)
cat("written:", OUT, "\n"); print(res[, .(Amino_acid, median = round(median,3),
    fold = round(fold), CV_pct = round(CV_pct))])

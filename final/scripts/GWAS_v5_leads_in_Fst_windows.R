## Goodman-Buckler amino acid GWAS leads (lifted AGPv4 -> v5) inside Indian Chief / Jarvis
## top-5% F_ST windows, vs 10,000 circular shifts of the window track. Leads within
## 25 kb are merged into one locus first, so a locus counts once.
## The null shifts all loci together by the same random offset around the genome
## (concatenated chromosomes), which keeps the clustering of both the loci and the
## windows. Tests: each population separately, and either population ("IC or J").
## Run from the repo root (or set GN_REPO). Outputs, in final/results/ic_j_fst/:
##   GWAS_v5_loci_in_Fst_windows.csv          loci inside a top-5% window
##   GWAS_v5_loci_in_Fst_windows_test.csv     observed vs circular-shift expectation
suppressPackageStartupMessages(library(data.table)); set.seed(1)
gv  <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", ".")
LDD  <- gv("GN_LDDIR", file.path(REPO, "final/results/goodman_buckler_amino_acid_gwas/LD"))
TAB  <- gv("GN_TABDIR", file.path(REPO, "final/supp_tables"))
OUT  <- gv("GN_OUTDIR", file.path(REPO, "final/results/ic_j_fst"))
lv <- fread(file.path(LDD, "amino_leads_v5.tsv"))[!is.na(pos_v5)]
tr <- fread(file.path(LDD, "amino_lead_traits.tsv"))
lv <- merge(lv, tr, by.x = "lead", by.y = "SNP")
lv[, asp := vapply(strsplit(traits, ","), function(v) any(grepl("^(D|N|T|M|K|I)(\\.|$)", v) | v == "IMNTDK"), logical(1))]
setorder(lv, chr, pos_v5)
lv[, locus := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(pos_v5) > 25000))]
LOC <- lv[, .(chr = chr[1], pos = pos_v5[which.min(minP)], minP = min(minP), asp = any(asp),
              traits = paste(sort(unique(unlist(strsplit(traits, ",")))), collapse = ","),
              leads = paste(lead, collapse = ";")), by = locus]
cat("leads", nrow(lv), "-> loci", nrow(LOC), " Asp-family loci", LOC[, sum(asp)], "\n")
run <- function(pop, f) {
  W <- fread(f)[, chr_num := as.integer(sub("chr", "", chr))]
  len <- W[, .(len = max(end_bp)), by = chr_num][order(chr_num)][, off := cumsum(c(0, head(len, -1)))]
  G <- sum(len$len); W <- merge(W, len[, .(chr_num, off)], by = "chr_num")
  top <- W[top5 == TRUE, .(s = start_bp + off, e = end_bp + off, chr, start_bp, end_bp, mean_fst)][order(s)]
  P <- merge(LOC, len[, .(chr = chr_num, off)], by = "chr")[, g := pos + off]
  inside <- function(g) { i <- findInterval(g, top$s); i > 0 & g <= top$e[pmax(i, 1)] }
  one <- function(sel, lab) { g <- P[sel, g]; o <- sum(inside(g))
    nul <- vapply(runif(10000, 0, G), function(k) sum(inside((g - k) %% G)), numeric(1))
    data.table(population = pop, loci = lab, n = length(g), observed = o, expected = round(mean(nul), 1),
               P = signif((sum(nul >= o) + 1) / 10001, 3)) }
  hits <- P[inside(g)][, i := findInterval(g, top$s)][, `:=`(window = sprintf("%s:%.2f-%.2f Mb", top$chr[i], top$start_bp[i]/1e6, top$end_bp[i]/1e6), window_fst = round(top$mean_fst[i], 3))]
  list(rbind(one(rep(TRUE, nrow(P)), "all"), one(P$asp, "Asp-family")), hits[, population := pop])
}
r <- list(run("Indian Chief", file.path(TAB, "SuppTable_S18_IC_Fst_windows.csv")),
          run("Jarvis",       file.path(TAB, "SuppTable_S19_J_Fst_windows.csv")))
## either population: a locus counts if it is inside an IC or a J top-5% window
both <- function() {
  rd <- function(f) { W <- fread(f)[, chr_num := as.integer(sub("chr", "", chr))]; W }
  WI <- rd(file.path(TAB, "SuppTable_S18_IC_Fst_windows.csv")); WJ <- rd(file.path(TAB, "SuppTable_S19_J_Fst_windows.csv"))
  len <- rbind(WI, WJ)[, .(len = max(end_bp)), by = chr_num][order(chr_num)][, off := cumsum(c(0, head(len, -1)))]
  G <- sum(len$len)
  tops <- lapply(list(WI, WJ), function(W) merge(W[top5 == TRUE], len[, .(chr_num, off)], by = "chr_num")[, .(s = start_bp + off, e = end_bp + off)][order(s)])
  inside1 <- function(g, top) { i <- findInterval(g, top$s); i > 0 & g <= top$e[pmax(i, 1)] }
  inside <- function(g) inside1(g, tops[[1]]) | inside1(g, tops[[2]])
  ## genome fraction covered by the union of the two tracks
  U <- rbind(tops[[1]], tops[[2]])[order(s)]; cov <- 0; cs <- U$s[1]; ce <- U$e[1]
  for (k in seq_len(nrow(U))[-1]) { if (U$s[k] <= ce) ce <- max(ce, U$e[k]) else { cov <- cov + ce - cs + 1; cs <- U$s[k]; ce <- U$e[k] } }
  cov <- cov + ce - cs + 1
  P <- merge(LOC, len[, .(chr = chr_num, off)], by = "chr")[, g := pos + off]
  one <- function(sel, lab) { g <- P[sel, g]; o <- sum(inside(g))
    nul <- vapply(runif(10000, 0, G), function(k) sum(inside((g - k) %% G)), numeric(1))
    data.table(population = "Indian Chief or Jarvis", loci = lab, n = length(g), observed = o,
               expected = round(mean(nul), 1), P = signif((sum(nul >= o) + 1) / 10001, 3),
               genome_fraction = round(cov / G, 4)) }
  rbind(one(rep(TRUE, nrow(P)), "all"), one(P$asp, "Asp-family"))
}
TEST <- rbindlist(c(lapply(r, `[[`, 1), list(both())), fill = TRUE)
print(TEST)
H <- rbindlist(lapply(r, `[[`, 2))
dir.create(OUT, FALSE, TRUE)
fwrite(H, file.path(OUT, "GWAS_v5_loci_in_Fst_windows.csv"))
fwrite(TEST, file.path(OUT, "GWAS_v5_loci_in_Fst_windows_test.csv"))
print(H[asp == TRUE][order(minP), .(population, window, window_fst, pos, minP = signif(minP, 3), traits = substr(traits, 1, 50))][1:15])

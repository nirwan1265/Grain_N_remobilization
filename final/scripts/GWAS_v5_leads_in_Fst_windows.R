## Goodman-Buckler amino acid GWAS leads (lifted AGPv4 -> v5) inside Indian Chief / Jarvis
## top-5% F_ST windows, vs 10,000 circular shifts of the window track. Leads within
## 25 kb are merged into one locus first, so a locus counts once.
suppressPackageStartupMessages(library(data.table)); set.seed(1)
lv <- fread("amino_leads_v5.tsv")[!is.na(pos_v5)]
tr <- fread("amino_lead_traits.tsv")
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
r <- list(run("Indian Chief", "SuppTable_S17_IC_Fst_windows.csv"), run("Jarvis", "SuppTable_S18_J_Fst_windows.csv"))
print(rbindlist(lapply(r, `[[`, 1)))
H <- rbindlist(lapply(r, `[[`, 2))
fwrite(H, "GWAS_v5_loci_in_Fst_windows.csv")
print(H[asp == TRUE][order(minP), .(population, window, window_fst, pos, minP = signif(minP, 3), traits = substr(traits, 1, 50))][1:15])

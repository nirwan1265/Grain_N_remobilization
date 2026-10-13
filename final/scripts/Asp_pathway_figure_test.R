## Inputs (local copies of final/results files): GO_universe_genomewide.txt, all_candidate_gene_sources.csv
## Asp / Asn pathway (the user's pathway figure) and Glu/Gln (GS/GOGAT/GDH) gene sets:
## how many fall in each candidate set vs chance (hypergeometric, genome background)
suppressPackageStartupMessages(library(data.table))
U <- unique(fread("/mnt/user-data/uploads/Grain_N_remobilization/final/results/GO/GO_universe_genomewide.txt", quote = "\""), by = "GeneID")
U[, nm := paste(sub(";.*", "", GeneName), "|", Family_Subfamily)]
core <- c(AspAT = "aspartate aminotransferase|glutamate.oxaloacetate transaminase",
          AK = "aspartate kinase|aspartokinase", ASD = "aspartate-semialdehyde dehydrogenase",
          HSD = "homoserine dehydrogenase", HSK = "homoserine kinase", TS = "threonine synthase",
          TD = "threonine dehydratase|threonine deaminase",
          Met = "cystathionine gamma-synthase|cystathionine beta-lyase|(^|[^l])methionine synthase|5-methyltetrahydropteroyltriglutamate--homocysteine",
          Lys = "dihydrodipicolinate|tetrahydrodipicolinate|diaminopimelate",
          AS = "asparagine synthetase", ASNase = "asparaginase", wAmid = "omega-amidase")
ext  <- c(AO = "L-aspartate oxidase", QS = "quinolinate synthase", ASS = "argininosuccinate synthase",
          ATC = "aspartate carbamoyltransferase", PAT = "prephenate aminotransferase|aspartate.prephenate",
          GH3 = "indole-3-acetic acid-amido synthetase")
glu  <- c(GS = "glutamine synthetase|GS catalytic", GOGAT = "(^|[^l])glutamate synthase", GDH = "glutamate dehydrogenase")
tag <- function(p) { out <- rep(NA_character_, nrow(U)); for (k in names(p)) out[is.na(out) & grepl(p[[k]], U$nm, ignore.case = TRUE)] <- k; out }
U[, core := tag(core)][, ext := tag(ext)][, glu := tag(glu)]
U[!is.na(core), ext := NA]
cat("genome: core", U[, sum(!is.na(core))], " extended", U[, sum(!is.na(core) | !is.na(ext))], " Glu/Gln", U[, sum(!is.na(glu))], "\n")
print(U[!is.na(core), .N, by = core]); print(U[!is.na(ext), .N, by = ext]); print(U[!is.na(glu), .N, by = glu])
S <- fread("all_candidate_gene_sources.csv")
sets <- list(`Amino acid GWAS (25 kb)` = S[aminoGWAS_25kb == 1, GeneID], `Soil-N GWAS (25 kb)` = S[soilN_25kb == 1, GeneID],
             `Indian Chief F_ST` = S[IC_Fst == 1, GeneID], `Jarvis F_ST` = S[J_Fst == 1, GeneID])
sets[["All four"]] <- unique(unlist(sets))
N <- nrow(U)
test <- function(g, flag) { g <- intersect(g, U$GeneID); K <- sum(flag); k <- sum(flag[U$GeneID %in% g])
  c(n = length(g), observed = k, expected = round(length(g) * K / N, 2), P = signif(phyper(k - 1, K, N - K, length(g), lower.tail = FALSE), 2)) }
R <- rbindlist(lapply(names(sets), function(s) rbind(
  data.table(set = s, pathway = "Asp/Asn core", t(test(sets[[s]], !is.na(U$core)))),
  data.table(set = s, pathway = "Asp/Asn + Asp-consuming (figure)", t(test(sets[[s]], !is.na(U$core) | !is.na(U$ext)))),
  data.table(set = s, pathway = "Glu/Gln (GS, GOGAT, GDH)", t(test(sets[[s]], !is.na(U$glu)))))))
print(R)
H <- rbindlist(lapply(names(sets)[1:4], function(s) U[GeneID %in% sets[[s]] & (!is.na(core) | !is.na(ext) | !is.na(glu)),
  .(set = s, GeneID, step = fcoalesce(core, ext, glu), name = substr(sub(";.*", "", GeneName), 1, 45))]))
print(H[order(set, step)])
fwrite(R, "Asp_pathway_figure_test.csv"); fwrite(H, "Asp_pathway_figure_genes.csv")

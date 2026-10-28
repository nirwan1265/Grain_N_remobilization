################################################################################
## FIGURE -- tissue expression of the aspartate pathway candidates and three
## transporter candidates (RNA-seq survey of seven tissues, flowering to two
## weeks after black layer; unpublished data from Alexa's lab, used with permission)
##
## A  Aspartate pathway genes      B  Transporter genes
## Fill = mean log(TPM + 1) per tissue (TMM-adjusted), same viridis "D" scale as
## the soil N map in Fig 2. The left strip marks the analysis that found each gene,
## in the population colours used elsewhere in the paper:
##   GB amino acid GWAS #C51B7D | soil N GWAS #1B7837 (Fig 2, 0-5 cm)
##   Indian Chief #2166AC | Jarvis #B35806 (Fig 4 and 5, Cycle 14)
##   B64 x CML14 Asp/Pro QTL interval #252525
##
## Usage: Rscript Fig_candidate_expression_tissues.R <pathway_csv> <transporter_csv> <out_dir>
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(viridis)
})
args <- commandArgs(trailingOnly = TRUE)
P_PATH <- if (length(args) >= 1) args[1] else "final/results/expression/Aspartate_Homoserine_Pathway_Summary_Table.csv"
P_TRAN <- if (length(args) >= 2) args[2] else "final/results/expression/Candidate_Transporters_Summary_Table.csv"
OUT    <- if (length(args) >= 3) args[3] else "final/main_figs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TISSUES <- c("Kernel", "Cob", "Husk", "Sheath", "Stem", "Leaf", "Root")
SRC_COL <- c(GB = "#C51B7D", SN = "#1B7837", IC = "#2166AC", J = "#B35806", QTL = "#252525")
SRC_LAB <- c(GB = "Amino acid GWAS", SN = "Soil N GWAS", IC = "Indian Chief", J = "Jarvis", QTL = "Asp/Pro QTL")

## genes, row order, labels and the analysis that found each one
GENES <- data.table(
  panel = c(rep("A", 6), rep("B", 3)),
  gene_id = paste0("Zm00001eb", c("152450", "094670", "220040", "214360", "405410", "294790",
                                  "277160", "132750", "429100")),
  label = c("Aspartate aminotransferase", "Aspartate kinase", "Aspartate kinase isozyme",
            "Asp-semialdehyde dehydrogenase", "Homoserine dehydrogenase", "Threonine synthase",
            "NPF4.6", "LHT-like", "AVT3A-like"),
  source = list("J", "J", "SN", "IC", "GB", "IC", "GB", c("J", "QTL"), "GB"))
GENES[, row := sprintf("%s (%s)", label, sub("Zm00001eb", "", gene_id))]
GENES[, ord := seq_len(.N), by = panel]

rd <- function(f) {
  d <- fread(f)
  m <- melt(d, id.vars = "gene_id", measure.vars = paste0("Mean_", TISSUES),
            variable.name = "tissue", value.name = "expr")
  m[, tissue := factor(sub("Mean_", "", tissue), levels = TISSUES)]
  m
}
E <- rbind(rd(P_PATH), rd(P_TRAN))
E <- merge(E, GENES, by = "gene_id")
stopifnot(E[, uniqueN(gene_id)] == nrow(GENES))
LIM <- c(0, ceiling(max(E$expr) * 2) / 2)

thm <- theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 13, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_blank(),
        plot.tag = element_text(size = 22, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        plot.margin = margin(8, 8, 8, 8))

heat <- function(p, tag, show_x = TRUE) {
  d <- E[panel == p]; g <- GENES[panel == p]
  lv <- rev(g[order(ord), row])
  d[, row := factor(row, levels = lv)]
  s <- copy(g)[, row := factor(row, levels = lv)]
  ## source strip, one column per analysis, filled where the gene was found
  st <- CJ(row = factor(lv, levels = lv), src = factor(names(SRC_COL), levels = names(SRC_COL)))
  hit <- s[, .(src = unlist(source)), by = row][, hit := TRUE]
  st <- merge(st, hit, by = c("row", "src"), all.x = TRUE)
  st[, fill := fifelse(!is.na(hit), SRC_COL[as.character(src)], "white")]
  st[, src := factor(as.character(src), levels = names(SRC_COL))]
  st[, row := factor(as.character(row), levels = lv)]
  strip <- ggplot(st, aes(src, row)) +
    geom_tile(aes(fill = I(fill)), colour = "grey70", linewidth = 0.4) +
    scale_x_discrete(position = "top", labels = SRC_LAB, drop = FALSE) +
    thm + theme(axis.text.x = element_text(size = 12, face = "bold", colour = "black",
                                           angle = 90, hjust = 0, vjust = 0.5))
  if (!show_x) strip <- strip + theme(axis.text.x = element_blank())
  hm <- ggplot(d, aes(tissue, row, fill = expr)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", expr),
                  colour = I(ifelse(expr > 0.55 * LIM[2], "black", "white"))), size = 4.2) +
    scale_fill_viridis_c(option = "D", limits = LIM, name = "log(TPM + 1)") +
    scale_x_discrete(position = "top") +
    thm + theme(axis.text.y = element_blank())
  if (!show_x) hm <- hm + theme(axis.text.x = element_blank())
  (strip + labs(tag = tag)) + hm + plot_layout(widths = c(5, 9.5))
}

pA <- heat("A", "A")
pB <- heat("B", "B", show_x = FALSE)
fig <- (pA / pB) + plot_layout(heights = c(6, 3), guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUT, "Fig_candidate_expression_tissues.png"), fig, width = 13.5, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_candidate_expression_tissues.pdf"), fig, width = 13.5, height = 7.5, bg = "white")
cat("wrote", file.path(OUT, "Fig_candidate_expression_tissues.{png,pdf}"), "\n")

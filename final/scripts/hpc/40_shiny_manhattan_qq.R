#!/usr/bin/env Rscript
################################################################################
## 40 -- MANHATTAN + QQ PLOTS FOR THE SHINY APP (one trait per call)
##
## For one trait's split_by_trait/<trait>.csv (all GAPIT models stacked):
##   * drops MLMM (it gave no significant SNP for any trait)
##   * lifts every SNP from AGPv4 to B73 RefGen_v5 with the UCSC chain file
##     (same rule as final/scripts/liftover_v4_v5.R; the lift is cached after
##     the first trait because every trait has the same markers)
##   * Manhattan: EVERY SNP is drawn, no thinning. One point per SNP, taking the
##     model with the lowest P. Grey bands = not significant. SNPs significant
##     at 5% Benjamini-Hochberg FDR (GAPIT's FDR column, per model) are coloured
##     by the best model among those calling them, as in Figure 2D. One dashed
##     line at the most lenient per-model FDR cutoff, so every significant SNP
##     sits on or above it. No title, no arrows, no gene labels.
##   * QQ: observed vs expected -log10(P), one series per model, lambda in the
##     legend.
##
## Output:  <out_dir>/manhattan/<trait>_manhattan.png
##          <out_dir>/qq/<trait>_qq.png
## These names are what shiny/amino_gwas_app/app.R looks for.
##
## Usage:   Rscript 40_shiny_manhattan_qq.R <trait> <in_dir> <out_dir> <chain_file>
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: Rscript 40_shiny_manhattan_qq.R <trait> <in_dir> <out_dir> <chain_file>")
TRAIT <- args[1]; IN_DIR <- args[2]; OUT_DIR <- args[3]; CHAIN <- args[4]
MODELS    <- c("MLM", "BLINK", "FarmCPU")
MODEL_COL <- c(MLM = "#0072B2", BLINK = "#009E73", FarmCPU = "#D55E00")
Q_FDR     <- 0.05
M_SNP     <- 4177796
dir.create(file.path(OUT_DIR, "manhattan"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "qq"),        recursive = TRUE, showWarnings = FALSE)

## png device: ragg if installed, otherwise cairo (no X11 on compute nodes)
DEV <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
  function(filename, ...) grDevices::png(filename, ..., type = "cairo")

## ------------------------------------------------------------ lift-over ----
read_chain_blocks <- function(chain_file) {
  L <- readLines(chain_file); out <- vector("list", 0); k <- 0L; i <- 1L; n <- length(L)
  while (i <= n) {
    if (startsWith(L[i], "chain")) {
      f <- strsplit(L[i], " ", fixed = TRUE)[[1]]
      score <- as.numeric(f[2]); tchr <- f[3]; tstart <- as.numeric(f[6])
      qchr <- f[8]; qstrand <- f[10]; qstart <- as.numeric(f[11])
      keep <- tchr == qchr && qstrand == "+" && grepl("^[0-9]+$", tchr)
      j <- i + 1L
      while (j <= n && nzchar(L[j]) && !startsWith(L[j], "chain")) j <- j + 1L
      if (keep) {
        m <- do.call(rbind, lapply(strsplit(L[(i + 1L):(j - 1L)], "[ \t]+"), function(x) as.numeric(c(x, 0, 0)[1:3])))
        size <- m[, 1]; dt <- m[, 2]; dq <- m[, 3]
        ts <- tstart + c(0, cumsum(size + dt)[-length(size)])
        qs <- qstart + c(0, cumsum(size + dq)[-length(size)])
        k <- k + 1L
        out[[k]] <- data.table(chr = as.integer(tchr), tS = ts, tE = ts + size, qS = qs, score = score)
      }
      i <- j
    } else i <- i + 1L
  }
  B <- rbindlist(out); setorder(B, chr, tS); B
}
lift_v4_v5 <- function(chr, pos, blocks) {
  res <- data.table(i = seq_along(pos), chr = as.integer(chr), x = pos - 1)
  res[, `:=`(pos_v5 = NA_real_)]
  for (c in unique(res$chr)) {
    B <- blocks[chr == c]; if (!nrow(B)) next
    idx <- res[chr == c, i]; x <- res$x[idx]
    j <- findInterval(x, B$tS)
    best <- rep(NA_real_, length(x)); bscore <- rep(-Inf, length(x))
    for (off in 0:30) {
      jj <- j - off; ok <- jj >= 1
      hit <- ok & B$tS[pmax(jj, 1)] <= x & B$tE[pmax(jj, 1)] > x & B$score[pmax(jj, 1)] > bscore
      best[hit] <- B$qS[jj[hit]] + (x[hit] - B$tS[jj[hit]]); bscore[hit] <- B$score[jj[hit]]
    }
    gap <- is.na(best)
    if (any(gap)) {
      jg <- pmax(j[gap], 1); jn <- pmin(jg + 1, nrow(B)); xg <- x[gap]
      dl <- abs(xg - (B$tE[jg] - 1)); dr <- abs(B$tS[jn] - xg)
      best[gap] <- ifelse(dr < dl, B$qS[jn], B$qS[jg] + (B$tE[jg] - B$tS[jg] - 1))
    }
    res[idx, pos_v5 := best + 1]
  }
  res$pos_v5
}

## ------------------------------------------------------------- read scan ---
f_in <- file.path(IN_DIR, paste0(TRAIT, ".csv"))
if (!file.exists(f_in)) stop("missing: ", f_in)
d <- fread(f_in, select = c("Chr", "Pos", "P.value", "model", "FDR"), showProgress = FALSE)
setnames(d, c("Chr", "Pos", "P", "Model", "FDR"))
d <- d[Model %in% MODELS & is.finite(P) & P > 0 & Chr >= 1 & Chr <= 10]
d[, Chr := as.integer(Chr)]
cat(sprintf("%s: %s rows, models %s\n", TRAIT, format(nrow(d), big.mark = ","),
            paste(sort(unique(d$Model)), collapse = ",")))

cache <- file.path(OUT_DIR, "lift_v4_v5_cache.rds")
key <- unique(d[, .(Chr, Pos)])
if (file.exists(cache)) {
  L <- readRDS(cache)
  miss <- key[!L, on = .(Chr, Pos)]
  if (nrow(miss)) {
    blocks <- read_chain_blocks(CHAIN)
    miss[, Pos5 := lift_v4_v5(Chr, Pos, blocks)]
    L <- rbind(L, miss); saveRDS(L, cache)
  }
} else {
  blocks <- read_chain_blocks(CHAIN)
  key[, Pos5 := lift_v4_v5(Chr, Pos, blocks)]
  L <- key; saveRDS(L, cache)
}
d <- L[d, on = .(Chr, Pos)]
d <- d[is.finite(Pos5)]

## ------------------------------------------------------------ thresholds ---
d[, sig := is.finite(FDR) & FDR <= Q_FDR]
cut <- d[sig == TRUE, .(cutoff = max(P)), by = Model]
fdr_line <- if (nrow(cut)) -log10(max(cut$cutoff)) else -log10(Q_FDR / M_SNP)
line_lab <- if (nrow(cut)) "FDR 5%" else "FDR 5% (no significant SNP)"

## one point per SNP: best significant model if any, else best model overall
setorder(d, Chr, Pos, P)
sig_best <- unique(d[sig == TRUE], by = c("Chr", "Pos"))
all_best <- unique(d, by = c("Chr", "Pos"))
best <- rbind(sig_best, all_best[!sig_best, on = .(Chr, Pos)])
best[, logp := -log10(P)]

chrlen <- best[, .(len = max(Pos5)), by = Chr][order(Chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
chrlen[, centre := offset + len / 2]
best <- chrlen[, .(Chr, offset)][best, on = "Chr"]
best[, x := Pos5 + offset]
bg <- best[sig == FALSE]; bg[, band := fifelse(Chr %% 2 == 0, "even", "odd")]
sg <- best[sig == TRUE]
cat(sprintf("  %s SNPs drawn, %d FDR-significant, line at -log10(p) = %.2f\n",
            format(nrow(best), big.mark = ","), nrow(sg), fdr_line))

## ------------------------------------------------------------ house theme --
plot_theme <- theme_minimal(base_size = 24) +
  theme(axis.title   = element_text(size = 18, face = "bold"),
        axis.text.x  = element_text(size = 16, colour = "black"),
        axis.text.y  = element_text(size = 16, face = "bold", colour = "black"),
        axis.line    = element_line(colour = "black"),
        panel.grid   = element_blank(),
        legend.title = element_blank(),
        legend.text  = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = "grey70", linewidth = 0.3),
        plot.margin  = margin(12, 15, 10, 12))
inside <- function(x, y, just) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just)
  else theme(legend.position = c(x, y), legend.justification = just)
}

## -------------------------------------------------------------- Manhattan --
ymax <- max(c(best$logp, fdr_line)) * 1.08
pm <- ggplot() +
  geom_point(data = bg, aes(x, logp, colour = band), size = 0.9, alpha = 0.6, shape = 16, show.legend = FALSE) +
  geom_hline(yintercept = fdr_line, linetype = "dashed", colour = "red3", linewidth = 0.9) +
  annotate("text", x = min(best$x), y = fdr_line, label = line_lab, hjust = 0, vjust = -0.6,
           size = 5, fontface = "bold", colour = "red3") +
  geom_point(data = sg, aes(x, logp, colour = Model), size = 2.6, alpha = 0.95, shape = 16) +
  scale_colour_manual(values = c(odd = "grey68", even = "grey45", MODEL_COL), breaks = MODELS, drop = FALSE) +
  scale_x_continuous(breaks = chrlen$centre, labels = chrlen$Chr, expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Chromosome", y = expression(bold(-log[10](bolditalic(p))))) +
  plot_theme + inside(0.005, 0.995, c(0, 1)) +
  theme(legend.direction = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))
ggsave(file.path(OUT_DIR, "manhattan", paste0(TRAIT, "_manhattan.png")), pm,
       width = 13, height = 5.5, dpi = 300, bg = "white", device = DEV)

## --------------------------------------------------------------------- QQ --
qq <- d[, {
  o <- sort(P); n <- length(o)
  lam <- median(qchisq(1 - P, 1)) / qchisq(0.5, 1)
  .(exp = -log10((seq_len(n) - 0.5) / n), obs = -log10(o), lambda = lam)
}, by = Model]
labs_l <- qq[, .(lambda = lambda[1]), by = Model]
labs_l <- labs_l[match(intersect(MODELS, Model), Model)]
lab_map <- parse(text = sprintf("%s~(lambda=='%.2f')", labs_l$Model, labs_l$lambda))
lim <- max(c(qq$exp, qq$obs))
pq <- ggplot(qq, aes(exp, obs, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1.1, shape = 16, alpha = 0.8) +
  scale_colour_manual(values = MODEL_COL, labels = lab_map, breaks = labs_l$Model) +
  coord_cartesian(xlim = c(0, max(qq$exp) * 1.02), ylim = c(0, lim * 1.02), expand = FALSE) +
  labs(x = expression(bold(Expected ~ -log[10](bolditalic(p)))),
       y = expression(bold(Observed ~ -log[10](bolditalic(p))))) +
  plot_theme + inside(0.02, 0.98, c(0, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))
ggsave(file.path(OUT_DIR, "qq", paste0(TRAIT, "_qq.png")), pq,
       width = 6.5, height = 6.5, dpi = 300, bg = "white", device = DEV)
cat("  written\n")

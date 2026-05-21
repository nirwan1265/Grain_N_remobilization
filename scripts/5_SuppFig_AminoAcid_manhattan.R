################################################################################
### LIBRARIES
################################################################################
library(dplyr)
library(data.table)
library(ggplot2)
library(vroom)
library(ggnewscale)
library(grid)

################################################################################
### CONFIGURATION
################################################################################

INPUT_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/Sarah_amino_acid/N_grain/Phenotypes_GWAS_Grain"
OUTPUT_DIR <- "Figs/Supp"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

AVAILABLE_TRAITS <- list.files(INPUT_DIR, pattern = "\\.csv$", full.names = FALSE)
AVAILABLE_TRAITS <- sub("\\.csv$", "", AVAILABLE_TRAITS)
AVAILABLE_TRAITS <- sort(AVAILABLE_TRAITS)

# Keep only the single amino-acid traits for the supplement.
TARGET_TRAITS <- c("D", "E", "N", "P", "Q", "Total_N", "Total_PBAA")
trait_files <- TARGET_TRAITS[TARGET_TRAITS %in% AVAILABLE_TRAITS]

N_TRAITS <- length(trait_files)
N_COLS <- 2
N_ROWS <- N_TRAITS

THR_SUGGESTIVE <- 5
THR_SIGNIFICANT <- 7
MANHATTAN_BG_MAX <- 150000
QQ_MAX_POINTS <- 120000

MODEL_COLORS_BASE <- c(
  "MLM" = "#0072B2",
  "MLMM" = "#E69F00",
  "BLINK" = "#009E73",
  "FarmCPU" = "#D55E00",
  "GLM" = "#CC79A7",
  "SUPER" = "#56B4E9"
)
MODEL_FALLBACK <- c("#999999", "#F0E442", "#000000", "#8C564B")

COMBINED_FILE <- file.path(OUTPUT_DIR, "Supp_AminoAcid_Selected_Manhattan_QQ.png")

cat("Found", N_TRAITS, "trait files\n")
print(trait_files)
cat("Layout:", N_ROWS, "rows x", N_COLS, "columns (Manhattan | QQ)\n")

################################################################################
### THEME
################################################################################

plot_theme <- theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6)
  )

################################################################################
### HELPERS
################################################################################

assign_model_colors <- function(models) {
  known <- intersect(models, names(MODEL_COLORS_BASE))
  colors <- MODEL_COLORS_BASE[known]
  unknown <- setdiff(models, names(MODEL_COLORS_BASE))

  if (length(unknown) > 0) {
    fallback <- rep(MODEL_FALLBACK, length.out = length(unknown))
    names(fallback) <- unknown
    colors <- c(colors, fallback)
  }

  colors[models]
}

thin_rows <- function(df, max_n) {
  if (is.null(df) || nrow(df) <= max_n) {
    return(df)
  }

  keep_idx <- unique(round(seq(1, nrow(df), length.out = max_n)))
  df[keep_idx, , drop = FALSE]
}

legend_models <- function() {
  base_order <- c("MLM", "MLMM", "BLINK", "FarmCPU")
  present <- base_order[base_order %in% names(MODEL_COLORS_BASE)]
  if (length(present) == 0) {
    present <- names(MODEL_COLORS_BASE)
  }
  present
}

make_shared_legend_plot <- function() {
  legend_df <- data.frame(
    Model = factor(legend_models(), levels = legend_models()),
    x = seq_along(legend_models()),
    y = 1
  )

  ggplot(legend_df, aes(x = x, y = y, color = Model)) +
    geom_point(size = 3) +
    scale_color_manual(values = MODEL_COLORS_BASE[legend_models()], drop = FALSE) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4))) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
}

load_trait_data <- function(trait_name) {
  file_path <- file.path(INPUT_DIR, paste0(trait_name, ".csv"))
  raw_dt <- data.table::fread(
    file_path,
    select = c("SNP", "Chr", "Pos", "P.value", "model"),
    showProgress = FALSE
  )

  setnames(raw_dt, c("Chr", "Pos", "P.value", "model"), c("CHR", "BP", "P", "Model"))
  raw_dt[, `:=`(
    CHR = as.integer(CHR),
    BP = as.numeric(BP),
    P = as.numeric(P)
  )]
  raw_dt <- raw_dt[!is.na(CHR) & !is.na(BP) & !is.na(P) & is.finite(P) & P > 0]
  raw_dt[, logp := -log10(P)]
  raw_dt <- raw_dt[is.finite(logp)]

  models <- unique(as.character(raw_dt$Model))
  model_colors <- assign_model_colors(models)

  setorder(raw_dt, SNP, CHR, BP, P)
  best_dt <- raw_dt[, .SD[1], by = .(SNP, CHR, BP)]
  setorder(best_dt, CHR, BP)
  setnames(best_dt, "Model", "best_model")

  chr_tbl <- best_dt[, .(chr_len = max(BP, na.rm = TRUE)), by = CHR]
  setorder(chr_tbl, CHR)
  chr_tbl[, offset := shift(cumsum(chr_len), fill = 0)]
  chr_tbl[, center := offset + chr_len / 2]

  plot_dt <- merge(
    best_dt[, .(SNP, CHR, BP, P, logp, best_model)],
    chr_tbl[, .(CHR, offset)],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  plot_dt[, pos_cum := BP + offset]
  setorder(plot_dt, CHR, BP)

  bg_df <- as.data.frame(plot_dt[logp < THR_SUGGESTIVE])
  sig_df <- as.data.frame(plot_dt[logp >= THR_SUGGESTIVE])
  bg_df <- thin_rows(bg_df, MANHATTAN_BG_MAX)

  keep_n <- min(QQ_MAX_POINTS, nrow(best_dt))
  keep_idx <- unique(round(seq(1, nrow(best_dt), length.out = keep_n)))
  best_sorted <- best_dt[order(P)]
  qq_df <- data.frame(
    expected = -log10(ppoints(nrow(best_sorted)))[keep_idx],
    observed = best_sorted$logp[keep_idx],
    best_model = best_sorted$best_model[keep_idx]
  )

  ci_n <- min(1000, nrow(best_sorted))
  ci_idx <- sort(unique(round(seq(1, nrow(best_sorted), length.out = ci_n))))
  conf_int <- data.frame(
    expected = -log10(ppoints(nrow(best_sorted)))[ci_idx],
    lower = -log10(qbeta(0.975, ci_idx, nrow(best_sorted) - ci_idx + 1)),
    upper = -log10(qbeta(0.025, ci_idx, nrow(best_sorted) - ci_idx + 1))
  )

  lambda <- median(qchisq(1 - best_dt$P, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

  list(
    bg_df = bg_df,
    sig_df = sig_df,
    qq_df = qq_df,
    conf_int = conf_int,
    chr_tbl = as.data.frame(chr_tbl),
    model_colors = model_colors,
    lambda = lambda
  )
}

make_manhattan_plot <- function(dat, trait_name, letter_label) {
  p <- ggplot() +
    geom_point(
      data = dat$bg_df,
      aes(x = pos_cum, y = logp, color = factor(CHR %% 2)),
      alpha = 0.35,
      size = 0.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("0" = "grey65", "1" = "grey40")) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = dat$sig_df,
      aes(x = pos_cum, y = logp, color = best_model),
      size = 1,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    scale_color_manual(values = dat$model_colors) +
    geom_hline(yintercept = THR_SUGGESTIVE, linetype = "dashed", color = "black", linewidth = 0.35) +
    geom_hline(yintercept = THR_SIGNIFICANT, linetype = "solid", color = "black", linewidth = 0.45) +
    scale_x_continuous(
      breaks = dat$chr_tbl$center,
      labels = dat$chr_tbl$CHR,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(p))),
      title = trait_name,
      tag = letter_label
    ) +
    plot_theme
  p
}

make_qq_plot <- function(dat, trait_name, letter_label) {
  p <- ggplot(dat$qq_df, aes(x = expected, y = observed)) +
    geom_ribbon(
      data = dat$conf_int,
      aes(x = expected, ymin = lower, ymax = upper),
      fill = "grey85",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.45) +
    geom_point(aes(color = best_model), alpha = 0.55, size = 0.7, show.legend = FALSE) +
    scale_color_manual(values = dat$model_colors) +
    annotate(
      "text",
      x = 0.3,
      y = max(dat$qq_df$observed, na.rm = TRUE) * 0.9,
      label = paste0("lambda = ", round(dat$lambda, 2)),
      hjust = 0,
      size = 2.7,
      fontface = "bold"
    ) +
    labs(
      x = expression(Expected~~-log[10](italic(p))),
      y = expression(Observed~~-log[10](italic(p))),
      title = trait_name,
      tag = letter_label
    ) +
    plot_theme
  p
}

init_plot_device <- function(output_file, width_in, height_in) {
  png(
    output_file,
    width = width_in,
    height = height_in,
    units = "in",
    res = 300,
    type = "cairo"
  )
  dev_id <- dev.cur()
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = N_ROWS + 1,
        ncol = N_COLS,
        heights = unit.c(unit(0.7, "in"), rep(unit(1, "null"), N_ROWS))
      )
    )
  )
  dev_id
}

render_combined_grid <- function() {
  shared_legend <- make_shared_legend_plot()

  combo_dev <- init_plot_device(
    output_file = COMBINED_FILE,
    width_in = 16,
    height_in = 3.0 * N_ROWS + 0.8
  )

  dev.set(combo_dev)
  print(shared_legend, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:N_COLS))

  for (i in seq_along(trait_files)) {
    trait_name <- trait_files[i]
    letter_label <- LETTERS[i]
    row_idx <- i + 1

    cat("Rendering [", letter_label, "]", trait_name, "\n")

    dat <- load_trait_data(trait_name)
    p_man <- make_manhattan_plot(dat, trait_name, letter_label)
    p_qq <- make_qq_plot(dat, trait_name, letter_label)

    print(
      p_man,
      vp = viewport(layout.pos.row = row_idx, layout.pos.col = 1)
    )
    cat("  Saved Manhattan panel for", trait_name, "\n")

    print(
      p_qq,
      vp = viewport(layout.pos.row = row_idx, layout.pos.col = 2)
    )
    cat("  Saved QQ panel for", trait_name, "\n")

    rm(dat, p_man, p_qq)
    gc()
  }

  dev.set(combo_dev)
  cat("Closing PNG device...\n")
  dev.off()
}

################################################################################
### SAVE COMPOSITE FIGURES
################################################################################

render_combined_grid()

################################################################################
### FINAL OUTPUT
################################################################################

cat("\nSaved combined amino-acid Manhattan + QQ figure to:\n")
cat("  ", COMBINED_FILE, "\n", sep = "")
cat("\nDone.\n")

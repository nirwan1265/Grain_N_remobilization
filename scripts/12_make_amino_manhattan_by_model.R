################################################################################
### AMINO-ACID GWAS MANHATTAN PLOTS BY MODEL
################################################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(vroom)

################################################################################
### CONFIGURATION
################################################################################

INPUT_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/Sarah_amino_acid/N_grain/Phenotypes_GWAS_Grain"
OUTPUT_DIR <- "Figs/Supplementary/amino_manhattan_by_model"

TARGET_TRAITS <- sort(sub("\\.csv$", "", list.files(INPUT_DIR, pattern = "\\.csv$")))
MODEL_ORDER <- c("MLM", "MLMM", "BLINK", "FarmCPU")
OVERLAY_SETS <- list(
  All = MODEL_ORDER
)

THR_SUGGESTIVE <- 5
THR_SIGNIFICANT <- 7
MAX_BG_POINTS <- 150000

MODEL_COLORS <- c(
  MLM = "#0072B2",
  MLMM = "#E69F00",
  BLINK = "#009E73",
  FarmCPU = "#D55E00"
)

TRAIT_LABELS <- c(
  D = "Aspartate",
  `D.IMNTDK` = "Aspartate relative to aspartate family",
  `D.Total` = "Aspartate relative to all free amino acids",
  E = "Glutamate",
  `E.EHPRQ` = "Glutamate relative to glutamate family",
  `E.Total` = "Glutamate relative to all free amino acids",
  EHPRQ = "Glutamate family",
  IMNTDK = "Aspartate family",
  N = "Asparagine",
  `N.E` = "Asparagine relative to glutamate",
  `N.IMNTDK` = "Asparagine relative to aspartate family",
  `N.Q` = "Asparagine relative to glutamine",
  `N.Total` = "Asparagine relative to all free amino acids",
  P = "Proline",
  `P.EHPRQ` = "Proline relative to glutamate family",
  `P.EHPRQ_NYC_` = "Proline relative to glutamate family (NYC)",
  `P.Total` = "Proline relative to all free amino acids",
  Q = "Glutamine",
  `Q.E` = "Glutamine relative to glutamate",
  `Q.Total` = "Glutamine relative to all amino acids",
  Total_N = "Total nitrogen",
  Total_PBAA = "Total protein-bound amino acids"
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "single_model"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "overlay"), showWarnings = FALSE, recursive = TRUE)

################################################################################
### THEME
################################################################################

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey35"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

################################################################################
### HELPERS
################################################################################

pretty_trait_name <- function(trait_name) {
  if (trait_name %in% names(TRAIT_LABELS)) {
    paste0(trait_name, " \u2014 ", TRAIT_LABELS[[trait_name]])
  } else {
    trait_name
  }
}

thin_rows <- function(df, max_n) {
  if (is.null(df) || nrow(df) <= max_n) {
    return(df)
  }

  keep_idx <- unique(round(seq(1, nrow(df), length.out = max_n)))
  df[keep_idx, , drop = FALSE]
}

load_trait_data <- function(trait_name) {
  file_path <- file.path(INPUT_DIR, paste0(trait_name, ".csv"))

  data.table::fread(
    file_path,
    select = c("SNP", "Chr", "Pos", "P.value", "model"),
    showProgress = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      Chr = as.integer(Chr),
      Pos = as.numeric(Pos),
      P.value = as.numeric(P.value),
      model = as.character(model),
      log10_P = -log10(P.value)
    ) %>%
    filter(
      !is.na(Chr),
      !is.na(Pos),
      !is.na(P.value),
      is.finite(P.value),
      P.value > 0,
      is.finite(log10_P)
    )
}

prepare_manhattan_data <- function(trait_df, selected_models) {
  df <- trait_df %>%
    filter(model %in% selected_models)

  if (nrow(df) == 0) {
    return(NULL)
  }

  best_df <- df %>%
    group_by(SNP, Chr, Pos) %>%
    slice_min(P.value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(best_model = model)

  chr_tbl <- best_df %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Pos, na.rm = TRUE), .groups = "drop") %>%
    arrange(Chr) %>%
    mutate(
      offset = lag(cumsum(chr_len), default = 0),
      center = offset + chr_len / 2
    )

  bg_df <- best_df %>%
    left_join(chr_tbl %>% select(Chr, offset), by = "Chr") %>%
    mutate(pos_cum = Pos + offset) %>%
    filter(log10_P < THR_SUGGESTIVE)

  sig_df <- df %>%
    left_join(chr_tbl %>% select(Chr, offset), by = "Chr") %>%
    mutate(pos_cum = Pos + offset) %>%
    filter(log10_P >= THR_SUGGESTIVE)

  list(
    bg_df = thin_rows(bg_df, MAX_BG_POINTS),
    sig_df = sig_df,
    chr_tbl = chr_tbl,
    ymax = max(c(best_df$log10_P, THR_SIGNIFICANT), na.rm = TRUE)
  )
}

make_manhattan_plot <- function(plot_dat, trait_name, selected_models, overlay_label) {
  if (is.null(plot_dat)) {
    return(NULL)
  }

  title_text <- pretty_trait_name(trait_name)
  subtitle_text <- if (length(selected_models) == 1) {
    paste0("Model: ", selected_models[[1]])
  } else {
    paste0("Overlay: ", overlay_label)
  }

  ggplot() +
    geom_point(
      data = plot_dat$bg_df,
      aes(x = pos_cum, y = log10_P, color = factor(Chr %% 2)),
      alpha = 0.4,
      size = 0.45,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("0" = "grey65", "1" = "grey40")) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = plot_dat$sig_df,
      aes(x = pos_cum, y = log10_P, color = model),
      alpha = 0.85,
      size = 1.1
    ) +
    scale_color_manual(
      values = MODEL_COLORS[names(MODEL_COLORS) %in% selected_models],
      breaks = selected_models,
      drop = FALSE
    ) +
    geom_hline(yintercept = THR_SUGGESTIVE, linetype = "dashed", color = "black", linewidth = 0.45) +
    geom_hline(yintercept = THR_SIGNIFICANT, linetype = "solid", color = "black", linewidth = 0.55) +
    annotate(
      "label",
      x = Inf,
      y = THR_SIGNIFICANT,
      label = paste0("-log10(p) = ", THR_SIGNIFICANT),
      hjust = 1.05,
      vjust = -0.25,
      size = 4,
      label.size = NA,
      fill = "white"
    ) +
    scale_x_continuous(
      breaks = plot_dat$chr_tbl$center,
      labels = plot_dat$chr_tbl$Chr,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, plot_dat$ymax * 1.05),
      expand = c(0, 0)
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Chromosome",
      y = expression(-log[10](italic(p))),
      color = NULL
    ) +
    plot_theme
}

################################################################################
### MAIN LOOP
################################################################################

run_all_manhattan <- function(traits = TARGET_TRAITS, output_dir = OUTPUT_DIR) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "single_model"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "overlay"), showWarnings = FALSE, recursive = TRUE)

  cat("Traits found:", length(traits), "\n")

  for (trait_name in traits) {
    cat("\n=== Trait:", trait_name, "===\n")
    trait_df <- load_trait_data(trait_name)
    present_models <- MODEL_ORDER[MODEL_ORDER %in% unique(trait_df$model)]

    cat("Rows:", nrow(trait_df), "\n")
    cat("Models:", paste(present_models, collapse = ", "), "\n")

    for (model_name in present_models) {
      plot_dat <- prepare_manhattan_data(trait_df, selected_models = model_name)
      p <- make_manhattan_plot(
        plot_dat = plot_dat,
        trait_name = trait_name,
        selected_models = model_name,
        overlay_label = model_name
      )

      out_file <- file.path(
        output_dir,
        "single_model",
        paste0(trait_name, "_", model_name, "_manhattan.png")
      )

      ggsave(out_file, p, width = 13, height = 5.5, dpi = 300, bg = "white")
      cat("  saved:", out_file, "\n")
    }

    for (overlay_name in names(OVERLAY_SETS)) {
      selected_models <- OVERLAY_SETS[[overlay_name]]
      selected_models <- selected_models[selected_models %in% present_models]

      if (length(selected_models) == 0) {
        next
      }

      plot_dat <- prepare_manhattan_data(trait_df, selected_models = selected_models)
      p <- make_manhattan_plot(
        plot_dat = plot_dat,
        trait_name = trait_name,
        selected_models = selected_models,
        overlay_label = overlay_name
      )

      out_file <- file.path(
        output_dir,
        "overlay",
        paste0(trait_name, "_", overlay_name, "_overlay_manhattan.png")
      )

      ggsave(out_file, p, width = 13, height = 5.5, dpi = 300, bg = "white")
      cat("  saved:", out_file, "\n")
    }
  }

  cat("\nDone.\n")
}

if (sys.nframe() == 0) {
  run_all_manhattan()
}

################################################################################
### AMINO-ACID GWAS MANHATTAN + QQ PLOTS BY MODEL
################################################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(vroom)

################################################################################
### CONFIGURATION
################################################################################

INPUT_DIR <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
OUTPUT_DIR <- "Figs"

TARGET_TRAITS <- sort(sub("\\.csv$", "", list.files(INPUT_DIR, pattern = "\\.csv$")))
MODEL_ORDER <- c("MLM", "MLMM", "BLINK", "FarmCPU")

THR_SUGGESTIVE <- 5
THR_SIGNIFICANT <- 7
MAX_BG_POINTS <- 150000
MAX_QQ_POINTS <- 120000

MODEL_COLORS <- c(
  MLM = "#0072B2",
  MLMM = "#E69F00",
  BLINK = "#009E73",
  FarmCPU = "#D55E00"
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "manhattan"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "qq"), showWarnings = FALSE, recursive = TRUE)

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

thin_rows <- function(df, max_n) {
  if (is.null(df) || nrow(df) <= max_n) {
    return(df)
  }

  keep_idx <- unique(round(seq(1, nrow(df), length.out = max_n)))
  df[keep_idx, , drop = FALSE]
}

load_trait_data <- function(trait_name) {
  file_path <- file.path(INPUT_DIR, paste0(trait_name, ".csv"))

  dt <- data.table::fread(
    file_path,
    select = c("SNP", "Chr", "Pos", "P.value", "model"),
    showProgress = FALSE
  )

  data.table::setDT(dt)
  dt[, Chr := as.integer(Chr)]
  dt[, Pos := as.numeric(Pos)]
  dt[, P.value := as.numeric(P.value)]
  dt[, model := trimws(as.character(model))]
  dt[toupper(model) == "MLM", model := "MLM"]
  dt[toupper(model) == "MLMM", model := "MLMM"]
  dt[toupper(model) == "BLINK", model := "BLINK"]
  dt[toupper(model) == "FARMCPU", model := "FarmCPU"]
  dt[, log10_P := -log10(P.value)]

  dt <- dt[
    !is.na(Chr) &
      !is.na(Pos) &
      !is.na(P.value) &
      is.finite(P.value) &
      P.value > 0 &
      is.finite(log10_P)
  ]

  dt[]
}

prepare_manhattan_data <- function(trait_df, selected_models) {
  df <- data.table::as.data.table(data.table::copy(trait_df))
  df <- df[model %chin% selected_models]

  if (nrow(df) == 0) {
    return(NULL)
  }

  data.table::setorder(df, SNP, Chr, Pos, P.value)
  best_df <- df[, .SD[1L], by = .(SNP, Chr, Pos)]
  data.table::setnames(best_df, "model", "best_model")

  chr_tbl <- best_df[, .(chr_len = max(Pos, na.rm = TRUE)), by = Chr]
  data.table::setorder(chr_tbl, Chr)
  chr_tbl[, offset := data.table::shift(cumsum(chr_len), fill = 0)]
  chr_tbl[, center := offset + chr_len / 2]

  bg_df <- merge(best_df, chr_tbl[, .(Chr, offset)], by = "Chr", sort = FALSE)
  bg_df[, pos_cum := Pos + offset]
  bg_df <- bg_df[log10_P < THR_SUGGESTIVE]

  sig_df <- merge(df, chr_tbl[, .(Chr, offset)], by = "Chr", sort = FALSE)
  sig_df[, pos_cum := Pos + offset]
  sig_df <- sig_df[log10_P >= THR_SUGGESTIVE]

  list(
    bg_df = thin_rows(bg_df, MAX_BG_POINTS),
    sig_df = sig_df,
    chr_tbl = chr_tbl,
    ymax = max(c(best_df$log10_P, THR_SIGNIFICANT), na.rm = TRUE)
  )
}

prepare_qq_data <- function(trait_df, selected_models) {
  df <- data.table::as.data.table(data.table::copy(trait_df))
  df <- df[model %chin% selected_models]

  if (nrow(df) == 0) {
    return(NULL)
  }

  qq_df <- bind_rows(lapply(selected_models, function(model_name) {
    model_df <- data.table::copy(df[model == model_name])
    data.table::setorder(model_df, P.value)

    if (nrow(model_df) == 0) {
      return(NULL)
    }

    keep_n <- min(MAX_QQ_POINTS, nrow(model_df))
    keep_idx <- unique(round(seq(1, nrow(model_df), length.out = keep_n)))

    tibble(
      expected = -log10(ppoints(nrow(model_df)))[keep_idx],
      observed = model_df$log10_P[keep_idx],
      model = model_name
    )
  }))

  if (nrow(qq_df) == 0) {
    return(NULL)
  }

  conf_int <- NULL
  lambda_df <- bind_rows(lapply(selected_models, function(model_name) {
    model_df <- df[model == model_name]

    if (nrow(model_df) == 0) {
      return(NULL)
    }

    lambda <- median(
      qchisq(1 - model_df$P.value, df = 1),
      na.rm = TRUE
    ) / qchisq(0.5, df = 1)

    tibble(
      model = model_name,
      lambda = lambda
    )
  }))

  if (length(selected_models) == 1) {
    model_n <- nrow(df[model == selected_models[[1]]])

    ci_n <- min(1000, model_n)
    ci_idx <- sort(unique(round(seq(1, model_n, length.out = ci_n))))

    conf_int <- tibble(
      expected = -log10(ppoints(model_n))[ci_idx],
      lower = -log10(qbeta(0.975, ci_idx, model_n - ci_idx + 1)),
      upper = -log10(qbeta(0.025, ci_idx, model_n - ci_idx + 1))
    )
  }

  list(
    qq_df = qq_df,
    conf_int = conf_int,
    lambda_df = lambda_df
  )
}

make_manhattan_plot <- function(plot_dat, trait_name, selected_models) {
  if (is.null(plot_dat)) {
    return(NULL)
  }

  title_text <- trait_name
  subtitle_text <- paste0("Models: ", paste(selected_models, collapse = ", "))

  bg_even <- plot_dat$bg_df[plot_dat$bg_df$Chr %% 2 == 0, , drop = FALSE]
  bg_odd <- plot_dat$bg_df[plot_dat$bg_df$Chr %% 2 == 1, , drop = FALSE]

  p <- ggplot() +
    geom_point(
      data = bg_even,
      aes(x = pos_cum, y = log10_P),
      color = "grey65",
      alpha = 0.4,
      size = 0.45,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = bg_odd,
      aes(x = pos_cum, y = log10_P),
      color = "grey40",
      alpha = 0.4,
      size = 0.45,
      inherit.aes = FALSE
    )

  if (nrow(plot_dat$sig_df) > 0) {
    p <- p +
      geom_point(
        data = plot_dat$sig_df,
        aes(x = pos_cum, y = log10_P, color = model),
        alpha = 0.85,
        size = 1.1,
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = MODEL_COLORS[selected_models],
        breaks = selected_models,
        drop = FALSE
      )
  }

  p +
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

make_qq_plot <- function(qq_dat, trait_name, selected_models) {
  if (is.null(qq_dat)) {
    return(NULL)
  }

  title_text <- trait_name
  subtitle_text <- paste0("Models: ", paste(selected_models, collapse = ", "))

  lambda_text <- paste(
    paste0(qq_dat$lambda_df$model, "=", sprintf("%.2f", qq_dat$lambda_df$lambda)),
    collapse = "  "
  )

  p <- ggplot(qq_dat$qq_df, aes(x = expected, y = observed)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "black",
      linetype = "dashed",
      linewidth = 0.45
    )

  if (!is.null(qq_dat$conf_int)) {
    p <- p +
      geom_ribbon(
        data = qq_dat$conf_int,
        aes(x = expected, ymin = lower, ymax = upper),
        fill = "grey85",
        alpha = 0.5,
        inherit.aes = FALSE
      )
  }

  p +
    geom_point(
      aes(color = model),
      alpha = 0.6,
      size = 0.75
    ) +
    scale_color_manual(
      values = MODEL_COLORS[selected_models],
      breaks = selected_models,
      drop = FALSE
    ) +
    annotate(
      "text",
      x = 0.3,
      y = max(qq_dat$qq_df$observed, na.rm = TRUE) * 0.92,
      label = paste0("lambda: ", lambda_text),
      hjust = 0,
      size = 3.3,
      fontface = "bold"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = expression(Expected~~-log[10](italic(p))),
      y = expression(Observed~~-log[10](italic(p))),
      color = NULL
    ) +
    plot_theme
}

################################################################################
### MAIN LOOP
################################################################################

run_all_gwas_plots <- function(traits = TARGET_TRAITS, output_dir = OUTPUT_DIR) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "manhattan"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "qq"), showWarnings = FALSE, recursive = TRUE)

  cat("Traits found:", length(traits), "\n")

  for (trait_name in traits) {
    cat("\n=== Trait:", trait_name, "===\n")
    trait_df <- load_trait_data(trait_name)
    present_models <- MODEL_ORDER[MODEL_ORDER %in% unique(trait_df[["model"]])]

    cat("Rows:", nrow(trait_df), "\n")
    cat("Models:", paste(present_models, collapse = ", "), "\n")

    if (length(present_models) == 0) {
      cat("  skipped: no recognized models present\n")
      next
    }

    plot_dat <- prepare_manhattan_data(trait_df, selected_models = present_models)
    qq_dat <- prepare_qq_data(trait_df, selected_models = present_models)

    p_manhattan <- make_manhattan_plot(
      plot_dat = plot_dat,
      trait_name = trait_name,
      selected_models = present_models
    )
    p_qq <- make_qq_plot(
      qq_dat = qq_dat,
      trait_name = trait_name,
      selected_models = present_models
    )

    out_file_manhattan <- file.path(
      output_dir,
      "manhattan",
      paste0(trait_name, "_manhattan.png")
    )
    out_file_qq <- file.path(
      output_dir,
      "qq",
      paste0(trait_name, "_qq.png")
    )

    ggsave(out_file_manhattan, p_manhattan, width = 13, height = 5.5, dpi = 300, bg = "white")
    ggsave(out_file_qq, p_qq, width = 6.8, height = 6.2, dpi = 300, bg = "white")
    cat("  saved:", out_file_manhattan, "\n")
    cat("  saved:", out_file_qq, "\n")
  }

  cat("\nDone.\n")
}

run_all_manhattan <- run_all_gwas_plots

if (sys.nframe() == 0) {
  run_all_gwas_plots()
}

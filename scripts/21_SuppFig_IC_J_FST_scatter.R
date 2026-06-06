################################################################################
### SUPPLEMENTARY FIGURE: GENE-ASSOCIATED FST OVERLAP, INDIAN CHIEF VS JARVIS
################################################################################

library(tidyverse)
library(ggplot2)

################################################################################
### CONFIGURATION
################################################################################

ic_windows_file <- "tables/supplementary/SuppTable_IC_GenWin_windows_Fst_Tajima_pi_genes.csv"
j_windows_file <- "tables/supplementary/SuppTable_J_GenWin_windows_Fst_Tajima_pi_genes.csv"

supp_fig_dir <- "Figs/Supplementary"
supp_table_dir <- "tables/supplementary"

supp_fig_file <- file.path(supp_fig_dir, "SuppFig_IC_J_gene_FST_scatter.png")
supp_table_file <- file.path(supp_table_dir, "SuppTable_IC_J_shared_FST_outlier_genes.csv")

dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_table_dir, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.caption = element_text(size = 12, hjust = 0.5, margin = margin(t = 10)),
    plot.margin = margin(15, 15, 15, 15)
  )

################################################################################
### HELPERS
################################################################################

normalize_flag <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }

  tolower(as.character(x)) %in% c("true", "yes", "1")
}

split_field <- function(x) {
  parts <- strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE)[[1]]
  trimws(parts[nzchar(trimws(parts))])
}

expand_windows_to_genes <- function(df, population_label) {
  expanded_list <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    gene_ids <- split_field(row$GeneIDs)

    if (length(gene_ids) == 0) {
      expanded_list[[i]] <- NULL
      next
    }

    gene_names <- split_field(row$GeneNames)
    gene_biotypes <- split_field(row$GeneBiotypes)

    if (length(gene_names) != length(gene_ids)) {
      gene_names <- gene_ids
    }

    if (length(gene_biotypes) != length(gene_ids)) {
      gene_biotypes <- rep(NA_character_, length(gene_ids))
    }

    expanded_list[[i]] <- tibble(
      population = population_label,
      gene_id = gene_ids,
      gene_name = gene_names,
      gene_biotype = gene_biotypes,
      fst = row$mean_fst,
      fst_outlier = normalize_flag(row$fst_outlier),
      chr = row$chr,
      window_start_bp = row$window_start_bp,
      window_end_bp = row$window_end_bp,
      snp_count = row$snp_count
    )
  }

  bind_rows(expanded_list)
}

collapse_to_gene_max <- function(gene_df) {
  gene_df %>%
    mutate(
      fst = as.numeric(fst),
      window_start_bp = as.numeric(window_start_bp),
      window_end_bp = as.numeric(window_end_bp),
      snp_count = as.numeric(snp_count)
    ) %>%
    arrange(desc(fst), desc(fst_outlier), chr, window_start_bp) %>%
    group_by(gene_id) %>%
    slice(1) %>%
    ungroup()
}

################################################################################
### LOAD DATA
################################################################################

required_files <- c(ic_windows_file, j_windows_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required window files:\n  ",
    paste(missing_files, collapse = "\n  ")
  )
}

ic_windows <- readr::read_csv(ic_windows_file, show_col_types = FALSE)
j_windows <- readr::read_csv(j_windows_file, show_col_types = FALSE)

ic_threshold <- quantile(ic_windows$mean_fst, 0.95, na.rm = TRUE)
j_threshold <- quantile(j_windows$mean_fst, 0.95, na.rm = TRUE)

ic_genes <- ic_windows %>%
  expand_windows_to_genes("Indian Chief") %>%
  collapse_to_gene_max() %>%
  rename_with(~ paste0("ic_", .x), -gene_id)

j_genes <- j_windows %>%
  expand_windows_to_genes("Jarvis") %>%
  collapse_to_gene_max() %>%
  rename_with(~ paste0("j_", .x), -gene_id)

scatter_df <- inner_join(ic_genes, j_genes, by = "gene_id") %>%
  mutate(
    gene_name = coalesce(ic_gene_name, j_gene_name, gene_id),
    gene_biotype = coalesce(ic_gene_biotype, j_gene_biotype),
    category = case_when(
      ic_fst_outlier & j_fst_outlier ~ "Both top 5%",
      ic_fst_outlier & !j_fst_outlier ~ "Indian Chief top 5% only",
      !ic_fst_outlier & j_fst_outlier ~ "Jarvis top 5% only",
      TRUE ~ "Neither"
    ),
    category = factor(
      category,
      levels = c(
        "Both top 5%",
        "Indian Chief top 5% only",
        "Jarvis top 5% only",
        "Neither"
      )
    )
  )

shared_gene_table <- scatter_df %>%
  filter(category == "Both top 5%") %>%
  transmute(
    GeneID = gene_id,
    GeneName = gene_name,
    GeneBiotype = gene_biotype,
    IC_FST = ic_fst,
    Jarvis_FST = j_fst,
    IC_chr = ic_chr,
    IC_window_start_bp = ic_window_start_bp,
    IC_window_end_bp = ic_window_end_bp,
    Jarvis_chr = j_chr,
    Jarvis_window_start_bp = j_window_start_bp,
    Jarvis_window_end_bp = j_window_end_bp,
    IC_snp_count = ic_snp_count,
    Jarvis_snp_count = j_snp_count
  ) %>%
  arrange(desc(IC_FST + Jarvis_FST), GeneID)

################################################################################
### PLOT
################################################################################

category_colors <- c(
  "Both top 5%" = "#7A4FA3",
  "Indian Chief top 5% only" = "#2166AC",
  "Jarvis top 5% only" = "#B35806",
  "Neither" = "#BFBFBF"
)

p <- ggplot(scatter_df, aes(x = ic_fst, y = j_fst, color = category)) +
  geom_point(
    data = scatter_df %>% filter(category != "Both top 5%"),
    size = 1.25,
    alpha = 0.68
  ) +
  geom_point(
    data = scatter_df %>% filter(category == "Both top 5%"),
    size = 2.6,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = ic_threshold,
    linetype = "dashed",
    color = "black",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = j_threshold,
    linetype = "dashed",
    color = "black",
    linewidth = 0.55
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4.8, alpha = 1))) +
  labs(
    title = "Gene-Associated FST: Indian Chief (X) vs Jarvis (Y)",
    x = "Indian Chief FST",
    y = "Jarvis FST",
    caption = paste0(
      "Dashed lines = 95th percentile window cutoffs ",
      "(IC = ", sprintf("%.4f", ic_threshold),
      ", Jarvis = ", sprintf("%.4f", j_threshold), ")"
    )
  ) +
  plot_theme

################################################################################
### WRITE OUTPUTS
################################################################################

readr::write_csv(shared_gene_table, supp_table_file)

ggsave(
  filename = supp_fig_file,
  plot = p,
  width = 11.5,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

cat("\nSaved figure:\n  ", supp_fig_file, "\n", sep = "")
cat("Saved table:\n  ", supp_table_file, "\n", sep = "")
cat("Shared top-5% genes:\n  ", nrow(shared_gene_table), "\n", sep = "")

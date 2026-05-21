################################################################################
### FIGURE 3: STRUCTURE, PCA, AND SNP SHARING IN INDIAN CHIEF AND JARVIS
################################################################################

library(SNPRelate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
### CONFIGURATION
################################################################################

data_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis"
vcf_file <- file.path(
  data_dir,
  "Variants.Filtered.PASS.renamed_chr_only.final.filtered.vcf.gz"
)
gds_file <- file.path(
  data_dir,
  "Variants.Filtered.PASS.renamed_chr_only.final.filtered.gds"
)

output_dir <- "Figs/main"
output_file <- file.path(output_dir, "Fig3.png")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

num_threads <- 4
ld_threshold <- 0.2
maf_threshold <- 0.005
missing_rate_threshold <- 1
chunk_size <- 10000
target_alt_freq_threshold <- 0.5

group_colors <- c(
  "Indian Chief G0" = "#92C5DE",
  "Indian Chief G14" = "#2166AC",
  "Jarvis G0" = "#FDB863",
  "Jarvis G14" = "#B35806"
)

summary_fill_colors <- c(
  "Indian Chief" = "#4A90C2",
  "Jarvis" = "#D9822B",
  "Generation 0" = "#B7C9D6",
  "Generation 14" = "#7B8FA1"
)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 15, 15, 15)
  )

################################################################################
### HELPERS
################################################################################

sample_metadata <- function(sample_ids) {
  tibble(sample_id = sample_ids) %>%
    mutate(
      line = case_when(
        grepl("^I", sample_id) ~ "Indian Chief",
        grepl("^J", sample_id) ~ "Jarvis",
        TRUE ~ "Other"
      ),
      generation = case_when(
        grepl("^I01|^J01", sample_id) ~ "G0",
        grepl("^I14|^J14", sample_id) ~ "G14",
        TRUE ~ "Other"
      ),
      group = case_when(
        line == "Indian Chief" & generation == "G0" ~ "Indian Chief G0",
        line == "Indian Chief" & generation == "G14" ~ "Indian Chief G14",
        line == "Jarvis" & generation == "G0" ~ "Jarvis G0",
        line == "Jarvis" & generation == "G14" ~ "Jarvis G14",
        TRUE ~ "Other"
      )
    ) %>%
    filter(group != "Other") %>%
    mutate(group = factor(group, levels = names(group_colors)))
}

calc_shared_prop <- function(genofile, sample_ids, group_indices, snp_ids, chunk_size, target_alt_freq_threshold) {
  n_samples <- length(sample_ids)
  n_snps <- length(snp_ids)

  shared_cross_gen <- numeric(n_samples)
  shared_cross_pop <- numeric(n_samples)
  alt_counts <- numeric(n_samples)

  cross_gen_target <- c(
    "Indian Chief G0" = "Indian Chief G14",
    "Indian Chief G14" = "Indian Chief G0",
    "Jarvis G0" = "Jarvis G14",
    "Jarvis G14" = "Jarvis G0"
  )

  cross_pop_target <- c(
    "Indian Chief G0" = "Jarvis G0",
    "Indian Chief G14" = "Jarvis G14",
    "Jarvis G0" = "Indian Chief G0",
    "Jarvis G14" = "Indian Chief G14"
  )

  chunk_starts <- seq(1, n_snps, by = chunk_size)

  for (start_idx in chunk_starts) {
    end_idx <- min(start_idx + chunk_size - 1, n_snps)
    chunk_snp_ids <- snp_ids[start_idx:end_idx]

    geno <- snpgdsGetGeno(
      genofile,
      sample.id = sample_ids,
      snp.id = chunk_snp_ids,
      with.id = FALSE,
      verbose = FALSE
    )

    alt_present <- !is.na(geno) & geno > 0
    alt_counts <- alt_counts + rowSums(alt_present)

    group_alt_presence <- lapply(group_indices, function(idx) {
      colMeans(alt_present[idx, , drop = FALSE], na.rm = TRUE) >= target_alt_freq_threshold
    })

    for (group_name in names(group_indices)) {
      idx <- group_indices[[group_name]]

      gen_target <- group_alt_presence[[cross_gen_target[[group_name]]]]
      pop_target <- group_alt_presence[[cross_pop_target[[group_name]]]]

      shared_cross_gen[idx] <- shared_cross_gen[idx] +
        rowSums(alt_present[idx, , drop = FALSE] &
          matrix(gen_target, nrow = length(idx), ncol = length(gen_target), byrow = TRUE))

      shared_cross_pop[idx] <- shared_cross_pop[idx] +
        rowSums(alt_present[idx, , drop = FALSE] &
          matrix(pop_target, nrow = length(idx), ncol = length(pop_target), byrow = TRUE))
    }
  }

  tibble(
    sample_id = sample_ids,
    shared_cross_generation = ifelse(alt_counts > 0, shared_cross_gen / alt_counts, NA_real_),
    shared_cross_population = ifelse(alt_counts > 0, shared_cross_pop / alt_counts, NA_real_)
  )
}

################################################################################
### BUILD OR LOAD GDS
################################################################################

if (!file.exists(gds_file)) {
  cat("Converting VCF to GDS:\n")
  cat("  ", vcf_file, "\n", sep = "")
  snpgdsVCF2GDS(
    vcf.fn = vcf_file,
    out.fn = gds_file,
    method = "copy.num.of.ref",
    verbose = TRUE
  )
}

genofile <- snpgdsOpen(gds_file)
on.exit(snpgdsClose(genofile), add = TRUE)

################################################################################
### LD PRUNING, IBS, PCA, AND ADMIXTURE-LIKE PROPORTIONS
################################################################################

set.seed(1000)

snpset <- snpgdsLDpruning(
  genofile,
  ld.threshold = ld_threshold,
  maf = maf_threshold,
  missing.rate = missing_rate_threshold,
  autosome.only = FALSE,
  num.thread = num_threads,
  verbose = TRUE
)

snp_ids <- unlist(snpset, use.names = FALSE)

if (length(snp_ids) == 0) {
  stop("LD pruning returned 0 SNPs. Please adjust the pruning filters.")
}

sample_ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
meta_df <- sample_metadata(sample_ids)
sample_ids <- meta_df$sample_id

group_sample_ids <- split(meta_df$sample_id, meta_df$group)
group_indices <- lapply(group_sample_ids, function(ids) match(ids, sample_ids))

ibs <- snpgdsIBS(
  genofile,
  sample.id = sample_ids,
  snp.id = snp_ids,
  autosome.only = FALSE,
  num.thread = num_threads,
  verbose = TRUE
)

hc <- hclust(as.dist(1 - ibs$ibs), method = "average")
sample_order <- ibs$sample.id[hc$order]

pca <- snpgdsPCA(
  genofile,
  sample.id = sample_ids,
  snp.id = snp_ids,
  num.thread = num_threads,
  autosome.only = FALSE,
  verbose = TRUE
)

admix_prop <- snpgdsAdmixProp(
  pca,
  groups = group_sample_ids,
  bound = TRUE
)

################################################################################
### SUMMARY SNP SHARING
################################################################################

sharing_df <- calc_shared_prop(
  genofile = genofile,
  sample_ids = sample_ids,
  group_indices = group_indices,
  snp_ids = snp_ids,
  chunk_size = chunk_size,
  target_alt_freq_threshold = target_alt_freq_threshold
)

plot_df <- meta_df %>%
  left_join(
    tibble(
      sample_id = pca$sample.id,
      PC1 = pca$eigenvect[, 1],
      PC2 = pca$eigenvect[, 2]
    ),
    by = "sample_id"
  ) %>%
  left_join(
    sharing_df,
    by = "sample_id"
  )

################################################################################
### PANEL A: STRUCTURE
################################################################################

structure_df <- as.data.frame(admix_prop) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(plot_df %>% select(sample_id, group, PC1), by = "sample_id") %>%
  mutate(sample_id = factor(sample_id, levels = sample_order)) %>%
  arrange(sample_id)

structure_long <- structure_df %>%
  pivot_longer(
    cols = all_of(names(group_colors)),
    names_to = "ancestry_group",
    values_to = "prop"
  ) %>%
  mutate(
    ancestry_group = factor(ancestry_group, levels = names(group_colors)),
    prop = pmin(pmax(prop, 0), 1)
  ) %>%
  filter(!is.na(prop))

structure_plot <- ggplot(structure_long, aes(x = sample_id, y = prop, fill = ancestry_group)) +
  geom_col(width = 1) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.25)) +
  labs(
    tag = "A",
    title = "Population Structure",
    x = NULL,
    y = "Ancestry Proportion"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  plot_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

################################################################################
### PANEL B: PCA
################################################################################

pca_plot <- ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.9, size = 3.2, stroke = 0) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  labs(
    tag = "B",
    title = "Principal Component Analysis",
    x = paste0("PC1 (", round(pca$varprop[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca$varprop[2] * 100, 2), "%)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  coord_fixed(ratio = 1) +
  plot_theme +
  theme(legend.position = "none")

################################################################################
### PANEL C: SHARED ACROSS GENERATIONS WITHIN POPULATION
################################################################################

cross_gen_df <- plot_df %>%
  mutate(population_label = factor(line, levels = c("Indian Chief", "Jarvis")))

p_cross_gen <- ggplot(cross_gen_df, aes(x = population_label, y = shared_cross_generation, fill = population_label)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = group), width = 0.12, size = 1.2, alpha = 0.65) +
  scale_fill_manual(values = summary_fill_colors[c("Indian Chief", "Jarvis")]) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    tag = "C",
    title = "Shared SNPs Across Generations",
    x = NULL,
    y = "Shared SNP Proportion"
  ) +
  plot_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

################################################################################
### PANEL D: SHARED WITH OPPOSITE POPULATION WITHIN GENERATION
################################################################################

cross_pop_df <- plot_df %>%
  mutate(generation_label = factor(
    ifelse(generation == "G0", "Generation 0", "Generation 14"),
    levels = c("Generation 0", "Generation 14")
  ))

p_cross_pop <- ggplot(cross_pop_df, aes(x = generation_label, y = shared_cross_population, fill = generation_label)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = group), width = 0.12, size = 1.2, alpha = 0.65) +
  scale_fill_manual(values = summary_fill_colors[c("Generation 0", "Generation 14")]) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    tag = "D",
    title = "Shared SNPs Across Populations",
    x = NULL,
    y = "Shared SNP Proportion"
  ) +
  plot_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

################################################################################
### COMBINE AND SAVE
################################################################################

layout_design <- "
AAA
BCD
"

final_plot <- structure_plot + pca_plot + p_cross_gen + p_cross_pop +
  plot_layout(design = layout_design, heights = c(1.15, 1))

ggsave(
  filename = output_file,
  plot = final_plot,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("\nSaved Figure 3 to:\n")
cat("  ", output_file, "\n", sep = "")

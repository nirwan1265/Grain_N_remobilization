#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- effective population size of Indian Chief and Jarvis
##
## Genome-wide mean F_ST between Cycle 0 and Cycle 14 expected under neutral
## Wright-Fisher drift for 14 generations at a grid of N_e, against the observed
## value in each population. N_e is read off by linear interpolation.
##
## The simulated values come from the calibration step of 04_drift_sim.py
## (thesis "Claude outputs"; run on the HPC from the ANGSD .mafs frequencies of
## each population-cycle group via 03_prep_freqs.sh): every site starts at its
## observed Cycle 0 frequency, drifts 14 generations at N_e, and is resampled at
## both ends with the number of individuals actually sequenced at that site;
## per-window Hudson F_ST (250-kb windows, ratio of sums) is averaged over windows.
## The dotted line marks N_e = 20, the value expected from the breeding design
## (20 selected parents intermated per cycle; Moll and Stuber 1971) when all
## parents contribute equally.
##
## Output: final/supp_figs/SuppFig_IC_J_Ne_drift.{png,pdf}
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
OUT  <- gv("GN_SUPPDIR", file.path(REPO, "final/supp_figs"))
dir.create(OUT, FALSE, TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = c(0.01, 0.99),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black"), panel.grid = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
        legend.title = element_blank(), legend.text = element_text(size = 13),
        plot.margin = margin(15, 15, 15, 15))

## calibration grid (simulated genome-wide mean F_ST; 04_drift_sim.py job output)
SIM <- data.table(Ne = c(20, 30, 50, 75, 100, 150, 250),
                  fst = c(0.14010, 0.09538, 0.05585, 0.03479, 0.02389, 0.01269, 0.00363))
OBS <- data.table(population = c("Indian Chief", "Jarvis"), fst = c(0.09479, 0.11371))
OBS[, Ne := approx(SIM$fst, SIM$Ne, xout = fst)$y]          # fst is monotone in Ne
print(OBS)
COL <- c("Indian Chief" = "#2166AC", "Jarvis" = "#B35806")

p <- ggplot(SIM, aes(Ne, fst)) +
  geom_vline(xintercept = 20, linetype = "dotted", colour = "grey30", linewidth = 0.7) +
  annotate("text", x = 22, y = 0.004, hjust = 0, vjust = 0, size = 4.6, colour = "grey30",
           label = "breeding design:\n20 parents per cycle") +
  geom_hline(data = OBS, aes(yintercept = fst, colour = population), linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = "grey25", linewidth = 0.9) +
  geom_point(shape = 21, fill = "white", colour = "grey25", size = 3.5, stroke = 1.2) +
  geom_point(data = OBS, aes(Ne, fst, colour = population), size = 5) +
  annotate("text", x = OBS$Ne + 12, y = OBS$fst + c(-0.007, 0.007), hjust = 0, size = 5,
           colour = COL[OBS$population],
           parse = TRUE,
           label = sprintf("'%s: observed %.3f, '*italic(N)[e] %%~~%% %.0f", OBS$population, OBS$fst, OBS$Ne)) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = c(20, 50, 100, 150, 200, 250), limits = c(15, 255)) +
  scale_y_continuous(limits = c(0, 0.15)) +
  labs(x = expression(bold(paste("Effective population size, ", italic(N)[e]))),
       y = expression(bold(paste("Mean ", italic(F)[ST], ", Cycle 0 vs Cycle 14")))) +
  plot_theme
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("SuppFig_IC_J_Ne_drift.", ext)), p, width = 10, height = 7, dpi = 300,
         bg = "white")

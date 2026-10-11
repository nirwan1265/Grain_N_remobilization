#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- soil-N association specificity against environmental
## control traits. Replaces the earlier genotype-correlation heatmap: this is
## built from the control GWAS themselves, not from dosage correlations.
##
## A  overlap of each trait's FDR-significant SNP set with the union of the four
##    geographic traits. Latitude and longitude are the positive controls; soil
##    N is the trait of interest and sits lowest.
## B  soil N against every other trait individually, observed overlap vs the
##    hypergeometric expectation.
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
D <- Sys.getenv("GN_ENV_DIR", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/soil_nitrogen_gwas/env_landrace")
OUT <- Sys.getenv("GN_OUT", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/supp_figs/SuppFig_soilN_specificity.png")
M <- 495755; MODELS <- c("MLM","BLINK","FarmCPU")
GEO <- c("latitude","abs_latitude","longitude","elevation")
LAB <- c(N_0_5="Soil N (0-5 cm)", N_5_15="Soil N (5-15 cm)", aridity="Aridity",
         pH="Soil pH", CEC="CEC", elevation="Elevation", latitude="Latitude",
         abs_latitude="|Latitude|", longitude="Longitude")

rdx <- function(tr,md) fread(file.path(D,sprintf("%s__%s.csv.gz",tr,md)),showProgress=FALSE)[
  fdr_sig==TRUE, .(trait=tr, SNP)]
traits <- unique(sub("__.*","",list.files(D,pattern="\\.csv\\.gz$")))
sig <- unique(rbindlist(lapply(traits,function(t) rbindlist(lapply(MODELS,function(m) rdx(t,m))))))
sets <- split(sig$SNP, sig$trait)

## ---- A: overlap with the geographic union, leave-self-out ------------------
A <- rbindlist(lapply(traits, function(t){
  s <- sets[[t]]; g <- unique(unlist(sets[setdiff(GEO,t)])); o <- length(intersect(s,g))
  data.table(trait=t, n=length(s), ov=o, pct=100*o/length(s),
             exp_pct=100*length(g)/M)}))
A[, lab := LAB[trait]][, grp := fifelse(trait=="N_0_5","focal",
                              fifelse(trait %in% GEO,"geography","other soil"))]
setorder(A, pct); A[, lab := factor(lab, levels=lab)]

pA <- ggplot(A, aes(pct, lab, fill=grp)) +
  geom_col(width=.7) +
  geom_point(aes(x=exp_pct), shape=124, size=4, colour="grey25") +
  geom_text(aes(label=sprintf("%.1f%%  (%d/%d)", pct, ov, n)), hjust=-0.08,
            size=3.3, fontface="bold") +
  scale_fill_manual(values=c(focal="#B35806", geography="#2166AC", `other soil`="#9E9E9E"),
                    name=NULL, breaks=c("focal","geography","other soil"),
                    labels=c("Soil N (focal)","Geographic trait","Other soil trait")) +
  scale_x_continuous(expand=expansion(mult=c(0,.30))) +
  labs(x="SNPs also significant for a geographic trait (%)", y=NULL, tag="A") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(colour="grey88", linewidth=.3),
        axis.line.x=element_line(colour="black"),
        axis.text=element_text(size=10, colour="black", face="bold"),
        axis.title=element_text(size=11, face="bold"),
        legend.position="top", legend.text=element_text(size=9),
        plot.tag=element_text(size=16, face="bold"))

## ---- B: soil N vs each trait, observed vs expected --------------------------
N <- sets[["N_0_5"]]
B <- rbindlist(lapply(setdiff(traits,"N_0_5"), function(t){
  o <- length(intersect(N, sets[[t]])); e <- length(N)*length(sets[[t]])/M
  data.table(trait=t, ov=o, exp=e,
    P=phyper(o-1, length(sets[[t]]), M-length(sets[[t]]), length(N), lower.tail=FALSE))}))
B[, lab := LAB[trait]]; setorder(B, ov); B[, lab := factor(lab, levels=lab)]
B[, star := fifelse(P<1e-20,"***", fifelse(P<1e-5,"**", fifelse(P<0.05,"*","n.s.")))]

pB <- ggplot(B, aes(ov, lab)) +
  geom_segment(aes(x=exp, xend=ov, yend=lab), colour="grey70", linewidth=.9) +
  geom_point(aes(x=exp), size=2.6, colour="grey45") +
  geom_point(size=3.8, colour="#B35806") +
  geom_text(aes(label=sprintf("%d  %s", ov, star)), hjust=-0.25, size=3.3, fontface="bold") +
  scale_x_continuous(expand=expansion(mult=c(.04,.26))) +
  labs(x="Soil N SNPs shared with trait (observed vs expected)", y=NULL, tag="B") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(colour="grey88", linewidth=.3),
        axis.line.x=element_line(colour="black"),
        axis.text=element_text(size=10, colour="black", face="bold"),
        axis.title=element_text(size=11, face="bold"),
        plot.tag=element_text(size=16, face="bold"))

ggsave(OUT, pA | pB, width=12.5, height=4.6, dpi=300, bg="white")
cat("Saved:", OUT, "\n"); print(A[, .(lab, n, ov, pct=round(pct,1))]); print(B[, .(lab, ov, exp=round(exp,2), star)])

#!/bin/bash
################################################################################
# Indian Chief / Jarvis: genotype-likelihood structure and PCA (Fig. 4A-B).
#
#   1. vcf_to_beagle_GL.py : PL field of the filtered GATK VCF -> Beagle GL file
#      (all 379 samples, same sites; >=50% samples with reads, MAF >= 0.05,
#      one site per 10 kb)                                   -> 82,061 sites
#   2. PCAngsd on the Beagle file                             -> IC_J_pcangsd.cov
#   3. NGSadmix K = 2-5, 3 seeds each                         -> *.qopt, *.log
#   4. Fig4_IC_J_structure_PCA_Fst.R draws the figure (also computes panel C,
#      between-population F_ST, from the winsfs 2D SFS in data/Between_pop_Fst_IC_J).
#
# Tools: pcangsd (pip install pcangsd); NGSadmix (github.com/aalbrechtsen/NGSadmix,
#        g++ -O3 -o NGSadmix NGSadmix.cpp -lpthread -lz). Runs in ~20 min on 2 cores.
################################################################################
set -euo pipefail
REPO=${REPO:-/Users/nirwantandukar/Documents/Github/Grain_N_remobilization}
VCF=${VCF:-/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/Variants.Filtered.PASS.renamed_chr_only.final.filtered.vcf.gz}
OUT=${OUT:-$REPO/final/results/ic_j_structure}
NGSADMIX=${NGSADMIX:-NGSadmix}
THREADS=${THREADS:-2}
mkdir -p "$OUT"
BGL=$OUT/IC_J_all_thin10kb.beagle.gz
python3 "$REPO/final/scripts/vcf_to_beagle_GL.py" "$VCF" "$BGL" "$OUT/IC_J_all_sample_order.txt"
pcangsd -b "$BGL" -t $THREADS -o "$OUT/IC_J_pcangsd"
for K in 2 3 4 5; do for S in 1 2 3; do
  $NGSADMIX -likes "$BGL" -K $K -seed $S -P $THREADS -minMaf 0.05 -o "$OUT/IC_J_K${K}_s${S}"
done; done
Rscript "$REPO/final/scripts/Fig4_IC_J_structure_PCA_Fst.R"

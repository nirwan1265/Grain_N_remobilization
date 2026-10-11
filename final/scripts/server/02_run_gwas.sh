#!/bin/bash
#BSUB -J gwas
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=64GB]"
#BSUB -W 2880
#BSUB -o logs/gwas.%J.out
#BSUB -e logs/gwas.%J.err

############################################################
##  CHANGE THIS LINE, then:   bsub < 02_run_gwas.sh
##
##  N_0_5  latitude  abs_latitude  longitude  elevation
##  N_5-15  pH  CEC  aridity
TRAIT=N_0_5
############################################################

BASE=/share/maize/ntanduk/GrainN/GWAS_env
ENV=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/env_GAPIT
RSCRIPT="$ENV/bin/Rscript"
[ -x "$RSCRIPT" ] || RSCRIPT=Rscript

mkdir -p "$BASE/logs"
cd "$BASE" || exit 1

export GN_GENO="$BASE/geno"
export GN_PHENO="$BASE/control_phenotypes_maize_landrace.csv"
export GN_PCA="$BASE/PCA_maize_romerro.csv"
export GN_OUTDIR="$BASE/out"

echo "=== $TRAIT on $(hostname) at $(date) ==="
"$RSCRIPT" 01_GWAS_models.R "$TRAIT"
echo "=== done $(date) ==="

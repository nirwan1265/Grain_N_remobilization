#!/bin/bash
#BSUB -J ld_phg
#BSUB -n 2
#BSUB -W 2:00
#BSUB -R "rusage[mem=16GB]"
#BSUB -o ld_phg.%J.out
#BSUB -e ld_phg.%J.err
################################################################################
## r^2 between the lead soil-N SNP and every PHG SNP within +/-1 Mb, for
##   NMX   (TATC, Zm00001eb399680) chr9  lead 151,733,953 (inside the gene)
##   ZmFd9 (Zm00001eb421930)       chr10 lead 116,706,616 (300 bp past the gene end)
## Leads are the strongest SNPs at each gene in GLM_SoilN_maize.txt (Pos column,
## B73 v5). SNP names in that file are wrong, so everything here keys on position:
## variant IDs are reset to CHROM:POS before PLINK.
##
## Set VCF to the PHG genotype VCF used for the soil-N GLM (bgzipped + indexed).
## If the PHG genotypes are in TASSEL HDF5 or HapMap, convert first, e.g.
##   run_pipeline.pl -Xmx32g -h5 phg.h5 -export phg -exportType VCF
##   bgzip phg.vcf && tabix -p vcf phg.vcf.gz
## KEEP (optional) = one sample ID per line, the accessions in the GLM.
##
## Output (copy these to final/results/soil_nitrogen_gwas/LD_PHG/ on the Mac):
##   NMX_r2.ld, ZmFd9_r2.ld   (PLINK --r2 table, lead SNP vs all SNPs)
################################################################################
set -euo pipefail
VCF=${VCF:-/path/to/PHG_genotypes.vcf.gz}   # <-- set this
KEEP=${KEEP:-}                              # optional sample list
FLANK=1000000
module load bcftools 2>/dev/null || true
module load plink/1.90 2>/dev/null || module load plink 2>/dev/null || true

# chromosome naming in the VCF: "9" or "chr9"
PFX=""; bcftools index -s "$VCF" | cut -f1 | grep -q '^chr' && PFX="chr"

while read NAME CHR LEAD; do
  S=$((LEAD-FLANK)); E=$((LEAD+FLANK)); C="${PFX}${CHR}"
  echo "== $NAME $C:$S-$E lead $LEAD"
  SAMP=(); [ -n "$KEEP" ] && SAMP=(-S "$KEEP" --force-samples)
  bcftools view -r "$C:$S-$E" "${SAMP[@]}" -m2 -M2 -v snps "$VCF" -Ou \
    | bcftools annotate --set-id '%CHROM:%POS' -Oz -o ${NAME}_region.vcf.gz
  plink --vcf ${NAME}_region.vcf.gz --double-id --allow-extra-chr \
        --make-bed --out ${NAME}_region
  if ! grep -qP "\t${C}:${LEAD}\t" ${NAME}_region.bim; then
    echo "lead SNP $C:$LEAD not in the VCF region, check positions/assembly"; exit 1
  fi
  plink --bfile ${NAME}_region --allow-extra-chr \
        --r2 --ld-snp ${C}:${LEAD} --ld-window 999999 --ld-window-kb 1100 --ld-window-r2 0 \
        --out ${NAME}_r2
  echo "$NAME: $(($(wc -l < ${NAME}_r2.ld)-1)) SNPs with r2"
done <<LOCI
NMX 9 151733953
ZmFd9 10 116706616
LOCI

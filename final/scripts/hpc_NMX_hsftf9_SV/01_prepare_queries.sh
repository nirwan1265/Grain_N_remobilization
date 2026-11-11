#!/bin/bash
################################################################################
## Step 1 (login node, ~5 min): build the query sequences from B73 v5
##  - genomic sequence of every B73 v5 gene in chr9:151.70-152.02 Mb
##    (NMX = Zm00001eb399680 ... hsftf9 = Zm00001eb399800 and everything between)
##  - 301 bp around four key SNPs (NMX lead, two 152 Mb SNPs, the HSF9 SNP)
##  - the pan-gene table rows for these genes
##  - the list of maize genomes to search (B73 v5 + NAM founders)
## Usage: bash 01_prepare_queries.sh
################################################################################
set -euo pipefail
DATA=/rsstu/users/r/rrellan/sara/maizegdb/data
OUT=${OUT:-/share/maize/ntanduk/NMX_hsftf9_SV}
mkdir -p $OUT && cd $OUT
eval "$(conda shell.bash hook)"; conda activate seqanal

REF=$(ls $DATA/All_assembly_sequence/Zm-B73-REFERENCE-NAM-5.0*.fa.gz | head -1)
GFF=$(ls $DATA/All_gene_model_GFF/Zm-B73-REFERENCE-NAM-5.0*Zm00001eb*.gff3.gz | head -1)
echo "B73 v5 assembly: $REF"; echo "B73 v5 GFF: $GFF"

# B73 v5 chr9 only, indexed
zcat $REF | awk '/^>/{p=($1==">chr9"||$1==">9")} p' > b73v5_chr9.fa
samtools faidx b73v5_chr9.fa
CHR=$(cut -f1 b73v5_chr9.fa.fai | head -1)

# genes in the window
zcat $GFF | awk -F'\t' -v c=$CHR '$1==c && $3=="gene" && $5>151700000 && $4<152020000' \
  | awk -F'\t' '{match($9,/ID=[^;]+/); print $1"\t"$4"\t"$5"\t"$7"\t"substr($9,RSTART+3,RLENGTH-3)}' > window_genes_b73v5.tsv
echo "genes in window:"; cat window_genes_b73v5.tsv
: > queries.fa
while read c s e str id; do
  samtools faidx b73v5_chr9.fa ${c}:${s}-${e} | sed "1s/.*/>${id}|${s}|${e}/" >> queries.fa
done < window_genes_b73v5.tsv
for P in 151742636 151983481 151984305 151992800; do
  samtools faidx b73v5_chr9.fa ${CHR}:$((P-150))-$((P+150)) | sed "1s/.*/>snp_${P}|$((P-150))|$((P+150))/" >> queries.fa
done
grep -c '^>' queries.fa

# pan-gene rows (format-independent grep)
zgrep -rh -E 'Zm00001eb399680|Zm00001eb399780|Zm00001eb399800' $DATA/Pan-genes 2>/dev/null > pangene_rows.txt || true
echo "pan-gene rows: $(wc -l < pangene_rows.txt)"

# genomes to search: B73 v5 + all NAM founder assemblies
ls $DATA/All_assembly_sequence/Zm-*-REFERENCE-NAM-*.fa.gz > genomes.txt
echo "genomes: $(wc -l < genomes.txt)"; cat genomes.txt
echo "next: bsub -J \"nam[1-$(wc -l < genomes.txt)]\" < 02_map_genomes.lsf"

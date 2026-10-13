#!/usr/bin/env python3
"""
LD-based candidate-gene assignment for GWAS hits (replaces the fixed +/-25 kb window).

For every significant SNP ("lead"), r^2 is computed between the lead and every SNP
within +/-WINDOW bp in the GWAS panel itself (PLINK .bed, only the accessions used in
the GWAS). Partners with r^2 >= R2 define the lead's LD block = [min, max] position of
the lead and its partners. Because markers are discrete, the causal variant can lie
anywhere up to the first flanking marker that is NOT in LD, so the LD interval used for
gene assignment runs from the nearest unlinked marker left of the block to the nearest
unlinked marker right of it (exclusive of those markers). A gene is assigned to the
lead if it overlaps this interval. The core block (partners only) and the +/-25 kb
assignment are written alongside for comparison.

r^2 = squared Pearson correlation of allele dosages over accessions genotyped at both
SNPs (missing genotypes dropped pairwise). Genes: B73 RefGen_v5 GFF3 gene models.
Coordinates of the .bed and the GFF must be the same assembly.

usage: LD_gene_assignment.py <bed_prefix> <keep_ids.txt|ALL> <sig_snps.tsv chr pos [id]>
                             <gff3> <out_prefix> [R2=0.4] [WINDOW=250000]
outputs: <out>_ld_blocks.tsv  (one row per lead)
         <out>_genes.tsv      (lead x gene, with in_ld_block and in_25kb flags)
"""
import sys, numpy as np
bed, keepf, sigf, gff, out = sys.argv[1:6]
R2 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.4
WIN = int(sys.argv[7]) if len(sys.argv) > 7 else 250_000
FLANK = 25_000

fam = [l.split()[1] for l in open(bed + ".fam")]
n = len(fam)
keep = np.arange(n) if keepf == "ALL" else \
    np.array([i for i, t in enumerate(fam) if t in set(open(keepf).read().split())])
bim = [l.split() for l in open(bed + ".bim")]
chrom = np.array([b[0] for b in bim]); pos = np.array([int(b[3]) for b in bim])
nb = (n + 3) // 4
raw = np.memmap(bed + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), nb)
lut = np.array([2, np.nan, 1, 0])                     # 00 hom A1, 01 missing, 10 het, 11 hom A2
def geno(idx):
    b = raw[idx]
    c = np.stack([(b >> (2 * k)) & 3 for k in range(4)], axis=-1).reshape(len(idx), nb * 4)[:, :n]
    return lut[c][:, keep]

genes = []
for l in open(gff):
    if "\tgene\t" not in l: continue
    f = l.split("\t")
    c = f[0].replace("chr", "")
    if not c.isdigit(): continue
    genes.append((c, int(f[3]), int(f[4]), f[8].split(";")[0].replace("ID=", "")))
G = {}
for c, s, e, g in genes: G.setdefault(c, []).append((s, e, g))
def overlap(c, a, b): return [g for s, e, g in G.get(c, []) if s <= b and e >= a]

sig = [l.split() for l in open(sigf) if l.strip() and not l.startswith("chr")]
fo = open(out + "_ld_blocks.tsv", "w"); fg = open(out + "_genes.tsv", "w")
fo.write("lead\tchr\tpos\tn_window_snps\tn_partners_r2\tmax_r2\tblock_start\tblock_end\tblock_kb\tinterval_start\tinterval_end\tinterval_kb\tn_genes_ld\tn_genes_25kb\n")
fg.write("lead\tchr\tpos\tGeneID\tin_ld_block\tin_25kb\n")
by_chr = {c: np.where(chrom == c)[0] for c in np.unique(chrom)}
for row in sig:
    c, p = row[0].replace("chr", ""), int(row[1]); lid = row[2] if len(row) > 2 else f"{c}_{p}"
    idx = by_chr[c]; ii = idx[(pos[idx] >= p - WIN) & (pos[idx] <= p + WIN)]
    j = np.where(pos[ii] == p)[0]
    if len(j) == 0:
        fo.write(f"{lid}\t{c}\t{p}\t{len(ii)}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"); continue
    X = geno(ii); x = X[j[0]]
    ok = ~np.isnan(X) & ~np.isnan(x)
    xm = np.where(ok, x, 0); Xm = np.where(ok, X, 0); m = ok.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, mX = xm.sum(1) / m, Xm.sum(1) / m
        cov = (xm * Xm).sum(1) / m - mx * mX
        vx = (xm * xm).sum(1) / m - mx ** 2; vX = (Xm * Xm).sum(1) / m - mX ** 2
        r2 = cov ** 2 / (vx * vX)
    r2[j[0]] = np.nan
    part = np.where(r2 >= R2)[0]
    ps = np.r_[pos[ii][part], p]; a, b = ps.min(), ps.max()
    wp = pos[ii]
    left = wp[wp < a]; right = wp[wp > b]
    ia = (left.max() + 1) if len(left) else max(1, p - WIN)
    ib = (right.min() - 1) if len(right) else p + WIN
    gl = overlap(c, ia, ib); g25 = overlap(c, p - FLANK, p + FLANK)
    fo.write(f"{lid}\t{c}\t{p}\t{len(ii)}\t{len(part)}\t{np.nanmax(r2) if np.isfinite(r2).any() else 'NA'}\t"
             f"{a}\t{b}\t{(b - a) / 1000:.1f}\t{ia}\t{ib}\t{(ib - ia) / 1000:.1f}\t{len(gl)}\t{len(g25)}\n")
    for g in sorted(set(gl) | set(g25)):
        fg.write(f"{lid}\t{c}\t{p}\t{g}\t{g in gl}\t{g in g25}\n")
print(f"{len(sig)} leads | samples {len(keep)} | r2 >= {R2} within +/-{WIN // 1000} kb -> {out}_*")

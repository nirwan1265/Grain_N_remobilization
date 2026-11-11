#!/usr/bin/env python3
"""Step 3 (login node, seconds): summarise the PAFs.
Per genome: where NMX (Zm00001eb399680) and hsftf9 (Zm00001eb399800) land, whether on the
same chromosome/contig, the distance between them, how many of the B73 genes between them
are found between them, and how many places each SNP flank aligns to (>= 90% of its length).
Writes NMX_hsftf9_per_genome.tsv and snp_flank_hits.tsv in OUT.
usage: python3 03_summarize.py [OUT]
"""
import sys, os, glob, collections
OUT = sys.argv[1] if len(sys.argv) > 1 else "/share/maize/ntanduk/NMX_hsftf9_SV"
NMX, HSF = "Zm00001eb399680", "Zm00001eb399800"
genes = [l.split()[4] for l in open(os.path.join(OUT, "window_genes_b73v5.tsv"))]
between = genes[genes.index(NMX) + 1: genes.index(HSF)] if NMX in genes and HSF in genes else []
rows, snp = [], []
for f in sorted(glob.glob(os.path.join(OUT, "paf", "*.paf"))):
    g = os.path.basename(f)[:-4]; hits = collections.defaultdict(list)
    for l in open(f):
        x = l.split("\t"); q = x[0].split("|")[0]; qlen = int(x[1])
        cov = (int(x[3]) - int(x[2])) / qlen; ident = int(x[9]) / max(int(x[10]), 1)
        hits[q].append((int(x[9]), cov, ident, x[5], int(x[7]), int(x[8]), x[4]))
    best = {q: max(h) for q, h in hits.items()}
    def loc(q): b = best.get(q); return (b[3], b[4], round(b[1], 2), round(b[2], 3)) if b else (None, None, 0, 0)
    n, h = loc(NMX), loc(HSF)
    same = n[0] is not None and n[0] == h[0]
    dist = abs(h[1] - n[1]) if same else None
    lo, hi = (min(n[1], h[1]), max(n[1], h[1])) if same else (None, None)
    found_between = sum(1 for q in between if q in best and best[q][3] == n[0] and lo <= best[q][4] <= hi and best[q][1] >= 0.8) if same else None
    rows.append([g, n[0], n[1], n[2], n[3], h[0], h[1], h[2], h[3], same, dist, f"{found_between}/{len(between)}" if same else "NA"])
    for q in [k for k in hits if k.startswith("snp_")]:
        good = [x for x in hits[q] if x[1] >= 0.9 and x[2] >= 0.9]
        snp.append([g, q, len(good), ";".join(f"{x[3]}:{x[4]}({x[2]:.3f})" for x in sorted(good, reverse=True)[:5])])
hdr = ["genome", "NMX_contig", "NMX_pos", "NMX_cov", "NMX_ident", "HSF9_contig", "HSF9_pos", "HSF9_cov", "HSF9_ident",
       "same_contig", "distance_bp", "B73_genes_between_found"]
with open(os.path.join(OUT, "NMX_hsftf9_per_genome.tsv"), "w") as o:
    o.write("\t".join(hdr) + "\n"); [o.write("\t".join(map(str, r)) + "\n") for r in rows]
with open(os.path.join(OUT, "snp_flank_hits.tsv"), "w") as o:
    o.write("genome\tsnp\tn_hits_ge90pct\ttop_hits\n"); [o.write("\t".join(map(str, r)) + "\n") for r in snp]
print("\t".join(hdr)); [print("\t".join(map(str, r))) for r in rows]
print("\nSNP flanks with more than one >=90% hit:")
[print("\t".join(map(str, r))) for r in snp if r[2] > 1]

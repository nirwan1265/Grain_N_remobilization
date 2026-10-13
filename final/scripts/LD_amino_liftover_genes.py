#!/usr/bin/env python3
"""
Lift the Goodman-Buckler GWAS lead SNPs and their LD intervals from B73 AGPv4 to
Zm-B73-REFERENCE-NAM-5.0 with the MaizeGDB/Ensembl chain file, then assign v5 genes.

Point lift-over: the aligned chain block (same chromosome, highest chain score)
containing the v4 position; a position in an alignment gap is moved to the nearest
aligned base of the best chain spanning it (flagged 'gap'). Interval end points are
lifted separately; if they land on different chromosomes or reverse order the
interval is replaced by the lifted lead position (flagged).

Genes per lead (B73 v5 GFF3 gene models):
  in_ld      gene overlaps the lifted LD interval (always includes a gene containing
             the lead)
  nearest    if no gene is in the interval, the closest gene and its distance
  in_25kb    gene within +/-25 kb of the lifted lead (the corrected window method)

usage: LD_amino_liftover_genes.py <ld_intervals_v4.tsv> <chain> <gff3> <out_prefix>
outputs: <out>_leads_v5.tsv, <out>_lead_genes.tsv
"""
import sys, bisect
ldf, chainf, gff, out = sys.argv[1:5]

# ---- chain: aligned blocks per v4 chromosome, same-chromosome chains only ----------
blocks = {}                                   # tchr -> list of (tS, tE, qS, score, chain_id)
cid = 0
with open(chainf) as fh:
    for line in fh:
        if line.startswith("chain"):
            f = line.split(); score = float(f[1]); tchr, tstart = f[2], int(f[5])
            qchr, qstrand, qstart = f[7], f[9], int(f[10])
            ok = (tchr == qchr and qstrand == "+" and tchr.isdigit()); t, q = tstart, qstart; cid += 1
            continue
        f = line.split()
        if not f: continue
        size = int(f[0])
        if ok: blocks.setdefault(tchr, []).append((t, t + size, q, score, cid))
        if len(f) == 3: t += size + int(f[1]); q += size + int(f[2])
for c in blocks: blocks[c].sort()
starts = {c: [b[0] for b in blocks[c]] for c in blocks}

def lift(c, p):
    """0-based chain coordinates; p is 1-based"""
    x = p - 1; B = blocks.get(c)
    if not B: return None, "no_chain"
    i = bisect.bisect_right(starts[c], x) - 1
    cand = []
    for k in range(max(0, i - 50), min(len(B), i + 50)):
        tS, tE, qS, sc, ch = B[k]
        if tS <= x < tE: cand.append((sc, qS + (x - tS)))
    if cand: return max(cand)[1] + 1, "ok"
    near = min(B[max(0, i - 2):i + 3], key=lambda b: min(abs(x - b[0]), abs(x - (b[1] - 1))))
    tS, tE, qS, sc, ch = near
    y = qS if x < tS else qS + (tE - tS - 1)
    return y + 1, "gap"

genes = {}
for l in open(gff):
    if "\tgene\t" not in l: continue
    f = l.split("\t"); c = f[0].replace("chr", "")
    if c.isdigit(): genes.setdefault(c, []).append((int(f[3]), int(f[4]), f[8].split(";")[0][3:]))
for c in genes: genes[c].sort()
def overlap(c, a, b): return [g for s, e, g in genes.get(c, []) if s <= b and e >= a]
def nearest(c, p):
    best = min(genes[c], key=lambda g: 0 if g[0] <= p <= g[1] else min(abs(p - g[0]), abs(p - g[1])))
    d = 0 if best[0] <= p <= best[1] else min(abs(p - best[0]), abs(p - best[1]))
    return best[2], d

rows = [l.rstrip("\n").split("\t") for l in open(ldf)]; H = rows[0]; rows = rows[1:]
ix = {h: i for i, h in enumerate(H)}
fo = open(out + "_leads_v5.tsv", "w"); fg = open(out + "_lead_genes.tsv", "w")
fo.write("lead\tchr\tpos_v4\tpos_v5\tlift\tld_found\tinterval_v5_start\tinterval_v5_end\tinterval_kb\tn_genes_ld\tnearest_gene\tnearest_dist_bp\tn_genes_25kb\n")
fg.write("lead\tchr\tpos_v5\tGeneID\tin_ld\tin_25kb\n")
stat = {"ok": 0, "gap": 0, "no_chain": 0}
for r in rows:
    lid, c, p4 = r[ix["lead"]], r[ix["chr"]], int(r[ix["pos_v4"]])
    p5, st = lift(c, p4); stat[st] += 1
    if p5 is None:
        fo.write(f"{lid}\t{c}\t{p4}\tNA\t{st}\t{r[ix['found']]}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"); continue
    if r[ix["found"]] == "TRUE":
        a5, _ = lift(c, int(r[ix["interval_start"]])); b5, _ = lift(c, int(r[ix["interval_end"]]))
        if a5 is None or b5 is None or not (a5 <= p5 <= b5) or b5 - a5 > 2e6: a5 = b5 = p5
    else:
        a5 = b5 = p5
    gl = overlap(c, a5, b5); g25 = overlap(c, p5 - 25000, p5 + 25000)
    ng, nd = (gl[0], 0) if gl else nearest(c, p5)
    fo.write(f"{lid}\t{c}\t{p4}\t{p5}\t{st}\t{r[ix['found']]}\t{a5}\t{b5}\t{(b5 - a5) / 1000:.1f}\t{len(gl)}\t{ng}\t{nd}\t{len(g25)}\n")
    for g in sorted(set(gl) | set(g25)): fg.write(f"{lid}\t{c}\t{p5}\t{g}\t{g in gl}\t{g in g25}\n")
print("lift:", stat)

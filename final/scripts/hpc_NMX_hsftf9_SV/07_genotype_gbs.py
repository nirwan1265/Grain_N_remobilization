#!/usr/bin/env python3
"""Step 7 (one genome per call): read each assembly's base at the 1,198 GBS SNPs.
For every SNP flank the best alignment (most matching bases) is used; the SNP is the
flank's centre base (query name gbs_<pos>|<start>|<end>, centre = pos - start), placed on
the genome by walking the CIGAR (cg tag) and read from the assembly. Bases are reported
on the B73 forward strand. NA = no alignment, or the SNP falls in an indel.
Also reports how many alignments cover >= 90% of the flank at >= 90% identity (n_good),
so multi-copy flanks can be dropped later.
usage: python3 07_genotype_gbs.py <genome.fa.gz> <OUT>   -> OUT/gbs_alleles/<genome>.tsv
"""
import sys, os, re, gzip
G, OUT = sys.argv[1], sys.argv[2]
name = os.path.basename(G).replace(".fa.gz", "")
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
best, ngood = {}, {}
for l in open(os.path.join(OUT, "paf_gbs", name + ".paf")):
    x = l.rstrip("\n").split("\t")
    q = x[0].split("|"); pos, start = int(q[0][4:]), int(q[1]); m, alen, qlen = int(x[9]), int(x[10]), int(x[1])
    if (int(x[3]) - int(x[2])) / qlen >= 0.9 and m / alen >= 0.9: ngood[pos] = ngood.get(pos, 0) + 1
    if pos not in best or m > best[pos][0]: best[pos] = (m, x, pos - start)
want = {}
for pos, (m, x, ctr) in best.items():
    qlen, qs, qe, strand, tname, ts = int(x[1]), int(x[2]), int(x[3]), x[4], x[5], int(x[7])
    cg = next(f[5:] for f in x[12:] if f.startswith("cg:Z:"))
    mapq = int(x[11]); s = ctr if strand == "+" else qlen - 1 - ctr
    q = qs if strand == "+" else qlen - qe; t = ts; tp = None
    for n, op in re.findall(r"(\d+)([MIDNX=])", cg):
        n = int(n)
        if op in "M=X":
            if q <= s < q + n: tp = t + (s - q); break
            q += n; t += n
        elif op == "I":
            if q <= s < q + n: break
            q += n
        elif op in "DN": t += n
    if tp is not None: want.setdefault(tname, []).append((pos, tp, strand, round(m / int(x[10]), 3), round((qe - qs) / qlen, 3), mapq))
res = {}
if want:
    need = sum(len(v) for v in want.values()); cur, seq = None, []
    def flush():
        if cur in want:
            s = "".join(seq)
            for pos, tp, strand, idt, cov, mq in want[cur]:
                b = s[tp].upper(); res[pos] = (b.translate(COMP) if strand == "-" else b, cur, tp + 1, strand, idt, cov, mq)
    with gzip.open(G, "rt") as f:
        for line in f:
            if line.startswith(">"):
                flush(); cur = line[1:].split()[0]; seq = []
                if len(res) == need: break
            elif cur in want: seq.append(line.strip())
        else: flush()
os.makedirs(os.path.join(OUT, "gbs_alleles"), exist_ok=True)
with open(os.path.join(OUT, "gbs_alleles", name + ".tsv"), "w") as o:
    o.write("genome\tpos\tallele\tcontig\ttpos\tstrand\tidentity\tqcov\tmapq\tn_good\n")
    for pos in sorted(best):
        r = res.get(pos)
        o.write("\t".join(map(str, [name, pos] + (list(r) if r else ["NA"] * 7) + [ngood.get(pos, 0)])) + "\n")
print(name, len(res), "of", len(best), "SNPs genotyped")

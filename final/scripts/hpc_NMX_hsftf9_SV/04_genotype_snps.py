#!/usr/bin/env python3
"""Step 4 (one genome per call): read the allele each genome carries at the four key SNPs.
Uses the SNP-flank alignments already in paf/<genome>.paf (step 2). For each SNP the best
alignment is taken, the SNP (centre base of the 301 bp flank) is placed on the genome by
walking the CIGAR (cg tag), and the base is read from the assembly. Bases are reported on the
B73 forward strand. A SNP inside an indel, or with no alignment, is reported as NA.
usage: python3 04_genotype_snps.py <genome.fa.gz> <OUT>
writes OUT/snp_alleles/<genome>.tsv
"""
import sys, os, re, gzip
G, OUT = sys.argv[1], sys.argv[2]
name = os.path.basename(G).replace(".fa.gz", "")
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
best = {}
for l in open(os.path.join(OUT, "paf", name + ".paf")):
    x = l.rstrip("\n").split("\t")
    if not x[0].startswith("snp_"): continue
    snp = x[0].split("|")[0]; m = int(x[9])
    if snp not in best or m > best[snp][0]: best[snp] = (m, x)
want = {}                                   # contig -> list of (snp, target pos 0-based, strand, identity)
for snp, (m, x) in best.items():
    qlen, qs, qe, strand, tname, ts = int(x[1]), int(x[2]), int(x[3]), x[4], x[5], int(x[7])
    cg = next(f[5:] for f in x[12:] if f.startswith("cg:Z:"))
    s = 150 if strand == "+" else qlen - 1 - 150      # SNP index on the strand the CIGAR walks
    q = qs if strand == "+" else qlen - qe; t = ts; pos = None
    for n, op in re.findall(r"(\d+)([MIDNX=])", cg):
        n = int(n)
        if op in "M=X":
            if q <= s < q + n: pos = t + (s - q); break
            q += n; t += n
        elif op == "I":
            if q <= s < q + n: break
            q += n
        elif op in "DN": t += n
    if pos is not None: want.setdefault(tname, []).append((snp, pos, strand, m / int(x[10])))
res = {}
if want:
    seq, cur = [], None
    def flush():
        if cur in want:
            s = "".join(seq)
            for snp, pos, strand, idt in want[cur]:
                b = s[pos].upper(); res[snp] = (b.translate(COMP) if strand == "-" else b, cur, pos + 1, strand, round(idt, 3))
    with gzip.open(G, "rt") as f:
        for line in f:
            if line.startswith(">"):
                flush(); cur = line[1:].split()[0]; seq = []
                if len(res) == sum(len(v) for v in want.values()): break
            elif cur in want: seq.append(line.strip())
        else: flush()
os.makedirs(os.path.join(OUT, "snp_alleles"), exist_ok=True)
with open(os.path.join(OUT, "snp_alleles", name + ".tsv"), "w") as o:
    o.write("genome\tsnp\tallele\tcontig\tpos\tstrand\tflank_identity\n")
    for snp in ["snp_151742636", "snp_151983481", "snp_151984305", "snp_151992800"]:
        r = res.get(snp)
        o.write("\t".join(map(str, [name, snp] + (list(r) if r else ["NA"] * 5))) + "\n")
print(name, {k: v[0] for k, v in res.items()})

#!/usr/bin/env python3
"""
all_maize2.vcf.gz -> a raw int8 genotype matrix plus the marker map.

Layout: one SNP's n bytes written consecutively, i.e. column-major for a
(taxa x SNP) matrix.  Markers stay in VCF order, so each chromosome is a
contiguous block of bytes and 01_GWAS_models.R can seek straight to it
without holding the whole genome in memory.

    python3 00_vcf_to_matrix.py <vcf.gz> <outdir>

Produces  GD.int8.bin  GM.csv  taxa.txt  dims.txt
Needs only numpy.  ~2 minutes for 3,539 x 495,755.
"""
import sys, gzip, os
import numpy as np

if len(sys.argv) != 3:
    sys.exit("usage: 00_vcf_to_matrix.py <vcf.gz> <outdir>")
vcf, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

fh = gzip.open(vcf, "rb")
for line in fh:
    if line.startswith(b"#CHROM"):
        taxa = line.rstrip().split(b"\t")[9:]
        break
n = len(taxa)
open(os.path.join(out, "taxa.txt"), "wb").write(b"\n".join(taxa) + b"\n")
print(f"taxa: {n}", flush=True)

gm = open(os.path.join(out, "GM.csv"), "w")
gm.write("SNP,Chromosome,Position\n")
k = nmiss = 0
with open(os.path.join(out, "GD.int8.bin"), "wb") as fb:
    for line in fh:
        f = line.split(b"\t", 9)
        if len(f) < 10:
            continue
        b = np.frombuffer(f[9].rstrip(), dtype=np.uint8)
        if b.size < 4 * n - 1:          # not a plain GT-only record
            continue
        a1 = b[0::4][:n].astype(np.int16)
        a2 = b[2::4][:n].astype(np.int16)
        miss = (a1 == 46) | (a2 == 46)  # '.'
        d = (a1 - 48) + (a2 - 48)
        if miss.any():
            obs = ~miss
            d[miss] = int(round(float(d[obs].mean()))) if obs.any() else 1
            nmiss += int(miss.sum())
        fb.write(d.astype(np.int8).tobytes())
        gm.write(f"{f[2].decode()},{f[0].decode()},{f[1].decode()}\n")
        k += 1
        if k % 100000 == 0:
            print(f"  {k} markers", file=sys.stderr, flush=True)
gm.close()
open(os.path.join(out, "dims.txt"), "w").write(f"{n}\t{k}\n")
print(f"markers: {k}   imputed calls: {nmiss} ({100*nmiss/(n*k):.5f}%)")
print(f"wrote {out}/GD.int8.bin ({n*k/1e9:.2f} GB)")

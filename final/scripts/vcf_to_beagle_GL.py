#!/usr/bin/env python3
"""
Convert the filtered Indian Chief / Jarvis VCF (GATK PL field) to a Beagle
genotype-likelihood file for NGSadmix and PCAngsd -- all 379 samples, same sites.

  GL = 10^(-PL/10), normalised to sum 1; a sample with no reads (DP = 0, ./.,
  or PL 0,0,0) gets the uninformative 1/3, 1/3, 1/3.
Site filters: biallelic SNP; >= 50% of samples with reads; MAF >= 0.05
(estimated from the genotype-likelihood-weighted dosage); then thinned to
one site per 10 kb so that linked SNPs do not dominate (NGSadmix assumes
independent sites).

usage: vcf_to_beagle_GL.py in.vcf.gz out.beagle.gz out_sample_order.txt
"""
import gzip, sys
import numpy as np
vcf, out, order = sys.argv[1:4]
MINFRAC, MINMAF, THIN = 0.5, 0.05, 10_000
fin = gzip.open(vcf, "rt")
for line in fin:
    if line.startswith("#CHROM"):
        samples = line.rstrip("\n").split("\t")[9:]; break
n = len(samples)
open(order, "w").write("\n".join(samples) + "\n")
fo = gzip.open(out, "wt")
fo.write("marker\tallele1\tallele2\t" + "\t".join(f"{s}\t{s}\t{s}" for s in samples) + "\n")
code = {"A": 0, "C": 1, "G": 2, "T": 3}
last_chr, last_pos = None, -10**12
seen = kept = drop_cov = drop_maf = 0
for line in fin:
    f = line.rstrip("\n").split("\t")
    seen += 1
    ref, alt = f[3], f[4]
    if len(ref) != 1 or len(alt) != 1 or ref not in code or alt not in code: continue
    chrom, pos = f[0], int(f[1])
    if chrom == last_chr and pos - last_pos < THIN: continue
    fmt = f[8].split(":"); ipl = fmt.index("PL"); idp = fmt.index("DP") if "DP" in fmt else None
    gl = np.full((n, 3), 1 / 3)
    has = np.zeros(n, bool)
    for i, s in enumerate(f[9:]):
        x = s.split(":")
        if len(x) <= ipl or x[ipl] in (".", ""): continue
        if idp is not None and len(x) > idp and x[idp] not in (".", "") and int(x[idp]) == 0: continue
        pl = x[ipl].split(",")
        if len(pl) != 3 or "." in pl: continue
        pl = np.array(pl, float)
        if pl.sum() == 0: continue
        l = 10 ** (-pl / 10); gl[i] = l / l.sum(); has[i] = True
    if has.mean() < MINFRAC: drop_cov += 1; continue
    dos = gl[has] @ np.array([0, 1, 2])           # expected alt dosage
    p = dos.mean() / 2
    if min(p, 1 - p) < MINMAF: drop_maf += 1; continue
    fo.write(f"{chrom}_{pos}\t{code[ref]}\t{code[alt]}\t" + "\t".join(f"{v:.6f}" for v in gl.ravel()) + "\n")
    kept += 1; last_chr, last_pos = chrom, pos
fo.close()
print(f"samples {n} | VCF records {seen:,} | kept {kept:,} (one per {THIN//1000} kb) | "
      f"dropped <{int(MINFRAC*100)}% coverage {drop_cov:,}, MAF<{MINMAF} {drop_maf:,}")

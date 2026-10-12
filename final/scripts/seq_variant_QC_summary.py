#!/usr/bin/env python3
"""Variant-level QC summaries for the 512,667 filtered SNPs used in all IC/J analyses.
Writes: site_depth.csv (mean depth per sample at each site = INFO DP / n samples, binned),
        maf_spectrum.csv (minor allele frequency, 0.025-wide bins, from called genotypes),
        titv_by_chr.csv (transition/transversion ratio per chromosome).
usage: seq_variant_QC_summary.py filtered.vcf.gz plink_prefix outdir"""
import gzip, sys, os, numpy as np
vcf, pre, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
ts = {("A","G"),("G","A"),("C","T"),("T","C")}
dp, titv = [], {}
n = None
for l in gzip.open(vcf, "rt"):
    if l.startswith("#CHROM"): n = len(l.split("\t")) - 9; continue
    if l[0] == "#": continue
    f = l.split("\t", 8)
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    dp.append(float(info.get("DP", "nan")) / n)
    c = titv.setdefault(f[0], [0, 0]); c[0 if (f[3], f[4]) in ts else 1] += 1
dp = np.array(dp)
edges = np.arange(0, 3.01, 0.1)
h, _ = np.histogram(np.clip(dp, 0, 2.999), edges)
np.savetxt(os.path.join(out, "site_depth.csv"), np.column_stack([edges[:-1], h]), delimiter=",",
           header="mean_depth_bin_start,n_sites", comments="", fmt=["%.1f", "%d"])
fam = sum(1 for _ in open(pre + ".fam")); m = sum(1 for _ in open(pre + ".bim")); nb = (fam + 3) // 4
raw = np.fromfile(pre + ".bed", np.uint8)[3:].reshape(m, nb)
alt = np.zeros(m); called = np.zeros(m)
for s in range(0, m, 50000):
    c = np.stack([(raw[s:s+50000] >> (2*k)) & 3 for k in range(4)], 2).reshape(-1, nb*4)[:, :fam]
    g = np.select([c == 0, c == 2, c == 3], [2, 1, 0], -1)
    called[s:s+50000] = (g >= 0).sum(1); alt[s:s+50000] = np.where(g >= 0, g, 0).sum(1)
p = alt / np.maximum(2 * called, 1); maf = np.minimum(p, 1 - p)
e2 = np.arange(0, 0.5001, 0.025); h2, _ = np.histogram(maf, e2)
np.savetxt(os.path.join(out, "maf_spectrum.csv"), np.column_stack([e2[:-1], h2]), delimiter=",",
           header="maf_bin_start,n_sites", comments="", fmt=["%.3f", "%d"])
with open(os.path.join(out, "titv_by_chr.csv"), "w") as fo:
    fo.write("chromosome,transitions,transversions,titv\n")
    for k in sorted(titv, key=lambda x: int(x) if x.isdigit() else 99):
        a, b = titv[k]; fo.write(f"{k},{a},{b},{a/b:.4f}\n")
tot = np.array(list(titv.values())).sum(0)
print(f"sites {len(dp):,} | mean site depth per sample: median {np.nanmedian(dp):.2f}x, IQR {np.nanpercentile(dp,25):.2f}-{np.nanpercentile(dp,75):.2f}x")
print(f"MAF: median {np.median(maf):.3f}; MAF<0.05 {100*np.mean(maf<0.05):.1f}%")
print(f"Ti/Tv genome-wide {tot[0]/tot[1]:.3f}; by chr " + ", ".join(f"{k}:{a/b:.2f}" for k,(a,b) in sorted(titv.items(), key=lambda x:int(x[0]))))

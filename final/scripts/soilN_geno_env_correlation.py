#!/usr/bin/env python3
"""
Genotype-environment correlation test for the soil-N GWAS SNPs.

For every one of the 495,755 landrace SNPs: Pearson r between allele dosage and
each of 9 collection-site variables (soil N 0-5 / 5-15 cm, latitude, |latitude|,
longitude, elevation, pH, CEC, aridity), over the 3,107 accessions with all
variables present. Each SNP's PERCENTILE is the share of the genome with an
equal or larger |r| for that variable (0.1 = top 0.1 % of the genome).
No model, no kinship, no PCs: a mixed model would absorb exactly the
geographic structure this test is meant to expose.

usage:  soilN_geno_env_correlation.py BED_PREFIX PHENO_CSV OUT_CSV_GZ [chunk]
"""
import sys, os, gzip, numpy as np

bed_prefix, pheno_csv, out_csv = sys.argv[1:4]
CHUNK = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
VARS = ["N_0_5", "N_5_15", "latitude", "abs_latitude", "longitude",
        "elevation", "pH", "CEC", "aridity"]

fam = [l.split()[1] for l in open(bed_prefix + ".fam")]
snps = [l.split()[1] for l in open(bed_prefix + ".bim")]
n, M = len(fam), len(snps)

hdr = open(pheno_csv).readline().strip().split(",")
ph = {}
for l in open(pheno_csv).readlines()[1:]:
    f = l.strip().split(",")
    try:
        ph[f[0]] = [float(f[hdr.index(v)]) for v in VARS]
    except ValueError:
        continue
keep = np.array([i for i, t in enumerate(fam)
                 if t in ph and all(np.isfinite(ph[t]))])
E = np.array([ph[fam[i]] for i in keep])
E = (E - E.mean(0)) / E.std(0)
print(f"{n} genotyped, {len(keep)} with all 9 variables, {M} SNPs", flush=True)

# PLINK 2-bit codes -> A2 dosage (01 = missing)
lut = np.array([0, np.nan, 1, 2], dtype=np.float32)
nb = (n + 3) // 4
idx = np.arange(256, dtype=np.uint8)
codes = np.stack([(idx >> (2 * k)) & 3 for k in range(4)], 1)      # 256 x 4
dec = lut[codes]                                                    # 256 x 4

R = np.empty((M, len(VARS)), dtype=np.float32)
with open(bed_prefix + ".bed", "rb") as fh:
    assert fh.read(3) == b"\x6c\x1b\x01", "not a SNP-major PLINK bed"
    for s in range(0, M, CHUNK):
        k = min(CHUNK, M - s)
        raw = np.frombuffer(fh.read(nb * k), dtype=np.uint8).reshape(k, nb)
        G = dec[raw].reshape(k, nb * 4)[:, :n][:, keep]
        mu = np.nanmean(G, 1, keepdims=True)
        G = np.where(np.isnan(G), mu, G) - mu
        sd = np.sqrt((G * G).mean(1, keepdims=True)); sd[sd == 0] = np.nan
        R[s:s + k] = (G / sd) @ E / len(keep)
        print(f"  {s + k}/{M}", flush=True)

# percentile: % of genome with |r| >= this SNP's |r|  (smaller = more extreme)
A = np.abs(np.nan_to_num(R, nan=0.0))
P = np.empty_like(A)
for j in range(A.shape[1]):
    o = np.argsort(-A[:, j], kind="stable")
    rk = np.empty(M); rk[o] = np.arange(1, M + 1)
    P[:, j] = 100.0 * rk / M

with gzip.open(out_csv, "wt") as out:
    out.write("SNP," + ",".join("r_" + v for v in VARS) + "," +
              ",".join("pct_" + v for v in VARS) + "\n")
    for i in range(M):
        out.write(snps[i] + "," + ",".join(f"{x:.5f}" for x in R[i]) + "," +
                  ",".join(f"{x:.5g}" for x in P[i]) + "\n")
print("written:", out_csv)

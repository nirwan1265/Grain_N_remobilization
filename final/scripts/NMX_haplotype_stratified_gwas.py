#!/usr/bin/env python3
"""Genome-wide GWAS within each NMX-hsftf9 haplotype group (sanity analysis, not for the PlantPhys paper).

Groups are landraces homozygous at both tag SNPs of the interval
  short (B73 type) : lead chr9:151,742,636 A/A and 152 Mb SNP chr9:151,983,481 T/T
  long             : lead G/G and 152 Mb SNP C/C
plus the full soil-N panel (3,107 landraces) run the same way for comparison.

Model  y ~ 1 + SNP + g,  g ~ N(0, s2 K). K = VanRaden kinship from 20,000 random SNPs
(MAF >= 0.05 within the group), leave-one-chromosome-out. delta by ML on the null model per
trait and chromosome (P3D), then a GLS Wald t test per SNP. SNPs with MAF < 0.01 in the group
are skipped. Traits (z-scored within group): soil N 0-5 cm, aridity, pH.

usage: python3 NMX_haplotype_stratified_gwas.py <repo> <bed_prefix> <out_dir> [groups, default short,long,full] [chromosomes, default 1-10]
(the full panel is slow, so it can be run a few chromosomes at a time and the gwas_full_chr*.csv.gz files joined)
"""
import sys, os, time, numpy as np, pandas as pd
from scipy import stats, optimize
repo, bed, OUT = sys.argv[1:4]
GROUPS = (sys.argv[4] if len(sys.argv) > 4 else "short,long,full").split(",")
CHROMS = [int(c) for c in sys.argv[5].split(",")] if len(sys.argv) > 5 else list(range(1, 11))
os.makedirs(OUT, exist_ok=True)
SN = os.path.join(repo, "final/results/soil_nitrogen_gwas")
fam = [l.split()[1] for l in open(bed + ".fam")]; n = len(fam)
bim = pd.read_csv(bed + ".bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
raw = np.memmap(bed + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), (n + 3) // 4)
lut = np.array([2, np.nan, 1, 0], dtype=np.float32)                       # dosage of a1
env = pd.read_csv(os.path.join(SN, "control_phenotypes_maize_landrace.csv")).drop_duplicates("Taxa").set_index("Taxa")
keepset = set(open(os.path.join(SN, "LD/landrace_N_keep.txt")).read().split())
keep = np.array([i for i, t in enumerate(fam) if t in keepset and t in env.index]); taxa = [fam[i] for i in keep]
E = env.loc[taxa]
TRAITS = ["N_0_5", "aridity", "pH"]

def G(idx, cols):
    b = np.asarray(raw[np.atleast_1d(idx)])
    g = np.stack([(b >> s) & 3 for s in (0, 2, 4, 6)], axis=-1).reshape(len(b), -1)[:, :n]
    return lut[g[:, cols]]
def dose(pos, allele):
    i = bim.index[(bim.chr == 9) & (bim.pos == pos)][0]; x = G(i, keep)[0]
    return x if bim.a1[i] == allele else (2 - x if bim.a2[i] == allele else None)
leadA, t152 = dose(151742636, "A"), dose(151983481, "T")
sel_all = {"short": (leadA == 2) & (t152 == 2), "long": (leadA == 0) & (t152 == 0), "full": np.ones(len(taxa), bool)}
print({k: int(v.sum()) for k, v in sel_all.items()}, flush=True)

for grp in GROUPS:
    t0 = time.time(); sel = sel_all[grp]; cols = keep[sel]; m = len(cols)
    Y = {t: ((lambda v: (v - v.mean()) / v.std())(E[t].values[sel].astype(float))) for t in TRAITS}
    # kinship SNPs
    rng = np.random.default_rng(1); cand = np.sort(rng.choice(len(bim), 60000, replace=False))
    M = G(cand, cols); f = np.nanmean(M, 1) / 2; ok = (f >= 0.05) & (f <= 0.95)
    kidx = cand[ok][:20000]; M = M[ok][:20000]; f = f[ok][:20000]; kchr = bim.chr.values[kidx]
    Z = M - 2 * f[:, None]; Z[np.isnan(Z)] = 0; del M
    ZtZ = Z.T.astype(np.float64) @ Z.astype(np.float64); h = 2 * f * (1 - f)
    res = []
    for c in CHROMS:
        zc = Z[kchr == c].astype(np.float64)
        K = (ZtZ - zc.T @ zc) / h[kchr != c].sum(); del zc
        lam, U = np.linalg.eigh(K); lam = np.clip(lam, 0, None); del K
        one = U.T @ np.ones(m); W, YT = {}, {}
        for t in TRAITS:
            yt = U.T @ Y[t]
            def nll(ld):
                d = lam + np.exp(ld); w = 1 / d; b = (one * w * yt).sum() / (one * w * one).sum(); r = yt - one * b
                return 0.5 * (m * np.log((r * r * w).sum() / m) + np.log(d).sum())
            W[t] = 1 / (lam + np.exp(optimize.minimize_scalar(nll, bounds=(-8, 8), method="bounded").x)); YT[t] = yt
        idx = bim.index[bim.chr == c].values
        for s in range(0, len(idx), 4000):
            bi = idx[s:s + 4000]; X = G(bi, cols).astype(np.float64)
            mu = np.nanmean(X, 1); X = np.where(np.isnan(X), mu[:, None], X); maf = np.minimum(mu / 2, 1 - mu / 2)
            okk = maf >= 0.01
            if not okk.any(): continue
            bi, X, maf = bi[okk], X[okk], maf[okk]
            XT = U.T @ X.T                                                   # m x snps
            out = {"chr": c, "pos": bim.pos.values[bi], "snp": bim.id.values[bi], "maf": np.round(maf, 4)}
            for t in TRAITS:
                w = W[t]; yt = YT[t]; cw = (one * w * one).sum()
                xr = XT - np.outer(one, (one * w) @ XT / cw); yr = yt - one * ((one * w * yt).sum() / cw)
                sxx = (w[:, None] * xr * xr).sum(0); sxy = (w * yr) @ xr; syy = (w * yr * yr).sum()
                b = sxy / sxx; s2 = (syy - b * sxy) / (m - 2); tt = b / np.sqrt(s2 / sxx)
                out[f"beta_{t}"] = np.round(b, 4); out[f"P_{t}"] = 2 * stats.t.sf(np.abs(tt), m - 2)
            res.append(pd.DataFrame(out))
        print(grp, "chr", c, "done", round(time.time() - t0), "s", flush=True)
    tag = grp if len(CHROMS) == 10 else grp + "_chr" + "-".join(map(str, CHROMS))
    R = pd.concat(res); R.to_csv(os.path.join(OUT, f"gwas_{tag}.csv.gz"), index=False, float_format="%.4g")
    print(grp, "n", m, "SNPs", len(R), "time", round(time.time() - t0), "s", flush=True)

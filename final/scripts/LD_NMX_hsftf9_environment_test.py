#!/usr/bin/env python3
"""Is the long-range LD between the NMX lead (chr9:151,742,636) and the hsftf9-region SNP
(chr9:151,983,481) explained by both loci tracking the environment or ancestry?
r^2 between the two SNPs (allele dosage, 3,107 landraces, GBS all_maize2) is compared
  (1) overall and within narrow elevation / aridity / latitude bins, within Mexico and overall,
  (2) after regressing both dosages on site variables (lat, long, elevation, pH, CEC, aridity, soil N),
  (3) after regressing both dosages on k genotype PCs (k = 3, 10, 20, 50; PCs from 20,000 random
      SNPs with MAF >= 0.05 on chromosomes other than 9),
  (4) against a genome-wide null of SNP pairs 200-280 kb apart with MAF >= 0.3.
Output: final/results/soil_nitrogen_gwas/LD_GBS_regional/NMX_hsftf9_LD_environment_test.csv
usage: python3 LD_NMX_hsftf9_environment_test.py <repo> <bed_prefix>
"""
import sys, os, numpy as np, pandas as pd
repo, bed = sys.argv[1:3]; rng = np.random.default_rng(1)
SN = os.path.join(repo, "final/results/soil_nitrogen_gwas")
fam = [l.split()[1] for l in open(bed + ".fam")]; n = len(fam)
env = pd.read_csv(os.path.join(SN, "control_phenotypes_maize_landrace.csv")).drop_duplicates("Taxa").set_index("Taxa")
geo = pd.read_csv(os.path.join(SN, "focal_SNP_genotypes_geography.csv")).drop_duplicates("Taxa").set_index("Taxa")[["region", "country"]]
keep = np.array([i for i, t in enumerate(fam) if t in env.index]); taxa = [fam[i] for i in keep]
E = env.loc[taxa].join(geo, how="left"); print("accessions", len(keep))
bim = pd.read_csv(bed + ".bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
nb = (n + 3) // 4; raw = np.memmap(bed + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), nb)
lut = np.array([2, np.nan, 1, 0])
def G(idx):
    b = np.asarray(raw[np.atleast_1d(idx)])
    g = np.stack([(b >> s) & 3 for s in (0, 2, 4, 6)], axis=-1).reshape(len(np.atleast_1d(idx)), -1)[:, :n]
    return lut[g][:, keep]
ia = bim.index[(bim.chr == 9) & (bim.pos == 151742636)][0]; ib = bim.index[(bim.chr == 9) & (bim.pos == 151983481)][0]
A, B = G(ia)[0], G(ib)[0]
def r2(x, y, m=None):
    m = (~np.isnan(x) & ~np.isnan(y)) if m is None else (m & ~np.isnan(x) & ~np.isnan(y))
    return np.corrcoef(x[m], y[m])[0, 1] ** 2 if m.sum() > 20 and x[m].std() > 0 and y[m].std() > 0 else np.nan
def resid(y, X):
    m = ~np.isnan(y) & ~np.isnan(X).any(1); out = np.full_like(y, np.nan)
    Xm = np.column_stack([np.ones(m.sum()), X[m]]); out[m] = y[m] - Xm @ np.linalg.lstsq(Xm, y[m], rcond=None)[0]; return out
rows = [dict(test="overall", group="all", n=len(A), r2=r2(A, B))]
# (1) strata
for scope, sel in [("all", np.ones(len(A), bool)), ("Mexico", (E.country == "MEX").values)]:
    for v in ["elevation", "aridity", "latitude"]:
        q = pd.qcut(E[v].where(sel), 5, labels=False, duplicates="drop").values
        for k in range(5):
            m = sel & (q == k)
            rows.append(dict(test=f"within {v} quintile", group=f"{scope} Q{k+1} ({E[v][m].min():.0f}-{E[v][m].max():.0f})", n=int(m.sum()), r2=r2(A, B, m)))
# (2) environment-adjusted
X = E[["latitude", "longitude", "elevation", "pH", "CEC", "aridity", "N_0_5"]].values.astype(float)
rows.append(dict(test="adjusted for site variables", group="all", n=len(A), r2=r2(resid(A, X), resid(B, X))))
# (3) PC-adjusted
cand = bim.index[(bim.chr != 9)].values; pick = np.sort(rng.choice(cand, 40000, replace=False))
Gp = G(pick); f = np.nanmean(Gp, 1) / 2; Gp = Gp[(f >= 0.05) & (f <= 0.95)][:20000]
Gp = Gp - np.nanmean(Gp, 1, keepdims=True); Gp[np.isnan(Gp)] = 0
U, S, _ = np.linalg.svd(Gp.T, full_matrices=False)
for k in (3, 10, 20, 50):
    rows.append(dict(test=f"adjusted for {k} genotype PCs", group="all", n=len(A), r2=r2(resid(A, U[:, :k]), resid(B, U[:, :k]))))
# (4) genome-wide null of distant common pairs
f_all = None; null = []
for c in range(1, 11):
    w = bim[bim.chr == c]; idx = rng.choice(w.index.values, 150, replace=False)
    for i in idx:
        p = bim.pos[i]; j = w.index[(w.pos > p + 200000) & (w.pos < p + 280000)].values
        if len(j) == 0: continue
        j = rng.choice(j, min(5, len(j)), replace=False); g = G(np.r_[i, j]); fr = np.nanmean(g, 1) / 2
        if min(fr[0], 1 - fr[0]) < 0.3: continue
        for t in range(1, len(g)):
            if min(fr[t], 1 - fr[t]) >= 0.3: null.append(r2(g[0], g[t]))
null = np.array([x for x in null if not np.isnan(x)])
obs = rows[0]["r2"]
rows.append(dict(test="genome-wide null, pairs 200-280 kb, MAF>=0.3", group=f"{len(null)} pairs", n=len(null),
                 r2=float(np.median(null)), note=f"95th pct {np.percentile(null,95):.3f}, 99th {np.percentile(null,99):.3f}, max {null.max():.3f}, frac >= observed {np.mean(null>=obs):.4f}"))
out = pd.DataFrame(rows); out["r2"] = out.r2.round(3)
out.to_csv(os.path.join(SN, "LD_GBS_regional/NMX_hsftf9_LD_environment_test.csv"), index=False)
pd.set_option("display.width", 200); print(out.to_string(index=False))

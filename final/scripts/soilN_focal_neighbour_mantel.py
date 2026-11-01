#!/usr/bin/env python3
"""Distance-based checks of the two focal soil-N SNPs (S1_16128632 omega-amidase,
S5_22338979 aspartate kinase isozyme) against geography.

1. Nearest-neighbour test: for each minor-allele carrier, compare its soil N with
   the mean soil N of its k nearest non-carriers (great-circle distance). One-sided
   Wilcoxon signed-rank test that carriers sit on richer soil than their own
   neighbours. Run for k = 5 and k = 10, for neighbours restricted to <= 50 km and
   <= 100 km, and at the site level (one value per collection site, to avoid
   counting accessions from the same site twice).
2. Mantel and partial Mantel tests (accession level): genotype distance
   |dosage_i - dosage_j| vs soil N distance |N_i - N_j|, controlling for
   great-circle distance; 999 permutations of the genotype matrix.

Input : focal_SNP_genotypes_geography.csv (SNP, Taxa, alt_dosage, N_0_5, lat, lon, ...)
Output: focal_SNP_neighbour_test.csv, focal_SNP_mantel_test.csv
Usage : python3 soilN_focal_neighbour_mantel.py <input.csv> <out_dir>
"""
import sys, os
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

inp = sys.argv[1] if len(sys.argv) > 1 else "final/results/soil_nitrogen_gwas/focal_SNP_genotypes_geography.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "final/results/soil_nitrogen_gwas"
rng = np.random.default_rng(20261029)
NPERM = 999
GENE = {"S1_16128632": "omega-amidase (Zm00001eb005710)",
        "S5_22338979": "aspartate kinase isozyme (Zm00001eb220040)"}

def hav(lat1, lon1, lat2, lon2):
    """great-circle distance in km; lat2/lon2 may be arrays"""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))

def neighbour_test(car, non, k, maxkm, label):
    diffs, dists = [], []
    for _, c in car.iterrows():
        d = hav(c.lat, c.lon, non.lat.values, non.lon.values)
        o = np.argsort(d, kind="stable")[:k]
        o = o[d[o] <= maxkm] if maxkm else o
        if len(o) == 0:
            continue
        diffs.append(c.N_0_5 - non.N_0_5.values[o].mean()); dists.append(np.median(d[o]))
    diffs = np.array(diffs)
    nz = diffs[diffs != 0]
    p = wilcoxon(nz, alternative="greater").pvalue if len(nz) >= 5 else np.nan
    return dict(unit=label, k=k, max_km=maxkm if maxkm else "none", n_carriers_tested=len(diffs),
                n_dropped_no_neighbour=len(car) - len(diffs),
                median_neighbour_km=round(float(np.median(dists)), 1) if dists else np.nan,
                median_diff_cg_kg=round(float(np.median(diffs)), 1) if len(diffs) else np.nan,
                mean_diff_cg_kg=round(float(diffs.mean()), 1) if len(diffs) else np.nan,
                frac_carrier_higher=round(float((diffs > 0).mean()), 3) if len(diffs) else np.nan,
                P_signed_rank_one_sided=p)

d0 = pd.read_csv(inp)
NB, MT = [], []
for snp, d in d0.groupby("SNP"):
    d = d.reset_index(drop=True)
    car, non = d[d.alt_dosage > 0], d[d.alt_dosage == 0]
    # accession level
    for k in (5, 10):
        for mk in (None, 100, 50):
            NB.append({"SNP": snp, "gene": GENE.get(snp, ""), **neighbour_test(car, non, k, mk, "accession")})
    # site level: carrier sites = sites with >= 1 carrier, compared with sites without carriers
    s = d.groupby(["lat", "lon"], as_index=False).agg(N_0_5=("N_0_5", "first"), anyc=("alt_dosage", lambda x: (x > 0).any()))
    for k in (5, 10):
        for mk in (None, 100, 50):
            NB.append({"SNP": snp, "gene": GENE.get(snp, ""),
                       **neighbour_test(s[s.anyc], s[~s.anyc], k, mk, "site")})

    # Mantel / partial Mantel, accession level
    n = len(d); iu = np.triu_indices(n, 1)
    g = d.alt_dosage.values.astype(float); N = d.N_0_5.values.astype(float)
    lat, lon = d.lat.values, d.lon.values
    D = np.vstack([hav(lat[i], lon[i], lat, lon) for i in range(n)])[iu]
    Nv = np.abs(N[:, None] - N[None, :])[iu]
    def zs(v): v = v - v.mean(); return v / np.sqrt((v * v).sum())
    Dz, Nz = zs(D), zs(Nv)
    r_ND = float(Dz @ Nz)
    def stats(gg):
        G = zs(np.abs(gg[:, None] - gg[None, :])[iu])
        r_GN, r_GD = float(G @ Nz), float(G @ Dz)
        part = (r_GN - r_GD * r_ND) / np.sqrt((1 - r_GD ** 2) * (1 - r_ND ** 2))
        return r_GN, r_GD, part
    obs = stats(g)
    null = np.array([stats(g[rng.permutation(n)]) for _ in range(NPERM)])
    pv = [(1 + (null[:, j] >= obs[j]).sum()) / (NPERM + 1) for j in range(3)]
    MT.append(dict(SNP=snp, gene=GENE.get(snp, ""), n=n, r_genotype_soilN=round(obs[0], 4), P_mantel_soilN=pv[0],
                   r_genotype_geodist=round(obs[1], 4), P_mantel_geodist=pv[1],
                   r_soilN_geodist=round(r_ND, 4),
                   partial_r_genotype_soilN_given_geodist=round(obs[2], 4), P_partial_mantel=pv[2],
                   permutations=NPERM))
    print(snp, "done", flush=True)

nb, mt = pd.DataFrame(NB), pd.DataFrame(MT)
nb.to_csv(os.path.join(out, "focal_SNP_neighbour_test.csv"), index=False)
mt.to_csv(os.path.join(out, "focal_SNP_mantel_test.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(nb.to_string(index=False)); print(mt.to_string(index=False))

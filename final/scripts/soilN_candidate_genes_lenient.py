#!/usr/bin/env python3
"""Soil N (0-5 cm) environmental GWAS: lenient candidate-gene list for overlap
with other resources (e.g. cis-eQTL). Any SNP with P < P_CUT in MLM, BLINK or
FarmCPU (3 PCs) -> genes within WIN bp on B73 RefGen_v5, with PANTHER/GO.

Output: final/results/soil_nitrogen_gwas/soilN_candidate_genes_P1e-4.csv (one row
per gene x SNP x model) and ..._by_gene.csv (one row per gene, best P).
"""
import os, bisect
import numpy as np, pandas as pd

REPO = os.environ.get("GN_REPO", os.path.expanduser("~/mnt/Grain_N_remobilization"))
D = f"{REPO}/final/results/soil_nitrogen_gwas"
P_CUT, WIN = 1e-4, 25_000
MODELS = ["MLM", "BLINK", "FarmCPU"]

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

rows = []
for m in MODELS:
    g = pd.read_csv(f"{D}/N_0_5_full/raw_GWAS_{m}_3PC_N_0_5.txt", sep="\t",
                    usecols=["SNP", "Chr", "Pos", "P.value", "MAF", "Effect"])
    g = g.dropna(subset=["P.value"])
    g["FDR"] = bh(g["P.value"].values)          # BH over the full scan, per model
    g = g[g["P.value"] < P_CUT].copy(); g["model"] = m; rows.append(g)
S = pd.concat(rows, ignore_index=True)
S["sig_FDR5"] = S["FDR"] <= 0.05
S["n_models_P1e-4"] = S.groupby("SNP")["model"].transform("nunique")

# genes (v5)
gff = pd.read_csv(f"{REPO}/final/results/annotation/Zm-B73-REFERENCE-NAM-5.0_genes_only.gff3",
                  sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8])
gff.columns = ["chr", "type", "start", "end", "attr"]
gff = gff[gff["type"] == "gene"]
gff["GeneID"] = gff["attr"].str.extract(r"ID=([^;]+)")
gff["Chr"] = pd.to_numeric(gff["chr"].str.replace("chr", ""), errors="coerce")
gff = gff.dropna(subset=["Chr"]); gff["Chr"] = gff["Chr"].astype(int)

out = []
for c, gc in gff.groupby("Chr"):
    gc = gc.sort_values("start"); st = gc["start"].values; en = gc["end"].values; ids = gc["GeneID"].values
    for _, s in S[S["Chr"] == c].iterrows():
        p = s["Pos"]; lo = bisect.bisect_left(st, p - WIN - 5_000_000)  # genes are < 5 Mb
        for i in range(lo, bisect.bisect_right(st, p + WIN)):
            if en[i] < p - WIN: continue
            rel = "within" if st[i] <= p <= en[i] else ("upstream_coord" if p < st[i] else "downstream_coord")
            dist = 0 if rel == "within" else (st[i] - p if p < st[i] else p - en[i])
            out.append({**s.to_dict(), "GeneID": ids[i], "GeneStart": st[i], "GeneEnd": en[i],
                        "Relation": rel, "Distance_bp": int(dist)})
G = pd.DataFrame(out)

ann = pd.read_csv(f"{REPO}/final/results/GO/GO_universe_genomewide.txt", sep="\t", dtype=str)
ann["GeneName"] = ann["GeneName"].str.split(";").str[0]
G = G.merge(ann[["GeneID", "GeneName", "Family_Subfamily", "GO_BP", "GO_MF"]], on="GeneID", how="left")
cols = ["GeneID", "GeneName", "SNP", "Chr", "Pos", "model", "P.value", "FDR", "sig_FDR5",
        "n_models_P1e-4", "MAF", "Effect", "Relation", "Distance_bp", "GeneStart", "GeneEnd",
        "Family_Subfamily", "GO_BP", "GO_MF"]
G = G[cols].sort_values(["P.value", "GeneID"])
G.to_csv(f"{D}/soilN_candidate_genes_P1e-4.csv", index=False)

best = G.sort_values("P.value").groupby("GeneID", as_index=False).agg(
    GeneName=("GeneName", "first"), best_P=("P.value", "first"), best_model=("model", "first"),
    best_SNP=("SNP", "first"), Chr=("Chr", "first"), min_distance_bp=("Distance_bp", "min"),
    models=("model", lambda x: "+".join(sorted(set(x)))), any_FDR5=("sig_FDR5", "any"),
    n_SNPs=("SNP", "nunique"), GO_BP=("GO_BP", "first"), GO_MF=("GO_MF", "first"),
    Family_Subfamily=("Family_Subfamily", "first"))
best = best.sort_values("best_P")
best.to_csv(f"{D}/soilN_candidate_genes_P1e-4_by_gene.csv", index=False)
print(f"SNP x model rows P<{P_CUT}: {len(S)}; unique SNPs {S.SNP.nunique()}; "
      f"genes within {WIN//1000} kb: {best.GeneID.nunique()}; of which any 5% FDR: {best.any_FDR5.sum()}")
print(S.groupby("model").size())

# ---- geography/environment call for SNPs that pass 5% FDR (Supplemental Table S12)
s12 = pd.read_csv(f"{REPO}/final/supp_tables/SuppTable_S12_soilN_geno_env_correlation.csv", usecols=["SNP", "call"])
G = G.merge(s12.rename(columns={"call": "geo_env_call_FDR_SNPs"}), on="SNP", how="left")
best = best.merge(s12.rename(columns={"SNP": "best_SNP", "call": "geo_env_call_best_SNP"}), on="best_SNP", how="left")
G.to_csv(f"{D}/soilN_candidate_genes_P1e-4.csv", index=False)
best.to_csv(f"{D}/soilN_candidate_genes_P1e-4_by_gene.csv", index=False)

readme = pd.DataFrame({"Item": [
  "Trait", "Panel", "Models", "Inclusion", "Gene window", "FDR", "geo/env call", "Annotation", "Sheets"],
 "Description": [
  "Soil total N at 0-5 cm at the collection site (SoilGrids 2.0, 250 m, cg/kg)",
  "3,107 georeferenced CIMMYT landraces, 495,755 SNPs, B73 RefGen_v5 coordinates",
  "GAPIT MLM, BLINK and FarmCPU with 3 genotype PCs",
  f"Lenient: any SNP with P < {P_CUT:g} in at least one model (not a significance threshold)",
  f"RefGen_v5 genes within {WIN//1000} kb of the SNP (Relation and Distance_bp given)",
  "Benjamini-Hochberg over the full scan, per model. sig_FDR5 = the SNPs used in the paper",
  "Only for SNPs passing 5% FDR: linked to environment / not linked / inconclusive (genotype-environment correlation test, paper Table S12)",
  "PANTHER family and GO BP/MF terms (PANTHER export)",
  "by_gene = one row per gene with its best P; gene_SNP_model = every gene x SNP x model row"]})
with pd.ExcelWriter(f"{D}/soilN_candidate_genes_P1e-4.xlsx") as xw:
    readme.to_excel(xw, sheet_name="README", index=False)
    best.to_excel(xw, sheet_name="by_gene", index=False)
    G.to_excel(xw, sheet_name="gene_SNP_model", index=False)
print("xlsx written; calls among best SNPs:", best["geo_env_call_best_SNP"].value_counts(dropna=False).to_dict())

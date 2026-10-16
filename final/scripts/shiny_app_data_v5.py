#!/usr/bin/env python3
"""Rebuild the data files of the amino acid GWAS Shiny app from the current (v5) results.

The deployed app used SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv, which
assigned genes from AGPv4 SNP positions read as if they were v5, kept every SNP with
P < 1e-3, and listed 114 traits (36 of them pairwise ratios that are not in the paper).
This script replaces it with the paper's tables:
  hits        <- final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv
                 (SNPs significant at 5% BH FDR per phenotype x model, lifted to
                  B73 RefGen_v5, genes within 25 kb)
  annotation  <- final/results/GO/GO_universe_genomewide.txt (v5 gene annotation)
  phenotypes  <- final/supp_tables/SuppTable_S1_FAA_phenotypes.csv (78 traits)
  definitions <- final/supp_tables/SuppTable_S2_phenotype_definitions.csv
Output: shiny/amino_gwas_app/data/supplementary/{amino_gwas_hits_v5,gene_annotation_v5,
        FAA_phenotypes,trait_definitions}.csv
Run from the repo root.
"""
import csv, math, os
OUT = "shiny/amino_gwas_app/data/supplementary"
os.makedirs(OUT, exist_ok=True)

ann = {}
for r in csv.DictReader(open("final/results/GO/GO_universe_genomewide.txt"), delimiter="\t"):
    ann[r["GeneID"]] = {"GeneID": r["GeneID"], "GeneSymbol": (r["GeneName"] or "").split(";")[0],
                        "Family_Subfamily": r.get("Family_Subfamily", ""), "Protein_Class": r.get("Protein_Class", ""),
                        "GO_MF": r.get("GO_MF", ""), "GO_BO": r.get("GO_BP", ""), "GO_CC": r.get("GO_CC", "")}

hits, used = [], set()
for r in csv.DictReader(open("final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv")):
    if not r["GeneID"]:
        continue
    p = float(r["min_P"]); used.add(r["GeneID"])
    hits.append({"GeneID": r["GeneID"], "GeneSymbol": ann.get(r["GeneID"], {}).get("GeneSymbol", r["GeneID"]),
                 "GeneChr": f"chr{r['Chr']}", "GeneStart": r["GeneStart"], "GeneEnd": r["GeneEnd"],
                 "Phenotype": r["Phenotype"], "Model": r["models_sig"].replace("+", ";"), "SNP": r["SNP"],
                 "Chr": r["Chr"], "SNP_Pos": r["Pos_v5"], "P.value": r["min_P"], "FDR": r["min_FDR"],
                 "log10_P": f"{-math.log10(p):.4f}", "Relation": r["Relation"],
                 "Distance_to_Gene_bp": r["Distance_to_Gene_bp"]})
with open(f"{OUT}/amino_gwas_hits_v5.csv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(hits[0])); w.writeheader(); w.writerows(hits)
with open(f"{OUT}/gene_annotation_v5.csv", "w", newline="") as fo:
    cols = ["GeneID", "GeneSymbol", "Family_Subfamily", "Protein_Class", "GO_MF", "GO_BO", "GO_CC"]
    w = csv.DictWriter(fo, fieldnames=cols); w.writeheader()
    w.writerows([ann.get(g, {"GeneID": g, **{c: "" for c in cols[1:]}}) for g in sorted(used)])
for src, dst in (("SuppTable_S1_FAA_phenotypes.csv", "FAA_phenotypes.csv"),
                 ("SuppTable_S2_phenotype_definitions.csv", "trait_definitions.csv")):
    open(f"{OUT}/{dst}", "w").write(open(f"final/supp_tables/{src}", encoding="utf-8-sig").read())
print(f"hits {len(hits)} rows, {len({h['Phenotype'] for h in hits})} phenotypes, {len(used)} genes, "
      f"{len({(h['Phenotype'], h['SNP']) for h in hits})} phenotype x SNP")

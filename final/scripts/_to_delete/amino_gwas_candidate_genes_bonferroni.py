#!/usr/bin/env python3
"""
Results paragraph 3 -- amino-acid GWAS candidate genes (Goodman-Buckler panel).

Recomputes every number quoted in Results paragraph 3 from the annotated GWAS
hit table, replacing the previous fixed -log10(p) >= 7 cutoff with an explicit
Bonferroni threshold derived from the number of markers tested.

    alpha = 0.05,  m = 4,117,796 SNPs
    threshold = alpha / m = 1.214e-08   (-log10 = 7.916)

Input
    shiny/amino_gwas_app/data/supplementary/
        SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv
    (one row per Gene x Phenotype, carrying the single best SNP across the four
     GAPIT models; produced by scripts/12.2_make_amino_gene_summary.R)

Output
    final/supp_tables/SuppTable2_amino_gwas_candidate_genes_bonferroni.csv
    final/supp_tables/SuppTable2b_amino_gwas_genes_per_trait.csv

NOTE ON THE "Model" COLUMN
    The input table keeps only the best SNP per Gene x Phenotype across models,
    so the Model column records which model happened to give that best p-value.
    It is NOT a model-overlap statistic and must not be used to claim that
    associations are "dominated by" any model.  A genuine multi-model consensus
    requires the full per-model scans (see consensus_four_model_snps.R).

Usage:  python3 final/scripts/amino_gwas_candidate_genes_bonferroni.py [repo_root]
"""

import csv
import sys
import math
import os
import collections

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."

IN_FILE = os.path.join(
    ROOT, "shiny", "amino_gwas_app", "data", "supplementary",
    "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv",
)
OUT_DIR = os.path.join(ROOT, "final", "supp_tables")
os.makedirs(OUT_DIR, exist_ok=True)

N_SNP = 4_117_796
ALPHA = 0.05
P_BONF = ALPHA / N_SNP

PROXY_TRAITS = ("Total_N", "Total_PBAA")


def main():
    with open(IN_FILE, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        rows = list(reader)

    n_pheno_tested = len({r["Phenotype"] for r in rows})
    sig = [r for r in rows if float(r["P.value"]) <= P_BONF]

    pheno_by_gene = collections.defaultdict(set)
    genes_by_pheno = collections.defaultdict(set)
    model_tally = collections.Counter()
    for r in sig:
        pheno_by_gene[r["GeneID"]].add(r["Phenotype"])
        genes_by_pheno[r["Phenotype"]].add(r["GeneID"])
        model_tally[r["Model"]] += 1

    n_single = sum(1 for v in pheno_by_gene.values() if len(v) == 1)
    n_multi = sum(1 for v in pheno_by_gene.values() if len(v) > 1)
    per_trait = sorted(
        ((p, len(g)) for p, g in genes_by_pheno.items()),
        key=lambda x: (-x[1], x[0]),
    )

    print(f"threshold                       : P <= {P_BONF:.4e}  (-log10 = {-math.log10(P_BONF):.4f})")
    print(f"unique candidate genes          : {len(pheno_by_gene)}")
    print(f"traits with >=1 candidate gene  : {len(genes_by_pheno)} of {n_pheno_tested}")
    print(f"genes in exactly one phenotype  : {n_single}")
    print(f"genes in more than one phenotype: {n_multi}")
    print("top traits                      : "
          + ", ".join(f"{p}={n}" for p, n in per_trait[:8]))
    for t in PROXY_TRAITS:
        print(f"  {t:<11}                   : {len(genes_by_pheno.get(t, ()))}")
    print(f"best-hit model tally (see note) : {dict(model_tally)}")

    sig.sort(key=lambda r: float(r["P.value"]))
    out1 = os.path.join(OUT_DIR, "SuppTable2_amino_gwas_candidate_genes_bonferroni.csv")
    with open(out1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sig)

    out2 = os.path.join(OUT_DIR, "SuppTable2b_amino_gwas_genes_per_trait.csv")
    with open(out2, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Phenotype", "N_candidate_genes"])
        w.writerows(per_trait)

    print(f"\nwrote {out1}\nwrote {out2}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Main Table: GWAS results for the focal free amino acids of the Goodman-Buckler panel.

Rows are the FAAs singled out in the Results: the two most abundant (Asn, Pro),
the most heritable (Lys, Asn, Asp, Pro by line-mean H2), the amide Gln (second
most variable, genomic h2 0.85) and the most variable (Met).
For each FAA:
  median, CV, H2, genomic h2   -> Supplementary Table S3
  loci   -> Supplementary Table S7 loci with >= 1 phenotype reporting on this
            amino acid (absolute trait and its ratios: phenotype name starts
            with the amino acid code); the number for the absolute trait alone
            is given in parentheses
  best P -> smallest P among those phenotypes (amino_gwas_v5_fdr_snps.csv: every
            FDR-significant phenotype x SNP, written by amino_gwas_tables_v5.R)
  genes  -> trait-matched candidates (amino_gwas_pathway_matched_genes.py),
            listed when hit by a phenotype that reports on this amino acid;
            the P, phenotype and SNP-gene distance are for that phenotype.
Output : final/main_tables/MainTable_amino_gwas_focus_traits.tex
"""
import csv, os
R = os.environ.get("GN_REPO", os.path.expanduser("~/Documents/Github/Grain_N_remobilization"))
S3 = {r["AA"]: r for r in csv.DictReader(open(f"{R}/final/supp_tables/SuppTable_S3_FAA_abundance_CV.csv"))}
S7 = list(csv.DictReader(open(f"{R}/final/supp_tables/SuppTable_S7_amino_gwas_loci.csv")))
S6 = list(csv.DictReader(open(f"{R}/final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv")))
FDR = list(csv.DictReader(open(f"{R}/final/results/goodman_buckler_amino_acid_gwas/amino_gwas_v5_fdr_snps.csv")))

ROWS = ["N", "P", "K", "D", "Q", "M"]
## curated display labels for the trait-matched genes (annotation in brackets)
GENES = {
 "N": [("Zm00001eb368680", "amino acid transporter AVT6C"),
       ("Zm00001eb329710", "glutamine amidotransferase type 2"),
       ("Zm00001eb405410", "homoserine dehydrogenase"),
       ("Zm00001eb036060", "asparagine synthetase domain"),
       ("Zm00001eb036070", "asparagine synthetase domain"),
       ("Zm00001eb389290", "dihydrodipicolinate synthase"),
       ("Zm00001eb176170", "dihydrodipicolinate reductase")],
 "P": [("Zm00001eb217220", "glutamine amidotransferase")],
 "K": [("Zm00001eb277160", "NRT1/PTR family 4.6"),
       ("Zm00001eb360480", "glutamate synthase (NADH)"),
       ("Zm00001eb352400", "WAT1-related (UmamiT)")],
 "D": [("Zm00001eb207630", "WAT1-related (UmamiT)")],
 "Q": [("Zm00001eb379220", "WAT1-related (UmamiT)"),
       ("Zm00001eb261470", "amino acid transporter"),
       ("Zm00001eb261480", "amino acid transporter"),
       ("Zm00001eb261490", "amino acid transporter")],
 "M": []}

def own(ph, aa): return ph.split(".")[0] == aa
def sci(p):
    m, e = f"{p:.1e}".split("e"); return f"${m}\\times10^{{{int(e)}}}$"
def dist(r):
    d = int(float(r["Distance_to_Gene_bp"]))
    return "in gene" if d == 0 else f"{d/1000:.1f} kb"
def tt(s): return "\\texttt{" + s.replace("_", "\\_") + "}"

lines = []
for aa in ROWS:
    s = S3[aa]
    loci = [l for l in S7 if any(own(p, aa) for p in l["phenotypes"].split(", "))]
    loci_abs = [l for l in S7 if aa in l["phenotypes"].split(", ")]
    bestp = min(float(r["min_P"]) for r in FDR if own(r["Phenotype"], aa))   # all FDR SNPs, with or without a gene
    cells = []
    for g, lab in GENES[aa]:
        hits = [r for r in S6 if r["GeneID"] == g and own(r["Phenotype"], aa)]
        b = min(hits, key=lambda r: float(r["min_P"]))
        two = "$^{\\dagger}$" if b["replicated"].upper() == "TRUE" else ""
        cells.append(f"\\textit{{{g}}} ({lab}; {tt(b['Phenotype'])}, {sci(float(b['min_P']))}{two}, {dist(b)})")
    lines.append(" & ".join([
        f"{s['Amino_acid']} ({aa})", f"{float(s['median']):.3f}", f"{float(s['CV_pct']):.0f}",
        f"{float(s['H2_line_mean']):.2f}", f"{float(s['h2_genomic']):.2f}",
        f"{len(loci)} ({len(loci_abs)})", sci(bestp), "; ".join(cells) if cells else "none"]) + " \\\\")

tex = r"""\begin{table*}[t]
\centering
\caption{Association mapping of the focal free amino acids (FAAs) in the Goodman--Buckler panel: the two most abundant (asparagine, proline), the most heritable (lysine, asparagine, aspartate, proline), the amide glutamine and the most variable (methionine). Median, coefficient of variation (CV), line-mean broad-sense heritability ($H^2$) and genomic heritability ($h^2_g$) are from Supplementary Table~S3. Loci: loci (Supplementary Table~S7) significant at 5\% Benjamini--Hochberg FDR for at least one phenotype reporting on the amino acid (its absolute abundance and its ratios); in parentheses, loci for the absolute trait alone. Best $P$: smallest $P$ among those phenotypes. Trait-matched candidate genes: B73 RefGen\_v5 genes within 25~kb of a significant SNP for one of those phenotypes whose annotation places them in the metabolism of the amino acid or its biosynthetic family, in amino acid or nitrogen transport, or in primary nitrogen assimilation (Methods); for each, the phenotype with the smallest $P$, that $P$, and the distance from the SNP to the gene are given. $^{\dagger}$Significant in both BLINK and FarmCPU.}
\label{tab:amino_gwas_focus}
\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabularx}{\textwidth}{@{}l r r r r r r >{\raggedright\arraybackslash}X@{}}
\hline
FAA & Median & CV (\%) & $H^2$ & $h^2_g$ & Loci & Best $P$ & Trait-matched candidate genes (annotation; phenotype, $P$, distance) \\
\hline
""" + "\n".join(lines) + r"""
\hline
\end{tabularx}
\end{table*}
"""
out = f"{R}/final/main_tables/MainTable_amino_gwas_focus_traits.tex"
open(out, "w").write(tex); print(tex)

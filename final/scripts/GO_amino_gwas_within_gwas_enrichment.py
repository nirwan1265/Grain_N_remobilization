#!/usr/bin/env python3
"""
GO over-representation WITHIN the amino acid GWAS: are the loci of one biosynthetic
family (or one focal amino acid) richer in a GO term than amino acid GWAS loci in
general?  Companion to GO_amino_gwas_locus_enrichment.py (genome background).

Units are loci, not genes, so tandem gene clusters count once:
  * all FDR-significant SNPs (any phenotype) are merged into loci when within 25 kb
    (as in Supplementary Table S7); each locus gets the annotated v5 genes within
    25 kb of its SNPs and the set of phenotypes significant there;
  * universe = the GWAS loci that carry >= 1 GO-annotated gene;
  * a set (family / focal amino acid) = universe loci with >= 1 of its phenotypes;
  * for each GO term (BP + MF, PANTHER, direct annotations) a locus "has" the term
    if any of its genes does.
Test: one-sided hypergeometric (= one-sided Fisher exact test of set vs rest) on
locus counts, for terms present at >= 3 loci of the set; BH within each set.
Caveat: set loci may carry more genes per locus than the rest; the number of genes
per locus is reported for each set so this can be checked.

Inputs/annotation as in GO_amino_gwas_locus_enrichment.py.
Output : results/GO/GO_amino_gwas_within_gwas_enrichment.csv
"""
import csv, os, re, bisect, random, math, collections, sys
gv = lambda k, d: os.environ.get(k, d)
REPO = gv("GN_REPO", os.path.expanduser("~/Documents/Github/Grain_N_remobilization"))
FDR  = gv("GN_FDRSNP", f"{REPO}/final/results/goodman_buckler_amino_acid_gwas/amino_gwas_v5_fdr_snps.csv")
GOU  = gv("GN_GOU", f"{REPO}/final/results/GO/GO_universe_genomewide.txt")
GFF  = gv("GN_GFF", os.path.expanduser("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"))
OUT  = gv("GN_OUT", f"{REPO}/final/results/GO/GO_amino_gwas_locus_enrichment.csv")
WIN, N_PERM, MIN_LOCI = 25000, int(gv("GN_NPERM", "1000")), 3
random.seed(1)

# ---- genes and genome coordinates
chrlen, genes = {}, []
for line in open(GFF):
    if line.startswith("#"): continue
    x = line.rstrip("\n").split("\t")
    if len(x) < 9 or not re.fullmatch(r"chr(\d+)", x[0]): continue
    c = int(x[0][3:])
    if not 1 <= c <= 10: continue
    if x[2] == "chromosome": chrlen[c] = int(x[4])
    elif x[2] == "gene": genes.append((c, int(x[3]), int(x[4]), re.search(r"ID=([^;]+)", x[8]).group(1)))
off, o = {}, 0
for c in sorted(chrlen): off[c] = o; o += chrlen[c]
G = o
# ---- GO annotation
term_name, gene_terms = {}, {}
for r in csv.DictReader(open(GOU), delimiter="\t"):
    ts = set()
    for col in ("GO_BP", "GO_MF"):
        for it in (r.get(col) or "").split(";"):
            m = re.match(r"(.*)\((GO:\d{7})\)\s*$", it.strip())
            if m: ts.add(m.group(2)); term_name[m.group(2)] = (m.group(1).strip(), col[3:])
    if ts: gene_terms[r["GeneID"]] = ts
G_ann = [(off[c] + s, off[c] + e, g) for c, s, e, g in genes if g in gene_terms]
G_ann.sort(); starts = [a for a, _, _ in G_ann]
MAXLEN = max(b - a for a, b, _ in G_ann)
N_bg = len(G_ann); K_bg = collections.Counter(t for _, _, g in G_ann for t in gene_terms[g])
def genes_in(a, b):
    i0 = bisect.bisect_left(starts, a - MAXLEN); i1 = bisect.bisect_right(starts, b)
    return [G_ann[i][2] for i in range(i0, i1) if G_ann[i][1] >= a]
# ---- phenotype sets
AA_FAM = dict(D="Aspartate", N="Aspartate", T="Aspartate", M="Aspartate", K="Aspartate", I="Aspartate",
              E="Glutamate", Q="Glutamate", P="Glutamate", R="Glutamate", H="Glutamate",
              S="Serine", G="Serine", C="Serine", F="Aromatic", Y="Aromatic", W="Aromatic",
              A="Pyruvate", V="Pyruvate", L="Pyruvate")
COMP = dict(IMNTDK="Aspartate", EHPRQ="Glutamate", SGC="Serine", WFY="Aromatic", LAV="Pyruvate", IVL="Pyruvate")
def fam_of(p):
    if p.startswith("Total"): return "Whole pool"
    if "." in p:
        q = p.split(".")
        for k in (q[1], q[0]):
            if k in COMP: return COMP[k]
        return AA_FAM.get(q[0], "Whole pool")
    return COMP.get(p, AA_FAM.get(p, "Whole pool"))
snps = list(csv.DictReader(open(FDR)))
SETS = collections.OrderedDict()
SETS["All phenotypes"] = lambda p: True
for f in ("Aspartate", "Glutamate", "Serine", "Pyruvate", "Aromatic"): SETS[f"Family: {f}"] = (lambda f: lambda p: fam_of(p) == f)(f)
for aa, nm in (("N", "Asn"), ("P", "Pro"), ("Q", "Gln"), ("K", "Lys"), ("D", "Asp")):
    SETS[f"Amino acid: {nm}"] = (lambda aa: lambda p: p.split(".")[0] == aa)(aa)

def loci_for(sel):
    pts = sorted({(int(r["Chr"]), int(float(r["Pos_v5"]))) for r in snps if sel(r["Phenotype"])})
    L = []
    for c, p in pts:
        if L and L[-1][0] == c and p - L[-1][2] <= WIN: L[-1][2] = p
        else: L.append([c, p, p])
    return [(off[c] + s - WIN, off[c] + e + WIN) for c, s, e in L]


OUT2 = gv("GN_OUT2", f"{REPO}/final/results/GO/GO_amino_gwas_within_gwas_enrichment.csv")
pts = collections.defaultdict(set)
for r in snps: pts[(int(r["Chr"]), int(float(r["Pos_v5"])))].add(r["Phenotype"])
L = []
for (c, p) in sorted(pts):
    if L and L[-1]["c"] == c and p - L[-1]["e"] <= WIN: L[-1]["e"] = p; L[-1]["ph"] |= pts[(c, p)]
    else: L.append(dict(c=c, s=p, e=p, ph=set(pts[(c, p)])))
for l in L:
    l["genes"] = genes_in(off[l["c"]] + l["s"] - WIN, off[l["c"]] + l["e"] + WIN)
    l["terms"] = set().union(*[gene_terms[g] for g in l["genes"]]) if l["genes"] else set()
U = [l for l in L if l["genes"]]
N = len(U); Kt = collections.Counter(t for l in U for t in l["terms"])
def lc(a, b): return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
def hyper_sf(k, K, n, N):
    tot = lc(N, n); return min(1.0, sum(math.exp(lc(K, x) + lc(N - K, n - x) - tot) for x in range(k, min(K, n) + 1)))
rows = []
print(f"GWAS loci: {len(L)}; with annotated genes (universe): {N}; mean genes/locus {sum(len(l['genes']) for l in U)/N:.2f}")
for name, sel in SETS.items():
    if name == "All phenotypes": continue
    S = [l for l in U if any(sel(p) for p in l["ph"])]
    n = len(S); kt = collections.Counter(t for l in S for t in l["terms"])
    gpl = sum(len(l["genes"]) for l in S) / n
    gpl_rest = sum(len(l["genes"]) for l in U if l not in S) / max(1, N - n)
    res = []
    for t, k in kt.items():
        if k < MIN_LOCI: continue
        exp = n * Kt[t] / N
        res.append(dict(gene_set=name, set_loci=n, universe_loci=N, genes_per_locus_set=round(gpl, 2),
                        genes_per_locus_rest=round(gpl_rest, 2), GO_id=t, term=term_name[t][0], ontology=term_name[t][1],
                        loci_in_set=k, loci_in_universe=Kt[t], loci_expected=round(exp, 2), fold=round(k / exp, 2),
                        P=hyper_sf(k, Kt[t], n, N),
                        genes=";".join(sorted({g for l in S for g in l["genes"] if t in gene_terms[g]}))))
    o = sorted(range(len(res)), key=lambda j: res[j]["P"]); m = len(res); prev = 1.0
    for rank in range(m - 1, -1, -1):
        j = o[rank]; prev = min(prev, res[j]["P"] * m / (rank + 1)); res[j]["FDR"] = prev
    res.sort(key=lambda r: r["P"]); rows += res
    print(f"{name}: {n} loci (genes/locus {gpl:.2f} vs rest {gpl_rest:.2f}), {m} terms tested, FDR<0.05: {sum(r['FDR'] < 0.05 for r in res)}, FDR<0.10: {sum(r['FDR'] < 0.10 for r in res)}")
    for r in res[:3]:
        print(f"    {r['term'][:55]:55s} {r['loci_in_set']}/{n} vs {r['loci_in_universe']}/{N} (exp {r['loci_expected']}) P={r['P']:.2g} FDR={r['FDR']:.2g}")
w = csv.DictWriter(open(OUT2, "w", newline=""), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("wrote", OUT2, len(rows), "rows")

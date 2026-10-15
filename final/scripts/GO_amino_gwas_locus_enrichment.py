#!/usr/bin/env python3
"""
GO over-representation for the Goodman-Buckler amino acid GWAS candidate genes,
counted per locus and tested against a circular-shift null.

Why per locus: many gene families (UGTs, peptidases, ...) sit in tandem arrays, so
one association within 25 kb of a cluster adds several genes with the same GO term
and inflates gene-level enrichment. Here each term is counted at most once per locus.

Gene sets (from every FDR-significant phenotype x SNP, amino_gwas_v5_fdr_snps.csv):
  * biosynthetic families (phenotype -> family as in amino_gwas_tables_v5.R):
    Aspartate, Glutamate, Serine, Pyruvate, Aromatic; and All phenotypes
  * focal amino acids (phenotype name starts with the amino acid code):
    Asn (N), Pro (P), Gln (Q), Lys (K), Asp (D)
For each set, its significant SNPs are merged into loci when within 25 kb, and the
locus window is [first SNP - 25 kb, last SNP + 25 kb] (B73 RefGen_v5). Candidate
genes are v5 genes overlapping the window, as in Supplementary Table S6.

Annotation: PANTHER GO molecular function and biological process terms for all v5
genes (GO_universe_genomewide.txt; direct annotations, not propagated to parent
terms). Background: all annotated v5 genes on chromosomes 1-10.

Test per term (terms annotated to >= MIN_LOCI loci of the set):
  observed = number of the set's loci with >= 1 gene carrying the term
  null     = gene-number-matched blocks: each locus that carries m annotated genes is
             replaced by a block of m consecutive annotated genes starting at a random
             gene (N_PERM draws). This keeps the number of loci, the number of genes per
             locus and the clustering of neighbouring (often tandem) genes. A simple
             circular shift of the windows is NOT used: GWAS loci lie in gene-rich
             regions, and shifted windows fall mostly in gene-poor repeat space, which
             makes every common term look enriched.
  P        = (#null >= observed + 1) / (N_PERM + 1); BH within each set.
The conventional gene-level hypergeometric test against the same background is
reported alongside for comparison only (it ignores gene clustering).

Inputs : results/goodman_buckler_amino_acid_gwas/amino_gwas_v5_fdr_snps.csv,
         results/GO/GO_universe_genomewide.txt, v5 GFF3 (GN_GFF)
Output : results/GO/GO_amino_gwas_locus_enrichment.csv (all tested terms)
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

rows = []
for name, sel in SETS.items():
    W = loci_for(sel)
    obs_loci_genes = [genes_in(a, b) for a, b in W]
    obs = collections.Counter(t for gl in obs_loci_genes for t in set().union(*[gene_terms[g] for g in gl]) if gl)
    set_genes = set(g for gl in obs_loci_genes for g in gl)
    gcount = collections.Counter(t for g in set_genes for t in gene_terms[g])
    tested = [t for t, k in obs.items() if k >= MIN_LOCI]
    idx = {t: i for i, t in enumerate(tested)}
    gidx = {}
    null_ge = [0] * len(tested); null_sum = [0] * len(tested)
    obs_v = [obs[t] for t in tested]
    sizes = [len(gl) for gl in obs_loci_genes if gl]
    nA = len(G_ann)
    for it in range(N_PERM):
        cnt = [0] * len(tested)
        for m in sizes:
            i0 = random.randrange(0, nA - m + 1)
            hit = set()
            for j in range(i0, i0 + m):
                g = G_ann[j][2]
                if g not in gidx: gidx[g] = [idx[t] for t in gene_terms[g] if t in idx]
                hit.update(gidx[g])
            for i in hit: cnt[i] += 1
        for i, v in enumerate(cnt):
            null_sum[i] += v
            if v >= obs_v[i]: null_ge[i] += 1
    n_set = len(set_genes)
    def hyper_sf(k, K, n, N):  # P(X >= k), hypergeometric, via log-binomials
        lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
        tot = lc(N, n); return min(1.0, sum(math.exp(lc(K, x) + lc(N - K, n - x) - tot) for x in range(k, min(K, n) + 1)))
    res = []
    for i, t in enumerate(tested):
        exp = null_sum[i] / N_PERM
        res.append(dict(gene_set=name, n_loci=len(W), n_loci_with_genes=len(sizes), n_genes=n_set, GO_id=t, term=term_name[t][0], ontology=term_name[t][1],
                        loci_observed=obs_v[i], loci_expected=round(exp, 2), fold=round(obs_v[i] / exp, 2) if exp else "",
                        P_locus_perm=(null_ge[i] + 1) / (N_PERM + 1),
                        genes_observed=gcount[t], genes_background=K_bg[t],
                        P_gene_hypergeom=hyper_sf(gcount[t], K_bg[t], n_set, N_bg),
                        genes=";".join(sorted(g for g in set_genes if t in gene_terms[g]))))
    for key, fk in (("P_locus_perm", "FDR_locus"), ("P_gene_hypergeom", "FDR_gene")):
        o = sorted(range(len(res)), key=lambda j: res[j][key]); m = len(res); prev = 1.0
        for rank in range(m - 1, -1, -1):
            j = o[rank]; prev = min(prev, res[j][key] * m / (rank + 1)); res[j][fk] = prev
    rows += sorted(res, key=lambda r: r["P_locus_perm"])
    sig = [r for r in res if r["FDR_locus"] < 0.05]
    print(f"{name}: {len(W)} loci ({len(sizes)} with genes), {n_set} annotated genes, {len(tested)} terms tested, FDR_locus<0.05: {len(sig)}; "
          f"gene-level FDR<0.05: {sum(r['FDR_gene'] < 0.05 for r in res)}", flush=True)
    for r in sorted(res, key=lambda r: r["P_locus_perm"])[:4]:
        print(f"    {r['term'][:55]:55s} loci {r['loci_observed']} vs {r['loci_expected']}  P={r['P_locus_perm']:.4f} FDR={r['FDR_locus']:.3f} | genes {r['genes_observed']} geneFDR={r['FDR_gene']:.2g}", flush=True)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
w = csv.DictWriter(open(OUT, "w", newline=""), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("wrote", OUT, len(rows), "rows; background", N_bg, "annotated genes")

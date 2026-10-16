#!/usr/bin/env python3
"""Does any CornCyc pathway recur across the four candidate gene sets more than chance?

Gene sets (B73 RefGen_v5 gene IDs):
  GB    Goodman-Buckler amino acid GWAS genes (S6; FDR 5%; 'Bonf' version keeps
        genes whose best SNP has P < 0.05 / 4,177,796)
  soilN landrace soil-N GWAS genes (S8; FDR 5%; 'Bonf' keeps min_P < 0.05 / 495,755)
  IC    Indian Chief top-5% F_ST window genes (S20)
  J     Jarvis top-5% F_ST window genes (S21)

For every CornCyc pathway (523 pathways) and every set we compute
  hit      = number of pathway genes in the set
  P_any    = chance that a random gene set of the same size contains >= 1 pathway gene
             (hypergeometric; background = all annotated v5 genes + all CornCyc genes)
A pathway "recurs" if it is hit in every set. Its chance probability is
  P_all = product of P_any over the sets (sets are from independent analyses).
Because the pathway of interest was chosen after seeing the results, P_all is judged
against all 523 pathways: we report its rank, how many pathways recur at least as
improbably, and a Bonferroni/BH correction over pathways.

Two schemes: 4 sets (GB, soilN, IC, J) and 3 populations (GB, soilN, IC or J).
Two thresholds: FDR gene lists and Bonferroni gene lists for the two GWAS.

usage: pathway_recurrence_test.py CORNCYC_FILE  (run from repo root)
out:   final/results/pathway_recurrence/pathway_recurrence_{FDR,Bonf}.csv
"""
import csv, sys, os, math
from collections import defaultdict

cyc = sys.argv[1]
R = "final/supp_tables"
OUT = "final/results/pathway_recurrence"; os.makedirs(OUT, exist_ok=True)

def col(path, key="GeneID", filt=None):
    s = set()
    for r in csv.DictReader(open(path)):
        if r.get(key) and r[key].startswith("Zm") and (filt is None or filt(r)):
            s.add(r[key])
    return s

B_GB, B_SN = 0.05 / 4177796, 0.05 / 495755
sets = {
  "FDR":  {"GB": col(f"{R}/SuppTable_S6_amino_gwas_candidate_genes.csv"),
           "soilN": col(f"{R}/SuppTable_S9_soilN_candidate_genes.csv")},
  "Bonf": {"GB": col(f"{R}/SuppTable_S6_amino_gwas_candidate_genes.csv", filt=lambda r: float(r["min_P"]) < B_GB),
           "soilN": col(f"{R}/SuppTable_S9_soilN_candidate_genes.csv", filt=lambda r: float(r["min_P"]) < B_SN)},
}
IC = col(f"{R}/SuppTable_S20_IC_Fst_top5_genes.csv"); J = col(f"{R}/SuppTable_S21_J_Fst_top5_genes.csv")
for k in sets: sets[k].update({"IC": IC, "J": J})

pw_genes, pw_name = defaultdict(set), {}
for r in csv.DictReader(open(cyc), delimiter="\t"):
    g = (r.get("Gene-id") or "").strip()
    if g.upper().startswith("ZM00001EB"):          # file mixes Zm00001eb and ZM00001EB
        g = "Zm00001eb" + g[9:].split("_")[0]
        pw_genes[r["Pathway-id"]].add(g); pw_name[r["Pathway-id"]] = r["Pathway-name"]

universe = set(r["GeneID"] for r in csv.DictReader(open("final/results/GO/GO_universe_genomewide.txt"), delimiter="\t"))
for s in pw_genes.values(): universe |= s
for v in sets.values():
    for s in v.values(): universe |= s
N = len(universe)

lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
def p_any(K, n):          # P(>=1 of K pathway genes in a random set of n genes)
    return 1.0 - math.exp(lc(N - K, n) - lc(N, n)) if K and n else 0.0

for thr, S in sets.items():
    S4 = {k: S[k] & universe for k in ("GB", "soilN", "IC", "J")}
    S3 = {"GB": S4["GB"], "soilN": S4["soilN"], "ICorJ": S4["IC"] | S4["J"]}
    rows = []
    for pid, genes in pw_genes.items():
        K = len(genes); row = {"pathway_id": pid, "pathway": pw_name[pid], "n_genes": K}
        for scheme, SS in (("4set", S4), ("3pop", S3)):
            pall, allhit = 1.0, True
            for k, s in SS.items():
                h = sorted(genes & s); pa = p_any(K, len(s))
                if scheme == "4set" or k == "ICorJ":
                    row[f"{k}_hits"] = len(h); row[f"{k}_genes"] = ";".join(h); row[f"{k}_P_any"] = round(pa, 4)
                pall *= pa; allhit &= bool(h)
            row[f"{scheme}_all_hit"] = allhit; row[f"{scheme}_P_all"] = pall
        rows.append(row)
    for scheme in ("4set", "3pop"):
        hit = [r for r in rows if r[f"{scheme}_all_hit"]]
        hit.sort(key=lambda r: r[f"{scheme}_P_all"])
        m = len(rows)
        for i, r in enumerate(hit, 1):
            r[f"{scheme}_rank"] = i
            r[f"{scheme}_P_bonf"] = min(1.0, r[f"{scheme}_P_all"] * m)
        # BH over all pathways (non-recurring pathways get P = 1)
        ps = sorted([(r[f"{scheme}_P_all"] if r[f"{scheme}_all_hit"] else 1.0, j) for j, r in enumerate(rows)])
        q = [0.0] * m; prev = 1.0
        for i in range(m - 1, -1, -1):
            p, j = ps[i]; prev = min(prev, p * m / (i + 1)); q[j] = prev
        for j, r in enumerate(rows): r[f"{scheme}_P_BH"] = q[j]
    rows.sort(key=lambda r: (not r["4set_all_hit"], r["4set_P_all"]))
    keys = list(rows[0].keys()) + [k for k in ("4set_rank", "4set_P_bonf", "3pop_rank", "3pop_P_bonf") if k not in rows[0]]
    with open(f"{OUT}/pathway_recurrence_{thr}.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    print(f"\n=== {thr} gene lists | N={N} | set sizes: " + ", ".join(f"{k} {len(v)}" for k, v in S4.items()) + f" | pathways {len(rows)}")
    for scheme in ("4set", "3pop"):
        hit = sorted([r for r in rows if r[f"{scheme}_all_hit"]], key=lambda r: r[f"{scheme}_P_all"])
        print(f"  [{scheme}] pathways hit in every set: {len(hit)}")
        for r in hit[:15]:
            print(f"    {r[scheme+'_rank']:>3}  P_all={r[scheme+'_P_all']:.4f}  BH={r[scheme+'_P_BH']:.2f}  K={r['n_genes']:>3}  {r['pathway'][:70]}")

# ---- multiple-testing-aware tail test (Tarone-style) -------------------------
# Most pathways cannot recur improbably (big pathways recur easily), so instead of a
# Bonferroni over all pathways we ask: summed over every pathway whose chance of
# recurring is <= that of the observed best pathway, how many such rare recurrences
# would be expected? P(>=1) = 1 - exp(-expected).
# ---- custom Asp/Asn core set (post hoc; reported for completeness) ------------
core_ids = ["ASPARTATESYN-PWY", "ASPARAGINE-BIOSYNTHESIS", "ASPARAGINE-DEG1-PWY", "HOMOSERSYN-PWY"]
core = set().union(*[pw_genes[i] for i in core_ids if i in pw_genes])
with open(f"{OUT}/pathway_recurrence_summary.csv", "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["thresholds", "scheme", "n_pathways", "observed_recurring", "expected_recurring",
                                    "best_pathway", "best_P_all", "expected_as_rare_or_rarer", "P_tail"])
    for thr, S in sets.items():
        S4 = {k: S[k] & universe for k in ("GB", "soilN", "IC", "J")}
        S3 = {"GB": S4["GB"], "soilN": S4["soilN"], "ICorJ": S4["IC"] | S4["J"]}
        rows = list(csv.DictReader(open(f"{OUT}/pathway_recurrence_{thr}.csv")))
        for sc in ("4set", "3pop"):
            P = [float(r[f"{sc}_P_all"]) for r in rows]
            hit = [r for r in rows if r[f"{sc}_all_hit"] == "True"]
            best = min(hit, key=lambda r: float(r[f"{sc}_P_all"]))
            pb = float(best[f"{sc}_P_all"]); e = sum(p for p in P if p <= pb)
            w.writerow([thr, sc, len(rows), len(hit), round(sum(P), 2), best["pathway"], f"{pb:.4g}", round(e, 4), round(1 - math.exp(-e), 3)])
            print(f"{thr} {sc}: best '{best['pathway']}' P_all={pb:.4f}; expected rare recurrences={e:.3f}; P_tail={1-math.exp(-e):.3f}")
        SS = S4
        pc = 1.0; hits = {}
        for k, s in SS.items(): pc *= p_any(len(core), len(s)); hits[k] = sorted(core & s)
        print(f"  {thr} Asp/Asn core (K={len(core)}): P_all={pc:.4f}; hits " + "; ".join(f"{k}:{','.join(v) or '-'}" for k, v in hits.items()))

#!/usr/bin/env python3
"""Genome-aware permutation test: does a CornCyc pathway recur across the four candidate
sets more often than chance, after accounting for clustered intervals and for having
examined every pathway?

Sets are rebuilt as genomic intervals (B73 RefGen_v5), exactly as the candidate genes were
defined:
  GB    amino acid GWAS SNPs +/- 25 kb   (Bonf: min_P < 0.05/4,177,796; FDR: all)
  soilN soil-N GWAS SNPs +/- 25 kb       (Bonf: min_P < 0.05/495,755;   FDR: all 234)
  IC/J  top-5% F_ST windows (S17 / S18)
A gene belongs to a set if it overlaps an interval. A pathway is 'hit' by a set if any of
its genes is in the set; R_p = number of sets hitting pathway p (0-4).

Null: each set is shifted independently, chromosome by chromosome, by a random circular
offset (all intervals on a chromosome move together, keeping their number, sizes,
spacing and clustering). B permutations.

Reported for each threshold:
  P_specific  pathway-specific empirical P of R_p = 4
  P_global    P(number of pathways with R = 4 >= observed)
  P_any       P(at least one pathway with R = 4)
  P_FWER      min-P correction over all pathways (INRICH-style): fraction of permutations
              whose most extreme pathway-specific P is <= the observed one
usage: pathway_recurrence_permutation.py CORNCYC GFF3 [B] [seed]
"""
import csv, sys, os, numpy as np
from collections import defaultdict
cyc, gff = sys.argv[1], sys.argv[2]
B = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
T = "final/supp_tables"; OUT = "final/results/pathway_recurrence"; os.makedirs(OUT, exist_ok=True)
FL = 25000

# ---- genes ---------------------------------------------------------------------
L, genes = {}, defaultdict(list)                   # chr -> [(start, end, id)]
for line in open(gff):
    if line.startswith("##sequence-region"):
        _, c, a, b = line.split();  L[c] = int(b)
    elif not line.startswith("#"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "gene" and f[0] in [f"chr{i}" for i in range(1, 11)]:
            gid = [x[3:] for x in f[8].split(";") if x.startswith("ID=")][0]
            genes[f[0]].append((int(f[3]), int(f[4]), gid))
CH = [f"chr{i}" for i in range(1, 11)]
for c in CH:
    genes[c].sort(); L.setdefault(c, max(e for _, e, _ in genes[c]))

# ---- pathways ------------------------------------------------------------------
pw = defaultdict(set); pname = {}
for r in csv.DictReader(open(cyc), delimiter="\t"):
    g = (r.get("Gene-id") or "").strip()
    if g.upper().startswith("ZM00001EB"):
        g = "Zm00001eb" + g[9:].split("_")[0]; pw[r["Pathway-id"]].add(g); pname[r["Pathway-id"]] = r["Pathway-name"]
pids = sorted(pw); P = len(pids); pidx = {p: i for i, p in enumerate(pids)}
g2p = defaultdict(list)
for p, gs in pw.items():
    for g in gs: g2p[g].append(pidx[p])
# per chromosome: arrays for pathway genes only (all that matters for R_p)
G = {}
for c in CH:
    rows = [(s, e, gid) for s, e, gid in genes[c] if gid in g2p]
    s = np.array([x[0] for x in rows]); e = np.array([x[1] for x in rows])
    M = np.zeros((len(rows), P), bool)
    for i, (_, _, gid) in enumerate(rows): M[i, g2p[gid]] = True
    G[c] = (s, e, M, [x[2] for x in rows])
inpw = set(g2p)
print(f"pathways {P}; pathway genes on chr1-10 {sum(len(G[c][0]) for c in CH)}")

# ---- intervals -----------------------------------------------------------------
def snp_iv(rows):
    d = defaultdict(set)
    for c, p in rows: d[f"chr{int(c)}"].add(int(p))
    return {c: np.array(sorted((max(1, p - FL), p + FL) for p in d[c])).reshape(-1, 2) for c in CH}
gb = list(csv.DictReader(open("final/results/goodman_buckler_amino_acid_gwas/amino_gwas_v5_fdr_snps.csv")))
sn = list(csv.DictReader(open(f"{T}/SuppTable_S11_soilN_SNPs_env_relations.csv")))
def win(f):
    d = defaultdict(list)
    for r in csv.DictReader(open(f)):
        if r["top5"].upper() == "TRUE": d[r["chr"]].append((int(r["start_bp"]), int(r["end_bp"])))
    return {c: np.array(sorted(d[c])).reshape(-1, 2) for c in CH}
IC, J = win(f"{T}/SuppTable_S18_IC_Fst_windows.csv"), win(f"{T}/SuppTable_S19_J_Fst_windows.csv")
SETS = {
 "Bonf": {"GB": snp_iv([(r["Chr"], r["Pos_v5"]) for r in gb if float(r["min_P"]) < 0.05 / 4177796]),
          "soilN": snp_iv([(r["Chr"], r["Pos"]) for r in sn if float(r["min_P"]) < 0.05 / 495755]), "IC": IC, "J": J},
 "FDR":  {"GB": snp_iv([(r["Chr"], r["Pos_v5"]) for r in gb]),
          "soilN": snp_iv([(r["Chr"], r["Pos"]) for r in sn]), "IC": IC, "J": J},
}

def hits(iv, c, genes_only=False):
    """boolean pathway vector (or gene ids) for genes on chr c overlapping intervals iv"""
    s, e, M, ids = G[c]
    if len(iv) == 0 or len(s) == 0: return (np.zeros(P, bool) if not genes_only else [])
    o = np.argsort(iv[:, 0]); a = iv[o, 0]; b = np.maximum.accumulate(iv[o, 1])
    k = np.searchsorted(a, e, side="right") - 1          # last interval starting <= gene end
    ok = (k >= 0) & (b[np.clip(k, 0, None)] >= s)
    return [ids[i] for i in np.where(ok)[0]] if genes_only else M[ok].any(0)

def shift(iv, c):
    if len(iv) == 0: return iv
    Lc = L[c]; off = rng.integers(0, Lc); a = (iv[:, 0] + off) % Lc; w = iv[:, 1] - iv[:, 0]
    b = a + w; over = b > Lc
    out = np.vstack([np.c_[a, np.minimum(b, Lc)], np.c_[np.zeros(over.sum(), int), (b - Lc)[over]]])
    return out

def R_of(S, perm):
    R = np.zeros(P, int)
    for k, d in S.items():
        h = np.zeros(P, bool)
        for c in CH: h |= hits(shift(d[c], c) if perm else d[c], c)
        R += h
    return R

# ---- all-pathway scoring ---------------------------------------------------------
# For every pathway with R_obs >= 2 (hit by at least two sets):
#   P_specific = P(R_null >= R_obs)                       (all B permutations)
#   P_FWER     = min-P correction over all pathways, with independent halves: half 1
#                gives each pathway's tail probabilities, half 2 the distribution of the
#                smallest pathway P per permutation
#   P_BH       = Benjamini-Hochberg over all pathways (pathways hit by < 2 sets count as P = 1)
thr_list = sys.argv[5].split(",") if len(sys.argv) > 5 else list(SETS)
for thr in thr_list:
    S = SETS[thr]
    nint = {k: sum(len(v[c]) for c in CH) for k, v in S.items()}
    # observed, per set
    obs = {}
    for k, d in S.items():
        h = np.zeros(P, bool)
        for c in CH: h |= hits(d[c], c)
        obs[k] = h
    R0 = sum(obs[k].astype(int) for k in S)
    print(f"\n== {thr}: intervals {nint}; pathways with R=4: {int((R0==4).sum())}, R>=3: {int((R0>=3).sum())}")
    Rn = np.zeros((B, P), np.int8)
    for b in range(B):
        Rn[b] = R_of(S, True)
        if (b + 1) % 5000 == 0: print(f"   {b+1}/{B}", flush=True)
    ge = (Rn >= R0[None, :]).sum(0); pspec = (1 + ge) / (B + 1)
    h1 = B // 2; A, Bh = Rn[:h1], Rn[h1:]
    # tail table from half 1: T[p, r] = P(R >= r), r = 0..4
    T = np.ones((P, 5))
    for r in range(1, 5): T[:, r] = (1 + (A >= r).sum(0)) / (h1 + 1)
    p1_obs = T[np.arange(P), R0]
    tested = R0 >= 2
    pb = T[np.arange(P)[None, :], Bh]                       # (B/2, P) P value of each permuted R
    pb[Bh < 2] = 1.0              # family = every pathway; a permuted pathway counts once hit by >= 2 sets
    minp = pb.min(1)
    fwer = np.array([(1 + (minp <= p1_obs[i]).sum()) / (len(minp) + 1) for i in range(P)])
    # narrower family: four-way recurrence only (the question 'is any pathway hit by all four?')
    pb4 = np.where(Bh == 4, T[:, 4][None, :], 1.0); minp4 = pb4.min(1)
    fwer4 = np.array([(1 + (minp4 <= T[i, 4]).sum()) / (len(minp4) + 1) if R0[i] == 4 else np.nan for i in range(P)])
    idx = np.where(tested)[0]; m = P; order = idx[np.argsort(pspec[idx])]   # BH over all pathways (untested = 1)
    bh = np.ones(P); prev = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]; prev = min(prev, pspec[i] * m / rank); bh[i] = prev
    with open(f"{OUT}/pathway_permutation_all_{thr}.csv", "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["pathway_id", "pathway", "n_genes", "R_obs", "GB", "soilN", "IC", "J",
                    "mean_R_null", "P_specific", "P_FWER", "P_FWER_fourway", "P_BH", "GB_genes", "soilN_genes", "IC_genes", "J_genes"])
        allg = {k: set() for k in S}
        for k, d in S.items():
            for c in CH: allg[k] |= set(hits(d[c], c, genes_only=True))
        for i in sorted(idx, key=lambda i: pspec[i]):
            gs = pw[pids[i]]
            w.writerow([pids[i], pname[pids[i]], len(gs), int(R0[i])] + [int(obs[k][i]) for k in ("GB", "soilN", "IC", "J")] +
                       [round(float(Rn[:, i].mean()), 3), f"{pspec[i]:.3g}", f"{fwer[i]:.3g}", ("" if np.isnan(fwer4[i]) else f"{fwer4[i]:.3g}"), f"{bh[i]:.3g}"] +
                       [";".join(sorted(gs & allg[k])) for k in ("GB", "soilN", "IC", "J")])
    print(f"   tested pathways (R_obs >= 2): {len(idx)} of {P}")
    for i in sorted(idx, key=lambda i: pspec[i])[:10]:
        print(f"   R={R0[i]} K={len(pw[pids[i]]):>3} P={pspec[i]:.4f} FWER={fwer[i]:.3f} FWER4={fwer4[i]:.3f} BH={bh[i]:.3f}  {pname[pids[i]][:60]}")

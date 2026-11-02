#!/usr/bin/env python3
"""LD (r^2) around NMX (TATC, Zm00001eb399680) and ZmFd9 (Zm00001eb421930) from the
GBS landrace genotypes (all_maize2.bed, 495,755 SNPs, B73 v5), because the PHG
genotypes are not available locally.

Lead SNP per locus = the GBS SNP with the lowest soil-N (0-5 cm) MLM P within +/-25 kb
of the gene. r^2 = squared Pearson correlation of allele dosages over the 3,107
landraces used in the soil-N GWAS (missing genotypes dropped pairwise).

Writes to final/results/soil_nitrogen_gwas/LD_GBS_regional/:
  <locus>_gbs_r2.csv     pos, r2 with the lead, MLM/BLINK/FarmCPU P (GBS soil-N GWAS)
  <locus>_phg_glm.csv    pos, P from the PHG soil-N GLM, r2 where the position is a GBS SNP
  leads.csv
usage: python3 LD_GBS_regional_NMX_ZmFd9.py <repo> <bed_prefix> <phg_glm.txt>
"""
import sys, os, numpy as np, pandas as pd
repo, bed, phg = sys.argv[1:4]
OUT = os.path.join(repo, "final/results/soil_nitrogen_gwas/LD_GBS_regional"); os.makedirs(OUT, exist_ok=True)
SN = os.path.join(repo, "final/results/soil_nitrogen_gwas")
FLANK = 1_000_000
LOCI = [("NMX", "Zm00001eb399680", 9, 151731980, 151737594),
        ("ZmFd9", "Zm00001eb421930", 10, 116703178, 116706315)]

fam = [l.split()[1] for l in open(bed + ".fam")]; n = len(fam)
keepset = set(open(os.path.join(SN, "LD/landrace_N_keep.txt")).read().split())
keep = np.array([i for i, t in enumerate(fam) if t in keepset]); print("accessions", len(keep))
bim = pd.read_csv(bed + ".bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
nb = (n + 3) // 4
raw = np.memmap(bed + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), nb)
lut = np.array([2, np.nan, 1, 0])
def geno(idx):
    b = np.asarray(raw[idx])
    g = np.stack([(b >> s) & 3 for s in (0, 2, 4, 6)], axis=-1).reshape(len(idx), -1)[:, :n]
    return lut[g][:, keep]
def r2(lead, G):
    out = np.empty(G.shape[0])
    for i, x in enumerate(G):
        m = ~np.isnan(lead) & ~np.isnan(x)
        a, b = lead[m], x[m]
        out[i] = np.nan if a.std() == 0 or b.std() == 0 else np.corrcoef(a, b)[0, 1] ** 2
    return out

scans = {}
for mod in ("MLM", "BLINK", "FarmCPU"):
    s = pd.read_csv(os.path.join(SN, f"N_0_5_full/raw_GWAS_{mod}_3PC_N_0_5.txt"), sep="\t", usecols=["Chr", "Pos", "P.value"])
    scans[mod] = s.rename(columns={"P.value": f"P_{mod}"})
phgdf = pd.read_csv(phg, sep="\t", usecols=["Chr", "Pos", "p"])
phgdf["chr"] = phgdf.Chr.str.replace("chr", "").astype(int)

leads = []
for name, gene, c, gs, ge in LOCI:
    mlm = scans["MLM"]; near = mlm[(mlm.Chr == c) & (mlm.Pos >= gs - 25000) & (mlm.Pos <= ge + 25000)]
    lead = int(near.loc[near.P_MLM.idxmin(), "Pos"])
    w = bim[(bim.chr == c) & (bim.pos >= lead - FLANK) & (bim.pos <= lead + FLANK)]
    G = geno(w.index.values); li = int(np.where(w.pos.values == lead)[0][0])
    d = pd.DataFrame({"pos": w.pos.values, "r2": r2(G[li], G)})
    for mod, s in scans.items():
        d = d.merge(s[s.Chr == c][["Pos", f"P_{mod}"]].rename(columns={"Pos": "pos"}), on="pos", how="left")
    d.to_csv(os.path.join(OUT, f"{name}_gbs_r2.csv"), index=False)
    p = phgdf[(phgdf.chr == c) & (phgdf.Pos >= lead - FLANK) & (phgdf.Pos <= lead + FLANK)]
    p = p.groupby("Pos", as_index=False).p.min().rename(columns={"Pos": "pos"}).merge(d[["pos", "r2"]], on="pos", how="left")
    p.to_csv(os.path.join(OUT, f"{name}_phg_glm.csv"), index=False)
    phg_lead = int(p.loc[p[(p.pos >= gs - 25000) & (p.pos <= ge + 25000)].p.idxmin(), "pos"])
    leads.append(dict(locus=name, gene=gene, chr=c, gene_start=gs, gene_end=ge, gbs_lead=lead,
                      gbs_lead_P_MLM=float(near.P_MLM.min()), n_gbs=len(d),
                      phg_lead=phg_lead, n_phg=len(p), n_phg_with_r2=int(p.r2.notna().sum())))
    L = leads[-1]; L["n_gbs_r2_ge08"] = int((d.r2 >= 0.8).sum())
    print(name, "GBS lead", lead, "P_MLM %.2e" % L["gbs_lead_P_MLM"], "GBS SNPs", len(d), "r2>=0.8:", L["n_gbs_r2_ge08"],
          "| PHG SNPs", len(p), "with r2", L["n_phg_with_r2"], "PHG lead", phg_lead)
pd.DataFrame(leads).to_csv(os.path.join(OUT, "leads.csv"), index=False)

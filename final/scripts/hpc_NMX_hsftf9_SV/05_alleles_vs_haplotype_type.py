#!/usr/bin/env python3
"""Step 5 (login node, seconds): does the allele at the NMX lead SNP and the 152 Mb SNP follow
the structural type of the NMX-hsftf9 region?
Type per genome from NMX_hsftf9_per_genome.tsv (step 3):
  B73-like   >= 8 of 9 B73 genes between, distance < 300 kb
  ~250kb     5 of 9 genes, 240-260 kb (PDJ, Il14H, Ms71, P39 ...)
  ~310kb     <= 4 of 9 genes, 280-360 kb (PT, CML459, CML69, HP301 ...; also teosintes)
  long       maize, distance >= 700 kb (Mo17, tropical NAM lines, ZAP, CML530 ...)
  other      anything else; NA = split contigs, misassembly, or non-Zea genome
GBS alleles (B73 v5 forward strand): S9_151742636 G/A, S9_151983481 C/T, S9_151984305 A/C,
S9_151992800 G/C. In the 3,107 landraces, S9_151742636 GG sites have median soil N 313 and
elevation 388 m, AA sites 227 and 1,885 m; S9_151983481 CC 310 / 390 m, TT 227 / 1,911 m.
writes OUT/snp_alleles_by_type.tsv and prints a cross-table per SNP
usage: python3 05_alleles_vs_haplotype_type.py [OUT]
"""
import sys, os, glob, csv, collections
OUT = sys.argv[1] if len(sys.argv) > 1 else "/share/maize/ntanduk/NMX_hsftf9_SV"
ty = {}
for r in csv.DictReader(open(os.path.join(OUT, "NMX_hsftf9_per_genome.tsv")), delimiter="\t"):
    d = float(r["distance_bp"]) if r["distance_bp"] not in ("None", "NA", "") else None
    g = r["B73_genes_between_found"]
    k = int(g.split("/")[0]) if "/" in g else -1
    if d is None or not r["genome"].startswith(("Zm-", "Zv-", "Zx-", "Zh-", "Zl-", "Zd-", "Zn-")): t = "NA"
    elif d > 2_000_000: t = "NA"                                   # misassembly (e.g. Xu178)
    elif r["genome"].startswith("Zm-") and 700_000 <= d: t = "long"
    elif k >= 8 and d < 300_000: t = "B73-like"
    elif 280_000 <= d <= 360_000 and k <= 4: t = "~310kb"
    elif 240_000 <= d < 260_000 and k == 5: t = "~250kb"
    else: t = "other"
    ty[r["genome"]] = t
rows = []
for f in sorted(glob.glob(os.path.join(OUT, "snp_alleles", "*.tsv"))):
    for r in csv.DictReader(open(f), delimiter="\t"):
        r["type"] = ty.get(r["genome"], "NA"); rows.append(r)
with open(os.path.join(OUT, "snp_alleles_by_type.tsv"), "w") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
for snp in ["snp_151742636", "snp_151983481", "snp_151984305", "snp_151992800"]:
    t = collections.Counter((r["type"], r["allele"]) for r in rows if r["snp"] == snp)
    alle = sorted({a for _, a in t}); types = ["B73-like", "~250kb", "~310kb", "other", "long", "NA"]
    print(f"\n{snp}\n" + "type\t" + "\t".join(alle))
    for k in types: print(k + "\t" + "\t".join(str(t.get((k, a), 0)) for a in alle))
print("\nper genome (NMX lead / 152 Mb SNP):")
by = collections.defaultdict(dict)
for r in rows: by[r["genome"]][r["snp"]] = r["allele"]
for g in sorted(by, key=lambda g: (ty.get(g, ""), g)):
    print(f"{ty.get(g,'NA'):9s}\t{g}\t{by[g].get('snp_151742636')}\t{by[g].get('snp_151983481')}")

#!/usr/bin/env python3
"""
Year-to-year consistency of amino acid line values in the Goodman-Buckler panel.

For every trait, each line's value is averaged over the two blocks within a year,
and the 2017 and 2018 line values are compared with a Spearman rank correlation
(average ranks for ties) across the lines measured in both years. A high rho means
lines that were high in one year were also high in the other -- the direct check
behind "differences among lines were consistent between years"; the line-mean H2
in the same tables is the model-based summary of the same property.

Inputs : data/FAA_Goodman_Buckler.csv, data/PBAA_Goodman_Buckler.csv (plot records:
         taxa, year, Block, traits)
Outputs: adds two columns, in place, to
         final/supp_tables/SuppTable_S3_FAA_abundance_CV.csv   (20 FAAs, matched on AA code)
         final/supp_tables/SuppTable_S4_heritability.csv       (all 116 FAA/PBAA traits)
           rho_2017_vs_2018  Spearman rho of line values, 2017 vs 2018
           n_lines_both_years  lines measured in both years
         and prints per-pool summaries. Re-running overwrites the two columns.
"""
import csv, os, statistics as st, collections
R = os.environ.get("GN_REPO", os.path.expanduser("~/Documents/Github/Grain_N_remobilization"))

def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r

def year_rho(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    traits = [k for k in rows[0] if k not in ("taxa", "year", "Block")]
    years = sorted({r["year"] for r in rows}); assert len(years) == 2, years
    out = {}
    for t in traits:
        d = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            try: x = float(r[t])
            except (ValueError, TypeError): continue
            if x == x: d[r["taxa"]][r["year"]].append(x)
        pairs = [(st.mean(v[years[0]]), st.mean(v[years[1]])) for v in d.values() if years[0] in v and years[1] in v]
        if len(pairs) < 10: out[t] = ("", len(pairs)); continue
        a, b = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
        out[t] = (round(st.correlation(a, b), 3), len(pairs))
    return out

FAA  = year_rho(f"{R}/data/FAA_Goodman_Buckler.csv")
PBAA = year_rho(f"{R}/data/PBAA_Goodman_Buckler.csv")

def add_cols(path, key):
    rows = list(csv.DictReader(open(path))); fields = [f for f in rows[0] if f not in ("rho_2017_vs_2018", "n_lines_both_years")]
    for r in rows:
        rho, n = key(r)
        r["rho_2017_vs_2018"], r["n_lines_both_years"] = rho, n
    w = csv.DictWriter(open(path, "w", newline=""), fieldnames=fields + ["rho_2017_vs_2018", "n_lines_both_years"])
    w.writeheader(); w.writerows(rows); return rows

s3 = add_cols(f"{R}/final/supp_tables/SuppTable_S3_FAA_abundance_CV.csv", lambda r: FAA.get(r["AA"], ("", "")))
s4 = add_cols(f"{R}/final/supp_tables/SuppTable_S4_heritability.csv",
              lambda r: (FAA if r["pool"] == "FAA" else PBAA).get(r["trait"], ("", "")))

PB_ABS = ["Ala","Arg","Asx","Glx","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Tyr","Val"]
f = {r["AA"]: r["rho_2017_vs_2018"] for r in s3}
print("FAA (20): median rho %.2f, range %.2f-%.2f; n lines %s; Asn %.2f, Pro %.2f" % (
      st.median(f.values()), min(f.values()), max(f.values()), s3[0]["n_lines_both_years"], f["N"], f["P"]))
p = {t: PBAA[t][0] for t in PB_ABS}
print("PBAA (15): median rho %.2f, range %.2f-%.2f; n lines %s" % (st.median(p.values()), min(p.values()), max(p.values()), PBAA["Ala"][1]))
print("S4 rows without rho:", [r["trait"] for r in s4 if r["rho_2017_vs_2018"] == ""])

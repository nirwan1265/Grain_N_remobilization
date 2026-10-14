#!/usr/bin/env python3
"""
Trait-matched candidate genes from the Goodman-Buckler amino acid GWAS.

A candidate gene (Supplementary Table S6: B73 RefGen_v5 gene within 25 kb of a
5% FDR SNP) is called trait-matched when its annotation (gene name + PANTHER
family) places it in
  * the metabolism of the amino acid the phenotype reports on
    ("enzyme: same amino acid"; numerator amino acid, both amino acids of an
    amino acid:amino acid ratio, or all members of a family composite), or
  * the metabolism of another amino acid of the same biosynthetic family
    ("enzyme: same family"), or
  * amino acid / ammonium / nitrate transport (UmamiT/WAT1-related, AAP, LHT,
    CAT, AVT, NRT1/PTR (NPF), AMT; "transporter (amino acid / N)"), or
  * primary N assimilation (GS, GOGAT; "N assimilation").
Enzyme names are matched with the regular expressions below (curated per amino
acid); methyltransferases, PEPC kinases and PAL are excluded as false matches.
This is a descriptive annotation, not an enrichment test: about 9% of all
annotated maize genes carry an amino acid or N annotation.

Inputs : final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv
         final/results/GO/GO_universe_genomewide.txt (PANTHER, all v5 genes)
Output : final/results/goodman_buckler_amino_acid_gwas/amino_gwas_pathway_matched_genes.csv
"""
import csv, re, os, collections
R = os.environ.get("GN_REPO", os.path.expanduser("~/Documents/Github/Grain_N_remobilization"))
U = {r["GeneID"]: r for r in csv.DictReader(open(f"{R}/final/results/GO/GO_universe_genomewide.txt"), delimiter="\t")}
S6 = list(csv.DictReader(open(f"{R}/final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv")))
P={ # amino acid -> own metabolism (synthesis, catabolism) patterns
'N':r'asparagine synth|asparaginase|omega-amidase|nitrilase',
'D':r'aspartate aminotransferase|aspartate transaminase|aspartate kinase|aspartokinase|aspartate-semialdehyde|asparaginase|asparagine synth|aspartate carbamoyl|phosphoenolpyruvate carboxylase',
'K':r'aspartate kinase|aspartokinase|aspartate-semialdehyde|dihydrodipicolinate|4-hydroxy-tetrahydrodipicolinate|diaminopimelate|lysine-ketoglutarate|saccharopine|lysine decarboxylase|LL-diaminopimelate',
'M':r'cystathionine|homocysteine|methionine synthase|5-methyltetrahydropteroyltriglutamate|S-adenosylmethionine synth|methionine adenosyltransferase|homoserine|methionine gamma-lyase|S-methylmethionine',
'T':r'threonine synthase|homoserine|threonine aldolase|threonine dehydratase|threonine deaminase',
'I':r'threonine dehydratase|threonine deaminase|acetolactate synth|ketol-acid reductoisomerase|dihydroxy-acid dehydratase|branched-chain(-amino-acid)? aminotransferase|branched chain amino acid aminotransferase',
'L':r'isopropylmalate|acetolactate synth|ketol-acid reductoisomerase|dihydroxy-acid dehydratase|branched-chain(-amino-acid)? aminotransferase|isovaleryl-CoA|methylcrotonoyl',
'V':r'acetolactate synth|ketol-acid reductoisomerase|dihydroxy-acid dehydratase|branched-chain(-amino-acid)? aminotransferase|3-hydroxyisobutyr|methylmalonate',
'A':r'alanine aminotransferase|alanine transaminase|alanine-glyoxylate|alanine dehydrogenase',
'G':r'glycine hydroxymethyltransferase|serine hydroxymethyltransferase|glycine dehydrogenase|glycine cleavage|aminomethyltransferase|glyoxylate aminotransferase|alanine-glyoxylate|serine-glyoxylate',
'S':r'phosphoglycerate dehydrogenase|phosphoserine|serine hydroxymethyltransferase|glycine hydroxymethyltransferase|serine acetyltransferase|serine-glyoxylate|serine racemase|serine dehydratase',
'C':r'serine acetyltransferase|cysteine synthase|O-acetylserine|cysteine desulfurase|cystathionine|sulfite reductase|APS reductase|adenylyl-sulfate',
'E':r'glutamate dehydrogenase|glutamate synthase|glutamate decarboxylase|glutamine synthetase|aspartate aminotransferase|glutamate-glyoxylate|GS catalytic|glutamate carboxypeptidase',
'Q':r'glutamine synthetase|GS catalytic|glutamate synthase|glutaminase|asparagine synth|glutamine amidotransferase|carbamoyl-phosphate synth',
'H':r'histidinol|imidazoleglycerol|ATP phosphoribosyltransferase|phosphoribosyl-AMP|phosphoribosyl-ATP|histidine decarboxylase',
'P':r'pyrroline-5-carboxylate|proline dehydrogenase|proline oxidase|ornithine aminotransferase|ornithine-oxo-acid',
'R':r'acetylglutamate|acetylornithine|ornithine carbamoyl|argininosuccinate|arginase|arginine decarboxylase|carbamoyl-phosphate synth|ornithine decarboxylase|glutamate N-acetyltransferase',
'F':r'chorismate mutase|arogenate|prephenate|phenylalanine ammonia|DAHP|3-deoxy-7-phosphoheptulonate|shikimate|EPSP|3-dehydroquinate|chorismate synthase',
'Y':r'chorismate mutase|arogenate|prephenate|tyrosine aminotransferase|DAHP|3-deoxy-7-phosphoheptulonate|shikimate|EPSP|3-dehydroquinate|chorismate synthase|hydroxyphenylpyruvate',
'W':r'anthranilate|tryptophan synth|indole-3-glycerol|phosphoribosylanthranilate|DAHP|3-deoxy-7-phosphoheptulonate|shikimate|EPSP|3-dehydroquinate|chorismate synthase|tryptophan aminotransferase',
}
TR=r'amino acid permease|amino acid transporter|aminoacid transporter|lysine histidine transporter|cationic amino acid|proline transporter|NRT1|PTR FAMILY|NPF|WAT1-related|UmamiT|ammonium transporter|nitrate transporter|bidirectional amino acid|GABA transporter|amino acid export'
TR = (r"amino acid permease|amino acid transporter|aminoacid transporter|lysine histidine transporter|"
      r"cationic amino acid|proline transporter|NRT1|PTR FAMILY|NPF|WAT1-related|UmamiT|ammonium transporter|"
      r"nitrate transporter|bidirectional amino acid|GABA transporter|amino acid export")
FAM = {"Asp": set("NDKMTI"), "Glu": set("EQPRH"), "Ser": set("SGC"), "Aro": set("WFY"), "Pyr": set("ALVI")}
COMPS = {"IMNTDK": set("IMNTDK"), "EHPRQ": set("EHPRQ"), "SGC": set("SGC"), "WFY": set("WFY"), "LAV": set("LAV"), "IVL": set("IVL")}
def num(ph):
    p = ph.split(".")
    if p[0] in COMPS: return COMPS[p[0]]
    if len(p) == 2 and len(p[0]) == 1 and len(p[1]) == 1: return {p[0], p[1]}
    if len(p[0]) == 1: return {p[0]}
    return set()
def fams(a): return {f for f, s in FAM.items() if a & s}
out = []
for r in S6:
    g = r["GeneID"]; u = U.get(g, {}); t = u.get("GeneName", "") + " " + u.get("Family_Subfamily", "")
    n = num(r["Phenotype"]); ga = set()
    if not re.search(r"methyltransferase|carboxylase kinase|ammonia-lyase", t, re.I):
        ga = {a for a, pat in P.items() if re.search(pat, t, re.I)}
    if ga & n: cls = "enzyme: same amino acid"
    elif ga and fams(ga) & fams(n): cls = "enzyme: same family"
    elif re.search(TR, t, re.I): cls = "transporter (amino acid / N)"
    elif re.search(r"glutamine synthetase|GS catalytic|glutamate synthase", t, re.I): cls = "N assimilation"
    else: continue
    out.append({"class": cls, "GeneID": g, "GeneName": u.get("GeneName", "").split(";")[0],
                "PANTHER": u.get("Family_Subfamily", ""), "Phenotype": r["Phenotype"], "SNP": r["SNP"],
                "Chr": r["Chr"], "Pos_v5": r["Pos_v5"], "Relation": r["Relation"],
                "Distance_bp": r["Distance_to_Gene_bp"], "min_P": r["min_P"], "models": r["models_sig"],
                "replicated": r["replicated"]})
out.sort(key=lambda x: (x["class"], float(x["min_P"])))
f = f"{R}/final/results/goodman_buckler_amino_acid_gwas/amino_gwas_pathway_matched_genes.csv"
w = csv.DictWriter(open(f, "w", newline=""), fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print(collections.Counter(x["class"] for x in out), "->", f)

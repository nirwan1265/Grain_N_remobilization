from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP_DIR = ROOT / "tables" / "supplementary"
ANNOT_DIR = ROOT / "tables" / "annotation"


CURATED_GENES = [
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb005710",
        "Interpretive_note": "Soil N omega-amidase/NIT2-like candidate linked to Asn/Gln-associated metabolism.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb220040",
        "Interpretive_note": "Soil N aspartate kinase candidate at the entry point to the Asp-family pathway.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb201770",
        "Interpretive_note": "Amino-acid GWAS transaminase candidate linking amino-acid balance to 2-oxoglutarate-dependent chemistry.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb257910",
        "Interpretive_note": "Methionine-associated GWAS candidate consistent with Asp-family transamination.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb214360",
        "Interpretive_note": "Indian Chief top-5% FST candidate for the step immediately after aspartate kinase.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb094670",
        "Interpretive_note": "Jarvis top-5% FST aspartate kinase candidate discussed as ask2-like.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb152450",
        "Interpretive_note": "Jarvis aspartate aminotransferase / glutamate-oxaloacetate transaminase candidate.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb231540",
        "Interpretive_note": "Jarvis L-aspartate oxidase candidate linking Asp metabolism to NAD biosynthesis.",
    },
    {
        "Functional_category": "Asp/Asn metabolism and C-N interconversion",
        "GeneID": "Zm00001eb294790",
        "Interpretive_note": "Indian Chief threonine synthase candidate in the downstream Asp-family branch.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb060520",
        "Interpretive_note": "Key bridge between Total_N GWAS and Indian Chief FST; glutamate dehydrogenase-related.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb214770",
        "Interpretive_note": "Indian Chief all-three candidate annotated as NADP-specific glutamate dehydrogenase.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb427810",
        "Interpretive_note": "Indian Chief top-5% FST glutamate dehydrogenase-related candidate.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb162300",
        "Interpretive_note": "Indian Chief top-5% FST glutamate dehydrogenase-related candidate.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb162900",
        "Interpretive_note": "Indian Chief top-5% FST glutamine synthetase candidate.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb217220",
        "Interpretive_note": "Indian Chief glutamine amidotransferase-domain candidate; broader glutamine-dependent biosynthesis.",
    },
    {
        "Functional_category": "Glu/Gln metabolism and glutamine-dependent reactions",
        "GeneID": "Zm00001eb232460",
        "Interpretive_note": "Jarvis glutamine amidotransferase-domain candidate interpreted cautiously.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb287980",
        "Interpretive_note": "Glutamine-associated NRT1/PTR-family transporter from the amino-acid GWAS.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb277160",
        "Interpretive_note": "Asp/Asn-related NRT1/PTR-family transporter associated with multiple amino-acid ratio traits.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb429100",
        "Interpretive_note": "AVT3A-like amino-acid transporter associated with glutamate and related traits.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb214890",
        "Interpretive_note": "Indian Chief all-three ammonium transporter candidate.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb287950",
        "Interpretive_note": "Indian Chief all-three NRT1/PTR-family transporter candidate.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb132750",
        "Interpretive_note": "Jarvis amino-acid transporter candidate with lysine/histidine transporter-like annotation.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb159580",
        "Interpretive_note": "Jarvis putative GABA transporter candidate.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb389560",
        "Interpretive_note": "Jarvis GLUTAMINE DUMPER 5 candidate consistent with intracellular amino-acid movement.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb393130",
        "Interpretive_note": "Soil N trafficking candidate annotated as vacuolar protein sorting-associated protein 52A.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb130210",
        "Interpretive_note": "Soil N regulatory/trafficking candidate annotated as a RING/FYVE/PHD zinc-finger ubiquitin ligase.",
    },
    {
        "Functional_category": "Transport and intracellular trafficking",
        "GeneID": "Zm00001eb393120",
        "Interpretive_note": "Soil N neighboring nitrilase-associated candidate in the VPS52A interval.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb372060",
        "Interpretive_note": "Indian Chief all-three autophagy-related gene and the main proline-overlap candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb287690",
        "Interpretive_note": "Indian Chief all-three autophagy-related candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb352300",
        "Interpretive_note": "Indian Chief all-three autophagy-related candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb214780",
        "Interpretive_note": "Indian Chief senescence-specific cysteine protease candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb082660",
        "Interpretive_note": "Amino-acid GWAS macroautophagy candidate from the broader Total trait.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb322660",
        "Interpretive_note": "Amino-acid GWAS senescence regulator candidate from the broader Total trait.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb293380",
        "Interpretive_note": "Aspartate GWAS RING-type E3 ubiquitin ligase candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb283510",
        "Interpretive_note": "Aspartate GWAS Derlin-1-like ERAD candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb429920",
        "Interpretive_note": "D.IMNTDK GWAS candidate annotated for autophagosome assembly.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb210360",
        "Interpretive_note": "Jarvis all-three OTU-domain deubiquitinating enzyme candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb146560",
        "Interpretive_note": "Jarvis all-three NPL4-like candidate involved in ubiquitin-dependent protein catabolism.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb175050",
        "Interpretive_note": "Jarvis all-three ubiquitin carboxyl-terminal hydrolase candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb210720",
        "Interpretive_note": "Jarvis all-three 26S proteasome regulatory subunit candidate.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb021640",
        "Interpretive_note": "Total_N / Total_PBAA glutathione S-transferase cluster in the proxy N interval.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb021650",
        "Interpretive_note": "Total_N / Total_PBAA glutathione S-transferase cluster in the proxy N interval.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb021690",
        "Interpretive_note": "Total_N / Total_PBAA glutathione S-transferase cluster in the proxy N interval.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb131830",
        "Interpretive_note": "Jarvis glutamine/glutathione-related candidate noted in the broader top-5% FST set.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb131840",
        "Interpretive_note": "Jarvis glutamine/glutathione-related candidate noted in the broader top-5% FST set.",
    },
    {
        "Functional_category": "Proteostasis, autophagy, senescence, and redox",
        "GeneID": "Zm00001eb210730",
        "Interpretive_note": "Jarvis glutamine/glutathione-related candidate noted in the broader top-5% FST set.",
    },
]


GO_COLS = [
    "SourceID",
    "GeneID",
    "Description",
    "Family_Subfamily",
    "Protein_Class",
    "GO_MF",
    "GO_BP",
    "GO_CC",
]

ANNOTATION_OVERRIDES = {
    "Zm00001eb201770": "Tyrosine aminotransferase",
    "Zm00001eb257910": "Aspartate aminotransferase / GOT-like protein",
    "Zm00001eb214770": "NADP-specific glutamate dehydrogenase",
    "Zm00001eb287950": "NRT1/PTR-family transporter",
    "Zm00001eb132750": "Lysine/histidine transporter-like protein",
    "Zm00001eb393120": "Nitrilase-associated protein",
    "Zm00001eb322660": "Putative senescence regulator",
    "Zm00001eb429920": "Autophagosome assembly-related protein",
}


def load_go_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=GO_COLS, dtype=str)
    for col in GO_COLS:
        df[col] = df[col].fillna("")
    return df


def first_nonempty(*values: str) -> str:
    for value in values:
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def short_desc(value: str) -> str:
    if not isinstance(value, str):
        return ""
    value = value.strip()
    if not value:
        return ""
    return value.split(";")[0].strip()


supp3 = pd.read_csv(SUPP_DIR / "SuppTable_3_GWAS_annotation_soilN_amino_with_GO.csv", dtype=str).fillna("")
amino_summary = pd.read_csv(SUPP_DIR / "SuppTable_amino_gwas_gene_summary_25kb.csv", dtype=str).fillna("")
ic_top5 = load_go_table(ANNOT_DIR / "IC_GO_Fst_5pc.txt")
j_top5 = load_go_table(ANNOT_DIR / "J_GO_Fst_5pc.txt")
ic_all3 = load_go_table(ANNOT_DIR / "IC_GO_all3_pop_test.txt")
j_all3 = load_go_table(ANNOT_DIR / "J_GO_all3_pop_test.txt")


supp3_by_gene = {gene: sub.iloc[0] for gene, sub in supp3.groupby("GeneID", sort=False)}
amino_by_gene = {gene: sub.iloc[0] for gene, sub in amino_summary.groupby("GeneID", sort=False)}
ic_top5_by_gene = {gene: sub.iloc[0] for gene, sub in ic_top5.groupby("GeneID", sort=False)}
j_top5_by_gene = {gene: sub.iloc[0] for gene, sub in j_top5.groupby("GeneID", sort=False)}
ic_all3_by_gene = {gene: sub.iloc[0] for gene, sub in ic_all3.groupby("GeneID", sort=False)}
j_all3_by_gene = {gene: sub.iloc[0] for gene, sub in j_all3.groupby("GeneID", sort=False)}


def phenotype_split(value: str) -> list[str]:
    return [x.strip() for x in value.split(";") if x.strip()]


def soiln_entry(gene_id: str) -> str:
    row = supp3_by_gene.get(gene_id)
    if row is None:
        return ""
    phenos = phenotype_split(row["Phenotypes"])
    if "SoilN" not in phenos:
        return ""
    pvals = phenotype_split(row["pvalues"])
    score = pvals[phenos.index("SoilN")] if len(pvals) == len(phenos) else ""
    return f"SoilN ({score})" if score else "SoilN"


def amino_entry(gene_id: str) -> str:
    row = amino_by_gene.get(gene_id)
    if row is not None:
        count = row["Phenotypes_count"]
        phenos = row["Phenotypes"]
        return f"{phenos} (n={count})" if phenos else ""
    row = supp3_by_gene.get(gene_id)
    if row is None:
        return ""
    phenos = [p for p in phenotype_split(row["Phenotypes"]) if p != "SoilN"]
    if not phenos:
        return ""
    return f"{';'.join(phenos)} (n={len(phenos)})"


def scan_entry(gene_id: str, top5_lookup: dict, all3_lookup: dict) -> str:
    if gene_id in all3_lookup:
        return "all_three"
    if gene_id in top5_lookup:
        return "top5_FST"
    return ""


def preferred_label(gene_id: str) -> str:
    row3 = supp3_by_gene.get(gene_id)
    rowa = amino_by_gene.get(gene_id)
    label = first_nonempty(
        row3["GeneSymbol"] if row3 is not None else "",
        rowa["GeneSymbol"] if rowa is not None else "",
    )
    if label == gene_id:
        label = ""
    return label


def annotation_text(gene_id: str) -> str:
    if gene_id in ANNOTATION_OVERRIDES:
        return ANNOTATION_OVERRIDES[gene_id]
    row3 = supp3_by_gene.get(gene_id)
    rowa = amino_by_gene.get(gene_id)
    row_ic_all3 = ic_all3_by_gene.get(gene_id)
    row_j_all3 = j_all3_by_gene.get(gene_id)
    row_ic_top5 = ic_top5_by_gene.get(gene_id)
    row_j_top5 = j_top5_by_gene.get(gene_id)
    return first_nonempty(
        short_desc(row3["Family_Subfamily"] if row3 is not None else ""),
        short_desc(row3["GeneSymbol"] if row3 is not None else ""),
        short_desc(row_ic_all3["Description"] if row_ic_all3 is not None else ""),
        short_desc(row_j_all3["Description"] if row_j_all3 is not None else ""),
        short_desc(row_ic_top5["Description"] if row_ic_top5 is not None else ""),
        short_desc(row_j_top5["Description"] if row_j_top5 is not None else ""),
        short_desc(rowa["GeneSymbol"] if rowa is not None else ""),
    )


rows = []
for item in CURATED_GENES:
    gene_id = item["GeneID"]
    rows.append(
        {
            "Functional_category": item["Functional_category"],
            "GeneID": gene_id,
            "Preferred_label": preferred_label(gene_id),
            "Annotation": annotation_text(gene_id),
            "SoilN_GWAS": soiln_entry(gene_id),
            "Amino_acid_GWAS": amino_entry(gene_id),
            "Indian_Chief_scan": scan_entry(gene_id, ic_top5_by_gene, ic_all3_by_gene),
            "Jarvis_scan": scan_entry(gene_id, j_top5_by_gene, j_all3_by_gene),
            "Interpretive_note": item["Interpretive_note"],
        }
    )


out_df = pd.DataFrame(rows)
out_df["Evidence_layers"] = out_df[
    ["SoilN_GWAS", "Amino_acid_GWAS", "Indian_Chief_scan", "Jarvis_scan"]
].apply(lambda s: sum(bool(str(v).strip()) for v in s), axis=1)
out_df = out_df[
    [
        "Functional_category",
        "GeneID",
        "Preferred_label",
        "Annotation",
        "SoilN_GWAS",
        "Amino_acid_GWAS",
        "Indian_Chief_scan",
        "Jarvis_scan",
        "Evidence_layers",
        "Interpretive_note",
    ]
]

out_path = SUPP_DIR / "SuppTable_11_integrative_candidate_gene_synthesis.csv"
out_df.to_csv(out_path, index=False)

print(f"Wrote {len(out_df)} rows to {out_path}")
print(out_df.groupby("Functional_category").size().to_string())

#!/usr/bin/env python3
"""Integrate external maize N-stress RNA-seq DEGs with manuscript candidate sets."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

RNA_WORKBOOK = ROOT / "data" / "List of significant DEG (including FPKM).xlsx"
AMINO_GWAS_FILE = ROOT / "tables" / "supplementary" / "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv"
AMINO_SUMMARY_FILE = ROOT / "tables" / "supplementary" / "SuppTable_amino_gwas_gene_summary_25kb.csv"
SOILN_ANNOTATION_FILE = ROOT / "tables" / "supplementary" / "SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv"
IC_FST_FILE = ROOT / "tables" / "supplementary" / "SuppTable_IC_GenWin_windows_Fst_Tajima_pi_genes.csv"
J_FST_FILE = ROOT / "tables" / "supplementary" / "SuppTable_J_GenWin_windows_Fst_Tajima_pi_genes.csv"
ALL_AMINO_ANNOTATION_FILE = ROOT / "tables" / "annotation" / "Grain_aminoacid_all.txt"
IC_FST_ANNOTATION_FILE = ROOT / "tables" / "annotation" / "IC_GO_Fst_5pc.txt"
J_FST_ANNOTATION_FILE = ROOT / "tables" / "annotation" / "J_GO_Fst_5pc.txt"

TABLE_DIR = ROOT / "tables" / "supplementary"
FIG_DIR = ROOT / "Figs" / "Supplementary"

FOCAL_AMINO_TRAITS = {
    "D",
    "E",
    "N",
    "P",
    "Q",
    "Total_N",
    "Total_PBAA",
    "IMNTDK",
    "D.IMNTDK",
    "N.E",
    "N.Total",
    "Q.E",
}

CONTRAST_ORDER = ["2022 V6", "2023 V6", "2022 V10", "2023 V10"]
PATHWAY_ORDER = [
    "GS/GOGAT + Glu/Gln",
    "Asp/Asn",
    "Proline Metabolism",
    "Amino Acid Transport",
    "Mineral N Transport",
    "Carbon + Photosynthesis",
    "Stress Response",
    "Senescence + Remobilization",
]
CORE_FIGURE_PATHWAY_ORDER = [
    "GS/GOGAT + Glu/Gln",
    "Asp/Asn",
    "Proline Metabolism",
    "Amino Acid Transport",
    "Mineral N Transport",
    "Senescence + Remobilization",
]
CORE_N_REMOBILIZATION_CATEGORIES = {
    "GS/GOGAT + Glu/Gln",
    "Asp/Asn",
    "Proline Metabolism",
    "Amino Acid Transport",
    "Mineral N Transport",
    "Senescence + Remobilization",
}
CANDIDATE_SOURCE_ORDER = [
    "RNA-seq DEGs",
    "All amino GWAS",
    "SoilN GWAS",
    "IC FST 5pc",
    "J FST 5pc",
]
GENERIC_ANNOTATION_COLUMNS = [
    "ExternalID",
    "GeneID",
    "GeneName",
    "Family_Subfamily",
    "Protein_Class",
    "GO_MF",
    "GO_BP",
    "GO_CC",
]

PATHWAY_KEYWORDS = {
    "GS/GOGAT + Glu/Gln": [
        (r"glutamine synthetase", "glutamine synthetase"),
        (r"glutamate synthase", "glutamate synthase"),
        (r"gogat", "GOGAT"),
        (r"glutamate dehydrogenase", "glutamate dehydrogenase"),
        (r"glutamine", "glutamine"),
        (r"glutamate", "glutamate"),
        (r"aminotransferase", "aminotransferase"),
        (r"omega-amidase", "omega-amidase"),
        (r"nit2", "NIT2"),
    ],
    "Asp/Asn": [
        (r"asparagine", "asparagine"),
        (r"aspartate", "aspartate"),
        (r"aspartyl", "aspartyl"),
        (r"aspartate kinase", "aspartate kinase"),
    ],
    "Proline Metabolism": [
        (r"p5cs", "P5CS"),
        (r"pyrroline-5-carboxylate", "pyrroline-5-carboxylate"),
        (r"delta 1-pyrroline-5-carboxylate synthase", "delta 1-pyrroline-5-carboxylate synthase"),
        (r"delta-1-pyrroline-5-carboxylate synthase", "delta-1-pyrroline-5-carboxylate synthase"),
        (r"proline dehydrogenase", "proline dehydrogenase"),
        (r"proline oxidase", "proline oxidase"),
        (r"ornithine aminotransferase", "ornithine aminotransferase"),
        (r"proline transporter", "proline transporter"),
        (r"p5c dehydrogenase", "P5C dehydrogenase"),
    ],
    "Amino Acid Transport": [
        (r"amino acid transporter", "amino acid transporter"),
        (r"amino acid permease", "amino acid permease"),
        (r"cationic amino acid transporter", "cationic amino acid transporter"),
        (r"bidirectional amino acid transporter", "bidirectional amino acid transporter"),
        (r"peptide transporter", "peptide transporter"),
        (r"oligopeptide transporter", "oligopeptide transporter"),
        (r"umamit", "UMAMIT"),
    ],
    "Mineral N Transport": [
        (r"nitrate", "nitrate"),
        (r"nitrite", "nitrite"),
        (r"ammonium transporter", "ammonium transporter"),
        (r"urea transporter", "urea transporter"),
    ],
    "Carbon + Photosynthesis": [
        (r"photosystem", "photosystem"),
        (r"photosynthesis", "photosynthesis"),
        (r"chlorophyll", "chlorophyll"),
        (r"carbonic anhydrase", "carbonic anhydrase"),
        (r"phosphoenolpyruvate", "phosphoenolpyruvate"),
        (r"pep carboxylase", "PEP carboxylase"),
        (r"sucrose", "sucrose"),
        (r"starch", "starch"),
        (r"calvin", "Calvin cycle"),
        (r"rubisco", "Rubisco"),
        (r"malate dehydrogenase", "malate dehydrogenase"),
    ],
    "Stress Response": [
        (r"stress", "stress"),
        (r"heat shock", "heat shock"),
        (r"jasmonate", "jasmonate"),
        (r"ethylene", "ethylene"),
        (r"abscisic", "abscisic"),
        (r"glutathione", "glutathione"),
    ],
    "Senescence + Remobilization": [
        (r"autophagy", "autophagy"),
        (r"ubiquitin", "ubiquitin"),
        (r"proteasome", "proteasome"),
        (r"protease", "protease"),
        (r"peptidase", "peptidase"),
        (r"senescence", "senescence"),
        (r"e3 ubiquitin", "E3 ubiquitin"),
        (r"deubiquitin", "deubiquitin"),
    ],
}


def first_nonempty(values: pd.Series) -> str:
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text.lower() != "nan":
            return text
    return ""


def clean_gene_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def split_semicolon(value: object) -> list[str]:
    if pd.isna(value):
        return []
    parts = [part.strip() for part in str(value).split(";")]
    return [part for part in parts if part]


def collapse_unique(values: list[str]) -> str:
    unique_values = sorted({value for value in values if value})
    return ";".join(unique_values)


def yes_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"yes", "true", "1"})


def format_float(value: float | int | None, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


def load_deg_workbook() -> pd.DataFrame:
    workbook = pd.ExcelFile(RNA_WORKBOOK)
    frames: list[pd.DataFrame] = []

    for sheet_name in workbook.sheet_names:
        df = workbook.parse(sheet_name)
        year, stage = sheet_name.split()

        description_col = "description" if "description" in df.columns else None
        if description_col is None and "description.x" in df.columns:
            description_col = "description.x"

        control_fpkm_cols = [col for col in df.columns if "control" in col.lower() and col.lower().endswith("_fpkm")]
        treatment_fpkm_cols = [col for col in df.columns if "urea" in col.lower() and col.lower().endswith("_fpkm")]

        frame = pd.DataFrame(
            {
                "GeneID": df["Geneid"].map(clean_gene_id),
                "contrast": sheet_name,
                "year": year,
                "stage": stage,
                "Regulated": df.get("Regulated", pd.Series(index=df.index, dtype=object)).astype(str).str.strip(),
                "padj": pd.to_numeric(df.get("padj"), errors="coerce"),
                "pvalue": pd.to_numeric(df.get("pvalue"), errors="coerce"),
                "log2FoldChange": pd.to_numeric(df.get("log2FoldChange"), errors="coerce"),
                "baseMean": pd.to_numeric(df.get("baseMean"), errors="coerce"),
                "description": df.get(description_col, pd.Series(index=df.index, dtype=object)),
                "function": df.get("function.", pd.Series(index=df.index, dtype=object)),
                "Family": df.get("Family", pd.Series(index=df.index, dtype=object)),
                "AGI": df.get("AGI", pd.Series(index=df.index, dtype=object)),
                "Control_Expressed": df.get("Control_Expressed", pd.Series(index=df.index, dtype=object)),
                "Urea_Expressed": df.get("Urea_Expressed", pd.Series(index=df.index, dtype=object)),
            }
        )

        if control_fpkm_cols:
            frame["control_mean_fpkm"] = df[control_fpkm_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        else:
            frame["control_mean_fpkm"] = np.nan

        if treatment_fpkm_cols:
            frame["treatment_mean_fpkm"] = df[treatment_fpkm_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        else:
            frame["treatment_mean_fpkm"] = np.nan

        frames.append(frame)

    deg_df = pd.concat(frames, ignore_index=True)
    deg_df = deg_df[deg_df["GeneID"].str.startswith("Zm")].copy()
    deg_df = deg_df.dropna(subset=["padj", "log2FoldChange"])
    deg_df["text_blob"] = (
        deg_df[["description", "function", "Family", "AGI"]]
        .fillna("")
        .agg(" ".join, axis=1)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
        .str.lower()
    )
    deg_df["abs_log2FoldChange"] = deg_df["log2FoldChange"].abs()
    return deg_df


def load_amino_metadata() -> tuple[dict[str, list[str]], dict[str, list[str]], set[str], set[str], set[str]]:
    amino_df = pd.read_csv(AMINO_GWAS_FILE)
    amino_df["GeneID"] = amino_df["GeneID"].map(clean_gene_id)
    amino_df["Phenotype"] = amino_df["Phenotype"].astype(str).str.strip()

    trait_map = (
        amino_df.groupby("GeneID")["Phenotype"]
        .apply(lambda values: sorted(set(values)))
        .to_dict()
    )

    focal_trait_map = {
        gene_id: [trait for trait in traits if trait in FOCAL_AMINO_TRAITS]
        for gene_id, traits in trait_map.items()
    }

    all_amino_genes = set(trait_map)
    focal_amino_genes = {gene_id for gene_id, traits in focal_trait_map.items() if traits}
    asparagine_genes = set(amino_df.loc[amino_df["Phenotype"] == "N", "GeneID"])

    return trait_map, focal_trait_map, all_amino_genes, focal_amino_genes, asparagine_genes


def standardize_annotation_frame(
    df: pd.DataFrame,
    *,
    gene_col: str = "GeneID",
    name_col: str = "GeneName",
    family_col: str = "Family_Subfamily",
    protein_col: str = "Protein_Class",
    mf_col: str = "GO_MF",
    bp_col: str = "GO_BP",
    cc_col: str = "GO_CC",
) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "GeneID": df.get(gene_col, pd.Series(index=df.index, dtype=object)).map(clean_gene_id),
            "GeneName": df.get(name_col, pd.Series(index=df.index, dtype=object)).fillna(""),
            "Family_Subfamily": df.get(family_col, pd.Series(index=df.index, dtype=object)).fillna(""),
            "Protein_Class": df.get(protein_col, pd.Series(index=df.index, dtype=object)).fillna(""),
            "GO_MF": df.get(mf_col, pd.Series(index=df.index, dtype=object)).fillna(""),
            "GO_BP": df.get(bp_col, pd.Series(index=df.index, dtype=object)).fillna(""),
            "GO_CC": df.get(cc_col, pd.Series(index=df.index, dtype=object)).fillna(""),
        }
    )
    out = out.loc[out["GeneID"].str.startswith("Zm")].copy()
    out["text_blob"] = (
        out[["GeneName", "Family_Subfamily", "Protein_Class", "GO_MF", "GO_BP", "GO_CC"]]
        .fillna("")
        .agg(" ".join, axis=1)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
        .str.lower()
    )
    return out.drop_duplicates("GeneID").reset_index(drop=True)


def load_fst_annotation_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        names=GENERIC_ANNOTATION_COLUMNS,
        header=None,
        skip_blank_lines=True,
    )
    return standardize_annotation_frame(df)


def parse_soiln_annotations() -> tuple[set[str], dict[str, list[str]], dict[str, str]]:
    soil_df = pd.read_csv(SOILN_ANNOTATION_FILE)
    soil_df["GeneID"] = soil_df["GeneID"].map(clean_gene_id)

    soil_genes: set[str] = set()
    phenotype_map: dict[str, list[str]] = {}
    soiln_log10p: dict[str, str] = {}

    for _, row in soil_df.iterrows():
        phenotypes = split_semicolon(row.get("Phenotypes"))
        if not phenotypes:
            continue

        phenotype_map[row["GeneID"]] = phenotypes
        if "SoilN" not in phenotypes:
            continue

        soil_genes.add(row["GeneID"])
        pvals = split_semicolon(row.get("pvalues"))
        pairs = list(zip(phenotypes, pvals))
        soiln_hits = [value for phenotype, value in pairs if phenotype == "SoilN"]
        soiln_log10p[row["GeneID"]] = ";".join(soiln_hits)

    return soil_genes, phenotype_map, soiln_log10p


def explode_selection_table(path: Path, population_code: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.loc[yes_mask(df["fst_outlier"])].copy()
    if df.empty:
        return pd.DataFrame(columns=["GeneID", "Population", "window_id", "mean_fst", "window_start_bp", "window_end_bp"])

    records: list[dict[str, object]] = []
    for _, row in df.iterrows():
        gene_ids = split_semicolon(row.get("GeneIDs"))
        for gene_id in gene_ids:
            records.append(
                {
                    "GeneID": gene_id,
                    "Population": population_code,
                    "window_id": row.get("window_id"),
                    "mean_fst": pd.to_numeric(row.get("mean_fst"), errors="coerce"),
                    "window_start_bp": pd.to_numeric(row.get("window_start_bp"), errors="coerce"),
                    "window_end_bp": pd.to_numeric(row.get("window_end_bp"), errors="coerce"),
                }
            )

    return pd.DataFrame.from_records(records)


def build_selection_summary(ic_hits: pd.DataFrame, j_hits: pd.DataFrame) -> tuple[set[str], dict[str, dict[str, object]]]:
    selection_df = pd.concat([ic_hits, j_hits], ignore_index=True)
    if selection_df.empty:
        return set(), {}

    summary: dict[str, dict[str, object]] = {}
    for gene_id, gene_df in selection_df.groupby("GeneID"):
        pop_counts = gene_df.groupby("Population")["window_id"].nunique().to_dict()
        pop_labels = sorted(pop_counts)
        summary[gene_id] = {
            "selection_populations": ";".join(pop_labels),
            "selection_window_count": int(gene_df["window_id"].nunique()),
            "selection_window_count_IC": int(pop_counts.get("IC", 0)),
            "selection_window_count_J": int(pop_counts.get("J", 0)),
            "selection_max_mean_fst": gene_df["mean_fst"].max(),
        }

    return set(summary), summary


def annotate_pathway_hits(gene_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    category_map: dict[str, list[str]] = defaultdict(list)
    keyword_map: dict[str, list[str]] = defaultdict(list)
    records: list[dict[str, str]] = []

    for row in gene_df.itertuples(index=False):
        text_blob = row.text_blob
        for category, pattern_items in PATHWAY_KEYWORDS.items():
            matched_labels = [label for pattern, label in pattern_items if re.search(pattern, text_blob, flags=re.IGNORECASE)]
            if matched_labels:
                category_map[row.GeneID].append(category)
                keyword_map[row.GeneID].extend(matched_labels)
                records.append(
                    {
                        "GeneID": row.GeneID,
                        "Pathway_category": category,
                        "Matched_keywords_by_category": collapse_unique(matched_labels),
                    }
                )

    category_df = pd.DataFrame({"GeneID": gene_df["GeneID"].unique()})
    category_df["Panel_categories"] = category_df["GeneID"].map(
        lambda gene_id: collapse_unique(category_map.get(gene_id, []))
    )
    category_df["Matched_keywords"] = category_df["GeneID"].map(
        lambda gene_id: collapse_unique(keyword_map.get(gene_id, []))
    )
    category_df["Annotation_panel_hit"] = category_df["Panel_categories"].ne("")
    long_df = pd.DataFrame.from_records(records).drop_duplicates()
    if not long_df.empty:
        long_df["Pathway_category"] = pd.Categorical(long_df["Pathway_category"], categories=PATHWAY_ORDER, ordered=True)
        long_df = long_df.sort_values(["Pathway_category", "GeneID"]).reset_index(drop=True)
    return category_df, long_df


def classify_annotation_source(annotation_df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    records: list[dict[str, str]] = []

    for row in annotation_df.itertuples(index=False):
        for category, pattern_items in PATHWAY_KEYWORDS.items():
            matched_labels = [label for pattern, label in pattern_items if re.search(pattern, row.text_blob, flags=re.IGNORECASE)]
            if not matched_labels:
                continue
            records.append(
                {
                    "Source": source_name,
                    "GeneID": row.GeneID,
                    "GeneName": row.GeneName,
                    "Family_Subfamily": row.Family_Subfamily,
                    "Protein_Class": row.Protein_Class,
                    "GO_MF": row.GO_MF,
                    "GO_BP": row.GO_BP,
                    "GO_CC": row.GO_CC,
                    "Pathway_category": category,
                    "Matched_keywords_by_category": collapse_unique(matched_labels),
                }
            )

    long_df = pd.DataFrame.from_records(records)
    if long_df.empty:
        return long_df

    long_df["Pathway_category"] = pd.Categorical(long_df["Pathway_category"], categories=CORE_FIGURE_PATHWAY_ORDER, ordered=True)
    long_df = long_df.sort_values(["Source", "Pathway_category", "GeneID"]).drop_duplicates(
        subset=["Source", "GeneID", "Pathway_category"]
    )
    return long_df.reset_index(drop=True)


def build_gene_level_deg(deg_df: pd.DataFrame) -> pd.DataFrame:
    gene_records: list[dict[str, object]] = []

    for gene_id, group in deg_df.groupby("GeneID"):
        ordered = group.copy()
        ordered["contrast"] = pd.Categorical(ordered["contrast"], categories=CONTRAST_ORDER, ordered=True)
        ordered = ordered.sort_values(["contrast", "padj", "abs_log2FoldChange"], ascending=[True, True, False])

        direction_values = ordered["Regulated"].dropna().astype(str).str.strip()
        direction_values = [value for value in direction_values if value and value.lower() != "nan"]
        direction_summary = collapse_unique(direction_values)

        details = [
            f"{row.contrast}:{row.Regulated} (log2FC={row.log2FoldChange:.2f}, padj={row.padj:.2e})"
            for row in ordered.itertuples(index=False)
        ]

        gene_records.append(
            {
                "GeneID": gene_id,
                "RNA_contrasts": ";".join(str(value) for value in ordered["contrast"].dropna().unique()),
                "RNA_n_contrasts": int(ordered["contrast"].nunique()),
                "RNA_years": ";".join(sorted(ordered["year"].dropna().astype(str).unique())),
                "RNA_stages": ";".join(sorted(ordered["stage"].dropna().astype(str).unique())),
                "RNA_direction": direction_summary,
                "RNA_min_padj": ordered["padj"].min(),
                "RNA_max_abs_log2FC": ordered["abs_log2FoldChange"].max(),
                "RNA_mean_log2FC": ordered["log2FoldChange"].mean(),
                "RNA_detail": "; ".join(details),
                "control_mean_fpkm": ordered["control_mean_fpkm"].mean(),
                "treatment_mean_fpkm": ordered["treatment_mean_fpkm"].mean(),
                "description": first_nonempty(ordered["description"]),
                "function": first_nonempty(ordered["function"]),
                "Family": first_nonempty(ordered["Family"]),
                "AGI": first_nonempty(ordered["AGI"]),
                "text_blob": first_nonempty(ordered["text_blob"]),
            }
        )

    gene_df = pd.DataFrame.from_records(gene_records)
    return gene_df.sort_values(["RNA_n_contrasts", "RNA_min_padj"], ascending=[False, True]).reset_index(drop=True)


def add_support_layers(
    gene_df: pd.DataFrame,
    amino_traits: dict[str, list[str]],
    focal_traits: dict[str, list[str]],
    all_amino_genes: set[str],
    focal_amino_genes: set[str],
    asparagine_genes: set[str],
    soil_genes: set[str],
    soiln_log10p: dict[str, str],
    selection_genes: set[str],
    selection_summary: dict[str, dict[str, object]],
    panel_df: pd.DataFrame,
) -> pd.DataFrame:
    panel_map = panel_df.set_index("GeneID")
    support_df = gene_df.copy()

    support_df["Amino_any_GWAS"] = support_df["GeneID"].isin(all_amino_genes)
    support_df["Amino_focal_N_traits"] = support_df["GeneID"].isin(focal_amino_genes)
    support_df["Asparagine_GWAS_N_trait"] = support_df["GeneID"].isin(asparagine_genes)
    support_df["SoilN_GWAS"] = support_df["GeneID"].isin(soil_genes)
    support_df["Selection_IC_or_J"] = support_df["GeneID"].isin(selection_genes)

    support_df["Amino_traits"] = support_df["GeneID"].map(lambda gene_id: collapse_unique(amino_traits.get(gene_id, [])))
    support_df["Amino_focal_traits_list"] = support_df["GeneID"].map(
        lambda gene_id: collapse_unique(focal_traits.get(gene_id, []))
    )
    support_df["SoilN_log10P"] = support_df["GeneID"].map(soiln_log10p).fillna("")

    support_df["Selection_populations"] = support_df["GeneID"].map(
        lambda gene_id: selection_summary.get(gene_id, {}).get("selection_populations", "")
    )
    support_df["Selection_window_count"] = support_df["GeneID"].map(
        lambda gene_id: selection_summary.get(gene_id, {}).get("selection_window_count", 0)
    )
    support_df["Selection_window_count_IC"] = support_df["GeneID"].map(
        lambda gene_id: selection_summary.get(gene_id, {}).get("selection_window_count_IC", 0)
    )
    support_df["Selection_window_count_J"] = support_df["GeneID"].map(
        lambda gene_id: selection_summary.get(gene_id, {}).get("selection_window_count_J", 0)
    )
    support_df["Selection_max_mean_fst"] = support_df["GeneID"].map(
        lambda gene_id: selection_summary.get(gene_id, {}).get("selection_max_mean_fst", np.nan)
    )

    support_df["Annotation_panel_hit"] = support_df["GeneID"].map(panel_map["Annotation_panel_hit"]).fillna(False)
    support_df["Panel_categories"] = support_df["GeneID"].map(panel_map["Panel_categories"]).fillna("")
    support_df["Matched_keywords"] = support_df["GeneID"].map(panel_map["Matched_keywords"]).fillna("")

    support_df["Core_panel_categories"] = support_df["Panel_categories"].map(
        lambda value: collapse_unique(
            [category for category in split_semicolon(value) if category in CORE_N_REMOBILIZATION_CATEGORIES]
        )
    )
    support_df["N_remobilization_annotation_hit"] = support_df["Core_panel_categories"].ne("")

    support_df["Genetic_support_count_broad"] = support_df[
        ["Amino_any_GWAS", "Asparagine_GWAS_N_trait", "SoilN_GWAS", "Selection_IC_or_J"]
    ].sum(axis=1)
    support_df["Genetic_support_count"] = support_df[
        ["Amino_focal_N_traits", "Asparagine_GWAS_N_trait", "SoilN_GWAS", "Selection_IC_or_J"]
    ].sum(axis=1)
    support_df["Support_score_broad"] = (
        1 + support_df["Genetic_support_count_broad"] + support_df["N_remobilization_annotation_hit"].astype(int)
    )
    support_df["Support_score"] = 1 + support_df["Genetic_support_count"] + support_df["N_remobilization_annotation_hit"].astype(int)

    support_df["Support_tier"] = np.select(
        [
            (support_df["Genetic_support_count"] >= 1) & support_df["N_remobilization_annotation_hit"],
            support_df["Genetic_support_count"] >= 1,
            support_df["N_remobilization_annotation_hit"],
        ],
        [
            "Tier 1: RNA-seq + genetic + annotation support",
            "Tier 2: RNA-seq + genetic support",
            "Tier 3: RNA-seq + annotation support",
        ],
        default="Tier 4: RNA-seq only",
    )

    return support_df


def build_overlap_summary(deg_df: pd.DataFrame, support_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for contrast in CONTRAST_ORDER + ["Any DEG"]:
        if contrast == "Any DEG":
            subset = support_df.copy()
        else:
            gene_ids = set(deg_df.loc[deg_df["contrast"] == contrast, "GeneID"])
            subset = support_df.loc[support_df["GeneID"].isin(gene_ids)].copy()

        rows.append(
            {
                "contrast": contrast,
                "n_deg_genes": int(subset["GeneID"].nunique()),
                "overlap_amino_any": int(subset.loc[subset["Amino_any_GWAS"], "GeneID"].nunique()),
                "overlap_amino_focal": int(subset.loc[subset["Amino_focal_N_traits"], "GeneID"].nunique()),
                "overlap_asparagine_N": int(subset.loc[subset["Asparagine_GWAS_N_trait"], "GeneID"].nunique()),
                "overlap_soilN": int(subset.loc[subset["SoilN_GWAS"], "GeneID"].nunique()),
                "overlap_selection": int(subset.loc[subset["Selection_IC_or_J"], "GeneID"].nunique()),
                "overlap_annotation_panel": int(subset.loc[subset["Annotation_panel_hit"], "GeneID"].nunique()),
                "overlap_n_remobilization_panel": int(
                    subset.loc[subset["N_remobilization_annotation_hit"], "GeneID"].nunique()
                ),
                "overlap_any_external_layer": int(
                    subset.loc[
                        subset[
                            [
                                "Amino_any_GWAS",
                                "Asparagine_GWAS_N_trait",
                                "SoilN_GWAS",
                                "Selection_IC_or_J",
                                "N_remobilization_annotation_hit",
                            ]
                        ].any(axis=1),
                        "GeneID",
                    ].nunique()
                ),
            }
        )

    return pd.DataFrame(rows)


def build_directionality_summary(deg_df: pd.DataFrame, panel_df: pd.DataFrame) -> pd.DataFrame:
    category_lookup = (
        panel_df.loc[panel_df["Annotation_panel_hit"], ["GeneID", "Panel_categories"]]
        .assign(Panel_categories=lambda df: df["Panel_categories"].str.split(";"))
        .explode("Panel_categories")
        .rename(columns={"Panel_categories": "Pathway_category"})
    )

    if category_lookup.empty:
        return pd.DataFrame(
            columns=["Pathway_category", "contrast", "n_unique_genes", "n_up", "n_down", "mean_log2FC", "median_log2FC", "direction_balance"]
        )

    category_deg = deg_df.merge(category_lookup, on="GeneID", how="inner")

    rows: list[dict[str, object]] = []
    for category in PATHWAY_ORDER:
        for contrast in CONTRAST_ORDER:
            year, stage = contrast.split()
            subset = category_deg.loc[
                (category_deg["Pathway_category"] == category) &
                (category_deg["contrast"] == contrast)
            ].copy()

            if subset.empty:
                rows.append(
                    {
                        "Pathway_category": category,
                        "contrast": contrast,
                        "year": year,
                        "stage": stage,
                        "condition_numerator": "urea",
                        "condition_denominator": "control",
                        "n_unique_genes": 0,
                        "n_up": 0,
                        "n_down": 0,
                        "mean_log2FC": np.nan,
                        "median_log2FC": np.nan,
                        "direction_balance": np.nan,
                        "mean_control_fpkm": np.nan,
                        "mean_urea_fpkm": np.nan,
                        "median_control_fpkm": np.nan,
                        "median_urea_fpkm": np.nan,
                    }
                )
                continue

            n_unique = subset["GeneID"].nunique()
            n_up = int((subset["log2FoldChange"] > 0).sum())
            n_down = int((subset["log2FoldChange"] < 0).sum())

            rows.append(
                {
                    "Pathway_category": category,
                    "contrast": contrast,
                    "year": year,
                    "stage": stage,
                    "condition_numerator": "urea",
                    "condition_denominator": "control",
                    "n_unique_genes": int(n_unique),
                    "n_up": n_up,
                    "n_down": n_down,
                    "mean_log2FC": subset["log2FoldChange"].mean(),
                    "median_log2FC": subset["log2FoldChange"].median(),
                    "direction_balance": (n_up - n_down) / n_unique if n_unique else np.nan,
                    "mean_control_fpkm": subset["control_mean_fpkm"].mean(),
                    "mean_urea_fpkm": subset["treatment_mean_fpkm"].mean(),
                    "median_control_fpkm": subset["control_mean_fpkm"].median(),
                    "median_urea_fpkm": subset["treatment_mean_fpkm"].median(),
                }
            )

    return pd.DataFrame(rows)


def build_pathway_candidate_map(support_df: pd.DataFrame, pathway_long_df: pd.DataFrame) -> pd.DataFrame:
    if pathway_long_df.empty:
        return pd.DataFrame()

    map_df = pathway_long_df.merge(
        support_df,
        on="GeneID",
        how="left",
        validate="many_to_one",
    )
    map_df["Pathway_category"] = pd.Categorical(map_df["Pathway_category"], categories=PATHWAY_ORDER, ordered=True)
    map_df = map_df.sort_values(
        [
            "Pathway_category",
            "Support_score",
            "RNA_n_contrasts",
            "Selection_IC_or_J",
            "Amino_focal_N_traits",
            "RNA_min_padj",
        ],
        ascending=[True, False, False, False, False, True],
    ).reset_index(drop=True)
    return map_df


def build_pathway_candidate_summary(pathway_map_df: pd.DataFrame) -> pd.DataFrame:
    if pathway_map_df.empty:
        return pd.DataFrame()

    summary_rows: list[dict[str, object]] = []
    for category in PATHWAY_ORDER:
        subset = pathway_map_df.loc[pathway_map_df["Pathway_category"] == category].copy()
        if subset.empty:
            summary_rows.append(
                {
                    "Pathway_category": category,
                    "n_candidate_genes": 0,
                    "n_recurrent_deg_genes": 0,
                    "n_amino_focal_genes": 0,
                    "n_asparagine_gwas_genes": 0,
                    "n_soilN_genes": 0,
                    "n_selection_genes": 0,
                    "Top_candidate_genes": "",
                }
            )
            continue

        unique_subset = subset.drop_duplicates("GeneID")
        top_genes = unique_subset.sort_values(
            ["Support_score", "RNA_n_contrasts", "RNA_min_padj"],
            ascending=[False, False, True],
        )["GeneID"].head(12).tolist()

        summary_rows.append(
            {
                "Pathway_category": category,
                "n_candidate_genes": int(unique_subset["GeneID"].nunique()),
                "n_recurrent_deg_genes": int(unique_subset.loc[unique_subset["RNA_n_contrasts"] >= 2, "GeneID"].nunique()),
                "n_amino_focal_genes": int(unique_subset.loc[unique_subset["Amino_focal_N_traits"], "GeneID"].nunique()),
                "n_asparagine_gwas_genes": int(unique_subset.loc[unique_subset["Asparagine_GWAS_N_trait"], "GeneID"].nunique()),
                "n_soilN_genes": int(unique_subset.loc[unique_subset["SoilN_GWAS"], "GeneID"].nunique()),
                "n_selection_genes": int(unique_subset.loc[unique_subset["Selection_IC_or_J"], "GeneID"].nunique()),
                "Top_candidate_genes": ";".join(top_genes),
            }
        )

    return pd.DataFrame(summary_rows)


def build_all_source_candidate_pathway_maps(
    support_df: pd.DataFrame,
    pathway_candidate_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_frames: list[pd.DataFrame] = []

    # RNA-seq directional candidates already carry the pathway assignments.
    rna_source = pathway_candidate_map.loc[
        pathway_candidate_map["Pathway_category"].isin(CORE_FIGURE_PATHWAY_ORDER),
        [
            "GeneID",
            "Pathway_category",
            "Matched_keywords_by_category",
            "description",
            "Family",
        ],
    ].copy()
    rna_source = rna_source.rename(columns={"description": "GeneName", "Family": "Family_Subfamily"})
    rna_source["Protein_Class"] = ""
    rna_source["GO_MF"] = ""
    rna_source["GO_BP"] = ""
    rna_source["GO_CC"] = ""
    rna_source["Source"] = "RNA-seq DEGs"
    source_frames.append(rna_source)

    # All amino-acid GWAS candidates.
    all_amino_annot = standardize_annotation_frame(pd.read_csv(ALL_AMINO_ANNOTATION_FILE, sep="\t"))
    source_frames.append(classify_annotation_source(all_amino_annot, "All amino GWAS"))

    # Soil N GWAS candidates.
    soil_df = pd.read_csv(SOILN_ANNOTATION_FILE)
    soil_df = soil_df.rename(columns={"GO_BO": "GO_BP"})
    soil_df = soil_df.loc[
        soil_df["Phenotypes"].fillna("").str.contains(r"(?:^|;)SoilN(?:;|$)", regex=True)
    ].copy()
    soil_annot = standardize_annotation_frame(
        soil_df,
        gene_col="GeneID",
        name_col="GeneSymbol",
        family_col="Family_Subfamily",
        protein_col="Protein_Class",
        mf_col="GO_MF",
        bp_col="GO_BP",
        cc_col="GO_CC",
    )
    source_frames.append(classify_annotation_source(soil_annot, "SoilN GWAS"))

    # IC/J selection candidate genes from the broader FST annotation sets.
    ic_annot = load_fst_annotation_file(IC_FST_ANNOTATION_FILE)
    j_annot = load_fst_annotation_file(J_FST_ANNOTATION_FILE)
    source_frames.append(classify_annotation_source(ic_annot, "IC FST 5pc"))
    source_frames.append(classify_annotation_source(j_annot, "J FST 5pc"))

    candidate_map = pd.concat(source_frames, ignore_index=True)
    candidate_map["Pathway_category"] = pd.Categorical(
        candidate_map["Pathway_category"], categories=CORE_FIGURE_PATHWAY_ORDER, ordered=True
    )
    candidate_map = candidate_map.sort_values(["Pathway_category", "Source", "GeneID"]).reset_index(drop=True)

    count_rows: list[dict[str, object]] = []
    for pathway in CORE_FIGURE_PATHWAY_ORDER:
        for source in CANDIDATE_SOURCE_ORDER:
            subset = candidate_map.loc[
                (candidate_map["Pathway_category"] == pathway) &
                (candidate_map["Source"] == source)
            ]
            count_rows.append(
                {
                    "Pathway_category": pathway,
                    "Source": source,
                    "n_genes": int(subset["GeneID"].nunique()),
                }
            )

    count_summary = pd.DataFrame(count_rows)
    return candidate_map, count_summary


def plot_candidate_source_counts(summary_df: pd.DataFrame, output_path: Path) -> None:
    heatmap = (
        summary_df.pivot(index="Pathway_category", columns="Source", values="n_genes")
        .reindex(index=CORE_FIGURE_PATHWAY_ORDER, columns=CANDIDATE_SOURCE_ORDER)
    )
    display_values = np.log10(heatmap.fillna(0).to_numpy(dtype=float) + 1.0)

    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    im = ax.imshow(display_values, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(CANDIDATE_SOURCE_ORDER)))
    ax.set_xticklabels(CANDIDATE_SOURCE_ORDER, rotation=20, ha="right")
    ax.set_yticks(range(len(CORE_FIGURE_PATHWAY_ORDER)))
    ax.set_yticklabels(CORE_FIGURE_PATHWAY_ORDER)
    ax.set_title("Pathway-mapped candidate genes across evidence layers")
    ax.set_xlabel("Evidence layer")

    for i, pathway in enumerate(CORE_FIGURE_PATHWAY_ORDER):
        for j, source in enumerate(CANDIDATE_SOURCE_ORDER):
            n_value = int(heatmap.iloc[i, j]) if not pd.isna(heatmap.iloc[i, j]) else 0
            color = "white" if display_values[i, j] >= np.nanmax(display_values) * 0.55 else "black"
            ax.text(j, i, f"n={n_value}", ha="center", va="center", fontsize=9, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("log10(candidate gene count + 1)")

    ax.set_xticks(np.arange(-0.5, len(CANDIDATE_SOURCE_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CORE_FIGURE_PATHWAY_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_directionality(
    summary_df: pd.DataFrame,
    output_path: Path,
    pathway_order: list[str],
    title: str,
    figsize: tuple[float, float],
) -> None:
    plot_df = summary_df.copy()
    heatmap = (
        plot_df.pivot(index="Pathway_category", columns="contrast", values="mean_log2FC")
        .reindex(index=pathway_order, columns=CONTRAST_ORDER)
    )

    annotation = (
        plot_df.assign(label=lambda df: np.where(
            df["n_unique_genes"] > 0,
            "n=" + df["n_unique_genes"].astype(str) + "\n" + df["n_up"].astype(str) + " up / " + df["n_down"].astype(str) + " down",
            "",
        ))
        .pivot(index="Pathway_category", columns="contrast", values="label")
        .reindex(index=pathway_order, columns=CONTRAST_ORDER)
    )

    finite_values = heatmap.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(finite_values)) if np.isfinite(finite_values).any() else 1.0
    vmax = max(vmax, 1.0)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(heatmap.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(CONTRAST_ORDER)))
    ax.set_xticklabels(["2022\nV6", "2023\nV6", "2022\nV10", "2023\nV10"], rotation=0)
    ax.set_yticks(range(len(pathway_order)))
    ax.set_yticklabels(pathway_order)
    ax.set_title(title)
    ax.set_xlabel("Year x vegetative stage")

    for i, pathway in enumerate(pathway_order):
        for j, contrast in enumerate(CONTRAST_ORDER):
            label = annotation.iloc[i, j]
            value = heatmap.iloc[i, j]
            if not label:
                continue
            text_color = "white" if pd.notna(value) and abs(value) >= (0.45 * vmax) else "black"
            ax.text(j, i, label, ha="center", va="center", fontsize=8.5, color=text_color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.92)
    cbar.set_label("Mean log2 fold change (urea / control)")

    ax.set_xticks(np.arange(-0.5, len(CONTRAST_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(pathway_order), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    deg_df = load_deg_workbook()
    gene_df = build_gene_level_deg(deg_df)

    amino_traits, focal_traits, all_amino_genes, focal_amino_genes, asparagine_genes = load_amino_metadata()
    soil_genes, _, soiln_log10p = parse_soiln_annotations()

    ic_selection = explode_selection_table(IC_FST_FILE, population_code="IC")
    j_selection = explode_selection_table(J_FST_FILE, population_code="J")
    selection_genes, selection_summary = build_selection_summary(ic_selection, j_selection)

    panel_df, pathway_long_df = annotate_pathway_hits(gene_df)
    support_df = add_support_layers(
        gene_df=gene_df,
        amino_traits=amino_traits,
        focal_traits=focal_traits,
        all_amino_genes=all_amino_genes,
        focal_amino_genes=focal_amino_genes,
        asparagine_genes=asparagine_genes,
        soil_genes=soil_genes,
        soiln_log10p=soiln_log10p,
        selection_genes=selection_genes,
        selection_summary=selection_summary,
        panel_df=panel_df,
    )

    overlap_gene_table = support_df.loc[
        support_df[
            [
                "Amino_any_GWAS",
                "Asparagine_GWAS_N_trait",
                "SoilN_GWAS",
                "Selection_IC_or_J",
            ]
        ].any(axis=1)
    ].copy()
    overlap_gene_table = overlap_gene_table.sort_values(
        ["Support_score", "RNA_n_contrasts", "RNA_min_padj", "Selection_window_count"],
        ascending=[False, False, True, False],
    )

    remobilization_panel = support_df.loc[support_df["N_remobilization_annotation_hit"]].copy()
    remobilization_panel = remobilization_panel.sort_values(
        ["Support_score", "N_remobilization_annotation_hit", "RNA_n_contrasts", "RNA_min_padj"],
        ascending=[False, False, False, True],
    )

    ranking_table = support_df.loc[
        support_df[
            [
                "Amino_focal_N_traits",
                "Asparagine_GWAS_N_trait",
                "SoilN_GWAS",
                "Selection_IC_or_J",
                "N_remobilization_annotation_hit",
            ]
        ].any(axis=1)
    ].copy()
    ranking_table = ranking_table.sort_values(
        [
            "Support_score",
            "N_remobilization_annotation_hit",
            "Asparagine_GWAS_N_trait",
            "SoilN_GWAS",
            "Amino_focal_N_traits",
            "Selection_IC_or_J",
            "RNA_n_contrasts",
            "RNA_min_padj",
            "RNA_max_abs_log2FC",
        ],
        ascending=[False, False, False, False, False, False, False, True, False],
    )

    overlap_summary = build_overlap_summary(deg_df, support_df)
    directionality_summary = build_directionality_summary(deg_df, panel_df)
    core_directionality_summary = directionality_summary.loc[
        directionality_summary["Pathway_category"].isin(CORE_FIGURE_PATHWAY_ORDER)
    ].copy()
    pathway_candidate_map = build_pathway_candidate_map(support_df, pathway_long_df)
    pathway_candidate_summary = build_pathway_candidate_summary(pathway_candidate_map)
    all_source_candidate_map, all_source_candidate_counts = build_all_source_candidate_pathway_maps(
        support_df=support_df,
        pathway_candidate_map=pathway_candidate_map,
    )

    overlap_summary.to_csv(TABLE_DIR / "SuppTable_RNAseq_DEG_overlap_summary.csv", index=False)
    overlap_gene_table.to_csv(TABLE_DIR / "SuppTable_RNAseq_DEG_overlap_genes.csv", index=False)
    remobilization_panel.to_csv(TABLE_DIR / "SuppTable_RNAseq_N_remobilization_panel.csv", index=False)
    ranking_table.to_csv(TABLE_DIR / "SuppTable_RNAseq_multilayer_support_ranking.csv", index=False)
    directionality_summary.to_csv(TABLE_DIR / "SuppTable_RNAseq_directional_pathway_summary.csv", index=False)
    core_directionality_summary.to_csv(TABLE_DIR / "SuppTable_RNAseq_directional_core_pathway_summary.csv", index=False)
    pathway_candidate_map.to_csv(TABLE_DIR / "SuppTable_RNAseq_pathway_candidate_map.csv", index=False)
    pathway_candidate_summary.to_csv(TABLE_DIR / "SuppTable_RNAseq_pathway_candidate_summary.csv", index=False)
    all_source_candidate_map.to_csv(TABLE_DIR / "SuppTable_candidate_pathway_source_genes.csv", index=False)
    all_source_candidate_counts.to_csv(TABLE_DIR / "SuppTable_candidate_pathway_source_counts.csv", index=False)

    plot_directionality(
        summary_df=directionality_summary,
        output_path=FIG_DIR / "SuppFig_RNAseq_directional_pathway_comparison.png",
        pathway_order=PATHWAY_ORDER,
        title="Published maize N-stress DEGs: year x stage pathway directionality",
        figsize=(10.8, 6.4),
    )
    plot_directionality(
        summary_df=core_directionality_summary,
        output_path=FIG_DIR / "SuppFig_RNAseq_directional_core_pathways.png",
        pathway_order=CORE_FIGURE_PATHWAY_ORDER,
        title="Maize N-stress RNA-seq: directional shifts in core N-remobilization pathways",
        figsize=(10.2, 5.5),
    )
    plot_candidate_source_counts(
        summary_df=all_source_candidate_counts,
        output_path=FIG_DIR / "SuppFig_candidate_pathway_counts_all_sources.png",
    )

    print("Figures written:")
    for fig_name in [
        "SuppFig_RNAseq_directional_pathway_comparison.png",
        "SuppFig_RNAseq_directional_core_pathways.png",
        "SuppFig_candidate_pathway_counts_all_sources.png",
    ]:
        print("  ", str(FIG_DIR / fig_name))

    print("RNA-seq DEG rows:", len(deg_df))
    print("RNA-seq unique genes:", gene_df['GeneID'].nunique())
    print("Overlap genes:", overlap_gene_table['GeneID'].nunique())
    print("Remobilization panel genes:", remobilization_panel['GeneID'].nunique())
    print("Pathway-mapped candidate genes:", pathway_candidate_map['GeneID'].nunique())
    print("Top support genes:")
    for row in ranking_table.head(10).itertuples(index=False):
        print(
            "  ",
            row.GeneID,
            "| score",
            row.Support_score,
            "| contrasts",
            row.RNA_contrasts,
            "| amino",
            row.Amino_traits or "None",
            "| soilN",
            row.SoilN_GWAS,
            "| selection",
            row.Selection_populations or "None",
            "| panel",
            row.Panel_categories or "None",
        )
    print("Pathway summary:")
    for row in pathway_candidate_summary.itertuples(index=False):
        print(
            "  ",
            row.Pathway_category,
            "| genes",
            row.n_candidate_genes,
            "| recurrent",
            row.n_recurrent_deg_genes,
            "| focal amino",
            row.n_amino_focal_genes,
            "| soilN",
            row.n_soilN_genes,
            "| selection",
            row.n_selection_genes,
        )
    print("All-source candidate pathway counts:")
    for row in all_source_candidate_counts.itertuples(index=False):
        if row.n_genes == 0:
            continue
        print("  ", row.Pathway_category, "|", row.Source, "| n", row.n_genes)


if __name__ == "__main__":
    main()

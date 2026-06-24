#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


BASE_DIR = Path(__file__).resolve().parent
INTERACTIONS_PATH = BASE_DIR / "string_interactions.tsv"
PROTEINS_PATH = BASE_DIR / "string_protein_annotations.tsv"
OUTPUT_DIR = BASE_DIR / "candidate_network_outputs"


@dataclass(frozen=True)
class ThemeSpec:
    name: str
    seeds: list[str]
    title: str


THEMES = [
    ThemeSpec(
        name="remobilization",
        title="Nitrogen Remobilization",
        seeds=[
            "A0A1D6HNF9",
            "Pco118382",
            "B4FFP2_MAIZE",
            "B4F9J5_MAIZE",
            "A0A1D6GZS8",
            "B4FA80_MAIZE",
            "A0A1D6DW09",
            "A0A1D6QM74",
        ],
    ),
    ThemeSpec(
        name="carbon_support",
        title="Carbon Skeleton Support",
        seeds=[
            "A0A1D6MNB9",
            "Pco104294",
            "Pco080711",
            "A0A1D6MRL4",
        ],
    ),
]

TRANSPORTER_CANDIDATES = [
    ("A0A096TZE2", "Ammonium transporter", "Ammonium transport"),
    ("A0A1D6MR48", "Peptide/nitrate transporter", "Nitrate / peptide transport"),
    ("A0A1D6MWS7", "Transmembrane amino acid transporter", "Amino acid transport"),
    ("A0A1D6MYM4", "Amino acid permease 6", "Amino acid transport"),
    ("A0A1D6MYN4", "Amino acid permease 2/3", "Amino acid transport"),
    ("Sut4", "Sucrose transporter 4", "Source-sink carbon transport"),
    ("A0A1D6MRZ5", "SWEET17a", "Source-sink carbon transport"),
    ("B4FTH1_MAIZE", "Bidirectional SWEET transporter", "Source-sink carbon transport"),
]

COLOR_MAP = {
    "seed": "#d1495b",
    "neighbor": "#b0b7c6",
    "amino_acid": "#00798c",
    "carbon": "#2a6f97",
    "nitrogen": "#3a7d44",
    "transport": "#edae49",
}

LABEL_OVERRIDES = {
    "A0A1D6HNF9": "Aspartate aminotransferase chloro",
    "Pco118382": "Aspartate aminotransferase cyto",
    "A0A1D6DW09": "L-asparaginase 2",
    "B4FFP2_MAIZE": "Glutamate dehydrogenase",
    "B4F9J5_MAIZE": "Glutamate dehydrogenase 2",
    "A0A1D6MRL4": "Quinolinate synthase chloro",
    "A0A1D6MNB9": "PEP carboxylase 4",
    "Pco104294": "Malate dehydrogenase",
    "Pco080711": "ATP-citrate synthase",
}


def read_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    interactions = pd.read_csv(INTERACTIONS_PATH, sep="\t")
    proteins = pd.read_csv(PROTEINS_PATH, sep="\t").rename(
        columns={
            "domain_summary_url": "protein_description",
            "annotation": "domain_summary_url",
        }
    )
    return interactions, proteins


def description_label(description: str, node_id: str) -> str:
    if node_id in LABEL_OVERRIDES:
        return LABEL_OVERRIDES[node_id]

    desc = str(description).split(";")[0].strip().rstrip(".")
    replacements = {
        "chloroplastic/mitochondrial": "chloro/mito",
        "chloroplastic": "chloro",
        "cytoplasmic": "cyto",
        "Putative ": "",
        "Uncharacterized protein ": "",
    }
    for old, new in replacements.items():
        desc = desc.replace(old, new)

    special = {
        "1,2-dihydroxy-3-keto-5-methylthiopentene dioxygenase": "ARD dioxygenase",
        "Aspartate aminotransferase chloroplastic": "Aspartate aminotransferase chloro",
        "Aspartate aminotransferase": "Aspartate aminotransferase",
        "L-aspartate oxidase": "L-aspartate oxidase",
        "Glutamate dehydrogenase": "Glutamate dehydrogenase",
        "Phosphoenolpyruvate carboxylase 4": "PEP carboxylase 4",
        "Malate dehydrogenase": "Malate dehydrogenase",
        "ATP-citrate synthase": "ATP-citrate synthase",
        "Quinolinate synthase chloroplastic": "Quinolinate synthase chloro",
        "Putative isoaspartyl peptidase/L-asparaginase 2": "L-asparaginase 2",
        "Tyrosine aminotransferase": "Tyrosine aminotransferase",
        "Zea CEFD homolog1": "Zea CEFD homolog1",
        "Aspartokinase": "Aspartokinase",
        "Enolase-phosphatase E1": "Enolase-phosphatase E1",
        "S-adenosylmethionine synthase": "SAM synthase",
        "Dihydrolipoamide acetyltransferase component of pyruvate dehydrogenase complex": "Pyruvate dehydrogenase E2",
    }
    for prefix, label in special.items():
        if desc.startswith(prefix):
            if label == "ARD dioxygenase":
                return f"{label} ({node_id[-4:]})"
            return label

    if len(desc) > 34:
        return textwrap.shorten(desc, width=34, placeholder="...")
    return desc or node_id


def wrap_label(label: str, width: int = 18) -> str:
    return textwrap.fill(label, width=width, break_long_words=False)


def build_neighborhood_edges(interactions: pd.DataFrame, seeds: set[str]) -> pd.DataFrame:
    subset = interactions[
        interactions["#node1"].isin(seeds) | interactions["node2"].isin(seeds)
    ].copy()
    subset["source"] = subset["#node1"]
    subset["target"] = subset["node2"]
    subset["edge_key"] = subset.apply(
        lambda row: tuple(sorted((row["source"], row["target"]))), axis=1
    )
    subset = (
        subset.sort_values("combined_score", ascending=False)
        .drop_duplicates("edge_key")
        .drop(columns=["#node1", "node2", "edge_key"])
    )
    return subset.reset_index(drop=True)


def component_layout(graph: nx.Graph, gap: float = 3.8) -> dict[str, tuple[float, float]]:
    if not graph.nodes:
        return {}

    positions: dict[str, tuple[float, float]] = {}
    x_offset = 0.0
    for component in sorted(nx.connected_components(graph), key=len, reverse=True):
        subgraph = graph.subgraph(component).copy()
        sub_pos = nx.spring_layout(subgraph, seed=11, k=1.5, iterations=400)
        xs = [coords[0] for coords in sub_pos.values()]
        ys = [coords[1] for coords in sub_pos.values()]
        width = max(xs) - min(xs) if xs else 0.0
        x_shift = x_offset - min(xs)
        y_shift = -(min(ys) + max(ys)) / 2 if ys else 0.0
        for node, coords in sub_pos.items():
            positions[node] = (coords[0] + x_shift, coords[1] + y_shift)
        x_offset += width + gap
    return positions


def build_graph(theme: ThemeSpec, interactions: pd.DataFrame, proteins: pd.DataFrame) -> tuple[nx.Graph, pd.DataFrame, pd.DataFrame]:
    seeds = set(theme.seeds)
    edge_table = build_neighborhood_edges(interactions, seeds)
    node_ids = set(edge_table["source"]) | set(edge_table["target"])

    protein_index = proteins.set_index("#node")
    node_records = []
    for node_id in sorted(node_ids):
        desc = protein_index.at[node_id, "protein_description"]
        node_records.append(
            {
                "node_id": node_id,
                "label": description_label(desc, node_id),
                "protein_description": desc,
                "node_type": "seed" if node_id in seeds else "neighbor",
            }
        )
    node_table = pd.DataFrame.from_records(node_records)

    graph = nx.Graph()
    for _, row in node_table.iterrows():
        graph.add_node(
            row["node_id"],
            label=row["label"],
            node_type=row["node_type"],
        )
    for _, row in edge_table.iterrows():
        graph.add_edge(
            row["source"],
            row["target"],
            combined_score=float(row["combined_score"]),
        )
    return graph, node_table, edge_table


def draw_network_figure(theme: ThemeSpec, graph: nx.Graph, node_table: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(13.5, 9.0), facecolor="white")
    ax.set_facecolor("white")

    positions = component_layout(graph)
    positions = {
        node: (coords[0] * 1.55, coords[1] * 1.18) for node, coords in positions.items()
    }

    color_lookup = node_table.set_index("node_id")["node_type"].map(
        {"seed": COLOR_MAP["seed"], "neighbor": COLOR_MAP["neighbor"]}
    )
    node_sizes = []
    for node_id in node_table["node_id"]:
        degree = graph.degree(node_id)
        node_sizes.append(1250 + degree * 240)

    nx.draw_networkx_edges(
        graph,
        pos=positions,
        ax=ax,
        width=[1.6 + 4.2 * graph[u][v]["combined_score"] for u, v in graph.edges()],
        edge_color="#8e99a8",
        alpha=0.65,
    )
    nx.draw_networkx_nodes(
        graph,
        pos=positions,
        ax=ax,
        node_size=node_sizes,
        node_color=color_lookup.tolist(),
        edgecolors="#25313f",
        linewidths=1.2,
    )
    nx.draw_networkx_labels(
        graph,
        pos=positions,
        labels={node: wrap_label(data["label"]) for node, data in graph.nodes(data=True)},
        ax=ax,
        font_size=8.5,
        font_weight="bold",
        font_family="DejaVu Sans",
    )

    xs = [coords[0] for coords in positions.values()]
    ys = [coords[1] for coords in positions.values()]
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(min(ys) - 0.6, max(ys) + 0.6)
    ax.set_axis_off()

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Theme seed gene",
            markerfacecolor=COLOR_MAP["seed"],
            markeredgecolor="#25313f",
            markersize=11,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="First-neighbor protein from STRING",
            markerfacecolor=COLOR_MAP["neighbor"],
            markeredgecolor="#25313f",
            markersize=11,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Larger circle = more STRING connections",
            markerfacecolor="#d9dde5",
            markeredgecolor="#25313f",
            markersize=18,
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.015),
        fontsize=10,
    )
    fig.tight_layout(rect=(0.01, 0.07, 0.99, 0.99))
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def draw_transporter_panel(proteins: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    fig, ax = plt.subplots(1, 1, figsize=(13.5, 7.2), facecolor="white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 8.2)
    ax.axis("off")

    protein_index = proteins.set_index("#node")
    records = []
    category_colors = {
        "Ammonium transport": "#3a7d44",
        "Nitrate / peptide transport": "#5c946e",
        "Amino acid transport": "#edae49",
        "Source-sink carbon transport": "#2a6f97",
    }

    xs = [1.0, 4.2, 7.4, 10.6]
    ys = [6.5, 4.4, 2.3]
    for idx, (node_id, short_label, category) in enumerate(TRANSPORTER_CANDIDATES):
        x = xs[idx % 4]
        y = ys[idx // 4]
        desc = protein_index.at[node_id, "protein_description"]
        box = FancyBboxPatch(
            (x - 0.85, y - 0.55),
            2.4,
            1.15,
            boxstyle="round,pad=0.18,rounding_size=0.12",
            linewidth=1.1,
            edgecolor="#25313f",
            facecolor=category_colors[category],
            alpha=0.94,
        )
        ax.add_patch(box)
        ax.text(
            x + 0.35,
            y + 0.17,
            wrap_label(short_label, width=18),
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="white",
        )
        ax.text(
            x + 0.35,
            y - 0.28,
            wrap_label(category, width=20),
            ha="center",
            va="center",
            fontsize=8.7,
            color="white",
        )
        records.append(
            {
                "node_id": node_id,
                "label": short_label,
                "category": category,
                "protein_description": desc,
                "string_interactions_in_export": 0,
            }
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            label=label,
            markerfacecolor=color,
            markeredgecolor="#25313f",
            markersize=12,
        )
        for label, color in category_colors.items()
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.03),
        fontsize=10,
    )
    ax.text(
        0.8,
        0.75,
        "Transporter candidates are shown as a panel because they carry no direct STRING interactions in this export.",
        fontsize=10,
        color="#25313f",
        ha="left",
        va="center",
    )
    fig.tight_layout(rect=(0.01, 0.08, 0.99, 0.99))
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame.from_records(records)


def write_summary(theme_tables: list[tuple[ThemeSpec, pd.DataFrame, pd.DataFrame]], transporter_table: pd.DataFrame) -> None:
    lines = [
        "# Candidate STRING Themes for Grain Nitrogen / Grain Filling",
        "",
    ]
    for theme, node_table, edge_table in theme_tables:
        seed_count = int((node_table["node_type"] == "seed").sum())
        neighbor_count = int((node_table["node_type"] == "neighbor").sum())
        lines.extend(
            [
                f"## {theme.title}",
                f"- Connected seed genes: {seed_count}",
                f"- First-neighbor proteins: {neighbor_count}",
                f"- Unique STRING interactions: {len(edge_table)}",
                "",
            ]
        )
    lines.extend(
        [
            "## Transporter Candidate Panel",
            f"- Candidate transport genes: {len(transporter_table)}",
            "- Direct STRING interactions in this export: 0",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines))


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    interactions, proteins = read_tables()

    theme_tables: list[tuple[ThemeSpec, pd.DataFrame, pd.DataFrame]] = []
    for theme in THEMES:
        graph, node_table, edge_table = build_graph(theme, interactions, proteins)
        draw_network_figure(theme, graph, node_table, OUTPUT_DIR / theme.name)
        node_table.sort_values(["node_type", "label"]).to_csv(
            OUTPUT_DIR / f"{theme.name}_nodes.tsv", sep="\t", index=False
        )
        edge_table.sort_values(
            ["combined_score", "source", "target"], ascending=[False, True, True]
        ).to_csv(OUTPUT_DIR / f"{theme.name}_edges.tsv", sep="\t", index=False)
        theme_tables.append((theme, node_table, edge_table))

    transporter_table = draw_transporter_panel(proteins, OUTPUT_DIR / "transport_candidates")
    transporter_table.to_csv(
        OUTPUT_DIR / "transport_candidates.tsv", sep="\t", index=False
    )
    write_summary(theme_tables, transporter_table)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import textwrap
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib.lines import Line2D


BASE_DIR = Path(__file__).resolve().parent
INTERACTIONS_PATH = BASE_DIR / "string_interactions.tsv"
PROTEINS_PATH = BASE_DIR / "string_protein_annotations.tsv"
FUNCTIONS_PATH = BASE_DIR / "string_functional_annotations.tsv"
OUTPUT_DIR = BASE_DIR / "network_outputs"


@dataclass(frozen=True)
class NetworkBuild:
    name: str
    title: str
    node_table: pd.DataFrame
    edge_table: pd.DataFrame
    note: str | None = None


def read_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    interactions = pd.read_csv(INTERACTIONS_PATH, sep="\t")
    proteins = pd.read_csv(PROTEINS_PATH, sep="\t").rename(
        columns={
            # STRING's header order is swapped in this export.
            "domain_summary_url": "protein_description",
            "annotation": "domain_summary_url",
        }
    )
    functions = pd.read_csv(FUNCTIONS_PATH, sep="\t")
    return interactions, proteins, functions


def find_nodes_by_keyword(
    proteins: pd.DataFrame,
    functions: pd.DataFrame,
    protein_pattern: str,
    function_pattern: str | None = None,
) -> tuple[set[str], set[str], set[str]]:
    function_pattern = function_pattern or protein_pattern

    protein_hits = proteins.loc[
        proteins["protein_description"].astype(str).str.contains(
            protein_pattern, case=False, regex=True, na=False
        ),
        "#node",
    ]
    function_hits = functions.loc[
        functions["term description"].astype(str).str.contains(
            function_pattern, case=False, regex=True, na=False
        ),
        "#node",
    ]

    protein_nodes = set(protein_hits)
    function_nodes = set(function_hits)
    return protein_nodes, function_nodes, protein_nodes | function_nodes


def build_node_table(
    node_ids: Iterable[str],
    proteins: pd.DataFrame,
    protein_nodes: set[str],
    function_nodes: set[str],
) -> pd.DataFrame:
    node_ids = sorted(set(node_ids))
    protein_index = proteins.set_index("#node")
    records = []

    for node_id in node_ids:
        row = protein_index.loc[node_id] if node_id in protein_index.index else None
        description = row["protein_description"] if row is not None else "Unannotated"
        aliases = row["other_names_and_aliases"] if row is not None else ""
        short_label = description_label(description, aliases, node_id)
        origin = []
        if node_id in protein_nodes:
            origin.append("protein_annotation")
        if node_id in function_nodes:
            origin.append("functional_annotation")
        records.append(
            {
                "node_id": node_id,
                "label": short_label,
                "protein_description": description,
                "match_origin": ",".join(origin) if origin else "none",
            }
        )

    node_table = pd.DataFrame.from_records(records)
    duplicated = node_table["label"].duplicated(keep=False)
    for idx in node_table.index[duplicated]:
        node_table.at[idx, "label"] = (
            f"{node_table.at[idx, 'label']} ({node_table.at[idx, 'node_id'][-4:]})"
        )
    return node_table


def alias_label(aliases: str, node_id: str) -> str:
    candidates = [part.strip() for part in str(aliases).split(",") if part.strip()]
    for candidate in candidates:
        upper = candidate.upper()
        if (
            candidate == node_id
            or candidate.startswith("4577.")
            or candidate.startswith("ZEAMMB73_")
            or candidate.startswith("NP_")
            or candidate.startswith("XP_")
            or candidate.startswith("XM_")
            or candidate.startswith("NM_")
            or candidate.startswith("ZMA:")
            or upper.endswith("_MAIZE")
            or upper.startswith("EC ")
            or candidate.replace(".", "", 1).isdigit()
        ):
            continue
        if len(candidate) <= 40:
            return candidate
    return node_id


def description_label(description: str, aliases: str, node_id: str) -> str:
    desc = str(description).split(";")[0].strip().rstrip(".")
    desc = desc.replace("chloroplastic/mitochondrial", "chloro/mito")
    desc = desc.replace("chloroplastic", "chloro")
    desc = desc.replace("cytoplasmic", "cyto")
    desc = desc.replace("Putative ", "")
    desc = desc.replace("Uncharacterized protein", "").strip(" -")

    if desc and desc.lower() not in {"unannotated", "protein"}:
        return desc
    return alias_label(aliases, node_id)


def wrap_label(label: str, width: int = 24) -> str:
    return textwrap.fill(label, width=width, break_long_words=False)


def component_layout(graph: nx.Graph) -> dict[str, tuple[float, float]]:
    if not graph.nodes:
        return {}

    positions: dict[str, tuple[float, float]] = {}
    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    x_offset = 0.0
    gap = 4.2

    for component in components:
        subgraph = graph.subgraph(component).copy()
        sub_pos = nx.spring_layout(subgraph, seed=7, k=1.7, iterations=400)
        xs = [coords[0] for coords in sub_pos.values()]
        ys = [coords[1] for coords in sub_pos.values()]
        width = max(xs) - min(xs) if xs else 0.0
        x_shift = x_offset - min(xs)
        y_shift = -(min(ys) + max(ys)) / 2 if ys else 0.0

        for node, coords in sub_pos.items():
            positions[node] = (coords[0] + x_shift, coords[1] + y_shift)

        x_offset += width + gap

    return positions


def deduplicate_edges(interactions: pd.DataFrame) -> pd.DataFrame:
    edges = interactions.copy()
    edges["source"] = edges["#node1"]
    edges["target"] = edges["node2"]
    edges["edge_key"] = edges.apply(
        lambda row: tuple(sorted((row["source"], row["target"]))), axis=1
    )
    edges = (
        edges.sort_values("combined_score", ascending=False)
        .drop_duplicates("edge_key")
        .drop(columns=["#node1", "node2", "edge_key"])
    )
    return edges[
        [
            "source",
            "target",
            "combined_score",
            "coexpression",
            "experimentally_determined_interaction",
            "database_annotated",
            "automated_textmining",
        ]
    ].reset_index(drop=True)


def subgraph_edges(interactions: pd.DataFrame, node_ids: set[str]) -> pd.DataFrame:
    subset = interactions[
        interactions["#node1"].isin(node_ids) & interactions["node2"].isin(node_ids)
    ]
    return deduplicate_edges(subset)


def neighborhood_edges(interactions: pd.DataFrame, seed_nodes: set[str]) -> pd.DataFrame:
    subset = interactions[
        interactions["#node1"].isin(seed_nodes) | interactions["node2"].isin(seed_nodes)
    ]
    return deduplicate_edges(subset)


def build_networks() -> list[NetworkBuild]:
    interactions, proteins, functions = read_tables()

    asp_protein, asp_function, asp_nodes = find_nodes_by_keyword(
        proteins=proteins,
        functions=functions,
        protein_pattern=r"aspartate",
        function_pattern=r"aspartate",
    )
    asp_edge_table = neighborhood_edges(interactions, asp_nodes)
    asp_path_nodes = set(asp_edge_table["source"]) | set(asp_edge_table["target"])
    asp_node_table = build_node_table(
        asp_path_nodes, proteins, asp_protein, asp_function
    )

    return [
        NetworkBuild(
            name="aspartate",
            title="Aspartate Seed Genes + STRING Neighborhood",
            node_table=asp_node_table,
            edge_table=asp_edge_table,
            note=(
                f"{len(asp_nodes & asp_path_nodes)} connected aspartate seed genes; "
                f"{len(asp_path_nodes - asp_nodes)} first-neighbor proteins; "
                f"{len(asp_edge_table)} unique STRING interactions"
            ),
        ),
    ]


def node_colors(node_table: pd.DataFrame) -> list[str]:
    colors = []
    for origin in node_table["match_origin"]:
        if origin == "protein_annotation,functional_annotation":
            colors.append("#d1495b")
        elif origin == "protein_annotation":
            colors.append("#edae49")
        elif origin == "functional_annotation":
            colors.append("#00798c")
        else:
            colors.append("#9ea3b0")
    return colors


def draw_network(ax: plt.Axes, network: NetworkBuild) -> None:
    graph = nx.Graph()
    node_table = network.node_table.copy()

    for _, row in node_table.iterrows():
        graph.add_node(
            row["node_id"],
            label=row["label"],
            match_origin=row["match_origin"],
        )

    for _, row in network.edge_table.iterrows():
        graph.add_edge(
            row["source"],
            row["target"],
            combined_score=float(row["combined_score"]),
        )

    positions = component_layout(graph)
    positions = {
        node: (coords[0] * 1.45, coords[1] * 1.20) for node, coords in positions.items()
    }

    node_sizes = []
    for node_id in node_table["node_id"]:
        degree = graph.degree(node_id) if node_id in graph else 0
        node_sizes.append(1300 + degree * 260)

    nx.draw_networkx_edges(
        graph,
        pos=positions,
        ax=ax,
        width=[1.5 + 4.5 * graph[u][v]["combined_score"] for u, v in graph.edges()],
        edge_color="#6c757d",
        alpha=0.55,
    )
    nx.draw_networkx_nodes(
        graph,
        pos=positions,
        ax=ax,
        node_size=node_sizes,
        node_color=node_colors(node_table),
        edgecolors="#1f2933",
        linewidths=1.1,
    )
    nx.draw_networkx_labels(
        graph,
        pos=positions,
        labels={node: wrap_label(data["label"], width=20) for node, data in graph.nodes(data=True)},
        ax=ax,
        font_size=8,
        font_weight="bold",
        font_family="DejaVu Sans",
    )

    if positions:
        xs = [coords[0] for coords in positions.values()]
        ys = [coords[1] for coords in positions.values()]
        ax.set_xlim(min(xs) - 0.60, max(xs) + 0.60)
        ax.set_ylim(min(ys) - 0.45, max(ys) + 0.45)
    ax.set_axis_off()


def write_tables(networks: list[NetworkBuild]) -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    obsolete = [
        OUTPUT_DIR / "amino_acid_permease_nodes.tsv",
        OUTPUT_DIR / "amino_acid_permease_edges.tsv",
    ]
    for path in obsolete:
        if path.exists():
            path.unlink()
    for network in networks:
        network.node_table.sort_values("label").to_csv(
            OUTPUT_DIR / f"{network.name}_nodes.tsv", sep="\t", index=False
        )
        network.edge_table.sort_values(
            ["combined_score", "source", "target"], ascending=[False, True, True]
        ).to_csv(OUTPUT_DIR / f"{network.name}_edges.tsv", sep="\t", index=False)


def write_summary(networks: list[NetworkBuild]) -> None:
    lines = [
        "# STRING Network Summary",
        "",
        "Built from `string_interactions.tsv`, `string_protein_annotations.tsv`, and `string_functional_annotations.tsv`.",
        "",
    ]
    for network in networks:
        lines.extend(
            [
                f"## {network.title}",
                f"- Nodes: {len(network.node_table)}",
                f"- Unique edges: {len(network.edge_table)}",
                f"- Note: {network.note or 'None'}",
                "",
            ]
        )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines))


def draw_figure(networks: list[NetworkBuild]) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(16, 10), facecolor="white")
    ax.set_facecolor("white")
    draw_network(ax, networks[0])

    color_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Matched in protein annotation + functional terms",
            markerfacecolor="#d1495b",
            markeredgecolor="#1f2933",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Matched in protein annotation only",
            markerfacecolor="#edae49",
            markeredgecolor="#1f2933",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Matched in functional annotation only",
            markerfacecolor="#00798c",
            markeredgecolor="#1f2933",
            markersize=10,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Neighbor added from STRING interactions",
            markerfacecolor="#9ea3b0",
            markeredgecolor="#1f2933",
            markersize=10,
        ),
    ]

    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Larger circle = more STRING connections",
            markerfacecolor="#d9dde5",
            markeredgecolor="#1f2933",
            markersize=18,
        )
    ]

    fig.legend(
        handles=color_handles + size_handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.03),
        fontsize=10,
    )
    fig.tight_layout(rect=(0.01, 0.09, 0.99, 0.99))
    fig.savefig(OUTPUT_DIR / "string_theme_networks.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    networks = build_networks()
    write_tables(networks)
    write_summary(networks)
    draw_figure(networks)


if __name__ == "__main__":
    main()

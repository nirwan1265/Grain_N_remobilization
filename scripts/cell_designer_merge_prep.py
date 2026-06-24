#!/usr/bin/env python3
"""Inventory and merge-prep helper for older SBML/CellDesigner XML models.

This script is intentionally tolerant of older BioModels exports that may not
parse cleanly with standard XML parsers. It uses lightweight regex extraction
to summarize models, identify candidate bridge metabolites, and generate a
first-pass merge plan.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from pathlib import Path


MODEL_NAME_RE = re.compile(r"<model\b[^>]*\bname=\"([^\"]+)\"")
MODEL_ID_RE = re.compile(r"<model\b[^>]*\bid=\"([^\"]+)\"")
SPECIES_RE = re.compile(r"<species\b([^>]*)/>|<species\b([^>]*)>")
REACTION_RE = re.compile(r"<reaction\b([^>]*)>")
PARAMETER_RE = re.compile(r"<parameter\b([^>]*)/>|<parameter\b([^>]*)>")
ATTR_RE = {
    "id": re.compile(r"\bid=\"([^\"]+)\""),
    "name": re.compile(r"\bname=\"([^\"]*)\""),
    "compartment": re.compile(r"\bcompartment=\"([^\"]+)\""),
    "initialConcentration": re.compile(r"\binitialConcentration=\"([^\"]+)\""),
    "value": re.compile(r"\bvalue=\"([^\"]+)\""),
}


BRIDGE_SYNONYMS = {
    "adp": "ADP",
    "atp": "ATP",
    "nadp": "NADP",
    "nadph": "NADPH",
    "asp": "Aspartate",
    "aspartate": "Aspartate",
    "aspp": "Aspartyl_phosphate",
    "aspartylp": "Aspartyl_phosphate",
    "aspartylphosphate": "Aspartyl_phosphate",
    "asadh": "Aspartate_semialdehyde",
    "aspartatesemialdehyde": "Aspartate_semialdehyde",
    "hser": "Homoserine",
    "homoserine": "Homoserine",
    "thr": "Threonine",
    "threonine": "Threonine",
    "lys": "Lysine",
    "lysine": "Lysine",
    "glu": "Glutamate",
    "glutamate": "Glutamate",
    "gln": "Glutamine",
    "glutamine": "Glutamine",
    "nh4": "Ammonium",
    "ammonium": "Ammonium",
    "kg": "2_oxoglutarate",
    "akg": "2_oxoglutarate",
    "alphaketoglutarate": "2_oxoglutarate",
    "oxoglutarate": "2_oxoglutarate",
    "oaa": "Oxaloacetate",
    "oxaloacetate": "Oxaloacetate",
    "pyruvate": "Pyruvate",
    "accoa": "Acetyl_CoA",
    "acetylcoa": "Acetyl_CoA",
    "coa": "CoA",
}


def clean_text(text: str) -> str:
    return " ".join(re.sub(r"<[^>]+>", " ", text).split())


def norm_name(name: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "", name.lower())
    return BRIDGE_SYNONYMS.get(key, name)


def parse_attrs(attr_blob: str) -> dict[str, str]:
    data: dict[str, str] = {}
    for key, pattern in ATTR_RE.items():
        match = pattern.search(attr_blob)
        if match:
            data[key] = html.unescape(match.group(1))
    return data


def infer_organism(text: str) -> str:
    lowered = text.lower()
    if "arabidopsis" in lowered:
        return "Arabidopsis"
    if "escherichia coli" in lowered or "e. coli" in lowered or "ecoli" in lowered:
        return "E. coli"
    if "lactococcus" in lowered:
        return "Lactococcus lactis"
    return "Unknown_or_unspecified"


def recommended_role(model_name: str, organism: str) -> str:
    if "Curien2009_Aspartate_Metabolism" in model_name:
        return "Backbone_for_plant_branch"
    if organism == "Arabidopsis":
        return "Plant_module_candidate"
    if organism == "E. coli":
        return "Bacterial_connector_only"
    if organism == "Lactococcus lactis":
        return "Hold_out_for_now"
    return "Needs_manual_curation"


def parse_model(xml_path: Path) -> dict:
    text = xml_path.read_text(errors="ignore")
    flattened_head = clean_text(text[:16000])

    model_name_match = MODEL_NAME_RE.search(text)
    model_id_match = MODEL_ID_RE.search(text)
    model_name = model_name_match.group(1) if model_name_match else xml_path.stem
    model_id = model_id_match.group(1) if model_id_match else ""
    organism = infer_organism(flattened_head)

    species = []
    for match in SPECIES_RE.finditer(text):
        attrs = parse_attrs(match.group(1) or match.group(2) or "")
        if attrs:
            species.append(attrs)

    reactions = []
    for match in REACTION_RE.finditer(text):
        attrs = parse_attrs(match.group(1))
        if attrs:
            reactions.append(attrs)

    parameters = []
    for match in PARAMETER_RE.finditer(text):
        attrs = parse_attrs(match.group(1) or match.group(2) or "")
        if attrs:
            parameters.append(attrs)

    bridge_metabolites = defaultdict(list)
    for item in species:
        if "name" not in item:
            continue
        normalized = norm_name(item["name"])
        if normalized in BRIDGE_SYNONYMS.values():
            bridge_metabolites[normalized].append(item["name"])

    return {
        "file": xml_path.name,
        "path": str(xml_path),
        "model_id": model_id,
        "model_name": model_name,
        "organism": organism,
        "species": species,
        "reactions": reactions,
        "parameters": parameters,
        "bridge_metabolites": bridge_metabolites,
        "recommended_role": recommended_role(model_name, organism),
    }


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_overlap(models: list[dict]) -> list[dict]:
    overlap = defaultdict(list)
    for model in models:
        seen = set()
        for species in model["species"]:
            name = species.get("name")
            if not name:
                continue
            normalized = norm_name(name)
            if normalized in BRIDGE_SYNONYMS.values() and normalized not in seen:
                overlap[normalized].append(model["model_name"])
                seen.add(normalized)

    rows = []
    for metabolite, model_names in sorted(
        overlap.items(), key=lambda item: (-len(item[1]), item[0])
    ):
        rows.append(
            {
                "bridge_metabolite": metabolite,
                "model_count": len(model_names),
                "models": "; ".join(model_names),
            }
        )
    return rows


def write_markdown_summary(path: Path, models: list[dict], overlap_rows: list[dict]) -> None:
    backbone = next(
        (m for m in models if "Curien2009_Aspartate_Metabolism" in m["model_name"]),
        None,
    )

    lines = [
        "# CellDesigner Merge Prep",
        "",
        "## What makes sense",
        "",
        "Yes: combining these into one working pathway makes sense as a curation project.",
        "The current files should not be pasted together directly for quantitative simulation, because they mix organisms and kinetic assumptions.",
        "",
        "## Model snapshot",
        "",
    ]

    for model in models:
        bridges = ", ".join(sorted(model["bridge_metabolites"])) or "none"
        lines.append(
            f"- `{model['model_name']}` ({model['organism']}): "
            f"{len(model['species'])} species, {len(model['reactions'])} reactions, "
            f"role `{model['recommended_role']}`, bridge metabolites: {bridges}"
        )

    lines.extend(
        [
            "",
            "## Recommended merge order",
            "",
            "1. Use `Curien2009_Aspartate_Metabolism` as the backbone for the amino-acid branch.",
            "2. Treat the ammonium assimilation and TCA models as donor modules for nitrogen and carbon supply logic, not as direct one-click kinetic inserts.",
            "3. Add bridge metabolites manually: `Aspartate`, `Glutamate`, `Glutamine`, `Ammonium`, `2_oxoglutarate`, and `Oxaloacetate`.",
            "4. Replace organism-specific enzyme pools and rate constants only after choosing a target species for the final model.",
            "",
            "## Best immediate simulation targets",
            "",
            "- Increase aspartate kinase capacity in the Curien model by raising enzyme species `AK1`, `AK2`, `AKI`, or `AKII`.",
            "- Increase aspartate kinase capacity in the Chassagnole model by raising `vm11` and/or `vm13`.",
            "- Reduce glutamate in the Bruggeman model by lowering species `GLU` or by reducing upstream `vgdh` / `vgog` capacity.",
            "- Impose nitrogen limitation in the Bruggeman model by lowering `NH4` and tracking the response of `GLU`, `GLN`, and assimilation fluxes.",
            "",
            "## Strong caution",
            "",
            "- The Curien backbone is Arabidopsis.",
            "- Bruggeman and Singh are E. coli.",
            "- Hoefnagel is Lactococcus lactis.",
            "- A single merged file can still be useful for hypothesis generation, but any quantitative claim will need species-specific re-parameterization.",
            "",
            "## Shared bridge metabolites detected",
            "",
        ]
    )

    for row in overlap_rows:
        if row["model_count"] < 2:
            continue
        lines.append(
            f"- `{row['bridge_metabolite']}` appears in {row['model_count']} models: {row['models']}"
        )

    if backbone is not None:
        lines.extend(
            [
                "",
                "## Why the Curien model is the best starting point",
                "",
                "It already contains the aspartate-derived branch logic you care about, including multiple aspartate kinase isoforms and downstream threonine / lysine competition.",
            ]
        )

    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        default="/Users/nirwantandukar/Documents/Research/data/Cell_designer/Grain_N",
        help="Directory containing SBML / CellDesigner XML files.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/cell_designer_merge",
        help="Directory for generated reports.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    models = [parse_model(path) for path in sorted(input_dir.glob("*.xml"))]

    inventory_rows = []
    for model in models:
        inventory_rows.append(
            {
                "file": model["file"],
                "model_name": model["model_name"],
                "organism": model["organism"],
                "species_count": len(model["species"]),
                "reaction_count": len(model["reactions"]),
                "parameter_count": len(model["parameters"]),
                "recommended_role": model["recommended_role"],
                "bridge_metabolites": ", ".join(sorted(model["bridge_metabolites"])),
                "source_path": model["path"],
            }
        )

    overlap_rows = summarize_overlap(models)

    write_tsv(
        output_dir / "model_inventory.tsv",
        inventory_rows,
        [
            "file",
            "model_name",
            "organism",
            "species_count",
            "reaction_count",
            "parameter_count",
            "recommended_role",
            "bridge_metabolites",
            "source_path",
        ],
    )
    write_tsv(
        output_dir / "bridge_metabolite_overlap.tsv",
        overlap_rows,
        ["bridge_metabolite", "model_count", "models"],
    )
    write_markdown_summary(output_dir / "merge_plan.md", models, overlap_rows)


if __name__ == "__main__":
    main()

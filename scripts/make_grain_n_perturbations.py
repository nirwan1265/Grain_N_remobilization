#!/usr/bin/env python3
"""Generate perturbation XML files for the Grain_N merged model."""

from __future__ import annotations

import argparse
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path


SBML_NS = "http://www.sbml.org/sbml/level2/version4"


def q(tag: str) -> str:
    return f"{{{SBML_NS}}}{tag}"


def find_species(model: ET.Element, species_id: str) -> ET.Element:
    for sp in model.find(q("listOfSpecies")) or []:
        if sp.get("id") == species_id:
            return sp
    raise KeyError(species_id)


def find_parameter_elems(model: ET.Element, param_id: str) -> list[ET.Element]:
    matches = []
    for parameter in model.iter(q("parameter")):
        if parameter.get("id") == param_id:
            matches.append(parameter)
    return matches


def set_initial(species: ET.Element, value: float) -> None:
    if "initialConcentration" in species.attrib:
        species.set("initialConcentration", f"{value:.15g}")
    elif "initialAmount" in species.attrib:
        species.set("initialAmount", f"{value:.15g}")
    else:
        species.set("initialConcentration", f"{value:.15g}")


def get_initial(species: ET.Element) -> float:
    if "initialConcentration" in species.attrib:
        return float(species.get("initialConcentration"))
    if "initialAmount" in species.attrib:
        return float(species.get("initialAmount"))
    raise ValueError(f"No initial value set for {species.get('id')}")


def get_parameter_value(model: ET.Element, param_id: str) -> float:
    matches = find_parameter_elems(model, param_id)
    if not matches:
        raise KeyError(param_id)
    return float(matches[0].get("value"))


def set_parameter_value(model: ET.Element, param_id: str, value: float) -> int:
    matches = find_parameter_elems(model, param_id)
    if not matches:
        raise KeyError(param_id)
    for parameter in matches:
        parameter.set("value", f"{value:.15g}")
    return len(matches)


def scale_species_initial(model: ET.Element, species_id: str, factor: float) -> None:
    sp = find_species(model, species_id)
    set_initial(sp, get_initial(sp) * factor)


def scale_parameter(model: ET.Element, param_id: str, factor: float) -> None:
    set_parameter_value(model, param_id, get_parameter_value(model, param_id) * factor)


SCENARIOS: dict[str, list[tuple[str, str, float]]] = {
    "asp_kinase_x2": [
        ("species", "asp_AK1", 2.0),
        ("species", "asp_AK2", 2.0),
        ("species", "asp_AKHSDHI", 2.0),
        ("species", "asp_AKHSDHII", 2.0),
    ],
    "asp_kinase_half": [
        ("species", "asp_AK1", 0.5),
        ("species", "asp_AK2", 0.5),
        ("species", "asp_AKHSDHI", 0.5),
        ("species", "asp_AKHSDHII", 0.5),
    ],
    "asadh_x2": [
        ("species", "asp_ASADH", 2.0),
    ],
    "asadh_half": [
        ("species", "asp_ASADH", 0.5),
    ],
    "aspat_capacity_x2": [
        ("parameter", "Vf_AspAT", 2.0),
        ("parameter", "Vr_AspAT", 2.0),
    ],
    "aspat_capacity_half": [
        ("parameter", "Vf_AspAT", 0.5),
        ("parameter", "Vr_AspAT", 0.5),
    ],
    "nh4_half": [
        ("species", "ext_NH4", 0.5),
    ],
    "nh4_x2": [
        ("species", "ext_NH4", 2.0),
    ],
    "nh4_import_half": [
        ("parameter", "nh4_Vmax_NH4_import", 0.5),
    ],
    "nh4_import_x2": [
        ("parameter", "nh4_Vmax_NH4_import", 2.0),
    ],
    "glu_half": [
        ("species", "nh4_GLU", 0.5),
    ],
    "glu_x2": [
        ("species", "nh4_GLU", 2.0),
    ],
    "gln_half": [
        ("species", "nh4_GLN", 0.5),
    ],
    "gln_x2": [
        ("species", "nh4_GLN", 2.0),
    ],
    "asns_half": [
        ("parameter", "Vmax_ASNS", 0.5),
    ],
    "asns_x2": [
        ("parameter", "Vmax_ASNS", 2.0),
    ],
    "omega_half": [
        ("parameter", "Vmax_OMEGA", 0.5),
    ],
    "omega_x2": [
        ("parameter", "Vmax_OMEGA", 2.0),
    ],
    "p5cs_half": [
        ("parameter", "Vmax_P5CS", 0.5),
    ],
    "p5cs_x2": [
        ("parameter", "Vmax_P5CS", 2.0),
    ],
    "p5cr_half": [
        ("parameter", "Vmax_P5CR", 0.5),
    ],
    "p5cr_x2": [
        ("parameter", "Vmax_P5CR", 2.0),
    ],
    "gs_half": [
        ("parameter", "nh4_Vgs", 0.5),
    ],
    "gs_x2": [
        ("parameter", "nh4_Vgs", 2.0),
    ],
    "gogat_half": [
        ("parameter", "nh4_Vgog", 0.5),
    ],
    "gogat_x2": [
        ("parameter", "nh4_Vgog", 2.0),
    ],
    "gs_gogat_half": [
        ("parameter", "nh4_Vgs", 0.5),
        ("parameter", "nh4_Vgog", 0.5),
    ],
    "gs_gogat_x2": [
        ("parameter", "nh4_Vgs", 2.0),
        ("parameter", "nh4_Vgog", 2.0),
    ],
    "gdh_half": [
        ("parameter", "nh4_Vgdh", 0.5),
    ],
    "gdh_x2": [
        ("parameter", "nh4_Vgdh", 2.0),
    ],
    "n_sink_glu_0p75": [
        ("parameter", "nh4_Vgludem", 0.75),
    ],
    "n_sink_gln_0p75": [
        ("parameter", "nh4_Vglndem", 0.75),
    ],
    "n_sink_both_0p75": [
        ("parameter", "nh4_Vgludem", 0.75),
        ("parameter", "nh4_Vglndem", 0.75),
    ],
    "n_sink_both_0p5": [
        ("parameter", "nh4_Vgludem", 0.5),
        ("parameter", "nh4_Vglndem", 0.5),
    ],
    "n_demand_x1p25": [
        ("parameter", "nh4_Vgludem", 1.25),
        ("parameter", "nh4_Vglndem", 1.25),
    ],
    "akg_pool_x2": [
        ("species", "tca_akg", 2.0),
    ],
    "akg_pool_half": [
        ("species", "tca_akg", 0.5),
    ],
    "kdh_capacity_half": [
        ("parameter", "tca_Vf_kdh", 0.5),
        ("parameter", "tca_Vr_kdh", 0.5),
    ],
    "kdh_capacity_x2": [
        ("parameter", "tca_Vf_kdh", 2.0),
        ("parameter", "tca_Vr_kdh", 2.0),
    ],
    "pc_bridge_half": [
        ("parameter", "Vmax_PC_bridge", 0.5),
    ],
    "pc_bridge_x2": [
        ("parameter", "Vmax_PC_bridge", 2.0),
    ],
    "pepck_bridge_half": [
        ("parameter", "Vf_PEPCK_bridge", 0.5),
        ("parameter", "Vr_PEPCK_bridge", 0.5),
    ],
    "pepck_bridge_x2": [
        ("parameter", "Vf_PEPCK_bridge", 2.0),
        ("parameter", "Vr_PEPCK_bridge", 2.0),
    ],
    "kbiosyn_half": [
        ("parameter", "tca_kbiosyn", 0.5),
    ],
    "kbiosyn_x2": [
        ("parameter", "tca_kbiosyn", 2.0),
    ],
    "asp_pool_x2": [
        ("species", "asp_Asp", 2.0),
    ],
    "asp_pool_half": [
        ("species", "asp_Asp", 0.5),
    ],
}


def apply_scenario(model: ET.Element, scenario: str) -> str:
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario: {scenario}")

    notes = []
    for kind, target, factor in SCENARIOS[scenario]:
        if kind == "species":
            scale_species_initial(model, target, factor)
            notes.append(f"{target} x{factor:g} (initial)")
        elif kind == "parameter":
            scale_parameter(model, target, factor)
            notes.append(f"{target} x{factor:g}")
        else:
            raise ValueError(kind)
    return "; ".join(notes)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_xml", help="Path to merged_nitrogen_connected.xml")
    parser.add_argument("output_dir", help="Directory to write scenario XML files")
    args = parser.parse_args()

    base_xml = Path(args.base_xml).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ET.register_namespace("", SBML_NS)
    ET.register_namespace("html", "http://www.w3.org/1999/xhtml")
    ET.register_namespace("math", "http://www.w3.org/1998/Math/MathML")
    ET.register_namespace("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    ET.register_namespace("dc", "http://purl.org/dc/elements/1.1/")
    ET.register_namespace("dcterms", "http://purl.org/dc/terms/")
    ET.register_namespace("vCard", "http://www.w3.org/2001/vcard-rdf/3.0#")
    ET.register_namespace("bqbiol", "http://biomodels.net/biology-qualifiers/")
    ET.register_namespace("bqmodel", "http://biomodels.net/model-qualifiers/")

    scenarios = list(SCENARIOS.keys())

    shutil.copy2(base_xml, output_dir / "baseline_merged_nitrogen_connected.xml")
    manifest_rows = [("baseline", "unmodified baseline model")]

    for scenario in scenarios:
        tree = ET.parse(base_xml)
        root = tree.getroot()
        model = root.find(q("model"))
        if model is None:
            raise RuntimeError("No SBML model element found")
        note = apply_scenario(model, scenario)
        model.set("name", f"{model.get('name')} [{scenario}]")
        outfile = output_dir / f"{scenario}_merged_nitrogen_connected.xml"
        tree.write(outfile, encoding="UTF-8", xml_declaration=True)
        print(f"Wrote {outfile.name}: {note}")
        manifest_rows.append((scenario, note))

    manifest_path = output_dir / "scenario_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8") as fh:
        fh.write("scenario\tdescription\n")
        for scenario, note in manifest_rows:
            fh.write(f"{scenario}\t{note}\n")
    print(f"Wrote {manifest_path.name}")


if __name__ == "__main__":
    main()

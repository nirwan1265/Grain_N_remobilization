#!/usr/bin/env python3
"""Directly patch the merged Grain_N SBML without libSBML.

This is a pragmatic fallback for environments where the builder scripts are
already updated but `libsbml` is not available to rebuild the merged file.
"""

from __future__ import annotations

import argparse
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path


SBML_NS = "http://www.sbml.org/sbml/level2/version4"
MATH_NS = "http://www.w3.org/1998/Math/MathML"
XHTML_NS = "http://www.w3.org/1999/xhtml"


def qname(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def fmt_number(value: float) -> str:
    return f"{value:.15g}"


def find_species(model: ET.Element, species_id: str) -> ET.Element | None:
    for sp in model.find(qname(SBML_NS, "listOfSpecies")) or []:
        if sp.get("id") == species_id:
            return sp
    return None


def find_reaction(model: ET.Element, reaction_id: str) -> ET.Element | None:
    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        return None
    for reaction in list_of_reactions:
        if reaction.get("id") == reaction_id:
            return reaction
    return None


def replace_ci_text(model: ET.Element, old: str, new: str) -> None:
    for ci in model.iter(qname(MATH_NS, "ci")):
        if ci.text is not None and ci.text.strip() == old:
            ci.text = new


def rename_species_refs(model: ET.Element, old: str, new: str) -> None:
    for elem in model.iter():
        if elem.get("species") == old:
            elem.set("species", new)
    replace_ci_text(model, old, new)


def remove_species(model: ET.Element, species_id: str) -> None:
    list_of_species = model.find(qname(SBML_NS, "listOfSpecies"))
    if list_of_species is None:
        return
    for sp in list(list_of_species):
        if sp.get("id") == species_id:
            list_of_species.remove(sp)
            break


def remove_reaction(model: ET.Element, reaction_id: str) -> None:
    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        return
    for reaction in list(list_of_reactions):
        if reaction.get("id") == reaction_id:
            list_of_reactions.remove(reaction)
            break


def ensure_species(
    model: ET.Element,
    species_id: str,
    name: str,
    compartment: str,
    initial_concentration: float,
    *,
    boundary_condition: bool = False,
    constant: bool = False,
    has_only_substance_units: bool = False,
) -> ET.Element:
    species = find_species(model, species_id)
    if species is None:
        list_of_species = model.find(qname(SBML_NS, "listOfSpecies"))
        if list_of_species is None:
            raise RuntimeError("No listOfSpecies found")
        species = ET.SubElement(
            list_of_species,
            qname(SBML_NS, "species"),
            {
                "id": species_id,
                "name": name,
                "compartment": compartment,
                "initialConcentration": fmt_number(initial_concentration),
                "boundaryCondition": str(boundary_condition).lower(),
                "constant": str(constant).lower(),
                "hasOnlySubstanceUnits": str(has_only_substance_units).lower(),
            },
        )
    else:
        species.set("name", name)
        species.set("compartment", compartment)
        species.attrib.pop("initialAmount", None)
        species.set("initialConcentration", fmt_number(initial_concentration))
        species.set("boundaryCondition", str(boundary_condition).lower())
        species.set("constant", str(constant).lower())
        species.set("hasOnlySubstanceUnits", str(has_only_substance_units).lower())
    return species


def update_unit_definitions(model: ET.Element) -> None:
    renames = {
        "asp_substance": "millimole",
        "asp_umole_per_litre": "mM",
        "asp_umole_per_litre_per_time": "mM_per_sec",
        "asp_umole2_per_litre2": "mM2",
        "asp_litre_per_umole_per_time": "litre_per_mmol_per_sec",
        "asp_umole_per_time": "mmole_per_sec",
    }
    list_of_unit_defs = model.find(qname(SBML_NS, "listOfUnitDefinitions"))
    if list_of_unit_defs is None:
        return

    for unit_def in list_of_unit_defs:
        unit_id = unit_def.get("id")
        if unit_id in renames:
            unit_def.set("name", renames[unit_id])
            list_of_units = unit_def.find(qname(SBML_NS, "listOfUnits"))
            if list_of_units is None:
                continue
            for unit in list_of_units:
                if unit.get("kind") == "mole":
                    unit.set("scale", "-3")

    if not any(unit_def.get("id") == "tca_permin" for unit_def in list_of_unit_defs):
        unit_def = ET.SubElement(
            list_of_unit_defs,
            qname(SBML_NS, "unitDefinition"),
            {"id": "tca_permin", "name": "per_min"},
        )
        list_of_units = ET.SubElement(unit_def, qname(SBML_NS, "listOfUnits"))
        ET.SubElement(
            list_of_units,
            qname(SBML_NS, "unit"),
            {"kind": "second", "exponent": "-1", "multiplier": "60"},
        )


def asp_units_are_already_rescaled(model: ET.Element) -> bool:
    list_of_unit_defs = model.find(qname(SBML_NS, "listOfUnitDefinitions"))
    if list_of_unit_defs is None:
        return False
    for unit_def in list_of_unit_defs:
        if unit_def.get("id") != "asp_umole_per_litre":
            continue
        for unit in unit_def.find(qname(SBML_NS, "listOfUnits")) or []:
            if unit.get("kind") == "mole" and unit.get("scale") == "-3":
                return True
    return False


def rescale_asp_quantities(model: ET.Element, already_rescaled: bool) -> None:
    if already_rescaled:
        return

    unit_scales = {
        "asp_umole_per_litre": 1.0 / 1000.0,
        "asp_umole_per_litre_per_time": 1.0 / 1000.0,
        "asp_umole2_per_litre2": 1.0 / 1000000.0,
        "asp_litre_per_umole_per_time": 1000.0,
        "asp_umole_per_time": 1.0 / 1000.0,
    }

    list_of_species = model.find(qname(SBML_NS, "listOfSpecies"))
    if list_of_species is not None:
        for sp in list_of_species:
            species_id = sp.get("id", "")
            if not species_id.startswith("asp_"):
                continue
            for attr in ("initialConcentration", "initialAmount"):
                if attr in sp.attrib:
                    sp.set(attr, fmt_number(float(sp.get(attr)) / 1000.0))

    for parameter in model.iter(qname(SBML_NS, "parameter")):
        units = parameter.get("units")
        value = parameter.get("value")
        if units not in unit_scales or value is None:
            continue
        parameter.set("value", fmt_number(float(value) * unit_scales[units]))


def patch_bridge_species(model: ET.Element) -> None:
    accoa = find_species(model, "pyr_AcCoA")
    if accoa is not None:
        accoa.attrib.pop("initialAmount", None)
        accoa.set("hasOnlySubstanceUnits", "false")
        # Millard2020 uses a substantial AcCoA baseline; seed the bridge so TCA
        # is not substrate-starved at t=0.
        accoa.set("initialConcentration", "1.46")

    oaa = find_species(model, "tca_oaa")
    if oaa is not None:
        oaa.set("boundaryCondition", "false")
        oaa.set("constant", "false")

    asp = find_species(model, "asp_Asp")
    if asp is not None:
        asp.set("boundaryCondition", "false")
        asp.set("constant", "false")


def fuse_coa_pool(model: ET.Element) -> None:
    rename_species_refs(model, "tca_coa", "pyr_CoA")
    remove_species(model, "tca_coa")


def replace_element(parent: ET.Element, old: ET.Element, new: ET.Element) -> None:
    children = list(parent)
    idx = children.index(old)
    parent.remove(old)
    parent.insert(idx, new)


def append_math(parent: ET.Element, xml_text: str) -> None:
    parent.append(ET.fromstring(xml_text))


def build_brenda_got_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "sboTerm": "SBO:0000176",
            "id": "GOT_AspAT",
            "name": "aspartate aminotransferase / GOT [OAA + Glu <-> Asp + 2-OG] (BRENDA-informed)",
            "reversible": "true",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Connector reaction absent from all source models. This replacement uses "
        "BRENDA-informed plant AspAT substrate affinities after rescaling the asp_ "
        "module onto an mM-style numeric scale. Km values were chosen from plant "
        "EC 2.6.1.1 records: oxaloacetate ~0.039 mM, glutamate ~8.2 mM, "
        "aspartate ~9.9 mM, 2-oxoglutarate ~0.25 mM. Forward/reverse Vmax values "
        "were tuned with the Arabidopsis AAT1 turnover-number ratio as a prior so "
        "the bridge is no longer arbitrary mass action."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "tca_oaa"})
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLU"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "asp_Asp"})
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "tca_akg"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <divide/>
    <apply>
      <minus/>
      <apply>
        <times/>
        <ci>Vf_AspAT</ci>
        <apply><divide/><ci>tca_oaa</ci><ci>Km_OAA_AspAT</ci></apply>
        <apply><divide/><ci>nh4_GLU</ci><ci>Km_GLU_AspAT</ci></apply>
      </apply>
      <apply>
        <times/>
        <ci>Vr_AspAT</ci>
        <apply><divide/><ci>asp_Asp</ci><ci>Km_Asp_AspAT</ci></apply>
        <apply><divide/><ci>tca_akg</ci><ci>Km_AKG_AspAT</ci></apply>
      </apply>
    </apply>
    <apply>
      <plus/>
      <cn type="integer">1</cn>
      <apply><divide/><ci>tca_oaa</ci><ci>Km_OAA_AspAT</ci></apply>
      <apply><divide/><ci>nh4_GLU</ci><ci>Km_GLU_AspAT</ci></apply>
      <apply><divide/><ci>asp_Asp</ci><ci>Km_Asp_AspAT</ci></apply>
      <apply><divide/><ci>tca_akg</ci><ci>Km_AKG_AspAT</ci></apply>
      <apply>
        <divide/>
        <apply><times/><ci>tca_oaa</ci><ci>nh4_GLU</ci></apply>
        <apply><times/><ci>Km_OAA_AspAT</ci><ci>Km_GLU_AspAT</ci></apply>
      </apply>
      <apply>
        <divide/>
        <apply><times/><ci>asp_Asp</ci><ci>tca_akg</ci></apply>
        <apply><times/><ci>Km_Asp_AspAT</ci><ci>Km_AKG_AspAT</ci></apply>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    list_of_parameters = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    params = [
        ("Vf_AspAT", "0.25", "tca_mmlmin"),
        ("Vr_AspAT", "0.39", "tca_mmlmin"),
        ("Km_OAA_AspAT", "0.039", "tca_mml"),
        ("Km_GLU_AspAT", "8.2", "tca_mml"),
        ("Km_Asp_AspAT", "9.9", "tca_mml"),
        ("Km_AKG_AspAT", "0.25", "tca_mml"),
    ]
    for pid, value, units in params:
        ET.SubElement(
            list_of_parameters,
            qname(SBML_NS, "parameter"),
            {"id": pid, "value": value, "units": units},
        )

    return reaction


def replace_got_reaction(model: ET.Element) -> None:
    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        return

    for reaction in list(list_of_reactions):
        if reaction.get("id") == "GOT_AspAT":
            list_of_reactions.remove(reaction)

    list_of_reactions.append(build_brenda_got_reaction())


def ensure_external_n_species(model: ET.Element) -> None:
    ext_nh4 = find_species(model, "ext_NH4")
    if ext_nh4 is None:
        list_of_species = model.find(qname(SBML_NS, "listOfSpecies"))
        if list_of_species is None:
            raise RuntimeError("No listOfSpecies found")
        ext_nh4 = ET.SubElement(
            list_of_species,
            qname(SBML_NS, "species"),
            {
                "id": "ext_NH4",
                "name": "external ammonium",
                "compartment": "nh4_compartment",
                "initialConcentration": "1",
                "boundaryCondition": "true",
                "constant": "true",
            },
        )
    else:
        ext_nh4.set("compartment", "nh4_compartment")
        ext_nh4.set("initialConcentration", ext_nh4.get("initialConcentration", "1"))
        ext_nh4.set("boundaryCondition", "true")
        ext_nh4.set("constant", "true")

    nh4 = find_species(model, "nh4_NH4")
    if nh4 is not None:
        nh4.set("boundaryCondition", "false")
        nh4.set("constant", "false")


def build_nh4_import_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "NH4_import",
            "name": "external ammonium import (Maeda-inspired reduced transport)",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Reduced active ammonium-import term inspired by the Maeda2019 "
        "ammonium transport and assimilation network. The original model "
        "tracks external, surface, and intracellular ammonium explicitly. "
        "Here that structure is compressed into a saturable import step so the "
        "merged scaffold can respond to low-N and high-N conditions without "
        "adding the full E. coli transport subnetwork."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "ext_NH4"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "nh4_NH4"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>nh4_Vmax_NH4_import</ci>
    <apply>
      <divide/>
      <ci>ext_NH4</ci>
      <apply>
        <plus/>
        <ci>nh4_Km_NH4_import</ci>
        <ci>ext_NH4</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    list_of_parameters = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        list_of_parameters,
        qname(SBML_NS, "parameter"),
        {"id": "nh4_Vmax_NH4_import", "value": "0.12", "units": "tca_permin"},
    )
    ET.SubElement(
        list_of_parameters,
        qname(SBML_NS, "parameter"),
        {"id": "nh4_Km_NH4_import", "value": "0.5", "units": "tca_mml"},
    )
    return reaction


def build_primary_assimilation_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "primary_N_assimilation",
            "name": "primary N assimilation support to glutamate",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Mild plant-facing assimilation support term added because the merged "
        "scaffold does not explicitly include nitrate reduction or the full "
        "photorespiratory ammonium-reassimilation context. The term is driven "
        "by external ammonium availability and modulated by 2-oxoglutarate so "
        "upstream N support remains coupled to the TCA carbon-skeleton pool."
    )

    modifiers = ET.SubElement(reaction, qname(SBML_NS, "listOfModifiers"))
    ET.SubElement(modifiers, qname(SBML_NS, "modifierSpeciesReference"), {"species": "ext_NH4"})
    ET.SubElement(modifiers, qname(SBML_NS, "modifierSpeciesReference"), {"species": "tca_akg"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLU"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>nh4_Vmax_primary_assim</ci>
    <apply>
      <divide/>
      <ci>ext_NH4</ci>
      <apply>
        <plus/>
        <ci>nh4_Km_extN_primary</ci>
        <ci>ext_NH4</ci>
      </apply>
    </apply>
    <apply>
      <divide/>
      <ci>tca_akg</ci>
      <apply>
        <plus/>
        <ci>nh4_Km_AKG_primary</ci>
        <ci>tca_akg</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    list_of_parameters = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        list_of_parameters,
        qname(SBML_NS, "parameter"),
        {"id": "nh4_Vmax_primary_assim", "value": "0.015", "units": "tca_permin"},
    )
    ET.SubElement(
        list_of_parameters,
        qname(SBML_NS, "parameter"),
        {"id": "nh4_Km_extN_primary", "value": "0.5", "units": "tca_mml"},
    )
    ET.SubElement(
        list_of_parameters,
        qname(SBML_NS, "parameter"),
        {"id": "nh4_Km_AKG_primary", "value": "0.05", "units": "tca_mml"},
    )
    return reaction


def ensure_round2_species(model: ET.Element) -> None:
    # Keep these in the nh4-facing compartment so the new manuscript-facing
    # branches couple to the Glu/Gln pools the user is already perturbing.
    ensure_species(
        model,
        "asn_c",
        "asparagine",
        "nh4_compartment",
        0.0,
    )
    ensure_species(
        model,
        "p5c_c",
        "1-pyrroline-5-carboxylate",
        "nh4_compartment",
        0.0,
    )
    ensure_species(
        model,
        "pro_c",
        "proline",
        "nh4_compartment",
        0.0,
    )


def build_asns_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "ASNS",
            "name": "asparagine synthetase [Asp + Gln -> Asn + Glu] (BRENDA-informed reduced form)",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Reduced-form asparagine synthetase branch added to represent the "
        "manuscript-facing Asn transport/storage node. The full plant reaction "
        "is ATP-dependent, but the merged model already mixes donor modules and "
        "energy bookkeeping, so this connector keeps only the Asp/Gln -> Asn/Glu "
        "nitrogen-transfer core. Km priors were chosen from maize ZmAsn1-4 "
        "records in BRENDA, with Asp around 0.91-0.98 mM and Gln around "
        "0.09-0.54 mM."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "asp_Asp"})
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLN"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "asn_c"})
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLU"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>Vmax_ASNS</ci>
    <apply>
      <divide/>
      <ci>asp_Asp</ci>
      <apply>
        <plus/>
        <ci>Km_Asp_ASNS</ci>
        <ci>asp_Asp</ci>
      </apply>
    </apply>
    <apply>
      <divide/>
      <ci>nh4_GLN</ci>
      <apply>
        <plus/>
        <ci>Km_Gln_ASNS</ci>
        <ci>nh4_GLN</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    for pid, value in [
        ("Vmax_ASNS", "1.0"),
        ("Km_Asp_ASNS", "0.95"),
        ("Km_Gln_ASNS", "0.25"),
    ]:
        ET.SubElement(
            params,
            qname(SBML_NS, "parameter"),
            {"id": pid, "value": value, "units": "tca_mml" if pid.startswith("Km_") else "tca_mmlmin"},
        )
    return reaction


def build_asnase_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "ASNASE",
            "name": "asparaginase [Asn -> Asp + NH4] (BRENDA-informed reduced form)",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Asparaginase branch added to represent asparagine remobilization back "
        "into Asp and ammonium. BRENDA contains broad kinetic coverage for "
        "L-asparaginase, including low-submillimolar substrate affinities. "
        "This reduced form now uses a concrete BRENDA-like Asn affinity prior "
        "instead of a purely scaffolded placeholder."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "asn_c"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "asp_Asp"})
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "nh4_NH4"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>Vmax_ASNASE</ci>
    <apply>
      <divide/>
      <ci>asn_c</ci>
      <apply>
        <plus/>
        <ci>Km_Asn_ASNASE</ci>
        <ci>asn_c</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Vmax_ASNASE", "value": "0.35", "units": "tca_mmlmin"},
    )
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Km_Asn_ASNASE", "value": "0.06", "units": "tca_mml"},
    )
    return reaction


def build_omega_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "OMEGA",
            "name": "collapsed glutamine transaminase + omega-amidase / NIT2 branch",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Collapsed glutamine-transaminase/omega-amidase recycling branch used "
        "to represent the manuscript-specific NIT2-like candidate without "
        "introducing a second underconstrained feeder step. BRENDA directly "
        "supports omega-amidase/NIT2 kinetics, including Arabidopsis and maize "
        "proteins, and the net branch converts glutamine-derived amide N back "
        "toward 2-oxoglutarate plus ammonium."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLN"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "tca_akg"})
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "nh4_NH4"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>Vmax_OMEGA</ci>
    <apply>
      <divide/>
      <ci>nh4_GLN</ci>
      <apply>
        <plus/>
        <ci>Km_Gln_OMEGA</ci>
        <ci>nh4_GLN</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Vmax_OMEGA", "value": "0.75", "units": "tca_mmlmin"},
    )
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Km_Gln_OMEGA", "value": "0.25", "units": "tca_mml"},
    )
    return reaction


def build_p5cs_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "P5CS",
            "name": "delta1-pyrroline-5-carboxylate synthase [Glu -> P5C] (BRENDA-informed reduced form)",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "Collapsed P5CS entry step for the proline axis. The true plant enzyme "
        "bundles glutamate 5-kinase plus glutamate-5-semialdehyde "
        "dehydrogenase chemistry, but the merged scaffold does not carry a "
        "fully consistent ATP/NADPH bookkeeping layer across modules. This "
        "reduced form therefore uses glutamate dependence only, with Km priors "
        "anchored on glutamate 5-kinase records."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "nh4_GLU"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "p5c_c"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>Vmax_P5CS</ci>
    <apply>
      <divide/>
      <ci>nh4_GLU</ci>
      <apply>
        <plus/>
        <ci>Km_Glu_P5CS</ci>
        <ci>nh4_GLU</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Vmax_P5CS", "value": "1.0", "units": "tca_mmlmin"},
    )
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Km_Glu_P5CS", "value": "10", "units": "tca_mml"},
    )
    return reaction


def build_p5cr_reaction() -> ET.Element:
    reaction = ET.Element(
        qname(SBML_NS, "reaction"),
        {
            "id": "P5CR",
            "name": "P5C reductase [P5C -> proline] (BRENDA-informed reduced form)",
            "reversible": "false",
            "fast": "false",
        },
    )

    notes = ET.SubElement(reaction, qname(SBML_NS, "notes"))
    body = ET.SubElement(notes, qname(XHTML_NS, "body"))
    p = ET.SubElement(body, qname(XHTML_NS, "p"))
    p.text = (
        "P5C reductase step for the proline axis. BRENDA provides strong plant "
        "P5CR support with low-submillimolar P5C affinity and NADPH dependence; "
        "the reduced form keeps the substrate dependence and omits explicit "
        "redox bookkeeping so the branch remains compatible with the stitched "
        "model."
    )

    reactants = ET.SubElement(reaction, qname(SBML_NS, "listOfReactants"))
    ET.SubElement(reactants, qname(SBML_NS, "speciesReference"), {"species": "p5c_c"})

    products = ET.SubElement(reaction, qname(SBML_NS, "listOfProducts"))
    ET.SubElement(products, qname(SBML_NS, "speciesReference"), {"species": "pro_c"})

    kinetic_law = ET.SubElement(reaction, qname(SBML_NS, "kineticLaw"))
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>nh4_compartment</ci>
    <ci>Vmax_P5CR</ci>
    <apply>
      <divide/>
      <ci>p5c_c</ci>
      <apply>
        <plus/>
        <ci>Km_P5C_P5CR</ci>
        <ci>p5c_c</ci>
      </apply>
    </apply>
  </apply>
</math>
""".strip(),
    )

    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Vmax_P5CR", "value": "0.6", "units": "tca_mmlmin"},
    )
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {"id": "Km_P5C_P5CR", "value": "0.18", "units": "tca_mml"},
    )
    return reaction


def add_round2_branches(model: ET.Element) -> None:
    ensure_round2_species(model)
    for reaction_id in ("ASNS", "ASNASE", "OMEGA", "P5CS", "P5CR"):
        remove_reaction(model, reaction_id)

    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        raise RuntimeError("No listOfReactions found")

    for reaction in (
        build_asns_reaction(),
        build_asnase_reaction(),
        build_omega_reaction(),
        build_p5cs_reaction(),
        build_p5cr_reaction(),
    ):
        list_of_reactions.append(reaction)


def add_external_n_support(model: ET.Element) -> None:
    ensure_external_n_species(model)
    remove_reaction(model, "NH4_import")
    remove_reaction(model, "GLU_source")
    remove_reaction(model, "primary_N_assimilation")

    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        raise RuntimeError("No listOfReactions found")
    list_of_reactions.append(build_nh4_import_reaction())
    list_of_reactions.append(build_primary_assimilation_reaction())


def build_tca_syn_kinetic_law() -> ET.Element:
    kinetic_law = ET.Element(qname(SBML_NS, "kineticLaw"), {"metaid": "tca__871723"})
    append_math(
        kinetic_law,
        f"""
<math xmlns="{MATH_NS}">
  <apply>
    <times/>
    <ci>tca_cell</ci>
    <ci>tca_kbiosyn</ci>
    <ci>tca_akg</ci>
  </apply>
</math>
""".strip(),
    )
    params = ET.SubElement(kinetic_law, qname(SBML_NS, "listOfParameters"))
    ET.SubElement(
        params,
        qname(SBML_NS, "parameter"),
        {
            "metaid": "tca__287462",
            "sboTerm": "SBO:0000320",
            "id": "tca_kbiosyn",
            "value": "0.005",
            "units": "tca_permin",
        },
    )
    return kinetic_law


def retune_tca_syn(model: ET.Element) -> None:
    list_of_reactions = model.find(qname(SBML_NS, "listOfReactions"))
    if list_of_reactions is None:
        return
    for reaction in list_of_reactions:
        if reaction.get("id") != "tca_SYN":
            continue
        reaction.set("name", "alpha-ketoglutarate biomass withdrawal (retuned)")
        notes = reaction.find(qname(SBML_NS, "notes"))
        if notes is not None:
            body = notes.find(qname(XHTML_NS, "body"))
            if body is not None:
                for child in list(body):
                    body.remove(child)
                p = ET.SubElement(body, qname(XHTML_NS, "p"))
                p.text = (
                    "Original Singh2006 SYN kinetics were copied from the ICD step and "
                    "pulled carbon out of alpha-ketoglutarate too aggressively for the "
                    "merged model. This retuned version keeps a small explicit biomass "
                    "withdrawal from tca_akg so the TCA module can run without the "
                    "borrowed ICD-like drain dominating the bridge."
                )
        modifiers = reaction.find(qname(SBML_NS, "listOfModifiers"))
        if modifiers is not None:
            reaction.remove(modifiers)
        kinetic_law = reaction.find(qname(SBML_NS, "kineticLaw"))
        if kinetic_law is not None:
            replace_element(reaction, kinetic_law, build_tca_syn_kinetic_law())
        else:
            reaction.append(build_tca_syn_kinetic_law())
        break


def set_parameter_value(model: ET.Element, param_id: str, value: float) -> None:
    for parameter in model.iter(qname(SBML_NS, "parameter")):
        if parameter.get("id") == param_id:
            parameter.set("value", fmt_number(value))
            return


def retune_nitrogen_pull(model: ET.Element) -> None:
    """Dampen the bacterial NH4 module so it behaves like a coarse plant connector.

    The Bruggeman module is useful structurally, but its original GS/GDH/GOGAT
    capacities are too aggressive for this stitched model and collapse the shared
    akg/GLU/GLN pools. Here we make it a gentler connector:
    - GS remains the dominant ammonium-assimilation route, but slower.
    - GDH is strongly downweighted because it is not the main assimilatory route
      in the plant-like use case.
    - GOGAT gets a less extreme affinity for akg, so it stops vacuuming the TCA
      pool at very low concentrations.
    """
    updates = {
        "nh4_Vgs": 30.0,
        "nh4_Vgdh": 10.0,
        "nh4_Vgog": 8.0,
        "nh4_Kgogkg": 0.2,
    }
    for param_id, value in updates.items():
        set_parameter_value(model, param_id, value)


def retune_n_demand_sinks(model: ET.Element) -> None:
    """Reduce abstract GLU/GLN demand sinks after integration.

    In the source Bruggeman model, vgludem/vglndem are lumped reduced-N demand
    terms. After merging in explicit downstream Asp-family and connected
    metabolism, leaving those sinks at their original values likely
    double-counts demand. We apply a conservative 0.75x reduction and treat it
    as the new default baseline.
    """
    updates = {
        "nh4_Vgludem": 90.0,
        "nh4_Vglndem": 52.5,
    }
    for param_id, value in updates.items():
        set_parameter_value(model, param_id, value)


def retune_thr_branch(model: ET.Element) -> None:
    """Use Chassagnole-style Thr-branch priors to relieve the TS1 bottleneck.

    PHser accumulates strongly while Thr barely rises, which points to the
    PHser -> Thr step as the main choke. The Chassagnole donor model has a
    simple threonine synthase step with Km(hsp) ~ 0.31 and no extreme phosphate
    inhibition. We borrow that directionally:
    - shrink the effective TS1 PHser pseudo-Km from 250 to 0.31
    - greatly weaken phosphate inhibition on TS1
    - modestly boost TS1 capacity
    """
    updates = {
        "asp_TS1_AdoMEt_Km_no_AdoMet_exp": 0.31,
        "asp_TS1_Phosphate_Ki_exp": 100.0,
        "asp_TS1_kcatmin_exp": 1.0,
        "asp_TS1_AdoMet_kcatmax_exp": 6.0,
    }
    for param_id, value in updates.items():
        set_parameter_value(model, param_id, value)


def patch_file(path: Path) -> Path:
    backup = path.with_suffix(path.suffix + ".pre_brenda_units.bak")
    if not backup.exists():
        shutil.copy2(path, backup)

    tree = ET.parse(path)
    root = tree.getroot()
    model = root.find(qname(SBML_NS, "model"))
    if model is None:
        raise RuntimeError("No SBML model element found")
    model.attrib.pop("volumeUnits", None)

    already_rescaled = asp_units_are_already_rescaled(model)
    rescale_asp_quantities(model, already_rescaled)
    update_unit_definitions(model)
    patch_bridge_species(model)
    fuse_coa_pool(model)
    retune_tca_syn(model)
    retune_nitrogen_pull(model)
    add_external_n_support(model)
    retune_n_demand_sinks(model)
    retune_thr_branch(model)
    replace_got_reaction(model)
    add_round2_branches(model)

    ET.indent(tree, space="  ")
    tree.write(path, encoding="UTF-8", xml_declaration=True)
    return backup


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("xml_path", help="Path to merged_nitrogen_connected.xml")
    args = parser.parse_args()

    for prefix, uri in [
        ("", SBML_NS),
        ("math", MATH_NS),
        ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#"),
        ("dc", "http://purl.org/dc/elements/1.1/"),
        ("dcterms", "http://purl.org/dc/terms/"),
        ("vCard", "http://www.w3.org/2001/vcard-rdf/3.0#"),
        ("bqbiol", "http://biomodels.net/biology-qualifiers/"),
        ("bqmodel", "http://biomodels.net/model-qualifiers/"),
    ]:
        ET.register_namespace(prefix, uri)

    xml_path = Path(args.xml_path).expanduser().resolve()
    backup = patch_file(xml_path)
    print(f"Patched {xml_path}")
    print(f"Backup saved to {backup}")


if __name__ == "__main__":
    main()

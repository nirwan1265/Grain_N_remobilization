#!/usr/bin/env python3
"""Patch the Grain_N builder scripts with conservative pyruvate/OAA bridges."""

from pathlib import Path
import sys


FUNCTION_BLOCK = '''
def add_pyruvate_oaa_bridges(model):
    """Add conservative pyruvate/OAA bridge reactions with BRENDA-backed PC affinities."""
    for rid in ["PYR_OAA_support", "PYR_PC_bridge", "OAA_PEPCK_bridge"]:
        rxn = model.getReaction(rid)
        if rxn is not None:
            model.removeReaction(rid)

    for sid in ["bridge_CO2", "bridge_Pi", "bridge_PEP"]:
        sp = model.getSpecies(sid)
        if sp is not None:
            model.removeSpecies(sid)

    tca_comp = model.getCompartment("tca_cell")
    if tca_comp is None:
        raise RuntimeError("missing tca_cell compartment for pyruvate/OAA bridge")
    tca_comp_id = tca_comp.getId()

    def add_boundary_species(sid, name, value):
        sp = model.createSpecies()
        sp.setId(sid)
        sp.setName(name)
        sp.setCompartment(tca_comp_id)
        sp.setInitialConcentration(value)
        sp.setBoundaryCondition(True)
        sp.setConstant(True)
        sp.setHasOnlySubstanceUnits(False)
        return sp

    add_boundary_species("bridge_CO2", "bridge CO2 support pool", 1.0)
    add_boundary_species("bridge_Pi", "bridge phosphate support pool", 1.0)

    pep = model.createSpecies()
    pep.setId("bridge_PEP")
    pep.setName("bridge PEP scaffold pool")
    pep.setCompartment(tca_comp_id)
    pep.setInitialConcentration(0.01)
    pep.setBoundaryCondition(False)
    pep.setConstant(False)
    pep.setHasOnlySubstanceUnits(False)

    r = model.createReaction()
    r.setId("PYR_OAA_support")
    r.setName("pyruvate to oxaloacetate support bridge")
    r.setReversible(False)
    r.setFast(False)
    rr = r.createReactant(); rr.setSpecies("pyr_pyruvate")
    pr = r.createProduct(); pr.setSpecies("tca_oaa")
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseFormula(
        "tca_cell * Vmax_PYR_OAA_support * pyr_pyruvate / (Km_PYR_OAA_support + pyr_pyruvate)"
    ))
    p = kl.createParameter(); p.setId("Vmax_PYR_OAA_support"); p.setValue(0.008); p.setUnits("tca_mmlmin")
    p = kl.createParameter(); p.setId("Km_PYR_OAA_support"); p.setValue(0.25); p.setUnits("tca_mml")

    r = model.createReaction()
    r.setId("PYR_PC_bridge")
    r.setName("pyruvate carboxylase-like bridge [pyruvate + CO2 + ATP -> OAA + ADP + Pi]")
    r.setReversible(False)
    r.setFast(False)
    r.setSBOTerm(176)
    rr = r.createReactant(); rr.setSpecies("bridge_CO2")
    rr = r.createReactant(); rr.setSpecies("pyr_pyruvate")
    rr = r.createReactant(); rr.setSpecies("pyr_ATP")
    pr = r.createProduct(); pr.setSpecies("bridge_Pi")
    pr = r.createProduct(); pr.setSpecies("tca_oaa")
    pr = r.createProduct(); pr.setSpecies("pyr_ADP")
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseFormula(
        "tca_cell * Vmax_PC_bridge * "
        "(bridge_CO2 / (Km_CO2_PC_bridge + bridge_CO2)) * "
        "(pyr_pyruvate / (Km_PYR_PC_bridge + pyr_pyruvate)) * "
        "(pyr_ATP / (Km_ATP_PC_bridge + pyr_ATP))"
    ))
    p = kl.createParameter(); p.setId("Vmax_PC_bridge"); p.setValue(0.018); p.setUnits("tca_mmlmin")
    p = kl.createParameter(); p.setId("Km_CO2_PC_bridge"); p.setValue(1.36); p.setUnits("tca_mml")
    p = kl.createParameter(); p.setId("Km_PYR_PC_bridge"); p.setValue(0.50); p.setUnits("tca_mml")
    p = kl.createParameter(); p.setId("Km_ATP_PC_bridge"); p.setValue(0.07); p.setUnits("tca_mml")

    r = model.createReaction()
    r.setId("OAA_PEPCK_bridge")
    r.setName("PEPCK-like scaffold [OAA <-> PEP]")
    r.setReversible(True)
    r.setFast(False)
    r.setSBOTerm(176)
    rr = r.createReactant(); rr.setSpecies("tca_oaa")
    pr = r.createProduct(); pr.setSpecies("bridge_PEP")
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseFormula(
        "tca_cell * (((Vf_PEPCK_bridge * (tca_oaa / Km_OAA_PEPCK_bridge)) - "
        "(Vr_PEPCK_bridge * (bridge_PEP / Km_PEP_PEPCK_bridge))) / "
        "(1 + tca_oaa / Km_OAA_PEPCK_bridge + bridge_PEP / Km_PEP_PEPCK_bridge)))"
    ))
    p = kl.createParameter(); p.setId("Vf_PEPCK_bridge"); p.setValue(0.004); p.setUnits("tca_mmlmin")
    p = kl.createParameter(); p.setId("Vr_PEPCK_bridge"); p.setValue(0.002); p.setUnits("tca_mmlmin")
    p = kl.createParameter(); p.setId("Km_OAA_PEPCK_bridge"); p.setValue(0.01); p.setUnits("tca_mml")
    p = kl.createParameter(); p.setId("Km_PEP_PEPCK_bridge"); p.setValue(0.05); p.setUnits("tca_mml")

    print("  added pyruvate/OAA bridge scaffold (support + BRENDA-backed PC-like + PEPCK-like)")
'''


def replace_once(text: str, old: str, new: str, path: Path) -> str:
    if old not in text:
        raise RuntimeError(f"expected block not found in {path}")
    return text.replace(old, new, 1)


def patch_builder(path: Path) -> None:
    text = path.read_text()

    if "def add_pyruvate_oaa_bridges(model):" not in text:
        text = replace_once(
            text,
            'def add_brenda_got_aspat(model):\n',
            FUNCTION_BLOCK + '\n\ndef add_brenda_got_aspat(model):\n',
            path,
        )

    if "add_pyruvate_oaa_bridges(M)" not in text:
        text = replace_once(
            text,
            "seed_bridge_state(M)\nretune_tca_syn(M)\n",
            "seed_bridge_state(M)\nadd_pyruvate_oaa_bridges(M)\nretune_tca_syn(M)\n",
            path,
        )

    path.write_text(text)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: add_pyruvate_oaa_bridges.py <builder.py> [<builder.py> ...]", file=sys.stderr)
        return 1
    for arg in argv[1:]:
        patch_builder(Path(arg))
        print(f"patched {arg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

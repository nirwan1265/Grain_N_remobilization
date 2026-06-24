#!/usr/bin/env python3
"""Patch the external Grain_N builder scripts with TCA bridge improvements."""

from __future__ import annotations

import argparse
from pathlib import Path


OLD_SEED = """def seed_bridge_state(model):
    \"\"\"Keep bridge species from starting substrate-starved.\"\"\"
    accoa = model.getSpecies(\"pyr_AcCoA\")
    if accoa is not None:
        accoa.setHasOnlySubstanceUnits(False)
        if hasattr(accoa, \"unsetInitialAmount\"):
            accoa.unsetInitialAmount()
        accoa.setInitialConcentration(0.5)
        print(\"  seeded pyr_AcCoA to 0.5 mM from the TCA source module\")

    oaa = model.getSpecies(\"tca_oaa\")
    if oaa is not None:
        oaa.setBoundaryCondition(False)
        oaa.setConstant(False)
        print(\"  un-clamped tca_oaa so the TCA cycle / AspAT bridge can move it\")

"""


NEW_SEED = """def seed_bridge_state(model):
    \"\"\"Keep bridge species from starting substrate-starved.\"\"\"
    accoa = model.getSpecies(\"pyr_AcCoA\")
    if accoa is not None:
        accoa.setHasOnlySubstanceUnits(False)
        if hasattr(accoa, \"unsetInitialAmount\"):
            accoa.unsetInitialAmount()
        accoa.setInitialConcentration(1.46)
        print(\"  seeded pyr_AcCoA to 1.46 mM using the Millard2020 carbon-entry baseline\")

    oaa = model.getSpecies(\"tca_oaa\")
    if oaa is not None:
        oaa.setBoundaryCondition(False)
        oaa.setConstant(False)
        print(\"  un-clamped tca_oaa so the TCA cycle / AspAT bridge can move it\")

"""


NEW_TCA_SYN = """
def ensure_tca_permin_unit(model):
    if model.getUnitDefinition("tca_permin") is not None:
        return
    ud = model.createUnitDefinition()
    ud.setId("tca_permin")
    ud.setName("per_min")
    u = ud.createUnit()
    u.setKind(libsbml.UNIT_KIND_SECOND)
    u.setExponent(-1)
    u.setMultiplier(60)
    print("  added tca_permin unit definition")

"""


NEW_TCA_SYN += """
def retune_tca_syn(model):
    \"\"\"Replace the copied ICD-like biomass drain with a mild explicit akg sink.\"\"\"
    rxn = model.getReaction(\"tca_SYN\")
    if rxn is None:
        return
    rxn.setName(\"alpha-ketoglutarate biomass withdrawal (retuned)\")
    rxn.unsetNotes()
    rxn.setNotes(
        '<body xmlns=\"http://www.w3.org/1999/xhtml\"><p>'
        'Original Singh2006 SYN kinetics were copied from the ICD step and pulled '
        'carbon out of alpha-ketoglutarate too aggressively for the merged model. '
        'This retuned version keeps a small explicit biomass withdrawal from '
        'tca_akg so the TCA module can run without the borrowed ICD-like drain '
        'dominating the bridge.</p></body>'
    )
    if rxn.getNumModifiers() > 0:
        for i in range(rxn.getNumModifiers() - 1, -1, -1):
            rxn.removeModifier(i)
    kl = rxn.getKineticLaw()
    if kl is None:
        kl = rxn.createKineticLaw()
    for i in range(kl.getNumParameters() - 1, -1, -1):
        kl.removeParameter(i)
    ensure_tca_permin_unit(model)
    kl.setMath(libsbml.parseFormula(\"tca_cell * tca_kbiosyn * tca_akg\"))
    p = kl.createParameter()
    p.setId(\"tca_kbiosyn\")
    p.setValue(0.005)
    p.setUnits(\"tca_permin\")
    print(\"  retuned tca_SYN to a mild explicit alpha-ketoglutarate sink\")

"""

NEW_TCA_SYN += """
def retune_nitrogen_pull(model):
    \"\"\"Dampen the bacterial NH4 module so it behaves like a coarse plant connector.\"\"\"
    updates = {
        \"nh4_Vgs\": 30.0,
        \"nh4_Vgdh\": 10.0,
        \"nh4_Vgog\": 8.0,
        \"nh4_Kgogkg\": 0.2,
    }
    for pid, value in updates.items():
        p = model.getParameter(pid)
        if p is not None:
            p.setValue(value)
    print(\"  retuned NH4 GS/GDH/GOGAT pull to a gentler integrated-model setting\")

"""

NITROGEN_HELPER = """
def retune_nitrogen_pull(model):
    \"\"\"Dampen the bacterial NH4 module so it behaves like a coarse plant connector.\"\"\"
    updates = {
        \"nh4_Vgs\": 30.0,
        \"nh4_Vgdh\": 10.0,
        \"nh4_Vgog\": 8.0,
        \"nh4_Kgogkg\": 0.2,
    }
    for pid, value in updates.items():
        p = model.getParameter(pid)
        if p is not None:
            p.setValue(value)
    print(\"  retuned NH4 GS/GDH/GOGAT pull to a gentler integrated-model setting\")

"""

THR_HELPER = """
def retune_thr_branch(model):
    \"\"\"Use Chassagnole-style Thr-branch priors to relieve the TS1 bottleneck.\"\"\"
    updates = {
        \"asp_TS1_AdoMEt_Km_no_AdoMet_exp\": 0.31,
        \"asp_TS1_Phosphate_Ki_exp\": 100.0,
        \"asp_TS1_kcatmin_exp\": 1.0,
        \"asp_TS1_AdoMet_kcatmax_exp\": 6.0,
    }
    for pid, value in updates.items():
        p = model.getParameter(pid)
        if p is not None:
            p.setValue(value)
    print(\"  retuned TS1 using Chassagnole-style Thr-branch priors\")

"""

EXT_N_HELPER = """
def add_external_n_support(model):
    \"\"\"Add a tunable external N feed and mild GLU support to the NH4 module.\"\"\"
    ext = model.getSpecies(\"ext_NH4\")
    if ext is None:
        ext = model.createSpecies()
        ext.setId(\"ext_NH4\")
        ext.setName(\"external ammonium\")
        ext.setCompartment(\"nh4_compartment\")
        ext.setInitialConcentration(1.0)
    ext.setBoundaryCondition(True)
    ext.setConstant(True)

    nh4 = model.getSpecies(\"nh4_NH4\")
    if nh4 is not None:
        nh4.setBoundaryCondition(False)
        nh4.setConstant(False)

    for rid in [\"NH4_import\", \"GLU_source\"]:
        rxn = model.getReaction(rid)
        if rxn is not None:
            model.removeReaction(rid)

    r = model.createReaction()
    r.setId(\"NH4_import\")
    r.setName(\"external ammonium import\")
    r.setReversible(False)
    r.setFast(False)
    rr = r.createReactant(); rr.setSpecies(\"ext_NH4\")
    pr = r.createProduct(); pr.setSpecies(\"nh4_NH4\")
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseFormula(\"nh4_compartment * nh4_kNH4_import * ext_NH4\"))
    p = kl.createParameter(); p.setId(\"nh4_kNH4_import\"); p.setValue(0.08); p.setUnits(\"tca_permin\")

    r = model.createReaction()
    r.setId(\"GLU_source\")
    r.setName(\"primary assimilation support to glutamate\")
    r.setReversible(False)
    r.setFast(False)
    md = r.createModifier()
    md.setSpecies(\"ext_NH4\")
    pr = r.createProduct(); pr.setSpecies(\"nh4_GLU\")
    kl = r.createKineticLaw()
    kl.setMath(libsbml.parseFormula(\"nh4_compartment * nh4_kGLU_source * ext_NH4\"))
    p = kl.createParameter(); p.setId(\"nh4_kGLU_source\"); p.setValue(0.01); p.setUnits(\"tca_permin\")
    print(\"  added external NH4 feed and mild GLU support term\")

"""


def patch_file(path: Path, label: str) -> None:
    text = path.read_text()

    text = text.replace('M.setVolumeUnits("litre")\n', "")

    if OLD_SEED in text:
        text = text.replace(OLD_SEED, NEW_SEED)

    marker = "def add_brenda_got_aspat(model):"
    if marker in text and "def retune_tca_syn(model):" not in text:
        text = text.replace(marker, NEW_TCA_SYN + "\n" + marker)
    elif marker in text and "def retune_nitrogen_pull(model):" not in text:
        text = text.replace(marker, NITROGEN_HELPER + "\n" + marker)
    elif marker in text and "def retune_thr_branch(model):" not in text:
        text = text.replace(marker, THR_HELPER + "\n" + marker)
    elif marker in text and "def add_external_n_support(model):" not in text:
        text = text.replace(marker, EXT_N_HELPER + "\n" + marker)

    old_connect = """fuse(M,\"tca_aca\",\"pyr_AcCoA\")   # acetyl-CoA: pyruvate metabolism -> TCA
fuse(M,\"nh4_KG\",\"tca_akg\")      # 2-oxoglutarate: TCA -> N assimilation
seed_bridge_state(M)
"""
    new_connect = """fuse(M,\"tca_aca\",\"pyr_AcCoA\")   # acetyl-CoA: pyruvate metabolism -> TCA
fuse(M,\"tca_coa\",\"pyr_CoA\")     # unify the CoA pool so pyruvate can replenish shared AcCoA
fuse(M,\"nh4_KG\",\"tca_akg\")      # 2-oxoglutarate: TCA -> N assimilation
seed_bridge_state(M)
retune_tca_syn(M)
retune_nitrogen_pull(M)
add_external_n_support(M)
retune_thr_branch(M)
"""
    if old_connect in text:
        text = text.replace(old_connect, new_connect)
    elif "retune_tca_syn(M)\n" in text and "retune_nitrogen_pull(M)\n" not in text:
        text = text.replace("retune_tca_syn(M)\n", "retune_tca_syn(M)\nretune_nitrogen_pull(M)\n")
    if "retune_nitrogen_pull(M)\n" in text and "add_external_n_support(M)\n" not in text:
        text = text.replace("retune_nitrogen_pull(M)\n", "retune_nitrogen_pull(M)\nadd_external_n_support(M)\n")
    if "retune_nitrogen_pull(M)\n" in text and "retune_thr_branch(M)\n" not in text:
        text = text.replace("retune_nitrogen_pull(M)\n", "retune_nitrogen_pull(M)\nretune_thr_branch(M)\n")

    old_connect_main = """fuse(M, \"tca_aca\", \"pyr_AcCoA\")    # acetyl-CoA bridge (keep dynamic pyruvate-side species)
fuse(M, \"nh4_KG\",  \"tca_akg\")      # 2-oxoglutarate bridge (keep dynamic TCA-side species)
seed_bridge_state(M)
"""
    new_connect_main = """fuse(M, \"tca_aca\", \"pyr_AcCoA\")    # acetyl-CoA bridge (keep dynamic pyruvate-side species)
fuse(M, \"tca_coa\", \"pyr_CoA\")      # unify the CoA pool so pyruvate can replenish shared AcCoA
fuse(M, \"nh4_KG\",  \"tca_akg\")      # 2-oxoglutarate bridge (keep dynamic TCA-side species)
seed_bridge_state(M)
retune_tca_syn(M)
retune_nitrogen_pull(M)
add_external_n_support(M)
retune_thr_branch(M)
"""
    if old_connect_main in text:
        text = text.replace(old_connect_main, new_connect_main)
    elif "retune_tca_syn(M)\n" in text and "retune_nitrogen_pull(M)\n" not in text:
        text = text.replace("retune_tca_syn(M)\n", "retune_tca_syn(M)\nretune_nitrogen_pull(M)\n")
    if "retune_nitrogen_pull(M)\n" in text and "add_external_n_support(M)\n" not in text:
        text = text.replace("retune_nitrogen_pull(M)\n", "retune_nitrogen_pull(M)\nadd_external_n_support(M)\n")
    if "retune_nitrogen_pull(M)\n" in text and "retune_thr_branch(M)\n" not in text:
        text = text.replace("retune_nitrogen_pull(M)\n", "retune_nitrogen_pull(M)\nretune_thr_branch(M)\n")

    path.write_text(text)
    print(f"Patched {label}: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("grain_n_dir", help="Folder containing build_connected*.py")
    args = parser.parse_args()

    grain_n_dir = Path(args.grain_n_dir).expanduser().resolve()
    patch_file(grain_n_dir / "build_connected_local.py", "local builder")
    patch_file(grain_n_dir / "build_connected.py", "main builder")


if __name__ == "__main__":
    main()

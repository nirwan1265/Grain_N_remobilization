#!/usr/bin/env python3
"""Patch the CellDesigner merge builders with unit harmonization + BRENDA AspAT.

This script targets the external Grain_N model-building folder and updates:
- build_connected_local.py
- build_connected.py

Then it can optionally rebuild merged_nitrogen_connected.xml.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


HELPERS = """
def rescale_asp_module_to_mM(model):
    \"\"\"Rescale Curien Asp module from uM-style numbers to mM-style numbers.\"\"\"

    def scale_quantity(obj):
        if obj is None or not obj.isSetUnits() or not obj.isSetValue():
            return
        unit_id = obj.getUnits()
        factor = None
        if unit_id == "asp_umole_per_litre":
            factor = 1.0 / 1000.0
        elif unit_id == "asp_umole_per_litre_per_time":
            factor = 1.0 / 1000.0
        elif unit_id == "asp_umole2_per_litre2":
            factor = 1.0 / 1000000.0
        elif unit_id == "asp_litre_per_umole_per_time":
            factor = 1000.0
        elif unit_id == "asp_umole_per_time":
            factor = 1.0 / 1000.0
        if factor is not None:
            obj.setValue(obj.getValue() * factor)

    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        if not sp.getId().startswith("asp_"):
            continue
        if sp.isSetInitialConcentration():
            sp.setInitialConcentration(sp.getInitialConcentration() / 1000.0)
        if sp.isSetInitialAmount():
            sp.setInitialAmount(sp.getInitialAmount() / 1000.0)

    for i in range(model.getNumParameters()):
        scale_quantity(model.getParameter(i))

    for i in range(model.getNumReactions()):
        kl = model.getReaction(i).getKineticLaw()
        if kl is None:
            continue
        for j in range(kl.getNumParameters()):
            scale_quantity(kl.getParameter(j))

    unit_names = {
        "asp_substance": "millimole",
        "asp_umole_per_litre": "mM",
        "asp_umole_per_litre_per_time": "mM_per_sec",
        "asp_umole2_per_litre2": "mM2",
        "asp_litre_per_umole_per_time": "litre_per_mmol_per_sec",
        "asp_umole_per_time": "mmole_per_sec",
    }
    for unit_id, name in unit_names.items():
        ud = model.getUnitDefinition(unit_id)
        if ud is None:
            continue
        if hasattr(ud, "setName"):
            ud.setName(name)
        for k in range(ud.getNumUnits()):
            u = ud.getUnit(k)
            if u.getKind() == libsbml.UNIT_KIND_MOLE:
                u.setScale(-3)

    print("  rescaled asp_ module from uM-style values to mM-style values")


def seed_bridge_state(model):
    \"\"\"Keep bridge species from starting substrate-starved.\"\"\"
    accoa = model.getSpecies("pyr_AcCoA")
    if accoa is not None:
        accoa.setHasOnlySubstanceUnits(False)
        if hasattr(accoa, "unsetInitialAmount"):
            accoa.unsetInitialAmount()
        accoa.setInitialConcentration(0.5)
        print("  seeded pyr_AcCoA to 0.5 mM from the TCA source module")

    oaa = model.getSpecies("tca_oaa")
    if oaa is not None:
        oaa.setBoundaryCondition(False)
        oaa.setConstant(False)
        print("  un-clamped tca_oaa so the TCA cycle / AspAT bridge can move it")


def add_brenda_got_aspat(model):
    \"\"\"Add a BRENDA-informed reversible AspAT connector in mM-scale units.\"\"\"
    r = model.createReaction()
    r.setId("GOT_AspAT")
    r.setName("aspartate aminotransferase / GOT [OAA + Glu <-> Asp + 2-OG] (BRENDA-informed)")
    r.setReversible(True)
    r.setFast(False)
    r.setSBOTerm(176)
    r.setNotes(
        '<body xmlns="http://www.w3.org/1999/xhtml"><p>'
        'Connector reaction absent from all source models. This replacement uses '
        'BRENDA-informed plant AspAT substrate affinities after rescaling the asp_ '
        'module onto an mM-style numeric scale. Km values were chosen from plant '
        'EC 2.6.1.1 records: oxaloacetate ~0.039 mM, glutamate ~8.2 mM, '
        'aspartate ~9.9 mM, 2-oxoglutarate ~0.25 mM. Forward/reverse Vmax values '
        'were tuned with the Arabidopsis AAT1 turnover-number ratio as a prior so '
        'the bridge is no longer arbitrary mass action.</p></body>'
    )
    for sp in ["tca_oaa", "nh4_GLU"]:
        sr = r.createReactant()
        sr.setSpecies(sp)
    for sp in ["asp_Asp", "tca_akg"]:
        sr = r.createProduct()
        sr.setSpecies(sp)

    kl = r.createKineticLaw()
    formula = (
        "((Vf_AspAT * (tca_oaa / Km_OAA_AspAT) * (nh4_GLU / Km_GLU_AspAT)) - "
        "(Vr_AspAT * (asp_Asp / Km_Asp_AspAT) * (tca_akg / Km_AKG_AspAT))) / "
        "(1 + tca_oaa / Km_OAA_AspAT + nh4_GLU / Km_GLU_AspAT + "
        "asp_Asp / Km_Asp_AspAT + tca_akg / Km_AKG_AspAT + "
        "(tca_oaa * nh4_GLU) / (Km_OAA_AspAT * Km_GLU_AspAT) + "
        "(asp_Asp * tca_akg) / (Km_Asp_AspAT * Km_AKG_AspAT))"
    )
    kl.setMath(libsbml.parseFormula(formula))

    params = [
        ("Vf_AspAT", 0.25, "tca_mmlmin"),
        ("Vr_AspAT", 0.39, "tca_mmlmin"),
        ("Km_OAA_AspAT", 0.039, "tca_mml"),
        ("Km_GLU_AspAT", 8.2, "tca_mml"),
        ("Km_Asp_AspAT", 9.9, "tca_mml"),
        ("Km_AKG_AspAT", 0.25, "tca_mml"),
    ]
    for pid, value, units in params:
        p = kl.createParameter()
        p.setId(pid)
        p.setValue(value)
        p.setUnits(units)

    print("  added BRENDA-informed GOT_AspAT connector")
"""


LOCAL_REPLACEMENT = """
print("connecting:")
rescale_asp_module_to_mM(M)
fuse(M,"tca_aca","pyr_AcCoA")   # acetyl-CoA: pyruvate metabolism -> TCA
fuse(M,"nh4_KG","tca_akg")      # 2-oxoglutarate: TCA -> N assimilation
seed_bridge_state(M)
asp=M.getSpecies("asp_Asp"); asp.setBoundaryCondition(False); asp.setConstant(False)
add_brenda_got_aspat(M)
"""


MAIN_REPLACEMENT = """
print("connecting modules:")
rescale_asp_module_to_mM(M)
fuse(M, "tca_aca", "pyr_AcCoA")    # acetyl-CoA bridge (keep dynamic pyruvate-side species)
fuse(M, "nh4_KG",  "tca_akg")      # 2-oxoglutarate bridge (keep dynamic TCA-side species)
seed_bridge_state(M)

asp = M.getSpecies("asp_Asp")
asp.setBoundaryCondition(False); asp.setConstant(False)
print("  un-clamped asp_Asp (now produced by GOT)")
add_brenda_got_aspat(M)
"""


def patch_file(path: Path, connect_start: str, connect_end: str, replacement: str) -> None:
    text = path.read_text()

    fuse_block = """def fuse(model, remove_id, keep_id):
    els=model.getListOfAllElements()
    for i in range(els.getSize()): els.get(i).renameSIdRefs(remove_id, keep_id)
    model.removeSpecies(remove_id); print("  fused %-12s -> %s"%(remove_id, keep_id))
"""
    if fuse_block in text and "def rescale_asp_module_to_mM" not in text:
        text = text.replace(fuse_block, fuse_block + "\n" + HELPERS.strip() + "\n")

    fuse_block_alt = """def fuse(model, remove_id, keep_id):
    els=model.getListOfAllElements()
    for i in range(els.getSize()): els.get(i).renameSIdRefs(remove_id, keep_id)
    model.removeSpecies(remove_id)
    print("  fused %-12s -> %s"%(remove_id, keep_id))
"""
    if fuse_block_alt in text and "def rescale_asp_module_to_mM" not in text:
        text = text.replace(fuse_block_alt, fuse_block_alt + "\n" + HELPERS.strip() + "\n")

    start_idx = text.index(connect_start)
    end_idx = text.index(connect_end)
    text = text[:start_idx] + replacement.strip() + "\n\n" + text[end_idx:]

    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("grain_n_dir", help="Folder containing build_connected*.py")
    parser.add_argument("--rebuild", action="store_true", help="Rebuild merged_nitrogen_connected.xml")
    args = parser.parse_args()

    grain_n_dir = Path(args.grain_n_dir).expanduser().resolve()
    local_path = grain_n_dir / "build_connected_local.py"
    main_path = grain_n_dir / "build_connected.py"

    patch_file(
        local_path,
        'print("connecting:")',
        'out=os.path.join(HERE,"merged_nitrogen_connected.xml")',
        LOCAL_REPLACEMENT,
    )
    patch_file(
        main_path,
        'print("connecting modules:")',
        '# ---- write + validate ----',
        MAIN_REPLACEMENT,
    )

    if args.rebuild:
        subprocess.run(
            ["python", str(local_path)],
            cwd=str(grain_n_dir),
            check=True,
        )

    print(f"Patched {local_path}")
    print(f"Patched {main_path}")
    if args.rebuild:
        print(f"Rebuilt {grain_n_dir / 'merged_nitrogen_connected.xml'}")


if __name__ == "__main__":
    main()

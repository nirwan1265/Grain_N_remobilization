#!/usr/bin/env python3
"""Patch merged_nitrogen_connected.xml with conservative pyruvate/OAA bridges."""

from pathlib import Path
import sys


SPECIES_SNIPPET = """
      <species id="bridge_CO2" name="bridge CO2 support pool" compartment="tca_cell" initialConcentration="1" boundaryCondition="true" constant="true" />
      <species id="bridge_Pi" name="bridge phosphate support pool" compartment="tca_cell" initialConcentration="1" boundaryCondition="true" constant="true" />
      <species id="bridge_PEP" name="bridge PEP scaffold pool" compartment="tca_cell" initialConcentration="0.01" boundaryCondition="false" constant="false" />
"""


REACTION_SNIPPET = """
      <reaction id="PYR_OAA_support" name="pyruvate to oxaloacetate support bridge" reversible="false" fast="false">
        <notes>
          <html:body>
            <html:p>Weak direct pyruvate to oxaloacetate fill step added to keep the C-N bridge from becoming substrate-starved. This is a conservative scaffold reaction rather than a claim of exact plant kinetics.</html:p>
          </html:body>
        </notes>
        <listOfReactants>
          <speciesReference species="pyr_pyruvate" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="tca_oaa" />
        </listOfProducts>
        <kineticLaw>
          <math:math>
            <math:apply>
              <math:times />
              <math:ci>tca_cell</math:ci>
              <math:ci>Vmax_PYR_OAA_support</math:ci>
              <math:apply>
                <math:divide />
                <math:ci>pyr_pyruvate</math:ci>
                <math:apply>
                  <math:plus />
                  <math:ci>Km_PYR_OAA_support</math:ci>
                  <math:ci>pyr_pyruvate</math:ci>
                </math:apply>
              </math:apply>
            </math:apply>
          </math:math>
          <listOfParameters>
            <parameter id="Vmax_PYR_OAA_support" value="0.008" units="tca_mmlmin" />
            <parameter id="Km_PYR_OAA_support" value="0.25" units="tca_mml" />
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction sboTerm="SBO:0000176" id="PYR_PC_bridge" name="pyruvate carboxylase-like bridge [pyruvate + CO2 + ATP -&gt; OAA + ADP + Pi]" reversible="false" fast="false">
        <notes>
          <html:body>
            <html:p>Pyruvate carboxylase topology added as a conservative bridge between the pyruvate donor module and the TCA/OAA node. The stoichiometry follows the BRENDA pyruvate carboxylase reaction. Km priors now follow wild-type Pyc1-like BRENDA values for pyruvate, ATP, and bicarbonate support, while the bridge Vmax is still kept mild so this step supplements rather than dominates the merged scaffold.</html:p>
          </html:body>
        </notes>
        <listOfReactants>
          <speciesReference species="bridge_CO2" />
          <speciesReference species="pyr_pyruvate" />
          <speciesReference species="pyr_ATP" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="bridge_Pi" />
          <speciesReference species="tca_oaa" />
          <speciesReference species="pyr_ADP" />
        </listOfProducts>
        <kineticLaw>
          <math:math>
            <math:apply>
              <math:times />
              <math:ci>tca_cell</math:ci>
              <math:ci>Vmax_PC_bridge</math:ci>
              <math:apply>
                <math:divide />
                <math:ci>bridge_CO2</math:ci>
                <math:apply>
                  <math:plus />
                  <math:ci>Km_CO2_PC_bridge</math:ci>
                  <math:ci>bridge_CO2</math:ci>
                </math:apply>
              </math:apply>
              <math:apply>
                <math:divide />
                <math:ci>pyr_pyruvate</math:ci>
                <math:apply>
                  <math:plus />
                  <math:ci>Km_PYR_PC_bridge</math:ci>
                  <math:ci>pyr_pyruvate</math:ci>
                </math:apply>
              </math:apply>
              <math:apply>
                <math:divide />
                <math:ci>pyr_ATP</math:ci>
                <math:apply>
                  <math:plus />
                  <math:ci>Km_ATP_PC_bridge</math:ci>
                  <math:ci>pyr_ATP</math:ci>
                </math:apply>
              </math:apply>
            </math:apply>
          </math:math>
          <listOfParameters>
            <parameter id="Vmax_PC_bridge" value="0.018" units="tca_mmlmin" />
            <parameter id="Km_CO2_PC_bridge" value="1.36" units="tca_mml" />
            <parameter id="Km_PYR_PC_bridge" value="0.50" units="tca_mml" />
            <parameter id="Km_ATP_PC_bridge" value="0.07" units="tca_mml" />
          </listOfParameters>
        </kineticLaw>
      </reaction>
      <reaction sboTerm="SBO:0000176" id="OAA_PEPCK_bridge" name="PEPCK-like scaffold [OAA &lt;-&gt; PEP]" reversible="true" fast="false">
        <notes>
          <html:body>
            <html:p>Low-capacity reversible PEPCK-like scaffold added as an optional follow-up carbon-routing step. This currently feeds a dedicated PEP scaffold pool and is intended as a bridge placeholder until a fuller plant-like PEP node is integrated.</html:p>
          </html:body>
        </notes>
        <listOfReactants>
          <speciesReference species="tca_oaa" />
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="bridge_PEP" />
        </listOfProducts>
        <kineticLaw>
          <math:math>
            <math:apply>
              <math:times />
              <math:ci>tca_cell</math:ci>
              <math:apply>
                <math:divide />
                <math:apply>
                  <math:minus />
                  <math:apply>
                    <math:times />
                    <math:ci>Vf_PEPCK_bridge</math:ci>
                    <math:apply>
                      <math:divide />
                      <math:ci>tca_oaa</math:ci>
                      <math:ci>Km_OAA_PEPCK_bridge</math:ci>
                    </math:apply>
                  </math:apply>
                  <math:apply>
                    <math:times />
                    <math:ci>Vr_PEPCK_bridge</math:ci>
                    <math:apply>
                      <math:divide />
                      <math:ci>bridge_PEP</math:ci>
                      <math:ci>Km_PEP_PEPCK_bridge</math:ci>
                    </math:apply>
                  </math:apply>
                </math:apply>
                <math:apply>
                  <math:plus />
                  <math:cn type="integer">1</math:cn>
                  <math:apply>
                    <math:divide />
                    <math:ci>tca_oaa</math:ci>
                    <math:ci>Km_OAA_PEPCK_bridge</math:ci>
                  </math:apply>
                  <math:apply>
                    <math:divide />
                    <math:ci>bridge_PEP</math:ci>
                    <math:ci>Km_PEP_PEPCK_bridge</math:ci>
                  </math:apply>
                </math:apply>
              </math:apply>
            </math:apply>
          </math:math>
          <listOfParameters>
            <parameter id="Vf_PEPCK_bridge" value="0.004" units="tca_mmlmin" />
            <parameter id="Vr_PEPCK_bridge" value="0.002" units="tca_mmlmin" />
            <parameter id="Km_OAA_PEPCK_bridge" value="0.01" units="tca_mml" />
            <parameter id="Km_PEP_PEPCK_bridge" value="0.05" units="tca_mml" />
          </listOfParameters>
        </kineticLaw>
      </reaction>
"""


def patch_xml(path: Path) -> None:
    text = path.read_text()

    if 'id="bridge_CO2"' not in text:
        if "</listOfSpecies>" not in text:
            raise RuntimeError("could not find </listOfSpecies> anchor")
        text = text.replace("</listOfSpecies>", SPECIES_SNIPPET + "    </listOfSpecies>", 1)

    if 'id="PYR_OAA_support"' not in text:
        anchor = '      <reaction sboTerm="SBO:0000176" id="GOT_AspAT"'
        if anchor not in text:
            raise RuntimeError("could not find GOT_AspAT reaction anchor")
        text = text.replace(anchor, REACTION_SNIPPET + anchor, 1)

    path.write_text(text)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: patch_merged_nitrogen_connected_with_pyruvate_oaa.py <merged.xml>", file=sys.stderr)
        return 1
    patch_xml(Path(argv[1]))
    print(f"patched {argv[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

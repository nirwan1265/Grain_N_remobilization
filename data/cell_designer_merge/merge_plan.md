# CellDesigner Merge Prep

## What makes sense

Yes: combining these into one working pathway makes sense as a curation project.
The current files should not be pasted together directly for quantitative simulation, because they mix organisms and kinetic assumptions.

## Model snapshot

- `Hoefnagel2002_PyruvateBranches` (Lactococcus lactis): 19 species, 14 reactions, role `Hold_out_for_now`, bridge metabolites: none
- `Chassagnole2001_Threonine Synthesis` (Unknown_or_unspecified): 11 species, 7 reactions, role `Needs_manual_curation`, bridge metabolites: ADP, ATP, Aspartate, Aspartyl_phosphate, Homoserine, NADP, NADPH, Threonine
- `Curien2009_Aspartate_Metabolism` (Arabidopsis): 31 species, 18 reactions, role `Backbone_for_plant_branch`, bridge metabolites: Aspartate, Aspartate_semialdehyde, Aspartyl_phosphate, Homoserine, Lysine, Threonine
- `Bruggeman2005_AmmoniumAssimilation` (E. coli): 22 species, 16 reactions, role `Bacterial_connector_only`, bridge metabolites: 2_oxoglutarate, ADP, ATP, Ammonium, Glutamate, Glutamine, NADP, NADPH
- `Singh2006_TCA_Ecoli_glucose` (E. coli): 12 species, 11 reactions, role `Bacterial_connector_only`, bridge metabolites: none

## Recommended merge order

1. Use `Curien2009_Aspartate_Metabolism` as the backbone for the amino-acid branch.
2. Treat the ammonium assimilation and TCA models as donor modules for nitrogen and carbon supply logic, not as direct one-click kinetic inserts.
3. Add bridge metabolites manually: `Aspartate`, `Glutamate`, `Glutamine`, `Ammonium`, `2_oxoglutarate`, and `Oxaloacetate`.
4. Replace organism-specific enzyme pools and rate constants only after choosing a target species for the final model.

## Best immediate simulation targets

- Increase aspartate kinase capacity in the Curien model by raising enzyme species `AK1`, `AK2`, `AKI`, or `AKII`.
- Increase aspartate kinase capacity in the Chassagnole model by raising `vm11` and/or `vm13`.
- Reduce glutamate in the Bruggeman model by lowering species `GLU` or by reducing upstream `vgdh` / `vgog` capacity.
- Impose nitrogen limitation in the Bruggeman model by lowering `NH4` and tracking the response of `GLU`, `GLN`, and assimilation fluxes.

## Strong caution

- The Curien backbone is Arabidopsis.
- Bruggeman and Singh are E. coli.
- Hoefnagel is Lactococcus lactis.
- A single merged file can still be useful for hypothesis generation, but any quantitative claim will need species-specific re-parameterization.

## Shared bridge metabolites detected

- `ADP` appears in 2 models: Chassagnole2001_Threonine Synthesis; Bruggeman2005_AmmoniumAssimilation
- `ATP` appears in 2 models: Chassagnole2001_Threonine Synthesis; Bruggeman2005_AmmoniumAssimilation
- `Aspartate` appears in 2 models: Chassagnole2001_Threonine Synthesis; Curien2009_Aspartate_Metabolism
- `Aspartyl_phosphate` appears in 2 models: Chassagnole2001_Threonine Synthesis; Curien2009_Aspartate_Metabolism
- `Homoserine` appears in 2 models: Chassagnole2001_Threonine Synthesis; Curien2009_Aspartate_Metabolism
- `NADP` appears in 2 models: Chassagnole2001_Threonine Synthesis; Bruggeman2005_AmmoniumAssimilation
- `NADPH` appears in 2 models: Chassagnole2001_Threonine Synthesis; Bruggeman2005_AmmoniumAssimilation
- `Threonine` appears in 2 models: Chassagnole2001_Threonine Synthesis; Curien2009_Aspartate_Metabolism

## Why the Curien model is the best starting point

It already contains the aspartate-derived branch logic you care about, including multiple aspartate kinase isoforms and downstream threonine / lysine competition.

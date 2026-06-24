# BRENDA Candidate Parameters For Grain_N Model

Source file:
- `/Users/nirwantandukar/Documents/Research/data/Cell_designer/Grain_N/brenda_2026_1.txt`

## Bottom line

Yes, BRENDA can help, but selectively.

Best uses:
- replace the placeholder `GOT_AspAT` connector with plant AspAT kinetic priors
- sanity-check / refine Asp kinase regulation and substrate affinity
- provide plant GS affinity ranges if we later rebuild the NH4 assimilation side

Less useful as direct drop-ins:
- the full bacterial GS/GDH/GOGAT module, because the current merged model still mixes organisms and unit systems
- any value without clear assay conditions, substrate definition, and organism context

## Strongest immediate candidates

### 1. Aspartate aminotransferase / GOT

EC:
- `2.6.1.1`

Why it matters:
- this is the current placeholder bridge `tca_oaa + nh4_GLU <-> asp_Asp + tca_akg`

Best records found:
- Arabidopsis thaliana, PR `#19`
  - protein summary includes isozymes `AAT1`, `AAT2`, `AAT5`
  - turnover numbers:
    - `TN #19 205 {2-oxoglutarate}` at pH 8.0, 25C
    - `TN #19 217 {2-oxoglutarate}` for AAT2
    - `TN #19 176 {2-oxoglutarate}` for AAT5
    - `TN #19 319 {oxaloacetate}` for AAT1
    - `TN #19 574 {oxaloacetate}` for AAT2
    - `TN #19 279 {oxaloacetate}` for AAT5
  - specific activity:
    - `SA #19 454.0` for recombinant AAT1
    - `SA #19 2934.0` for recombinant AAT2
    - `SA #19 165.0` for recombinant AAT5
- Arabidopsis thaliana, PR `#68`
  - recombinant wild-type `AtPAT` / bifunctional AspAT-prephenate aminotransferase
  - `KM #68 4.91 {L-aspartate}` at pH 7.5, 25C
  - `KM #68 9.24 {2-oxoglutarate}` at pH 7.5, 25C
  - `TN #68 9.33 {L-aspartate}` at pH 7.5, 25C
  - `TN #68 10.16 {2-oxoglutarate}` at pH 7.5, 25C

Interpretation:
- Arabidopsis PR `#19` is the better conceptual match for a plant AspAT bridge.
- PR `#68` has a cleaner paired `KM` and `TN`, but it is bifunctional and may not be the best biological analogue for maize GOT.

Recommendation:
- use Arabidopsis plant AspAT data as the first non-arbitrary replacement for the current placeholder mass-action `kf` / `kr`
- if we want a quick first pass, build reversible Michaelis-Menten style kinetics around plant-scale `Km(Asp)` and `Km(2-OG)` from BRENDA rather than keeping arbitrary `0.1` / `0.01`

### 2. Aspartate kinase

EC:
- `2.7.2.4`

Why it matters:
- the Asp branch is already running, but this is one of the main perturbation targets you wanted

Best records found:
- Zea mays, PR `#15`
  - explicit maize record
  - `KM #15 1.04 {L-aspartate}` at pH 8.0, 25C
  - `KM #15 0.43 {ATP}` and `KM #15 0.56 {ATP}` lines appear in merged record context
  - `IC50 #15 0.08 {L-lysine}` for maize ask variants
  - `SA #15 18.9 units/mg`
- Arabidopsis thaliana, PR `#17`
  - bifunctional AK-HSDH
  - `KM #17 11.6 {L-aspartate}` at pH 8.0, 30C
  - `KM #17 5.5 {ATP}` at pH 8.0, 30C
  - `KI #17 0.049 {L-threonine}`
  - `KI #17 0.091 {L-threonine}`
  - `SA #17 5.4` for forward reaction
- Arabidopsis thaliana, PR `#56` and `#57`
  - monofunctional AK1 / AK3
  - `KM #56 2.037 {L-aspartate}`
  - `KM #56 1.7 {ATP}`
  - `KM #57 1.095 {L-aspartate}`
  - `TN #56 23.4 {ATP}`
  - `TN #57 8.4 {ATP}`

Interpretation:
- maize values exist and are especially useful for inhibition / affinity sanity checks.
- Arabidopsis values are directly relevant to the Curien-style Asp branch already in the model.

Recommendation:
- use these to justify Asp kinase perturbation ranges and possibly to compare Curien embedded constants against independent BRENDA values
- do not blindly overwrite the Curien kinetic model unless we deliberately choose to re-parameterize the Asp module

### 3. Glutamine synthetase

EC:
- `6.3.1.2`

Why it matters:
- the NH4 assimilation side crashes early in the integrated run

Best records found:
- Arabidopsis thaliana, PR `#121` isoenzyme `GLN1;3`
  - `KM #121 0.85 {ATP}` at pH 7.8, 30C
  - `KM #121 3.9 {L-glutamate}` at pH 7.8, 30C
  - `KM #121 1.21 {NH4+}` at pH 7.8, 30C
  - additional mutant lines also present
- Arabidopsis thaliana, PR `#31` isoenzyme `GLN1;2 / GLN1;4`
  - `KM #31 1.1 {ATP}` at pH 7.8, 30C
  - `KM #31 3.8 {L-glutamate}` at pH 7.8, 30C
  - `KM #31 2.45 {NH4+}` at pH 7.8, 30C
  - `KM #31 0.12 {NH4+}` for GLN1;4 in another assay
- Arabidopsis thaliana, PR `#125` isoenzyme `GLN1;1`
  - `KM #125 0.01 {NH4+}` at pH 7.8, 30C
- Zea mays protein entries are present:
  - PR `#28`, `#130`, `#167`, `#175`
  - the extracted section clearly shows maize proteins and expression / isoform information
  - but the strongest readily visible `KM` lines in this pass are Arabidopsis rather than maize

Interpretation:
- these are useful if we decide to plant-ify the NH4 assimilation side.
- they are not a drop-in replacement for the current bacterial Bruggeman module.

Recommendation:
- use plant GS values as reference targets if we replace or simplify the bacterial NH4 assimilation module

## What I would use first

Priority order:
1. Replace `GOT_AspAT` placeholder kinetics using plant AspAT priors from EC `2.6.1.1`
2. Use maize / Arabidopsis Asp kinase values from EC `2.7.2.4` to justify perturbation ranges and inhibition settings
3. Keep GS values from EC `6.3.1.2` on deck for a later rebuild of the N assimilation side

## Round 2 manuscript-facing branches

### 4. Asparagine synthetase / asparaginase branch

Why it matters:
- the manuscript repeatedly frames asparagine as a dominant grain FAA and a key N transport / storage compound

Strong BRENDA support:
- glutamine-dependent asparagine synthetase includes maize isozymes `ZmAsn1-4`
- useful maize-like values:
  - `Km(Asp)` about `0.91-0.98`
  - `Km(Gln)` about `0.09-0.54`
  - `Km(ATP)` about `0.097-0.128`
- asparaginase also has abundant `Km(L-Asn)` records across systems

How I would use it in this merged model:
- add a reduced-form `Asp + Gln -> Asn + Glu` connector for the Asn branch
- add a conservative `Asn -> Asp + NH4` remobilization step
- keep the ATP bookkeeping implicit unless we later rebuild the energy layer consistently

### 5. Omega-amidase / NIT2-like branch

Why it matters:
- one of the cleaner soil-N manuscript candidates is `Zm00001eb005710`, annotated as omega-amidase / NIT2-like

Strong BRENDA support:
- BRENDA explicitly lists:
  - `Arabidopsis thaliana`
  - `Zea mays` protein entries
- direct omega-amidase kinetics for `2-oxoglutaramate` are available
  - `TN` roughly `50.8-72.7`
  - `Km` roughly `5.31-6.49`

How I would use it in this merged model:
- add a reduced glutamine-linked recycling branch representing the net glutamine-transaminase + omega-amidase pathway
- interpret it as manuscript-facing structural support for Gln-linked recycling into `2-OG + NH4`, not as a fully resolved plant transaminase subsystem

### 6. Proline axis

Why it matters:
- proline is one of the clearest Indian Chief FAA remodeling signals in the manuscript

Strongest BRENDA support:
- `pyrroline-5-carboxylate reductase` has strong plant-facing kinetics
  - common `Km(P5C)` values in the low sub-mM range
  - common `Km(NADPH)` values in the low sub-mM range
- `glutamate 5-kinase` has usable entry-step priors
  - `Km(ATP)` around `0.5-0.6` in some records
  - `Km(Glu)` around `10-12` in representative records
- `glutamate-5-semialdehyde dehydrogenase` also exists in BRENDA, but the whole plant `P5CS` axis is cleaner as a reduced-form scaffold than as a one-shot direct transplant

How I would use it in this merged model:
- add a reduced `Glu -> P5C` entry step using glutamate dependence only
- add `P5C -> Pro` using P5CR-like substrate affinity
- treat this branch as structurally meaningful now, but still quantitatively tentative because the integrated Glu pool remains under-scaled

### 7. Pyruvate carboxylase bridge

Why it matters:
- we added a pyruvate-to-OAA bridge so the pyruvate donor block can actually talk to the TCA/Asp node

Strong BRENDA support:
- `RN pyruvate carboxylase`
- directly matching reaction:
  - `ATP + pyruvate + HCO3- = ADP + phosphate + oxaloacetate`
- clean wild-type Pyc1-like affinity set visible in the local BRENDA dump:
  - `KM #12,56,57 0.50 {pyruvate}`
  - `KM #12,56,57 0.07 {ATP}`
  - `KM #12,56,57 1.36 {HCO3-}`
- additional wild-type set also exists:
  - `KM #61 0.15 {pyruvate}`
  - `KM #61 0.145 {ATP}`
  - `KM #61 10.8 {HCO3-}`

Interpretation:
- the bridge stoichiometry we added is biologically right
- the bridge affinity terms can be grounded in BRENDA
- the bridge `Vmax` should still be treated as a conservative integration setting unless we deliberately rebuild pyruvate/TCA entry quantitatively

Recommendation:
- use the coherent Pyc1-like set (`pyruvate 0.50`, `ATP 0.07`, `HCO3- 1.36`) as the default BRENDA-backed bridge affinities
- keep `Vmax_PC_bridge` explicitly labeled as provisional

## What I would not do yet

- I would not overwrite the Curien Asp model wholesale with BRENDA values
- I would not mix BRENDA assay values into the current bacterial NH4 module without reconciling units and model structure first
- I would not use expression-only or phenotype-only records as kinetic constants

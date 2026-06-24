# Grain N Parameter Provenance Audit

Updated: 2026-06-22

## Bottom line

Yes, some values in the current merged model are still scaffolded or integration-tuned.

The model currently uses a mix of:

- source Biomodel kinetics kept from donor modules
- BRENDA-informed plant priors
- conservative reduced-form branch constants
- integration retunes added so the merged network behaves sensibly enough for perturbation testing

So the right question is not "did we make things up at all", but "which values are still provisional and worth checking next".

## Buckets

### 1. Mostly donor-model values

These are primarily inherited from source models, even if some surrounding connectors were later adjusted.

- Asp-family core from Curien-style Asp model
- NH4 assimilation structure from Bruggeman ammonium-assimilation model
- TCA backbone from Singh/Millard-style donor logic
- Pyruvate core from donor pyruvate modules

These are not arbitrary, but they are also not necessarily maize-specific.

### 2. BRENDA-backed or BRENDA-anchored

These have direct plant-facing support in `brenda_2026_1.txt`, even if the final merged value is still a representative prior rather than a strict one-to-one transplant.

| Component | Current model values | Status | Notes |
| --- | --- | --- | --- |
| AspAT / GOT bridge | `Vf_AspAT=0.25`, `Vr_AspAT=0.39`, `Km_OAA=0.039`, `Km_GLU=8.2`, `Km_Asp=9.9`, `Km_AKG=0.25` | BRENDA-informed | Km side is clearly anchored to plant AspAT priors; Vmax side is still a merged-model choice |
| Asp kinase | donor Curien constants plus perturbation tests | donor + BRENDA sanity-check | BRENDA supports maize/Arabidopsis ranges; we mostly used it as validation, not wholesale replacement |
| ASNS branch | `Km_Asp_ASNS=0.95`, `Km_Gln_ASNS=0.25` | BRENDA-backed | Fits maize-like ASNS ranges from BRENDA |
| ASNASE branch | `Km_Asn_ASNASE=0.06` | BRENDA-backed | Uses a concrete BRENDA-like low-submillimolar asparaginase substrate affinity instead of a generic scaffold value |
| OMEGA/NIT2 branch | branch existence only | BRENDA-supported structurally | BRENDA supports omega-amidase in maize/Arabidopsis, but our collapsed `GLN -> AKG + NH4` rate law is still a scaffold |
| P5CS entry | `Km_Glu_P5CS=10` | BRENDA-backed at component level | Anchored to glutamate 5-kinase style Glu affinity, but the branch is reduced-form |
| P5CR step | `Km_P5C_P5CR=0.18` | BRENDA-backed | Fits plant-facing P5CR substrate affinity ranges |
| GS affinity ranges | not directly transplanted | BRENDA available | Plant GS records exist and are good candidates for the next pass |
| Pyruvate carboxylase bridge | `Km_CO2_PC_bridge=1.36`, `Km_PYR_PC_bridge=0.50`, `Km_ATP_PC_bridge=0.07` | BRENDA-backed affinities | Affinity terms now follow a coherent Pyc1-like BRENDA set; `Vmax_PC_bridge` remains provisional |

### 3. Scaffolded or integration-tuned

These are the main values we should treat as provisional.

| Component | Current model values | Why provisional |
| --- | --- | --- |
| NH4 import support | `nh4_Vmax_NH4_import=0.12` | added to keep external N testable in the merged model |
| Primary N support term | `nh4_Vmax_primary_assim=0.015`, `nh4_Km_extN_primary=0.5`, `nh4_Km_AKG_primary=0.05` | merged-model support term, not a direct donor reaction |
| Retuned NH4 module capacities | `nh4_Vgs=30`, `nh4_Vgdh=10`, `nh4_Vgog=8`, `nh4_Kgogkg=0.2` | post-merge stabilization retunes |
| Retuned reduced-N demand sinks | `nh4_Vgludem=90`, `nh4_Vglndem=52.5` | explicit integration retune of abstract Bruggeman demand terms |
| TCA biosynthetic sink | `tca_kbiosyn=0.005` | integration retune |
| Pyruvate carboxylase bridge | `Vmax_PC_bridge=0.018` | bridge capacity still provisional even though the affinity terms are now BRENDA-backed |
| PEPCK bridge | `Vf_PEPCK_bridge=0.004`, `Vr_PEPCK_bridge=0.002`, `Km_OAA=0.01`, `Km_PEP=0.05` | manually added bridge and currently low-impact |
| ASNS Vmax | `Vmax_ASNS=1.0` | reduced-form branch capacity chosen conservatively |
| ASNASE Vmax | `Vmax_ASNASE=0.35` | representative scaffold value |
| OMEGA/NIT2 collapsed branch | `Vmax_OMEGA=0.75`, `Km_Gln_OMEGA=0.25` | not directly measurable from omega-amidase BRENDA because the model collapses upstream chemistry |
| P5CS Vmax | `Vmax_P5CS=1.0` | reduced-form proline entry capacity |
| P5CR Vmax | `Vmax_P5CR=0.6` | reduced-form proline finishing capacity |
| Thr-side throttle | `asp_TS1_AdoMEt_Km_no_AdoMet_exp=0.31`, `asp_TS1_Phosphate_Ki_exp=100`, `asp_TS1_kcatmin_exp=1.0`, `asp_TS1_AdoMet_kcatmax_exp=6.0` | tuned after merge using Chassagnole-informed logic rather than direct one-shot maize kinetic transplant |

## Did we "make things up"?

Not in the bad sense, but yes in the modeling sense:

- we added reduced-form support reactions where the merged network otherwise did not carry flux well
- we introduced bridge reactions between modules that were not fully connected in the donor models
- we assigned conservative placeholder capacities to new manuscript-facing branches
- we retuned some abstract sink terms after integration

That is normal for a scaffold/integration model, but those parameters should be clearly labeled as provisional.

## Highest-priority BRENDA checks from here

If we want the next best return on effort, the order should be:

1. **GS / GOGAT / GDH / NH4 support side**
   - biggest biological payoff
   - most important for the manuscript nitrogen story

2. **Pyruvate carboxylase bridge**
   - now that pyruvate-to-OAA was added, this is the cleanest carbon-entry step to curate from BRENDA

3. **ASNASE**
   - easy cleanup because BRENDA has plenty of asparaginase kinetic records

4. **OMEGA/NIT2 branch**
   - worth revisiting if we want to replace the collapsed `GLN -> AKG + NH4` proxy with a more explicit two-step branch

5. **P5CS / P5CR Vmax values**
   - useful, but lower priority than the N support side

## What I would not spend time on right now

- PEPCK bridge, unless we decide we really need that route, because it currently has little effect
- wholesale replacement of the Curien Asp module
- blind substitution of random BRENDA values without matching organism, substrate definition, and assay context

## Practical interpretation

For the current model:

- **Asp-family directionality** is on relatively strong footing
- **AspAT / C-N bridge logic** is reasonably grounded
- **Asn / Pro / omega-amidase additions** are structurally manuscript-relevant, but still quantitatively provisional
- **NH4 support and reduced-N demand tuning** are the least source-grounded part of the current merged model and should be the first place we tighten with BRENDA

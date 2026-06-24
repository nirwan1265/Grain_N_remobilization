## Grain N modeling pipeline notes

This folder contains a paper-oriented bridge between `main.tex` candidate genes and the CellDesigner perturbation scenarios.

Core files:

- `grain_n_candidate_model_map.tsv`
  - maps manuscript candidates or pathway modules to scenario XML / export labels
- `scripts/plot_manuscript_perturbations.R`
  - scans a directory of exported CellDesigner `.txt` runs and writes manuscript-friendly summaries

Expected workflow:

1. Run scenario XML files in CellDesigner.
2. Export each simulation as a `.txt` file into one folder.
3. Keep the baseline export as `base_values.txt` or `baseline.txt`.
4. Keep perturbation exports named after the scenario, e.g. `asadh_x2.txt`, `gs_gogat_half.txt`, `akg_pool_x2.txt`.
5. Run:

```bash
Rscript scripts/plot_manuscript_perturbations.R \
  /path/to/exported_txt_runs \
  data/cell_designer_simulations/manuscript_panel
```

Outputs:

- `manuscript_endpoint_heatmap.png`
- `manuscript_nitrogen_proxy_heatmap.png`
- `manuscript_core_timecourses.png`
- `manuscript_endpoint_summary.tsv`
- `manuscript_nitrogen_proxy_summary.tsv`
- `scenario_file_manifest.tsv`

Important caveat:

- The nitrogen summaries are **model-side proxies**, not direct elemental nitrogen measurements.
- They are useful for comparing scenarios consistently within this reduced model.
- True nitrogen uptake, influx, utilization, or export rates would require explicit source/sink reactions in the SBML model.

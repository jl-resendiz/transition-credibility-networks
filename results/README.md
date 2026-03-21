# Results

Outputs produced by `src/` scripts. Do not edit manually — regenerate by running the pipeline.

## Structure

| Folder | Contents |
|---|---|
| `metrics/` | Regression coefficient tables (one `.md` per specification) |
| `summaries/` | Aggregated CSV and markdown summaries for key results |
| `tables/` | LaTeX table source files for the manuscript |
| `interconnectors/` | Cross-border interconnector event-study outputs |

## Metrics

Regression metrics are named by strategy and transform variant, e.g.:

- `strategy1_panel_metrics.md` — two-way FE panel regression
- `strategy2_panel_did_base.md` — spatial DiD, untransformed exposure
- `strategy3_phaseout_event_time_tier1_log1p.md` — phase-out event study, tier-1 binding laws, log(1+x) transform

## Tables

LaTeX tables imported directly into `manuscript/main.tex`:

| File | Table |
|---|---|
| `table_spec_progression.tex` | Specification progression (Table 2) |
| `table_channel_controls.tex` | Channel decomposition with controls |
| `table_channel_sensitivity.tex` | Robustness across transforms |
| `table_bandwidth_sensitivity.tex` | Kernel bandwidth sensitivity |
| `table_weight_correlations.tex` | Spatial weight matrix pairwise correlations |
| `table_vif_controls.tex` | Variance inflation factors |
| `table_placebo_shuffle.tex` | Placebo: shuffled fuel matrix |

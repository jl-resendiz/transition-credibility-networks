# Source Code

Python 3.8+ standard library only (no pandas, numpy, or scipy). The only non-stdlib dependency is `openpyxl` for reading GEM Excel files.

## Entry Point

```bash
python src/run_all.py              # full pipeline (build + analysis)
python src/run_all.py --analysis   # analysis only (skip build steps)
```

## Path Resolution

`_paths.py` resolves all data and results paths relative to the repository root:

```python
from _paths import raw_path, derived_path, results_path
```

## Return Data

Monthly equity returns use total return series exclusively:
- **US firms**: CRSP total returns (`trt1m`), the gold standard
- **Non-US firms**: Eikon TR.TotalReturn via `pull_eikon_returns.py`, which includes dividends

For non-US firms not covered by Eikon (~105 firms), deduplicated Compustat Global Security price returns are used as a fallback (Ince & Porter 2006). Compustat Global daily returns are used for all non-US firms (with deduplication) because Eikon does not provide daily total returns via the batch API.

Filters per Ince & Porter (2006): price >= $1, daily volume >= 1000, monthly returns capped at +/-100%, daily at +/-50%. Eikon returns cleaned to end-of-month only with illiquid firms removed.

## Build Scripts (Stage 1-5)

| Stage | Script | Output |
|---|---|---|
| 1 | `parse_gem.py` | `gem_parents_parsed.csv` |
| 2 | `match_gem_compustat.py` | `gem_compustat_matches.csv` |
| 3 | `build_fundamentals.py` | `firm_fundamentals.csv` |
| 3 | `build_ets_matrix.py` | `weight_matrix_W_regulatory.csv`, `firm_ets_membership.csv` |
| 3 | `build_weight_matrix.py` | `firm_centroids.csv`, `weight_matrix_W_geo.csv` |
| 3 | `build_fuel_matrix.py` | `weight_matrix_W_fuel.csv` |
| 3 | `build_time_varying_alpha.py` | `firm_alpha_panel.csv` |
| 3 | `compute_returns.py` | `daily_returns.csv`, `monthly_returns.csv` |
| 4 | `build_retirement_events.py` | `coal_retirement_events.csv` |
| 4 | `build_coal_phaseout_events.py` | `coal_phaseout_shocks_events.csv` |
| 4 | `build_eia860_announcement_events.py` | `eia860_announcement_events.csv` |
| 5 | `summary_statistics.py` | `results/summaries/summary_statistics.md` |

## Analysis Scripts (Stage 6)

Each script maps to a specific takeaway (T1-T10) from the paper.

| Script | Takeaway | Purpose |
|---|---|---|
| `strategy2_robust_inference.py` | T1, T8 | FM+NW, event-clustered, two-way clustered inference |
| `strategy2_difference_test_summary.py` | T8 | Difference tests: FM, sign, Wilcoxon, randomization |
| `strategy2_bandwidth_fmb.py` | T1 | Bandwidth sensitivity (250-1500km) |
| `strategy2_joint_tests.py` | T1, T8 | Joint F-test, pooled event-clustered coefficients |
| `strategy2_firm_level_test.py` | T1 | Two-way clustered (firm + event) |
| `strategy2_event_specific_geo.py` | T2 | Event-specific w_geo (minimum distance, not centroid) |
| `strategy2_geo_diversification.py` | T2 | Single-country subsample + diversification interaction |
| `strategy2_credibility_interaction.py` | T3, T4 | ETS binary and carbon price continuous interactions |
| `strategy2_esg_ets_fmb.py` | T5 | ESG horse race under FM+NW inference |
| `strategy2_esg_horse_race.py` | T5 | ESG horse race, event-clustered |
| `strategy2_learning_alternatives.py` | T6 | Learning order test (Bayesian vs cascading) |
| `strategy2_spatial_score.py` | T7 | STS in-sample and out-of-sample holdout |
| `strategy3_phaseout_wild_bootstrap.py` | T9 | Wild cluster bootstrap (14 clusters) |
| `strategy2_romano_wolf.py` | T10 | Romano-Wolf stepdown correction |
| `strategy2_referee_tables.py` | Appendix | Robustness tables (correlations, placebo, progression) |

## Data Pull Scripts (not in pipeline)

| Script | Purpose |
|---|---|
| `pull_eikon_returns.py` | Pull Eikon TR.TotalReturn for non-US firms (requires API key) |
| `pull_refinitiv_esg.py` | Pull LSEG Eikon ESG scores (requires API key) |
| `pull_refinitiv_panel.py` | Pull LSEG Eikon financial panel (requires API key) |
| `pull_refinitiv_extra.py` | Pull LSEG Eikon governance extras (requires API key) |

## R Scripts

| Script | Purpose |
|---|---|
| `robustness_conley_se.R` | Conley (1999) spatial SEs at 500/1000km (requires R + fixest) |

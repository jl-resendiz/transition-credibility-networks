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

9 scripts, each addressing a specific identification threat or result.

| Script | Purpose |
|---|---|
| `strategy2_robust_inference.py` | Main results: FM+NW, event-clustered, two-way clustered, window sensitivity |
| `strategy2_joint_tests.jl` | Joint F-test + difference test (Julia; falls back to .py) |
| `strategy2_esg_horse_race.py` | ESG vs fuel horse race on 153-firm subsample |
| `strategy2_bartik_shiftshare.py` | Shift-share causal diagnostics (GPS 2020, Oster bounds) |
| `strategy2_romano_wolf.py` | Romano-Wolf stepdown correction (3 hypotheses, hybrid Julia bootstrap) |
| `strategy2_geo_diversification.py` | Aggregation lemma test (HHI + multi-country interaction) |
| `strategy2_learning_alternatives.py` | Geographic heterogeneity (US vs non-US, calendar time) |
| `strategy2_referee_compute.py` | Heavy computation for appendix tables (exports JSON) |
| `strategy2_referee_tables.py` | LaTeX table formatting from JSON (< 1 second) |

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

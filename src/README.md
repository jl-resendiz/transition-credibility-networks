# Source Code

The main pipeline uses Python 3.8+ standard library only (no pandas, numpy, or scipy).

Econometric robustness checks use **R** where the R ecosystem is decisively stronger (Conley SEs, wild cluster bootstrap, CR2 few-cluster inference). See the R Scripts section below.

## Entry Point

```bash
python src/run_pipeline.py
```

Runs all build and strategy scripts in dependency order.

## Path Resolution

`_paths.py` resolves all data and results paths relative to the repository root.
Import it in any script with:

```python
from _paths import raw_path, derived_path, results_path
```

## Build Scripts

These construct derived datasets from raw inputs. Run in order before analysis.

| Script | Output |
|---|---|
| `build_fundamentals.py` | Firm financials panel from Compustat |
| `build_alpha_trajectories.py` | Fossil intensity trajectories |
| `build_time_varying_alpha.py` | Time-varying fossil intensity (alpha) |
| `build_weight_matrix.py` | Combined spatial weight matrix |
| `build_fuel_matrix.py` | Fuel-similarity weight matrix (cosine similarity) |
| `build_ets_matrix.py` | ETS membership binary matrix |
| `build_coal_phaseout_events.py` | Coal phase-out policy shocks |
| `build_eia860_announcement_events.py` | EIA-860 US coal retirement announcements |
| `build_retirement_events.py` | GEM plant retirement events |
| `compute_returns.py` | Market-adjusted daily/monthly returns |
| `parse_gem.py` | Parse GEM tracker XLSX files |
| `match_gem_compustat.py` | Match GEM utilities to Compustat identifiers |

## Analysis Scripts

Five research strategies corresponding to the empirical framework.

**Strategy 1: R² Channel Decomposition**

| Script | Purpose |
|---|---|
| `strategy1_r2_test.py` | R² split: fuel vs. geographic vs. regulatory channels |
| `strategy1_r2_kernel_sensitivity.py` | Kernel bandwidth sensitivity |
| `strategy1_r2_transform_experiments.py` | Robustness across transforms (winsor, log-asinh, z-score) |
| `strategy1_panel_regression.py` | Two-way fixed-effects panel regression |
| `strategy1_placebo.py` | Placebo: shuffle fuel similarity matrix |
| `strategy1_forward_r2_test.py` | Forward (lead-time) R² test |

**Strategy 2: Spatial Transmission in Returns**

| Script | Purpose |
|---|---|
| `strategy2_spatial_regression.py` | Core spatial event-study regression (CAR on network exposure) |
| `strategy2_spatial_regression_summarize.py` | Summarize and format spatial results |
| `strategy2_panel_did.py` | Difference-in-differences with spatial treatment |
| `strategy2_panel_did_two_part.py` | Two-part DiD (tier-1 binding vs. other) |
| `strategy2_panel_did_wild_bootstrap.py` | Wild bootstrap inference |
| `strategy2_event_time_fuel_plot.py` | Event-time figure by fuel-similarity channel |
| `strategy2_event_fmb_foreign.py` | Foreign-interconnector robustness (Fama-MacBeth) |
| `strategy2_eia860_event_time_did.py` | EIA-860 announcement event study |
| `strategy2_eia860_diagnostics.py` | EIA-860 fuel composition diagnostics |
| `strategy2_eia860_volatility_mediation.py` | Volatility mediation analysis |
| `strategy2_placebo_leaveout.py` | Leave-one-out placebo tests |

**Strategy 3: Policy Shocks and the Credibility Gap**

| Script | Purpose |
|---|---|
| `strategy3_policy_shocks.py` | Coal phase-out policy event study |
| `strategy3_phaseout_event_time_did.py` | Phase-out event-time DiD |
| `strategy3_phaseout_event_time_volatility.py` | Phase-out volatility response |
| `strategy3_phaseout_diagnostics.py` | Phase-out exposure diagnostics |
| `strategy3_phaseout_coalshare_event_time.py` | Coal-share response to phase-outs |
| `strategy3_phaseout_coalshare_panel_did.py` | Coal-share DiD |
| `strategy3_phaseout_coalshare_diagnostics.py` | Coal-share exposure diagnostics |
| `strategy3_phaseout_car_robustness.py` | Cumulative abnormal return robustness |

**Strategy 4: Quantile Regression**

| Script | Purpose |
|---|---|
| `strategy4_quantile_regression.py` | Quantile regression: heterogeneous network effects |

**Strategy 5: ESG vs. Physical Fundamentals**

| Script | Purpose |
|---|---|
| `strategy5_esg_forward_delivery.py` | ESG scores vs. alpha in forward fossil-share delivery |

## Diagnostic and Utility Scripts

| Script | Purpose |
|---|---|
| `diagnose_alpha_measurement.py` | Fossil intensity measurement quality |
| `validate_alpha_trajectory.py` | Trajectory cross-validation |
| `summary_statistics.py` | Summary statistics table |
| `strategy2_interconnector_diagnostics.py` | Cross-border interconnector exposure |
| `strategy2_referee_tables.py` | Referee-response robustness tables |
| `plot_maps.py` | Geographic mapping utility |

## Data Pull Scripts

| Script | Purpose |
|---|---|
| `pull_refinitiv_esg.py` | Download Refinitiv ESG panel |
| `pull_refinitiv_panel.py` | Download Refinitiv financial panel |
| `pull_refinitiv_extra.py` | Download Refinitiv governance extras |
| `scrape_gem_wiki_dates.py` | Scrape GEM wiki announcement dates |
| `update_event_dates.py` | Update event dates from external sources |
| `update_announcement_dates.py` | Update announcement dates |

## R Scripts (Econometric Robustness)

These scripts use R where the package ecosystem is stronger than Python. Requires R 4.5+ with `fixest`, `clubSandwich`, and `data.table`.

```bash
Rscript src/robustness_conley_se.R
Rscript src/robustness_wild_bootstrap.R
Rscript src/robustness_volatility_measures.R
```

| Script | Purpose | Key R Packages |
|---|---|---|
| `robustness_conley_se.R` | Conley (1999) SEs for Table 2 channel decomposition at 250/500/1000 km | `fixest::conley()` |
| `robustness_wild_bootstrap.R` | Wild cluster bootstrap + CR2 SEs for phase-out DiD (G=14) | `fixest`, `clubSandwich` |
| `robustness_volatility_measures.R` | Beaver (1968) \|AR\|-based and squared-return volatility measures | `fixest` |

### Why R?

The econometrics skill identifies clear cross-language gaps:

| Capability | R | Python | Used In |
|---|---|---|---|
| Conley SEs (integrated with FE) | `fixest::conley()` | No polished package | `robustness_conley_se.R` |
| Few-cluster inference (CR2) | `clubSandwich` | Not available | `robustness_wild_bootstrap.R` |
| Wild cluster bootstrap | Manual (fixest + base R) | No equivalent | `robustness_wild_bootstrap.R` |

### R Package Installation

```r
install.packages(c("fixest", "clubSandwich", "data.table"))
```

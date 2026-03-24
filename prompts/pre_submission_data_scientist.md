# Pre-Submission Audit: Data Scientist

You are a data scientist reviewing a research pipeline before journal submission. The paper studies how coal plant retirements propagate through technology networks to the stock prices of 703 listed power utilities worldwide.

## Your Task

Verify that every number in the manuscript (`manuscript/main.tex`) can be traced back to a specific pipeline output, and that no pipeline finding is misrepresented.

## Specific Checks

### 1. Number-to-Source Traceability

For EVERY number in the manuscript (coefficients, t-statistics, p-values, N, R², event counts, firm counts), verify:
- Which script in `src/` produces it
- Which output file in `results/metrics/` or `results/summaries/` contains it
- Whether the number in the manuscript EXACTLY matches the pipeline output (not rounded incorrectly, not from an older run)

Key numbers to verify:
- Abstract: 3.3 percentage points, 703 utilities, 80 countries, 175 events
- Table 1: All coefficients and SEs under three inference methods (FM, event-clustered, two-way)
- Table 2: Difference test statistics (t=5.703, sign test 82/117, Wilcoxon z=6.578, randomisation 0/999)
- Table 3: ETS interaction coefficients (-3.052, -4.264, +0.071)
- Table 4: ESG horse race (all four specifications)
- Section 4.2: Geographic channel t-stats across specifications (-0.41, -1.30, 1.68-1.80, 1.32, 0.69, 0.58, -1.04)
- Section 4.3: ETS portfolio sorts (-4.4% ETS vs +2.5% non-ETS, t=-5.11), geo x ETS placebo (t=-4.12)
- Section 4.5: Year tercile betas (+1.72, -3.58, -5.52), US vs non-US (0.11 vs -5.34)
- Section 4.6: Bandwidth range (-2.28 to -2.65), Conley t=-4.16, Romano-Wolf 0/9 with F=70.83, wild bootstrap p=0.091
- Identification section: Bartik t=-5.16, Oster δ*=20.8, pre-event balance t=-1.87, Rotemberg HHI=0.031

### 2. Pipeline Reproducibility

Run `python src/run_all.py` from scratch and verify:
- All scripts complete without errors
- Output files in `results/metrics/` match what the manuscript cites
- NOTE: `strategy2_bartik_shiftshare.py` is NOT in run_all.py — verify whether it should be added
- NOTE: `strategy2_conley_se.py` — check if it exists and whether it's in run_all.py

### 3. Figures

Verify that the three figures in `manuscript/figures/` match the data:
- `fig1_fuel_vs_geo.pdf`: Generated from `results/summaries/event_level_betas.csv` (117 rows, columns beta_fuel and beta_geo). Check that 82/117 points fall below the 45-degree line.
- `fig2_calendar_time.pdf`: Hardcoded from learning_alternatives.md tercile results. Verify the three beta values and SEs match.
- `fig3_world_map.pdf`: Generated from `data/derived/events/coal_retirement_events.csv` (filter is_first_mover=True) and `data/derived/networks/firm_centroids.csv`. Verify counts (344 events, 414 firms with GPS).

### 4. Stale Outputs

Check for stale outputs: files in `results/` that are no longer produced by any script in `src/`. These should be deleted before submission.

### 5. Data Integrity

Verify key data files:
- `data/derived/returns/monthly_returns.csv` and `daily_returns.csv`: row counts, date ranges, no obvious anomalies
- `data/derived/events/coal_retirement_events.csv`: 344+ first-mover events, lat/lon populated
- `data/derived/networks/weight_matrix_W.csv`: dimensions match 565 firms (if it exists)
- `data/derived/networks/firm_centroids.csv`: 414 firms with non-zero coordinates

## Output

Produce a table:

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|

Flag any mismatches, missing sources, or stale files.

## Constraints
- All scripts use stdlib-only Python (csv, math, collections, openpyxl). No pandas/numpy.
- Path resolution via `src/_paths.py`
- OLS implementation in `src/_ols.py`

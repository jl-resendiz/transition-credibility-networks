# Pre-Submission Data Architect Audit Report

**Date**: 2026-03-24
**Project**: When Coal Retires: The Propagation of Stranding Risk
**Pipeline root**: `c:\Users\jlres\vscode\TCB\transition-credibility-networks\`

---

## Deliverable 1: Logical Consistency Report (5 Claims)

### Claim 1: "Technology transmits, geography does not"

**Logical chain**: GEM plant records -> `parse_gem.py` -> `match_gem_compustat.py` -> `build_fuel_matrix.py` (fuel-mix similarity) + `build_weight_matrix.py` (geographic proximity) -> `strategy2_robust_inference.py` (FM betas) -> fuel significant, geo not.

**Verification**:
- The chain is complete and traceable. Weight matrices are pre-built and loaded from CSV files (`weight_matrix_W_fuel.csv`, `weight_matrix_W_geo.csv`, `weight_matrix_W_regulatory.csv`).
- FM + Newey-West: w_fuel t = -7.362, w_geo t = -1.297. Difference t = 5.703. The claim is well supported.
- The manuscript reports the result correctly: Table 1 shows FM fuel = -4.782 (SE 0.650), geo = -0.607 (SE 0.468). These match `strategy2_robust_inference.md` exactly.

**ISSUE (MEDIUM): Similarity metric mismatch with manuscript**. The manuscript states (Section 3.2, line 277): "I build a fuel-mix vector for each firm... and measure similarity as the cosine of the angle between vectors, following Hoberg and Phillips (2016)." However, `build_fuel_matrix.py` line 36-41 implements **L1 (Manhattan) similarity**: `sim = 1 - 0.5 * sum(|a-b|)`, NOT cosine similarity. These metrics are mathematically different. On normalised vectors (which fuel shares are not, unless total = 1), they produce different orderings. The manuscript must either be corrected to describe L1 similarity, or the code must be changed to implement cosine similarity.

### Claim 2: "ETS amplifies the fuel channel"

**Logical chain**: World Bank Carbon Pricing Dashboard -> `build_ets_matrix.py` -> `strategy2_credibility_interaction.py` -> ETS interaction on fuel.

**Verification**:
- The output file `strategy2_credibility_interaction.md` contains the full result including the placebo.
- **Honesty**: The manuscript is transparent. Section 4.3 (line 437) explicitly acknowledges three qualifications: (1) FM sign reversal, (2) geo x ETS placebo fails (t = -4.12, p = 0.0000), (3) extensive margin only.
- The placebo result IS in the output file (Specification 3), confirming pipeline honesty.
- The carbon price interaction has the WRONG SIGN (+0.071, t = 1.98), which the manuscript acknowledges.

**PASS**: The pipeline is honest about the fragility of this result.

### Claim 3: "ESG and fuel-mix measure different things"

**Logical chain**: LSEG Eikon ESG scores -> `strategy2_esg_horse_race.py` -> horse race regressions.

**Verification**:
- Output confirms N = 14,731 observations, 165 events, 153 firms (vs. 565 in full sample, 703 matched). The 153 firms are 21.8% of the full sample.
- **ISSUE (HIGH): No selection check**. The script does not test whether the 153 ESG-covered firms are representative of the 703-firm sample. There is no balance test on observables (fuel shares, region, size) between ESG-covered and non-covered firms. The manuscript (line 468) acknowledges the coverage difference but treats it as a practical point about data availability rather than a potential selection bias. A referee will ask whether the horse race result holds only because ESG-covered firms are systematically different.
- The manuscript correctly reports: ESG alone R2 = 0.012 vs. spatial alone R2 = 0.003; both survive FM horse race; under pooled OLS, fuel loses significance when ESG included.

### Claim 4: "The fuel signal strengthens over calendar time"

**Logical chain**: `strategy2_learning_alternatives.py` -> year tercile split, US vs non-US split.

**Verification**:
- The manuscript (line 473) reports: year tercile betas of +1.72, -3.58, -5.52 with t-stats 0.75, -3.02, -3.58. These match the output file.
- **Independence**: The year tercile and US/non-US splits use partially overlapping observations (non-US events span all terciles). However, these are presented as separate descriptive cuts, not as independent tests that multiply p-values, so this is acceptable practice.
- The continuous interaction w_fuel x year is reported (t = -2.34, p = 0.019), providing proper inference.
- The within-jurisdiction cascade test is honestly reported as weak (t = -1.57, p = 0.117).

**PASS**: Tests are honestly reported with appropriate caveats.

### Claim 5: "Fuel-mix shares are pre-determined (shift-share identification)"

**Logical chain**: GEM plant records (pre-2014 vintage) -> `strategy2_bartik_shiftshare.py` -> Bartik instrument.

**Verification**:
- The Bartik result is significant: pooled t = -5.16, FM t = -2.32. Rotemberg weights all non-negative (0/40 negative). Oster bound delta* = 20.8.
- Pre-event balance test: Bartik t = -1.87, p = 0.062. The manuscript reports this as "significant at 10 percent but not at 5 percent" (line 315). This is borderline and the manuscript is appropriately cautious.
- **ISSUE (MEDIUM): No cutoff sensitivity**. The PRE_CUTOFF = 2014 is hardcoded. The script does not report results for alternative cutoffs (2012, 2013, 2015). This is a researcher degree of freedom that a referee may challenge. The manuscript does not discuss sensitivity to this choice.
- **ISSUE (LOW): Very few FM events**. The Bartik FM specification uses only 3 valid events (line 53 of output). With T = 3, the Newey-West SE is barely identified (auto-lag selection yields max_lag = max(1, int(4*(3/100)^(2/9))) = 1). The t-statistics from FM with 3 events should be interpreted with extreme caution. The pooled specification (N = 24,070) is more reliable.

---

## Deliverable 2: Cross-Script Consistency Matrix

### CAR Computation

All strategy2_* scripts implement `compute_monthly_car()` with identical logic:
- Window: [-1, +POST_MONTHS] = [-1, +3] months
- Pre-event demeaning: 24 months of abnormal returns
- Market return: Fama-French value-weighted (MktRF + RF) from `F-F_Research_Data_Factors.csv`
- Minimum pre-event data: 12 months
- Event month identification: first month >= event_month in firm's return series

**EXCEPTION**: `strategy2_referee_tables.py` (line 78) requires `event_month in months` (exact match), while all other scripts use `m >= event_month` (first available month at or after). This is a subtle difference: if the firm has no return for the exact event month but has the following month, the referee_tables script drops the observation while other scripts keep it. This could create sample size discrepancies between the appendix tables and the main results.

### Weight Matrices

All scripts load from the same files:
- `derived/networks/weight_matrix_W_geo.csv`
- `derived/networks/weight_matrix_W_fuel.csv`
- `derived/networks/weight_matrix_W_regulatory.csv`

No script constructs weights internally (except `strategy2_bartik_shiftshare.py` which constructs pre-period weights from raw GEM data, and `strategy2_referee_tables.py` / `strategy2_bandwidth_fmb.py` which rebuild geo weights with alternative bandwidths).

### Sample Restrictions

| Script | MIN_OBS_PER_EVENT | Consistency |
|---|---|---|
| strategy2_robust_inference.py | 20 | OK |
| strategy2_credibility_interaction.py | (pooled, no FM filter) | Uses 175/179 events |
| strategy2_esg_horse_race.py | 20 (implied) | ESG subsample only |
| strategy2_learning_alternatives.py | 20 | OK |
| strategy2_bartik_shiftshare.py | 20 | OK |
| strategy2_difference_test_summary.py | 20 | OK |
| strategy2_joint_tests.py | (pooled) | Uses all 175 events |
| strategy2_romano_wolf.py | (pooled) | Uses all events |
| strategy2_bandwidth_fmb.py | 20 | OK |
| strategy2_event_specific_geo.py | 20 | OK |
| strategy2_geo_diversification.py | 20 | OK |
| strategy2_firm_level_test.py | 20 | OK |
| strategy2_esg_ets_fmb.py | 20 | OK |
| strategy2_spatial_score.py | (mixed) | OK |

The pooled scripts (joint_tests, romano_wolf, credibility_interaction) use all events with valid observations (175 of 179), while FM scripts filter to >= 20 firms per event (117 events). This is consistent and correct: pooled OLS does not require a minimum cross-section size per event, while FM does.

### Clustering Scale Factor

All scripts that implement clustered SEs use the same Cameron-Gelbach-Miller (2011) scale factor:
```
scale = (G / (G - 1)) * ((n - 1) / (n - k))
```
This is the standard small-sample correction for one-way clustering. The formula is consistent across all 12 scripts that implement it.

### Market Return

All scripts load the Fama-French value-weighted return as `(MktRF + RF) / 100` from the same file. Consistent.

---

## Deliverable 3: Sample Flow

```
GEM Plant Trackers (4 trackers, ~40,000 units)
    |
    v
parse_gem.py -> parsed plant records by technology
    |
    v
match_gem_compustat.py -> gem_compustat_matches.csv
    |                      703 firms matched (gvkeys)
    |
    +---> build_weight_matrix.py -> W_geo: ~565 firms with GPS centroids
    |                               (138 firms lack GPS data, dropped)
    |
    +---> build_fuel_matrix.py -> W_fuel: firms from W_geo with fuel shares
    |                             (~405 with fuel data, edges on W_geo structure)
    |
    +---> build_ets_matrix.py -> W_reg: ETS co-membership binary matrix
    |
    +---> build_retirement_events.py -> coal_retirement_events.csv
    |         ~1,844 total retirements
    |         -> 344 first-mover events (is_first_mover = True)
    |         -> 179 with matched_gvkeys (matched to Compustat firms)
    |         -> 175 with >= 20 firm observations per event
    |         -> 117 for FM (>= 20 firms AND sufficient cross-sectional variation)
    |
    v
compute_returns.py -> monthly_returns.csv, daily_returns.csv
    |                  703 firms, CRSP (US) + Eikon (non-US) + Compustat fallback
    |
    v
Analysis scripts:
    Full sample:    175 events x ~318 firms/event = 55,580 obs (pooled)
    FM sample:      117 events x ~244 firms/event = ~28,600 obs
    ESG subsample:  165 events x ~89 firms/event = 14,731 obs (153 unique firms)
    Bartik sample:  40 events with pre-period shares = 24,070 obs
    Bartik FM:      3 events (!)
```

**ISSUE (MEDIUM): The 344 -> 179 step is undocumented**. The manuscript says "175 first-mover retirement events" but does not explain how 344 first-movers become 179 matched events. The difference (165 events) represents events that could not be matched to any Compustat firm. The manuscript should document this step.

**ISSUE (HIGH): The Bartik FM uses only 3 events**. The manuscript cites "t = -2.32 Fama-MacBeth" (line 315) for the Bartik specification. With T = 3, the Newey-West standard error is unreliable. The output file confirms: "Valid events: 3". The manuscript does not disclose this. The pooled result (t = -5.16, N = 24,070) is reliable, but the FM claim requires a caveat.

---

## Deliverable 4: Architectural Issues (Ranked by Severity)

### HIGH SEVERITY

1. **Similarity metric mismatch**: The manuscript claims cosine similarity (Hoberg-Phillips 2016) but the code implements L1 (Manhattan) similarity. These are different metrics. A referee who checks the code will flag this inconsistency. **Fix**: Either change the manuscript text or change the code. If L1 is the correct choice, cite Czekanowski (1909) instead.

2. **No ESG selection check**: The 153-firm ESG subsample is not tested for representativeness. ESG coverage is endogenous to firm size, listing exchange, and geography. If ESG-covered firms are systematically larger or more geographically concentrated, the horse race is biased. **Fix**: Add a balance test comparing ESG-covered vs. non-covered firms on key observables (fuel mix, market cap, region).

3. **Bartik FM with T = 3**: The Fama-MacBeth Bartik t-statistic is based on 3 events. The manuscript cites this number without disclosing the event count. **Fix**: Either disclose the event count in the text or drop the FM Bartik result and rely on the pooled result only.

### MEDIUM SEVERITY

4. **No Bartik cutoff sensitivity**: The pre-2014 cutoff is a researcher degree of freedom. No sensitivity analysis for alternative cutoffs (2012, 2013, 2015) is reported. **Fix**: Run the Bartik specification with 2-3 alternative cutoffs and report in a footnote or appendix.

5. **CAR function inconsistency in referee_tables.py**: The `compute_monthly_car` function in `strategy2_referee_tables.py` requires exact month match (`event_month in months`), while all other scripts use first-available-month matching (`m >= event_month`). This could create small sample differences between appendix tables and main tables. **Fix**: Align the referee_tables function with the other scripts.

6. **344 -> 179 event funnel undocumented**: The manuscript does not explain how 344 first-mover events become 179 matched events. **Fix**: Add one sentence to the data section.

7. **Massive code duplication**: The `invert_matrix`, `mat_mul`, `ols_simple`, `_normal_cdf`, `newey_west_se`, `compute_monthly_car`, and `load_ff_factors_monthly` functions are copied verbatim across 15+ scripts. There is no shared `_ols.py` module in the repository (the `_paths.py` path helper exists, but no OLS module). Any bug fix in one copy must be manually propagated to all others. **Fix**: Extract common functions into a shared module.

8. **README.md references deleted strategies**: `results/README.md` line 18 still mentions `strategy1_panel_metrics.md`, which no longer exists. **Fix**: Update the README.

### LOW SEVERITY

9. **Appendix tables use different sample size**: The manuscript appendix note (line 559) states: "The appendix tables use the full event-firm panel (N ~ 72,600), which includes events with fewer than 20 firm observations." This is documented but creates a subtle inconsistency: the main tables (N = 55,580) and appendix tables (N ~ 72,600) use different samples. A referee may find this confusing.

10. **Romano-Wolf: 0/9 survive at 5%**: The manuscript (line 501) explains this carefully, but the headline result is that no individual channel survives family-wise correction. The joint F-test (F = 70.83) provides reassurance, but the multiple-testing section weakens the paper's statistical claims.

11. **No `_ols.py` shared module exists**: The prompt mentions `_ols.py` as a shared module, but it does not exist in the repository. Each script implements its own OLS. This is a documentation error in the project metadata.

12. **Git status shows many deleted files**: The working tree has ~90 deleted files from old strategies. These should be committed or restored before submission to avoid confusion in the replication package.

---

## Deliverable 5: Replication Package Recommendations

### Currently in place

| Component | Status |
|---|---|
| Master script (`run_all.py`) | EXISTS. Runs build + analysis in dependency order. |
| Path resolution (`_paths.py`) | EXISTS. Correctly anchored to repo root. |
| Raw data structure | EXISTS. Organised under `data/raw/` by source. |
| Derived data structure | EXISTS. Organised under `data/derived/` by type. |
| Output structure | EXISTS. `results/metrics/` (MD), `results/tables/` (LaTeX). |
| LaTeX compilation | EXISTS. Via `\input{../results/tables/...}` references. |
| Figure files | EXIST. All 3 referenced PDFs present in `manuscript/figures/`. |
| Table files | EXIST. All 5 `\input` tables present in `results/tables/`. |

### Missing

| Component | Recommendation |
|---|---|
| **README for replication** | Create a top-level `REPLICATION.md` with: (a) system requirements (Python 3.x, openpyxl); (b) data acquisition instructions for proprietary sources (Compustat, Eikon, CRSP); (c) `python src/run_all.py` as the single command; (d) expected runtime. |
| **Data manifest** | Create a file listing every input file, its source, and whether it is proprietary or public. The Fama-French factors are freely downloadable; GEM trackers require registration; Compustat/Eikon/CRSP require licenses. |
| **Version information** | Record Python version (3.x), openpyxl version. No other dependencies (stdlib-only). |
| **Estimated runtime** | The pipeline runs 15+ analysis scripts, each doing O(n^2) matrix operations on ~55,000 observations. Estimate and document. |
| **Proprietary data placeholders** | Provide empty CSV files with headers for Compustat, Eikon, and CRSP data, so a replicator can understand the expected format. |
| **Seed documentation** | Several scripts use different random seeds (42 in some, 20260222 in referee_tables, `hashlib.md5` in others). Document that results are seed-dependent and specify the seeds used. |
| **Clean git state** | Commit or discard the ~90 deleted/modified files shown in git status. The working tree should be clean before generating the replication package. |

### Manuscript-Pipeline Alignment Check

| Manuscript number | Source | Match? |
|---|---|---|
| Fuel FM beta = -4.782 | `strategy2_robust_inference.md` line 21 | YES |
| Fuel FM t = -7.362 | `strategy2_robust_inference.md` line 21 | YES |
| Pooled fuel = -5.474 | `strategy2_joint_tests.md` line 40 | YES |
| ETS interaction = -4.264 | `strategy2_credibility_interaction.md` line 21 | YES |
| ESG score beta = -0.114 | `strategy2_esg_horse_race.md` line 16 | YES |
| 117 FM events | `strategy2_robust_inference.md` line 13 | YES |
| 175 events | `strategy2_joint_tests.md` line 5 | YES |
| N = 55,580 | `strategy2_joint_tests.md` line 5 | YES |
| Difference t = 5.703 | `strategy2_difference_test_summary.md` line 28 | YES |
| Sign test 82/117 | `strategy2_difference_test_summary.md` line 35 | YES |
| Romano-Wolf 0/9 at 5% | `strategy2_romano_wolf.md` line 44 | YES |
| F = 70.83 | `strategy2_joint_tests.md` line 15 | YES |
| Bartik t = -5.16 | `strategy2_bartik_shiftshare.md` line 90 | YES |
| Oster delta* = 20.8 | `strategy2_bartik_shiftshare.md` line 176 | YES |
| STS Q5-Q1 out-of-sample = -2.76% | `strategy2_spatial_score.md` (to verify) | PARTIAL |
| Weight correlations < 0.16 | `table_weight_correlations.tex` max = 0.138 | YES |
| Conclusion: fuel FM t = -7.36 | `strategy2_robust_inference.md` (t = -7.362) | YES (rounding) |

**ISSUE**: The manuscript conclusion (line 542) says "t-statistic of -7.36" while the main results section (line 361) says "t = 2.75" for w_reg. The robust_inference output reports w_reg t = 2.749. This rounds to 2.75, consistent.

All hardcoded numbers checked match their pipeline sources. No orphaned references to deleted strategies were found in the manuscript.

---

## Summary of Findings

The pipeline is architecturally sound. The `run_all.py` orchestrator correctly orders dependencies. The CAR computation is consistent across all analysis scripts (with one minor exception in `referee_tables.py`). Weight matrices are loaded from shared files. Sample restrictions are consistently applied. The manuscript numbers match pipeline outputs.

Three issues require attention before submission:

1. **Fix the similarity metric description** (cosine vs. L1) in the manuscript or code.
2. **Add an ESG selection balance test** to address the 153/703 coverage gap.
3. **Disclose the T = 3 event count** for the Bartik FM result, or drop that particular statistic.

The replication package needs a README, data manifest, and clean git state. These are straightforward to produce from the existing structure.

# Pre-Submission Data Scientist Audit Report

**Project**: "When Coal Retires: The Propagation of Stranding Risk"
**Audit date**: 2026-03-24
**Auditor**: Automated pipeline verification (Claude Code)

---

## 1. Number-to-Source Traceability

### Abstract

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Abstract | 3.3 percentage points | `results/metrics/strategy2_joint_tests.md` (w_fuel = -5.474, x median w_fuel ~0.006) | -5.474 x 0.006 = -3.28 pp, rounded to 3.3 | YES |
| Abstract | 703 utilities | `data/derived/returns/monthly_returns.csv` (unique gvkeys) | 703 | YES |
| Abstract | 80 countries | `data/derived/fundamentals/firm_fundamentals.csv` (unique fic for 703 return-sample firms) | 80 | YES |
| Abstract | 175 events | Summary: 179 matched first-movers, 175 with >= 20 firms | 175 event clusters | MINOR: see note 1 |

### Table 1: Channel Decomposition

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Table 1 FM w_geo | -0.607 (0.468) | `strategy2_robust_inference.md` | -0.607185 (0.468053) | YES |
| Table 1 FM w_fuel | -4.782 (0.650) | `strategy2_robust_inference.md` | -4.782353 (0.649637) | YES |
| Table 1 FM w_reg | +2.643 (0.961) | `strategy2_robust_inference.md` | +2.642804 (0.961429) | YES |
| Table 1 FM same_sector | +0.021 (0.011) | `strategy2_robust_inference.md` | +0.021454 (0.011209) | YES |
| Table 1 EC w_geo | -0.023 (0.057) | `strategy2_joint_tests.md` | -0.022960 (0.056598) | YES |
| Table 1 EC w_fuel | -5.474 (0.730) | `strategy2_joint_tests.md` | -5.474254 (0.730054) | YES |
| Table 1 EC w_reg | +1.453 (1.052) | `strategy2_joint_tests.md` | +1.452522 (1.051538) | YES |
| Table 1 EC same_sector | +0.033 (0.009) | `strategy2_joint_tests.md` | +0.033201 (0.008872) | YES |
| Table 1 TW w_geo | -0.023 (0.119) | `strategy2_firm_level_test.md` Panel B | -0.022960 (0.119301) | YES |
| Table 1 TW w_fuel | -5.474 (1.271) | `strategy2_firm_level_test.md` Panel B | -5.474254 (1.271256) | YES |
| Table 1 TW w_reg | +1.453 (1.134) | `strategy2_firm_level_test.md` Panel B | +1.452522 (1.133554) | YES |
| Table 1 TW same_sector | +0.033 (0.011) | `strategy2_firm_level_test.md` Panel B | +0.033201 (0.011264) | YES |
| Table 1 FM diff | +4.175 (0.732) | `strategy2_robust_inference.md` | +4.175168 (0.732045) | YES |
| Table 1 EC diff | +5.451 (0.731) | `strategy2_joint_tests.md` | +5.451294 (0.730684) | YES |
| Table 1 TW diff | +5.451 (1.268) | `strategy2_firm_level_test.md` Panel B | +5.451294 (1.268) | YES |
| Table 1 FM Events/N | 117 / ~28,600 | `strategy2_robust_inference.md` | 117 events, avg 244.6 firms (= ~28,618) | YES |
| Table 1 EC Events/N | 175 / 55,580 | `strategy2_joint_tests.md` | 175 / 55580 | YES |
| Table 1 FM Avg R2 | 0.052 | `strategy2_robust_inference.md` | 0.0515 | YES (rounds to 0.052) |
| Table 1 EC R2 | 0.007 | `strategy2_joint_tests.md` | 0.0071 | YES |

### Table 2: Difference Test Battery

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Table 2 FM+NW t | 5.703 | `strategy2_difference_test_summary.md` | 5.703 | YES |
| Table 2 Sign test | 82/117 (70.1%) | `strategy2_difference_test_summary.md` | 82/117 (70.1%) | YES |
| Table 2 Wilcoxon z | 6.578 | `strategy2_difference_test_summary.md` | 6.578 | YES |
| Table 2 Randomisation | 0/999 | `strategy2_difference_test_summary.md` | 0/999 | YES |
| Table 2 p-values | 0.000, 0.000, 0.000, 0.002/0.001 | `strategy2_difference_test_summary.md` | 0.0000, 0.0000, 0.0000, 0.0020/0.0010 | YES |

### Table 3: ETS Interaction

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Table 3 w_fuel base | -3.052 (0.602) | `strategy2_credibility_interaction.md` Spec 1 | -3.051701 (0.601992) | YES |
| Table 3 w_fuel x ETS | -4.264 (1.405) | `strategy2_credibility_interaction.md` Spec 1 | -4.263692 (1.404737) | YES |
| Table 3 Carbon price | +0.071 (0.036) | `strategy2_credibility_interaction.md` Spec 2 | +0.071476 (0.036120) | YES |

### Table 4: ESG Horse Race

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Table 4 (1) ESG | -0.114 (0.014) | `strategy2_esg_horse_race.md` Model 1 | -0.114409 (0.013942) | YES |
| Table 4 (2) w_fuel | -1.559 (0.866) | `strategy2_esg_horse_race.md` Model 2 | -1.559102 (0.866465) | YES |
| Table 4 (3) ESG | -0.118 (0.014) | `strategy2_esg_horse_race.md` Model 3 | -0.117597 (0.014368) | YES |
| Table 4 (3) w_fuel | -1.135 (0.854) | `strategy2_esg_horse_race.md` Model 3 | -1.134982 (0.853913) | YES |
| Table 4 (4) ESG | -0.116 (0.015) | `strategy2_esg_horse_race.md` Model 4 | -0.116359 (0.014853) | YES |
| Table 4 (4) w_fuel | +0.397 (0.923) | `strategy2_esg_horse_race.md` Model 4 | +0.397032 (0.928586) | MINOR: SE 0.923 vs 0.929 |
| Table 4 (4) w_fuel x ETS | -3.489 (1.685) | `strategy2_esg_horse_race.md` Model 4 | -3.488802 (1.682291) | YES |
| Table 4 R2 values | 0.012, 0.003, 0.016, 0.017 | Source | 0.0119, 0.0034, 0.0159, 0.0173 | YES |
| Table 4 N | 14,731 | `strategy2_esg_horse_race.md` | 14731 | YES |
| Table 4 footnote FM+NW ESG t | -5.72 | `strategy2_esg_ets_fmb.md` Spec 3 | -5.72 | YES |
| Table 4 footnote FM+NW fuel t | -2.17 | `strategy2_esg_ets_fmb.md` Spec 3 | -2.17 | YES |

### Section 4.2: Geographic Channel

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 4.2 pooled t | -0.41 | `strategy2_geo_diversification.md` full sample OLS | -0.41 | YES |
| Sec 4.2 FM t | -1.30 | `strategy2_geo_diversification.md` full sample FM | -1.30 | YES |
| Sec 4.2 bandwidth FM t | 1.68 to 1.80 | `strategy2_bandwidth_fmb.md` | 1.680 to 1.798 | YES |
| Sec 4.2 event-specific pooled t | 1.32 | `strategy2_event_specific_geo.md` Spec B pooled | 1.320 | YES |
| Sec 4.2 event-specific FM t | 0.69 | `strategy2_event_specific_geo.md` Spec B FM | 0.685 | YES |
| Sec 4.2 single-country pooled t | 0.58 | `strategy2_geo_diversification.md` single-country OLS | 0.58 | YES |
| Sec 4.2 single-country FM t | -1.04 | `strategy2_geo_diversification.md` single-country FM | -1.04 | YES |

### Section 4.3: ETS Amplification

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 4.3 ETS penalty sum | -7.316 | Derived: -3.052 + -4.264 | -7.316 | YES |
| Sec 4.3 ETS Q5-Q1 | -4.4% | `strategy2_credibility_interaction.md` portfolio | -0.0438 = -4.38% | YES |
| Sec 4.3 non-ETS Q5-Q1 | +2.5% | `strategy2_credibility_interaction.md` portfolio | +0.0247 = +2.47% | YES |
| Sec 4.3 paired diff t | -5.11 | `strategy2_credibility_interaction.md` portfolio | -5.112 | YES |
| Sec 4.3 carbon price t | 1.98 | `strategy2_credibility_interaction.md` Spec 2 | 1.979 | YES |
| Sec 4.3 FM fuel t (ETS panel) | -2.60 | `strategy2_esg_ets_fmb.md` Panel B Spec 1 | -2.60 | YES |
| Sec 4.3 geo x ETS placebo t | -4.12 | `strategy2_credibility_interaction.md` Spec 3 | -4.118 | YES |

### Section 4.5: Calendar Time / Learning

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 4.5 year interaction t | -2.34 | `strategy2_learning_alternatives.md` D3 | -2.337 | YES |
| Sec 4.5 T1 beta | +1.72 (t=0.75) | `strategy2_learning_alternatives.md` D4 | +1.717125 (t=0.749) | YES |
| Sec 4.5 T2 beta | -3.58 (t=-3.02) | `strategy2_learning_alternatives.md` D4 | -3.577552 (t=-3.019) | YES |
| Sec 4.5 T3 beta | -5.52 (t=-3.58) | `strategy2_learning_alternatives.md` D4 | -5.518344 (t=-3.580) | YES |
| Sec 4.5 US t | 0.11 (91 events) | `strategy2_learning_alternatives.md` A1 | t=0.113 (N=91) | YES |
| Sec 4.5 non-US beta | -5.34 (t=-4.10) | `strategy2_learning_alternatives.md` A1 | -5.342905 (t=-4.098) | YES |
| Sec 4.5 non-US events | 84 | `strategy2_learning_alternatives.md` A1 | N=84 | YES |
| Sec 4.5 log_order t | -1.57 | `strategy2_learning_alternatives.md` summary | -1.569 | YES |
| Sec 4.5 Welch t | -2.48 | `strategy2_learning_alternatives.md` E1 | -2.482 | YES |

### Section 4.6: Robustness

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 4.6 bandwidth fuel t range | -2.28 to -2.65 | `strategy2_bandwidth_fmb.md` | -2.279 to -2.648 | YES |
| Sec 4.6 Conley 1000km fuel t | -4.16 | `strategy2_conley_se.md` | -4.16 | YES |
| Sec 4.6 Romano-Wolf 0/9 | 0/9 | `strategy2_romano_wolf.md` | 0/9 | YES |
| Sec 4.6 raw t=-7.50 | -7.50 | `strategy2_romano_wolf.md` | -7.498 | YES |
| Sec 4.6 F=70.83 | 70.83 | `strategy2_joint_tests.md` | 70.8335 | YES |
| Sec 4.6 wild bootstrap p | 0.091 | `strategy3_phaseout_wild_bootstrap.md` | 0.0910 | YES |

### Identification Section

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 3.5 Bartik pooled t | -5.16 | `strategy2_bartik_shiftshare.md` pooled | -5.160 | YES |
| Sec 3.5 Bartik FM t | -2.32 | `strategy2_bartik_shiftshare.md` FM | -2.317 | YES |
| Sec 3.5 Rotemberg neg weights | 0/40 | `strategy2_bartik_shiftshare.md` | 0/40 | YES |
| Sec 3.5 Rotemberg HHI | 0.031 | `strategy2_bartik_shiftshare.md` | 0.0311 | YES |
| Sec 3.5 pre-event t | -1.87 | `strategy2_bartik_shiftshare.md` | -1.870 | YES |
| Sec 3.5 Oster delta* | 20.8 | `strategy2_bartik_shiftshare.md` | 20.7993 | YES |
| Sec 3.5 standard Oster delta* | 61.6 | `strategy2_bartik_shiftshare.md` | 61.6174 | YES |

### Conclusion

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Conclusion FM t | -7.36 | `strategy2_robust_inference.md` | -7.362 | YES |
| Conclusion F | 70.83 | `strategy2_joint_tests.md` | 70.8335 | YES |
| Conclusion delta* | 20.8 | `strategy2_bartik_shiftshare.md` | 20.7993 | YES |
| Conclusion 105 of 703 | 105 | Not verified from pipeline output | -- | NOT TRACED |

### Other In-Text Numbers

| Manuscript Location | Number Cited | Source File | Source Value | Match? |
|---|---|---|---|---|
| Sec 3.1 565 firms | 565 | `strategy2_firm_level_test.md` | 565 firms | YES |
| Sec 3.1 153 ESG firms | 153 | `strategy2_esg_horse_race.md` | 153 unique firms | YES |
| Sec 3.1 117 FM events | 117 | `strategy2_robust_inference.md` | 117 | YES |
| Sec 3.1 55,580 obs | 55,580 | `strategy2_joint_tests.md` | 55580 | YES |
| Sec 3.1 35 countries | 35 | `results/summaries/summary_statistics.md` | 35 | YES |
| Sec 3.3 weight corr < 0.16 | < 0.16 | `results/tables/table_weight_correlations.tex` | max = 0.138 | YES |
| Sec 4.1 t=2.75 (w_reg FM) | 2.75 | `strategy2_robust_inference.md` | 2.749 | YES |
| Sec 4.1 pooled R2 = 0.007 | 0.007 | `strategy2_joint_tests.md` | 0.0071 | YES |
| Sec 4.1 FM avg R2 = 0.052 | 0.052 | `strategy2_robust_inference.md` | 0.0515 | YES |
| Sec 4.1 firm-level R2 = 0.026 | 0.026 | `strategy2_firm_level_test.md` | 0.025677 | YES |
| Sec 4.1 3.6x increase | 3.6x | `strategy2_firm_level_test.md` | 3.6x | YES |
| App STS Q5-Q1 = -2.76% (t=-2.98) | -2.76%, -2.98 | `strategy2_spatial_score.md` Post-2020 | -0.0276, t=-2.976 | YES |
| App STS in-sample t | 1.58 | `strategy2_spatial_score.md` | t=1.584 | YES |

---

## 2. Pipeline Reproducibility

### Scripts in `run_all.py`

The orchestrator at `src/run_all.py` includes 15 analysis scripts covering all main results. Two scripts that produce results cited in the manuscript are **NOT** in `run_all.py`:

| Script | Output | In run_all? | Cited in Manuscript? |
|---|---|---|---|
| `strategy2_bartik_shiftshare.py` | `strategy2_bartik_shiftshare.md` | **NO** | YES (Sec 3.5: Bartik t, Oster, Rotemberg, pre-event balance) |
| `robustness_conley_se.R` | `strategy2_conley_se.md` | **NO** | YES (Sec 4.6: Conley t=-4.16) |

**Finding**: `strategy2_bartik_shiftshare.py` exists as a Python script but is not called by `run_all.py`. `strategy2_conley_se.py` does **not exist** as a Python file; only `robustness_conley_se.R` exists. The Conley SE results were produced by an R script that cannot be run by the Python orchestrator.

**Recommendation**: Add `strategy2_bartik_shiftshare.py` to `run_all.py`. Document that the Conley SE result requires running the R script separately, or port it to Python.

### Pipeline Execution

Pipeline was not re-executed during this audit (would require Refinitiv API credentials and substantial runtime). Verification was performed by cross-referencing manuscript numbers against existing pipeline output files.

---

## 3. Figures

| Figure | Source Data | Verification | Status |
|---|---|---|---|
| `fig1_fuel_vs_geo.pdf` | `results/summaries/event_level_betas.csv` | 117 rows confirmed. 82/117 events have beta_fuel < beta_geo (below 45-degree line). Matches manuscript claim of 82/117 = 70%. | PASS |
| `fig2_calendar_time.pdf` | `strategy2_learning_alternatives.md` D4 tercile results | T1: +1.72, T2: -3.58, T3: -5.52, all match. | PASS |
| `fig3_world_map.pdf` | `data/derived/events/coal_retirement_events.csv` + `data/derived/networks/firm_centroids.csv` | 344 first-mover events confirmed. 414 firms with non-zero GPS coordinates confirmed. | PASS |

Note: Figure generation scripts are R files in `results/figures/` (`generate_fig1.R`, `generate_fig2.R`, `generate_fig3.R`). These are not called by the Python `run_all.py` orchestrator.

---

## 4. Stale Outputs

All files currently in `results/` are produced by scripts in `src/` or are supporting files:

| File | Producer | Status |
|---|---|---|
| `results/metrics/strategy2_bartik_shiftshare.md` | `strategy2_bartik_shiftshare.py` | Active but not in run_all.py |
| `results/metrics/strategy2_conley_se.md` | `robustness_conley_se.R` | Active but R-only |
| `results/README.md` | Manual | OK |
| `results/figure_specs.md` | Manual | OK |
| All other `results/metrics/*.md` | Corresponding `src/strategy2_*.py` or `src/strategy3_*.py` | Active |
| All `results/tables/*.tex` | `src/strategy2_referee_tables.py` | Active |
| `results/summaries/event_level_betas.csv` | `src/strategy2_difference_test_summary.py` | Active |
| `results/summaries/summary_statistics.md` | `src/summary_statistics.py` | Active |
| `results/figures/*.pdf` | R scripts in `results/figures/` | Active |

**No stale outputs detected.** All result files trace to an active script.

Note: The git status shows many **deleted** files in `results/` (e.g., `strategy2_panel_did_metrics_*.md`, `strategy3_policy_metrics_*.md`). These have already been removed from the working tree, which is correct -- they were from earlier pipeline versions.

---

## 5. Data Integrity

### Returns Files

| File | Rows | Date Range | Unique Firms | Status |
|---|---|---|---|---|
| `monthly_returns.csv` | 78,191 | 2009-06-26 to 2026-01-31 | 703 | OK |
| `daily_returns.csv` | 1,232,315 | 2010-01-03 to 2026-02-17 | -- | OK |

### Event File

| Check | Expected | Actual | Status |
|---|---|---|---|
| Total retirements | -- | 1,844 | OK |
| First-mover events | 344 | 344 | PASS |
| First-mover with lat/lon | 344 | 344 | PASS |
| First-mover matched to Compustat | 179 | 179 | PASS |
| Countries with first-mover events | 35 | 35 | PASS |

### Network Files

| Check | Expected | Actual | Status |
|---|---|---|---|
| `weight_matrix_W.csv` (single combined) | Prompt says 565 firms | File does not exist | N/A (superseded by separate layer files) |
| `weight_matrix_W_fuel.csv` | -- | 85,424 edges, sparse format | OK |
| `weight_matrix_W_geo.csv` | -- | 170,906 edges, sparse format | OK |
| `weight_matrix_W_regulatory.csv` | -- | 14,844 edges, sparse format | OK |
| Unique firms across all matrices | Prompt says 565 | 414 | SEE NOTE 2 |
| `firm_centroids.csv` | 414 firms with non-zero coords | 414 with non-zero coords | PASS |

---

## 6. Issues Found

### Note 1: Event Count Ambiguity (MINOR)

The manuscript states "175 first-mover retirement events" in Section 3.1. The actual count is 179 first-mover events matched to Compustat, of which 175 produce >= 20 firms for stable OLS estimation. The pooled regressions use all 179 events (N = 55,580 from 175 event clusters after removing 4 tiny events), while Fama-MacBeth uses 117 events. The manuscript should clarify: "179 first-mover retirement events matched to Compustat, of which 175 contain at least 20 peer firms."

### Note 2: Firm Count Discrepancy (MODERATE)

The manuscript claims "565 have complete spatial weight data for all three network layers." However, the weight matrix files contain only 414 unique gvkeys. The 565 figure comes from the event-firm panel (`strategy2_firm_level_test.md`: "565 firms"), which counts firms that appear in the regression sample. The discrepancy likely arises because the regression assigns zero weights to firms not in the weight matrices. This is not an error per se, but the claim "complete spatial weight data" is imprecise. 414 firms have non-trivial weight data; the remaining 151 have zero weights on at least some layers.

### Note 3: Table 4 SE Rounding (TRIVIAL)

In Table 4 specification (4), the manuscript reports w_fuel SE as (0.923). The source value is 0.928586, which rounds to 0.929 at 3 decimal places, not 0.923. This appears to be a transcription error of 6 basis points in the standard error.

### Note 4: Scripts Missing from run_all.py (MODERATE)

Two scripts producing manuscript-cited results are not in `run_all.py`:
- `strategy2_bartik_shiftshare.py` -- produces Bartik, Oster, Rotemberg, and pre-event balance results cited in the identification section
- `robustness_conley_se.R` -- produces Conley spatial SE result; exists only as R, not Python

A reviewer running `python src/run_all.py` will not reproduce these results.

### Note 5: 105 of 703 Firms with Price Returns (NOT TRACED)

The conclusion mentions "105 of 703" firms using price returns without dividend adjustment. This number is not found in any pipeline output file. It likely comes from `src/compute_returns.py` diagnostic output but is not persisted to a results file.

### Note 6: Figure Generation Not Automated (MINOR)

The three manuscript figures are generated by R scripts (`results/figures/generate_fig*.R`), not by the Python pipeline. Running `python src/run_all.py` will not regenerate figures.

---

## 7. Summary Assessment

| Category | Verdict |
|---|---|
| Number accuracy | 99%+ of numbers match exactly (within rounding). One trivial SE transcription error (Note 3). |
| Pipeline completeness | Two analysis scripts missing from orchestrator (Note 4). Figures require separate R execution (Note 6). |
| Data integrity | All key data files present with expected row counts, date ranges, and GPS coverage. |
| Stale outputs | None detected. Prior stale files already deleted from working tree. |
| Reproducibility | High for Python scripts in `run_all.py`. Bartik and Conley results require manual execution. |

### Recommended Actions Before Submission

1. **Fix Table 4 SE**: Change w_fuel SE in specification (4) from (0.923) to (0.929).
2. **Add `strategy2_bartik_shiftshare.py` to `run_all.py`**.
3. **Clarify event count**: State "179 first-mover events (175 with >= 20 peer firms)" in Section 3.1.
4. **Clarify firm count**: State "414 firms with plant-level GPS data" rather than "565 with complete spatial weight data," or explain the distinction.
5. **Document R dependencies**: Note that Conley SEs and figures require R scripts run separately.
6. **Persist the "105 of 703" diagnostic** to a pipeline output file so the number is traceable.

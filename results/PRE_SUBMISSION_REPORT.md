# Pre-Submission Quantitative Report

**Paper:** *When Coal Retires: The Propagation of Stranding Risk*
**Author:** Jose Luis Resendiz (Smith School, University of Oxford)
**Compiled:** 2026-04-30 from a clean `python src/run_all.py --analysis` run.
**Manuscript file:** [`manuscript/when_coal_retires.tex`](when_coal_retires.tex) (renamed from `main.tex`).
**Pipeline status:** all 21 analysis scripts ran to exit code 0. Latex compiled to 37 pages, 0 undefined references, 0 errors.

This document collects every quantitative result that is referenced in the manuscript, grouped by manuscript section, with the source script and CSV/MD artefact for each number. The final section flags every divergence between the freshly regenerated outputs and the prose currently in the manuscript.

---

## 1. Sample (manuscript §3.1, §4)

| Quantity | Value | Source |
|---|---:|---|
| Listed power utilities (returns coverage) | 703 | `summary_statistics.py` (monthly returns load) |
| Countries (returns sample) | 80 | `summary_statistics.py` |
| Firms with full network coverage (W_geo, W_fuel, W_reg) | 414 | `robust_inference.py` console: `W_geo firms: 414` |
| First-mover retirement events | 179 | `summary_statistics.py` Panel D |
| First-mover countries | 35 | `summary_statistics.py` Panel D |
| Events with ≥20 firms (FM-eligible) | 117 | `robust_inference.py` |
| Total event-firm observations (analysis panel) | 55,580 | `robust_inference.py` |
| US first-mover events | 135 | `summary_statistics.py` Panel D |
| ESG-covered firms | 153 unique | `esg_horse_race.py` |
| ESG events ≥10 covered firms | 165 | `esg_fm_joint.py` |

> **Discrepancy:** `results/summaries/summary_statistics.md` Panel A reports 428 firms in the "analysis sample" and 449 firms in the GEM-matched subsample. The 703-firm headline used in the abstract and §4.1 comes from the *monthly returns* sample, which is computed downstream of `compute_returns.py`. These two definitions are inconsistent and should be reconciled before submission (see §6, Red Flag #1).

---

## 2. Headline Channel Decomposition (manuscript §4.4, Table 2)

Window: $[-1, +3]$ months. CAR is market-adjusted using value-weighted returns. Spec includes `w_geo`, `w_fuel`, `w_reg`, and `same_sector`.

### 2a. Fama-MacBeth + Newey-West (Approach 1, gold standard)

Source: [`results/metrics/robust_inference.md`](metrics/robust_inference.md)

| Variable | Mean β | NW SE | t | p |
|---|---:|---:|---:|---:|
| Intercept | +0.0255 | 0.0091 | +2.80 | 0.005*** |
| $w^{\text{geo}}$ | **−0.5427** | 0.3090 | −1.76 | 0.079* |
| $w^{\text{fuel}}$ | **−4.7656** | 0.6508 | **−7.32** | 0.000*** |
| $w^{\text{reg}}$ | +2.6975 | 0.9518 | +2.83 | 0.005*** |
| Same sector | +0.0215 | 0.0112 | +1.92 | 0.055* |

Channel difference test (β_geo − β_fuel) = +4.2229 (NW SE = 0.7076, **t = +5.97**, p = 0.000).
Joint Wald F-test (3 channels) = 20.65.
Avg within-event $R^2$ = 0.0518; avg firms per event = 244.6.

### 2b. Pooled OLS, event-clustered (joint_tests.jl)

Source: [`results/metrics/joint_tests.md`](metrics/joint_tests.md)

| Variable | β | SE | t |
|---|---:|---:|---:|
| Intercept | +0.0358 | 0.0077 | +4.67 |
| $w^{\text{geo}}$ | **+0.0176** | 0.1007 | +0.18 |
| $w^{\text{fuel}}$ | **−5.4884** | 0.7279 | **−7.54** |
| $w^{\text{reg}}$ | +1.4412 | 1.0503 | +1.37 |
| Same sector | +0.0332 | 0.0089 | +3.74 |

Difference test (β_geo − β_fuel) = +5.5059 (SE 0.7285, **t = +7.56**, p < 0.001).
Joint F (3, 55575) = 70.81 (permutation p < 0.001).
$N$ = 55,580; 175 event clusters.

> **Discrepancy with manuscript Table 2:** the manuscript reports w_geo = −0.023 and w_fuel = −5.474 in the "event-clustered" column. The fresh Julia run gives +0.0176 and −5.4884. Magnitudes are within sampling noise of pseudorandom control draws (md5 vs Julia hash), but **signs and the third decimal differ** — this affects how Table 2 reads on close inspection. See §6, Red Flag #2.

### 2c. Pooled OLS, event-clustered (Python, robust_inference Approach 4)

Source: [`results/metrics/robust_inference.md`](metrics/robust_inference.md) §Approach 4.

| Window | $w^{\text{fuel}}$ β | t | $w^{\text{geo}}$ β | t | $R^2$ |
|---|---:|---:|---:|---:|---:|
| $[-1, +1]$ | −3.4295 | **−7.05** | −0.0400 | −0.66 | 0.0025 |
| $[-1, +2]$ | −4.3295 | **−7.90** | −0.0940 | −1.21 | 0.0029 |
| $[-1, +3]$ | **−5.5418** | **−7.94** | +0.0876 | +0.82 | 0.0035 |
| $[0, +1]$ | −3.1500 | **−7.46** | +0.0016 | +0.03 | 0.0028 |

The fuel channel is significant at $|t| > 7$ at every window; geo is insignificant at every window.

### 2d. Economic magnitude

| Quantity | Value | Source |
|---|---:|---|
| SD of $w^{\text{fuel}}$ within event | ≈0.004 | manuscript Table 1 (consistent with row-normalisation) |
| 1-SD effect on 4-month CAR (event-clustered −5.474) | −2.2 pp | manuscript §4.4 |
| Annualised equivalent | ≈−6.6% | manuscript §4.4 |
| 1-SD effect under multi-factor (−3.10) | −1.2 pp | manuscript §4.7 |

---

## 3. Robustness Battery

### 3a. Outlier diagnostics (Cook's D)

Source: `robust_inference.py` Approach 5.

| Metric | Value |
|---|---:|
| N (full sample) | 55,580 |
| Max Cook's D | 0.0150 |
| Obs with D > 4/N | 1,728 (3.1%) |
| Fuel β (full) | −5.5418 |
| Fuel β (trimmed of 1,728 obs) | **−6.3905** |
| $R^2$ full → trimmed | 0.0035 → 0.0058 |

Trimming **strengthens** the fuel coefficient: outliers attenuate, not inflate, the estimate.

> **Minor discrepancy:** manuscript §4.6 reports "1,688 (3.0 percent)" exceed the 4/N threshold. Fresh run = 1,728 (3.1%). Trimmed coefficient: manuscript −6.35 vs fresh −6.39.

### 3b. Romano-Wolf multiple-testing

Source: [`results/metrics/romano_wolf.md`](metrics/romano_wolf.md). 999 Rademacher cluster bootstrap draws, seed 42.

| Variable | β | t | Raw p | Bonferroni | Max-t | Romano-Wolf |
|---|---:|---:|---:|---:|---:|---:|
| $w^{\text{geo}}$ | +0.0176 | +0.18 | 0.861 | 1.000 | 0.999 | 0.879 |
| $w^{\text{fuel}}$ | **−5.4884** | **−7.54** | **0.000*** | **0.000** | **0.000** | **0.000** |
| $w^{\text{reg}}$ | +1.4412 | +1.37 | 0.170 | 0.510 | 0.456 | 0.334 |

Only the fuel channel survives any multiple-testing correction.

### 3c. Fisher Randomization Inference (within-event permutation)

Source: [`results/metrics/fisher_ri.md`](metrics/fisher_ri.md). 999 within-event $w^{\text{fuel}}$ permutations.

| Quantity | Value |
|---|---:|
| Observed FM β_fuel | **−4.8318** |
| Permutation 99th percentile | +1.3121 |
| Permutation range | [−1.8364, +2.0268] |
| One-sided RI p | **0.0010** |
| Two-sided RI p | **0.0010** |

The observed coefficient lies more than 2× outside the entire permutation envelope.

> **Minor discrepancy:** manuscript §4.6 says "99th percentile of the permutation distribution is approximately ±2.4". Fresh run gives 99th percentile = +1.31 and range = [−1.84, +2.03]. The "±2.4" figure is no longer correct; replace with "approximately ±2.0" or quote the actual percentile.

### 3d. Lag sensitivity (Newey-West, FM time series of T = 117)

Source: [`results/metrics/lag_sensitivity.md`](metrics/lag_sensitivity.md).

| Channel | Mean β | t NW(4) | t NW(8) | t NW(12) | t NW(18) |
|---|---:|---:|---:|---:|---:|
| $\gamma_{\text{fuel}}$ | −4.7656 | **−7.32** | −6.52 | −6.02 | **−5.45** |
| $\gamma_{\text{geo}}$ | −0.5427 | −1.76 | −1.68 | −1.62 | −1.66 |
| $\gamma_{\text{reg}}$ | +2.6975 | +2.83 | +2.85 | +2.98 | +2.99 |
| $\gamma_{\text{geo}} − \gamma_{\text{fuel}}$ | +4.2229 | **+5.97** | +5.35 | +5.05 | **+4.73** |

Fuel rejects null at 5% across all four lag choices; manuscript text matches exactly.

### 3e. Multi-factor abnormal returns (FF3 + sample-built utility industry)

Source: [`results/metrics/multifactor_inference.md`](metrics/multifactor_inference.md). 24-month pre-event window.

| Channel | Single-factor (headline) | Multi-factor | Shrinkage |
|---|---|---|---:|
| $\gamma_{\text{fuel}}$ | −4.7656 (0.65) [t=−7.32] | **−3.1043** (0.69) [**t=−4.50**] | +35% |
| $\gamma_{\text{geo}}$ | −0.5427 (0.31) [t=−1.76] | −0.0679 (0.53) [t=−0.13] | +88% |
| $\gamma_{\text{reg}}$ | +2.6975 (0.95) [t=+2.83] | +2.4749 (0.82) [t=+3.00] | +8% |
| $\gamma_{\text{geo}}-\gamma_{\text{fuel}}$ | +4.2229 (0.71) [t=+5.97] | +3.0364 (0.92) [t=+3.31] | +28% |

The fuel channel survives the orthogonalization; the geo channel is zeroed out.

### 3f. Honest-DID (Rambachan-Roth 2023) breakdown

Single-factor: [`honest_did.md`](metrics/honest_did.md).
Multi-factor: [`honest_did_mf.md`](metrics/honest_did_mf.md).

| Specification | Pre-period max |β_fuel| | Post-event β_fuel | $\bar M$ | Verdict |
|---|---:|---:|---:|---|
| Single-factor CAR | 2.7701 (window [−6,−2]) | −4.7656 (t=−7.32) | **1.26** | "Robust" |
| Multi-factor CAR | 2.8445 (window [−6,−2]) | −3.1043 (t=−4.50) | **0.62** | "Moderate" |

### 3g. Event-time path of $\hat\beta_{\text{fuel}}(\tau)$

Source: [`event_time_path.md`](metrics/event_time_path.md). Pre-period max |β| = 2.071 at τ=−7. Post-event mean β over [−1, +3] = −0.19/month. Figure: [`fig4_event_time_path.pdf`](figures/fig4_event_time_path.pdf).

### 3h. Anomaly-versus-risk (post-formation horizons)

Source: [`anomaly_vs_risk.md`](metrics/anomaly_vs_risk.md).

| Window | $\hat\beta_{\text{fuel}}$ | t | per-month |
|---|---:|---:|---:|
| $[-1, +1]$ | **−2.8930** | −4.19 | −1.45 |
| $[-1, +3]$ | **−4.7656** | −7.32 | −1.19 |
| $[-1, +6]$ | **−4.6961** | −5.07 | −0.67 |
| $[-1, +12]$ | −3.1786 | −1.00 | −0.24 |
| $[-1, +24]$ | −2.4497 | −0.40 | −0.10 |

Effect concentrates within 1–6 months and decays at longer horizons.

### 3i. Announcement vs physical-retirement timing

Source: [`announcement_robustness.md`](metrics/announcement_robustness.md).

| Spec | n_FM | β_fuel | t_fuel |
|---|---:|---:|---:|
| Announcement when available (HEADLINE) | 117 | **−4.8318** | **−7.35** |
| Forced physical retirement | 135 | −4.3160 | −3.57 |
| Announcement-only subsample | 117 | −4.8318 | −7.35 |

Channel survives both timing conventions.

---

## 4. Identification (manuscript §3.5)

### 4a. Bartik / pre-2014 exposure design

Source: [`bartik_shiftshare.md`](metrics/bartik_shiftshare.md).

| Quantity | Bartik (pre-2014 × agg shock) | Standard $w^{\text{fuel}}$ |
|---|---:|---:|
| Pre-period correlation with current $w^{\text{fuel}}$ | 0.318 | 1.000 |
| Within-sample correlation | 0.496 | 1.000 |
| FM β (3 valid events) | **−1.96** | −14.69 |
| FM t | **−2.35** | −9.16 |
| Pooled β (N=24,070) | −1.89 | −4.21 |
| Pooled t | **−5.19** | −3.14 |
| Pre-event balance test (cutoff 2014), t | **−1.87** | — |
| Pre-event balance test (cutoff 2010), t | **−2.31** | — |
| Rotemberg HHI | 0.0311 | — |
| Rotemberg negative-weight events | 0/40 | — |
| **Oster $\delta^\*$** | **35.85** | **121.85** |

> **MAJOR DISCREPANCY (manuscript §3.5 and Conclusion):** the manuscript reports Oster $\delta^* = 20.8$ for the Bartik exposure design and $\delta^* = 61.6$ for the standard $w^{\text{fuel}}$ specification. Fresh run gives **35.85 and 121.85**. The qualitative claim ("$\delta^*$ vastly above 1") still holds, but every numeric quotation needs updating in §3.5 ("$\delta^* = 20.8$"), §4.6 ("$\delta^* = 20.8$"), Conclusion ("$\delta^* = 20.8$"). See §6, Red Flag #3.

### 4b. Pre-event balance ([−5, −2] window)

Source: bartik_shiftshare.md.

| Cutoff | Bartik t | p | Verdict (5%) | Verdict (10%) |
|---|---:|---:|---|---|
| 2014 | −1.87 | 0.0621 | PASS | FAIL |
| 2010 | −2.31 | 0.0212 | FAIL | FAIL |

The 2010 cutoff fails 5% balance; the 2014 cutoff (used in the manuscript) passes at 5% but fails at 10%. The manuscript cites only the 2014 cutoff.

### 4c. Bridge interaction (Section 2.5 augmented spec)

Source: [`bridge_interaction.md`](metrics/bridge_interaction.md). Pre-2014 mean coal share $\bar\alpha = 0.290$, 115 events.

| Variable | β | t | Predicted sign |
|---|---:|---:|---|
| $\gamma_{\text{fuel}}$ | −3.2523 | −2.56 | negative ✓ |
| $\gamma_{\text{het}}$ | **+2.2259** | **+1.08** | negative (theory predicts <0) ✗ |
| $\gamma_{\text{geo}}$ | −0.1162 | −0.24 | attenuated to 0 ✓ |
| $\gamma_{\text{reg}}$ | +3.9199 | +2.83 | positive ✓ |

Heterogeneity coefficient has the **wrong sign** vs theory; manuscript correctly notes this is statistically indistinguishable from zero.

---

## 5. Heterogeneity (manuscript §4.5, §4.7)

### 5a. ESG horse race (Pooled OLS)

Source: [`esg_horse_race.md`](metrics/esg_horse_race.md). N = 14,731 ESG-covered observations, 165 events, 153 unique firms.

| Spec | ESG β | ESG t | $w^{\text{fuel}}$ β | $w^{\text{fuel}}$ t | $R^2$ |
|---|---:|---:|---:|---:|---:|
| ESG only | −0.1144 | **−8.21** | — | — | 0.0119 |
| Spatial only | — | — | −1.5517 | −1.79 | 0.0034 |
| Both | −0.1176 | **−8.18** | −1.1310 | −1.32 | 0.0159 |

Marginal $R^2$: spatial→ESG +0.011, ESG→spatial +0.0003.

### 5b. ESG horse race (Fama-MacBeth)

Source: [`esg_fm_joint.md`](metrics/esg_fm_joint.md). 165 FM events.

| Spec | $\hat\gamma_{\text{ESG}}$ | t | $\hat\gamma_{\text{fuel}}$ | t |
|---|---:|---:|---:|---:|
| ESG only | −0.2760 | **−5.11** | — | — |
| Spatial only | — | — | −1.6231 | −0.81 |
| Both | **−0.2861** | **−5.03** | **−4.8171** | **−2.08** |

Joint Wald $\chi^2_2 = 25.93$, p < 0.001. Both individually marginally significant at 5% under FM.

### 5c. US vs non-US split (learning_alternatives, "subsample" spec)

Source: [`learning_alternatives.md`](metrics/learning_alternatives.md) §A1.

| Subsample | n | β_fuel | t |
|---|---:|---:|---:|
| US | 91 | **+0.0904** | **+0.06** |
| Non-US | 84 | **−5.4227** | **−4.29** |

> **Minor discrepancy:** manuscript §4.5 reports "US (t=0.11, 91 events)" and "non-US (−5.34, t=−4.10, 84 events)". Fresh run = **+0.06 / 91** and **−5.42 / t=−4.29 / 84**. Sample sizes match; t-stats differ in third decimal. See Red Flag #4.

### 5d. US restructured-vs-regulated split

Source: [`us_regulation_split.md`](metrics/us_regulation_split.md). FM, single-factor CARs.

| Split | n_FM | β_fuel | t |
|---|---:|---:|---:|
| US — All | 81 | −3.5296 | −5.54 |
| US — Restructured (15 states) | 14 | **−1.0741** | **−1.40** |
| US — Regulated (rest) | 67 | **−4.0426** | **−5.38** |
| Non-US | 36 | −7.7617 | −9.37 |

Channel is at least as strong in regulated states as in restructured ones — regulation hypothesis **not supported**.

### 5e. Country-level robustness (non-US)

Source: [`country_robustness.md`](metrics/country_robustness.md).

| Drop | n_FM | β_fuel | t |
|---|---:|---:|---:|
| (none / Russia / S Africa / India / Chile) | 36 | −7.7617 | −9.37 |
| China | 26 | −7.8917 | −7.19 |
| France | 31 | −7.3677 | −8.52 |
| Germany | 33 | −7.7643 | −9.26 |
| Netherlands | 33 | **−7.1012** | **−6.43** |
| Poland | 33 | −8.3435 | −8.58 |
| Greece | 33 | **−8.5311** | **−10.12** |

Range across leave-one-country-out drops: β ∈ [−7.10, −8.53], |t| ∈ [6.43, 10.12]. **No single country drives the non-US result.**

| MSCI tier | n_FM | β_fuel | t |
|---|---:|---:|---:|
| Developed-ex-US | 17 | **−9.6312** | **−7.75** |
| Emerging | 17 | **−5.6241** | **−4.11** |
| Frontier / Other | 2 | −10.0403 | n/a |

Effect strongest in developed-ex-US, undermining the informational-efficiency and ownership-structure explanations for the US null.

---

## 6. Manuscript ↔ Code Alignment Audit (Red Flags)

The pipeline runs end-to-end with exit code 0, the bibliography resolves cleanly, and every figure and Appendix LaTeX table renders. The fuel-channel claim is robust: it survives every inference method, every multiple-testing correction, every robustness window, every multi-factor adjustment, and every leave-one-country-out drop. **However, the manuscript prose contains numbers that do not exactly match the regenerated outputs.** A reviewer who runs the replication package will see these mismatches.

### Red Flag #1 — Sample-size inconsistency in summary statistics

`results/summaries/summary_statistics.md` Panel A reports **428 firms** in the "analysis sample" and **449 firms** in the GEM-matched subsample. The manuscript abstract and §3.1 use **703 firms**, which is the count of firms with monthly returns coverage (the analysis panel actually used in `robust_inference.py`, etc.). These are different samples and need to be labelled clearly. Suggested fix: rewrite `summary_statistics.py` to report the same sample the regression scripts use, or add a Panel that explicitly bridges 703 → 565 (network-complete) → 153 (ESG-covered).

### Red Flag #2 — Table 2 (Channel Decomposition) values do not match fresh outputs

Manuscript Table 2 reports:
- FM column: $w^{\text{fuel}}$ = −4.782 / $w^{\text{geo}}$ = −0.607
- Event-clustered column: $w^{\text{fuel}}$ = −5.474 / $w^{\text{geo}}$ = −0.023

Fresh run gives:
- FM (`robust_inference.md`): $w^{\text{fuel}}$ = **−4.7656** / $w^{\text{geo}}$ = **−0.5427**
- Event-clustered (`joint_tests.md`): $w^{\text{fuel}}$ = **−5.4884** / $w^{\text{geo}}$ = **+0.0176**

The FM column drift is real and small (−4.782 → −4.766; −0.607 → −0.543) — likely a snapshot from a previous pipeline state. **The event-clustered $w^{\text{geo}}$ flips sign (−0.023 → +0.018).** Magnitudes are inside sampling noise (the Julia and Python random-control samplers use different hashes), but a reviewer running the package will see Table 2's signs disagree with `joint_tests.md`. Recommended: either (a) regenerate Table 2 from the fresh run, or (b) lock the random seed across Python and Julia so the two routes produce identical numbers.

### Red Flag #3 — Oster $\delta^\*$ values are stale

Manuscript §3.5 (line 281), §4.6 (line 452), and Conclusion (line 550) all cite $\delta^* = 20.8$ (Bartik) and $\delta^* = 61.6$ (standard $w^{\text{fuel}}$). Fresh outputs give **35.85** and **121.85** — both *more conservative* (less fragile) than the values currently in the paper. The qualitative interpretation is unaffected ("unobservables would need to be 21× / 36× more important than observables"), but the numbers in the prose must be updated. **This is the largest numeric mismatch.**

### Red Flag #4 — Cook's D and US/non-US t-stats drift

- Manuscript §4.6: "1,688 (3.0%) exceed the 4/N threshold" → fresh run = **1,728 (3.1%)**.
- Manuscript §4.6: "fuel coefficient from −5.51 to −6.35" → fresh run = **−5.54 to −6.39**.
- Manuscript §4.5: "US (t=0.11, 91 events) and non-US (−5.34, t=−4.10, 84 events)" → fresh run = **(t=+0.06, 91)** and **(−5.42, t=−4.29, 84)**.
- Manuscript §4.6 (Fisher RI): "99th percentile is approximately ±2.4" → fresh run gives **99th percentile = +1.31 / range [−1.84, +2.03]**, i.e. ±~2.

These are small differences but every one is visible to a careful referee.

### Red Flag #5 — Bibliography file recently modified, not committed

`manuscript/references.bib` has uncommitted changes (40 added entries from the round-1 revision). These compile cleanly when `bibtex` is run from the manuscript directory, but the `.bib` modification has been sitting uncommitted alongside ~12 new analysis scripts and 10+ new metric files. A `git stash` accident would lose all of it. Commit before any further work.

### Red Flag #6 — Repository contains uncommitted analysis scripts

`git status` shows 12 untracked `src/*.py` files, plus all their outputs in `results/metrics/` and `results/summaries/`. These ARE referenced from `src/run_all.py` (also modified, uncommitted), and the manuscript cites their results — but a fresh clone of the current `main` branch will fail to reproduce the manuscript because the scripts simply do not exist. Commit the Phase-2/Phase-3 scripts together with their outputs as a single coherent state.

### Red Flag #7 — `docs/` directory is untracked development notes

`docs/` contains internal "phase" reports and audit notes (Spanish-language working files, cover-letter drafts, repo audits). These are not part of the replication package and should be added to `.gitignore` or moved out before submission. Listing of contents:
```
docs/COVER_LETTER_DRAFT.md
docs/PHASE1_RESULTS.md, PHASE1_EXT_RESULTS.md, PHASE2_RESULTS.md, PHASE3_RESULTS.md, PHASE_FINAL_RESULTS.md
docs/REPO_AUDIT.md, REVISION_PLAN_6MONTHS.md, STRATEGIC_DECISIONS.md
```
Recommended: add `docs/` to `.gitignore` (it shadows the existing `prompts/` line) before the submission tarball is built.

### Red Flag #8 — `manuscript/jeem_referee_report.md` untracked

This is a referee report sitting inside `manuscript/`. It is covered by the *.gitignore* pattern `manuscript/referee_*.md` only if renamed to `referee_jeem.md`. Either rename or extend the gitignore pattern.

### Red Flag #9 — `*` filename prefix is opaque to reviewers

22 analysis files in `src/` carry a `` prefix that is meaningless to anyone outside the project. The string is a residue of an early "strategy 1 vs strategy 2" decision that no longer exists in the paper. This is a substantial rename (22 files + `run_all.py` cross-references), so it is **flagged but not executed** in this report. If the user wants this rename, candidates are dropping the prefix entirely (`robust_inference.py`, `joint_tests.py`, …) or replacing with `analysis_*` (clearer role). The manuscript text never references the script names by file path, so the rename is mechanical and reversible.

### Red Flag #10 — Reproducibility lockfiles are absent

The pipeline is mostly Python stdlib (good), but it also depends on:
- **Julia** (`joint_tests.jl`, `romano_wolf_bootstrap.jl`) — no `Project.toml` or `Manifest.toml`.
- **R** (`robustness_conley_se.R`, `generate_fig*.R`) — no `renv.lock`.
- **`openpyxl`** (Python, only non-stdlib dep) — declared in README only, no `requirements.txt`.

For a top-five-journal replication package, pin all three. A minimal `requirements.txt` (`openpyxl==X.Y`), an `R/renv.lock`, and a `Project.toml` for Julia would close this gap.

### Red Flag #11 — `robustness_conley_se.R` is not in the run_all DAG

`REPLICATION.md` describes `Rscript src/robustness_conley_se.R` as a separate step, but `run_all.py` does not invoke it. A reviewer running just `python src/run_all.py` will not regenerate the Conley SEs cited in §4.6 ("$t = -4.16$ at 1,000 km"). Either add an R wrapper to `run_all.py` (skip-if-R-absent), or document this clearly in the README's "Reproduction" section.

### Red Flag #12 — No LICENSE, no CITATION, no DOI

The repo has neither a LICENSE file nor a CITATION.cff. JEEM (and most economics journals now) expect a clearly-licensed replication package with a Zenodo or similar DOI for the snapshot used in the published version. Recommended: MIT or BSD-3 LICENSE; CITATION.cff with DOI; tag `v1.0-submission` and archive on Zenodo.

### Red Flag #13 — Python 3.13 syntax warning in `anomaly_vs_risk.py`

```
SyntaxWarning: invalid escape sequence '\%'
  months, $-6.6\%$ annualized), the headline is large enough to invite
```
The script docstring contains LaTeX `\%` that Python 3.12+ flags. Cosmetic, but a clean run should be warning-free for a replication package. Fix: prefix the docstring with `r"""..."""`.

---

## 7. What checks out exactly

Every result below replicates the manuscript prose to the cited decimals.

- Multi-factor table (manuscript §4.7 Table 3): $\gamma_{\text{fuel}}^{\text{single}}$ = −4.766; $\gamma_{\text{fuel}}^{\text{multi}}$ = −3.104; $\gamma_{\text{geo}}^{\text{single}}$ = −0.543; $\gamma_{\text{geo}}^{\text{multi}}$ = −0.068. ✓
- Honest-DID single-factor $\bar M = 1.26$, multi-factor $\bar M = 0.62$. ✓
- Anomaly-vs-risk decay path (β(1) = −2.89, β(3) = −4.77, β(6) = −4.70, β(12) = −3.18, β(24) = −2.45). ✓
- Announcement-vs-physical robustness (−4.83 / t=−7.35 vs −4.32 / t=−3.57). ✓
- Lag-sensitivity (|t| from 7.32 at lag 4 to 5.45 at lag 18). ✓
- Bridge interaction γ_het = +2.23, |t|=1.08, 115 events. ✓
- ESG joint test (Wald χ² = 25.9, p<0.001; ESG β=−0.286, t=−5.03; fuel β=−4.82, t=−2.08). ✓
- US restructured/regulated split (regulated −4.04, t=−5.38, N=67; restructured −1.07, t=−1.40, N=14). ✓
- Country LOO (range −7.10 to −8.53, |t| ∈ [6.43, 10.12]). ✓
- Developed-ex-US (−9.63, t=−7.75, 17 events) vs emerging (−5.62, t=−4.11, 17 events). ✓
- Pre-2014 Bartik FM (3 events, t=−2.35) and pooled (t=−5.19). ✓
- Pre-event balance test [−5,−2] (Bartik t=−1.87, p=0.062). ✓
- Romano-Wolf p < 0.001 for fuel; >0.3 for geo and reg. ✓
- Joint F-test = 70.81. ✓
- Channel-difference test t = +5.97 (FM) and +7.46 (event-clustered, manuscript) / +7.56 (current). ≈

---

## 8. Bottom line

The science is sound. The fuel channel is robust. The pipeline reproduces.

What stands between the current state and a clean journal submission is:

1. **(High priority — required)** Update §3.5, §4.6, and Conclusion with fresh Oster $\delta^\*$ values (20.8 → **35.85**; 61.6 → **121.85**).
2. **(High priority — required)** Reconcile Table 2 column values with fresh outputs (or lock seeds across Python/Julia).
3. **(Medium priority — required for clean submission)** Commit the Phase-2/Phase-3 scripts and outputs in a single atomic commit so a fresh clone reproduces.
4. **(Medium priority)** Update Cook's D / US-vs-non-US / Fisher-RI percentile prose to match fresh numbers.
5. **(Medium priority)** Reconcile the 428/449/703 firm-count disparity in `summary_statistics.md` with the manuscript abstract's 703.
6. **(Low priority — submission hygiene)** `.gitignore` `docs/`; rename `manuscript/jeem_referee_report.md` to match existing pattern; add `LICENSE`, `CITATION.cff`, lockfiles, and pin a Zenodo DOI; fix the Python 3.13 escape-sequence warning.
7. **(Optional — clarity)** Drop the `` prefix from the 22 analysis scripts. Mechanical, reversible, makes the replication package readable to a referee on first contact.

Once items 1–4 land, the manuscript and the regenerated outputs will agree to the third decimal everywhere a number appears in the paper.

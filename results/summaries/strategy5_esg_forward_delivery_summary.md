# Strategy 5B: ESG Forward-Delivery Horse Race (Clean Test)

**Script:** `strategy5_esg_forward_delivery.py`  
**Run date:** 2026-02-20

## Design
- Outcome: forward delivery change, `?delivery_{t?t+h} = (1-a_{t+h}) - (1-a_t)` using `firm_alpha_panel.csv`.
- Predictors: full theta (alpha, lambda, rho, kappa, delta) + ESG scores.
- Panel with firm clustering (HC1 robust, clustered by gvkey).
- Horizons: 3 years and 5 years.

## Horizon = 3 years
N = 1,568 firm-year pairs
- Theta only: R = 0.0127
- + ESG score: R = 0.0134 (dR = +0.0007), ESG t = 0.72
- + Env avg: R = 0.0151 (dR = +0.0024), env_avg t = 1.28
- + E subscores: R = 0.0225 (dR = +0.0098)
  - emissions t = 1.83, innovation t = -0.39, resource_use t = -1.25

## Horizon = 5 years
N = 1,321 firm-year pairs
- Theta only: R = 0.0223
- + ESG score: R = 0.0239 (dR = +0.0016), ESG t = 0.76
- + Env avg: R = 0.0272 (dR = +0.0049), env_avg t = 1.32
- + E subscores: R = 0.0389 (dR = +0.0165)
  - emissions t = 1.77, innovation t = -0.55, resource_use t = -0.98

## Density split (delivery ~ ESG score, latest-year cross-section)
Median density = 5
- Low density: beta = +0.0027, R = 0.0144 (N = 67)
- High density: beta = +0.0014, R = 0.0043 (N = 62)

## Interpretation
- The clean forward-looking test shows **minimal incremental power** from ESG scores once theta (including alpha) is included.
- ESG adds at most 12 percentage points in R under the richest spec, with weak t-stats.
- ESG predictive power is **lower in dense networks**, consistent with the regime-dependent information loss.

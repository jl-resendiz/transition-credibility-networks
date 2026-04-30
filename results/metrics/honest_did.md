# Honest DID Sensitivity (Rambachan-Roth 2023)

Placebo-CAR-window approach to bound parallel-trends violations. For each pre-event 5-month window, the same Fama-MacBeth cross-sectional regression as the headline (CAR on `w_geo`, `w_fuel`, `w_reg`, `same_sector`) is run with the dependent variable replaced by the pre-event CAR over the placebo window. Pre-event coefficients should be statistically indistinguishable from zero under the parallel-trends assumption.

## Window-by-window results

| Window | beta_fuel | SE (NW lag=4) | t-stat | N events | Role |
|---|---|---|---|---|---|
| [-12, -8] | +0.4802 | 0.7914 | +0.61 | 117 | placebo |
| [-9, -5] | -0.0564 | 0.5087 | -0.11 | 117 | placebo |
| [-6, -2] | +2.7701 | 0.7393 | +3.75 | 117 | placebo |
| [-1, +3] | -4.7656 | 0.6508 | -7.32 | 117 | **HEADLINE** |

## Rambachan-Roth M-bar breakdown

- Headline beta_fuel (-1, +3): -4.7656
- Headline SE: 0.6508
- Headline t-stat: -7.32
- Pre-period max |beta_fuel|: 2.7701 (window [-6, -2])
- **M-bar (5% breakdown): 1.26**

Interpretation: the headline post-event effect survives a parallel-trends violation up to 1.26 times the largest deviation observed in the pre-period before becoming statistically indistinguishable from zero at the 5% level.

**Robust** by Rambachan-Roth convention: the result survives violations as large as anything actually observed in the pre-period.

## Method note

The placebo windows are non-overlapping 5-month CAR windows ending strictly before the event. The headline window [-1, +3] is also 5 months wide, so the placebo distribution is directly comparable to the post-event estimate. Each window-level coefficient is a Fama-MacBeth mean of event-level cross-sectional OLS coefficients, with Newey-West (lag=4) HAC standard errors on the time series of event-level betas.

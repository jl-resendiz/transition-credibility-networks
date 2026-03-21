# Focused Hypothesis Tests: Channel Split

Window: [-1, +3] months (monthly CARs, vwretd)
Events: 179 first-mover coal retirements
N = 72398 observations, 175 event clusters
Standard errors: event-clustered

## Test 1: Joint F-test (H0: beta_geo = beta_fuel = beta_reg = 0)

Unrestricted: CAR = alpha + beta_geo * w^geo + beta_fuel * w^fuel + beta_reg * w^reg + beta_s * SameSector
Restricted:   CAR = alpha + beta_s * SameSector

SSR_restricted:   9211.130537
SSR_unrestricted: 9206.825727
F-statistic: 11.2829
df: (3, 72393)
p-value (permutation, B=999): 0.0010***
N: 72398

Interpretation: Reject H0 at 5%. The spatial network channels jointly predict CARs around coal retirement events.

## Test 2: Difference test (H0: beta_geo = beta_fuel)

beta_geo:  +0.354224 (SE 0.119191)
beta_fuel: -1.496868 (SE 0.473663)
Difference (beta_geo - beta_fuel): +1.851092
SE of difference: 0.507637
t-statistic: 3.646
p-value: 0.0003***
Cov(beta_geo, beta_fuel): -0.0095662284

Interpretation: The opposing-sign channel split is statistically significant in a single test (t = 3.646, p = 0.0003).

## Full regression coefficients

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.005606 | 0.005628 | 0.996 | 0.3192 |
| w_geo | +0.354224 | 0.119191 | 2.972 | 0.0030*** |
| w_fuel | -1.496868 | 0.473663 | -3.160 | 0.0016*** |
| w_reg | +1.126878 | 0.690894 | 1.631 | 0.1029 |
| same_sector | +0.021264 | 0.004969 | 4.279 | 0.0000*** |

## Comparison with individual tests

| Test | Hypotheses tested | Correction needed | Result |
|---|---|---|---|
| Individual t-tests | 9 | Romano-Wolf | 0/9 significant |
| Joint F-test | 1 | None | significant (p=0.0010) |
| Difference test | 1 | None | significant (p=0.0003) |

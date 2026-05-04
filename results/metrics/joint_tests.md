# Focused Hypothesis Tests: Channel Split

Window: [-1, +3] months (monthly CARs, vwretd)
Events: 179 first-mover-matched (175 used in pooled regression below; 117 with ≥20 firms qualify for FM)
N = 55580 observations, 175 event clusters
Standard errors: event-clustered

## Test 1: Joint F-test (H0: beta_geo = beta_fuel = beta_reg = 0)

Unrestricted: CAR = alpha + beta_geo * w^geo + beta_fuel * w^fuel + beta_reg * w^reg + beta_s * SameSector
Restricted:   CAR = alpha + beta_s * SameSector

SSR_restricted:   3946.450094
SSR_unrestricted: 3931.422831
F-statistic: 70.8090
df: (3, 55575)
p-value (permutation, B=999): 0.0000***
N: 55580

Reject H0 at 5%. The spatial network channels jointly predict CARs around coal retirement events.

## Test 2: Difference test (H0: beta_geo = beta_fuel)

beta_geo:  +0.017592 (SE 0.100693)
beta_fuel: -5.488354 (SE 0.727851)
Difference (beta_geo - beta_fuel): +5.505946
SE of difference: 0.728539
t-statistic: 7.558
p-value: 0.0000***
Cov(beta_geo, beta_fuel): +0.0045687476

Interpretation: The opposing-sign channel split is statistically significant in a single test (t = 7.558, p = 0.0000).

## Full regression coefficients

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035757 | 0.007665 | 4.665 | 0.0000*** |
| w_geo | +0.017592 | 0.100693 | 0.175 | 0.8613 |
| w_fuel | -5.488354 | 0.727851 | -7.540 | 0.0000*** |
| w_reg | +1.441175 | 1.050350 | 1.372 | 0.1700 |
| same_sector | +0.033207 | 0.008873 | 3.743 | 0.0002*** |

## Comparison with individual tests

| Test | Hypotheses tested | Correction needed | Result |
|---|---|---|---|
| Individual t-tests | 3 | Romano-Wolf | 1/3 significant (fuel) |
| Joint F-test | 1 | None | significant (p=0.0000) |
| Difference test | 1 | None | significant (p=0.0000) |

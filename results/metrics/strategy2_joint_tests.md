# Focused Hypothesis Tests: Channel Split

Window: [-1, +3] months (monthly CARs, vwretd)
Events: 179 first-mover coal retirements
N = 55580 observations, 175 event clusters
Standard errors: event-clustered

## Test 1: Joint F-test (H0: beta_geo = beta_fuel = beta_reg = 0)

Unrestricted: CAR = alpha + beta_geo * w^geo + beta_fuel * w^fuel + beta_reg * w^reg + beta_s * SameSector
Restricted:   CAR = alpha + beta_s * SameSector

SSR_restricted:   3946.450094
SSR_unrestricted: 3931.417645
F-statistic: 70.8335
df: (3, 55575)
p-value (permutation, B=999): 0.0000***
N: 55580

Reject H0 at 5%. The spatial network channels jointly predict CARs around coal retirement events.

## Test 2: Difference test (H0: beta_geo = beta_fuel)

beta_geo:  -0.022960 (SE 0.056598)
beta_fuel: -5.474254 (SE 0.730054)
Difference (beta_geo - beta_fuel): +5.451294
SE of difference: 0.730684
t-statistic: 7.461
p-value: 0.0000***
Cov(beta_geo, beta_fuel): +0.0011419292

Interpretation: The opposing-sign channel split is statistically significant in a single test (t = 7.461, p = 0.0000).

## Full regression coefficients

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035796 | 0.007657 | 4.675 | 0.0000*** |
| w_geo | -0.022960 | 0.056598 | -0.406 | 0.6850 |
| w_fuel | -5.474254 | 0.730054 | -7.498 | 0.0000*** |
| w_reg | +1.452522 | 1.051538 | 1.381 | 0.1672 |
| same_sector | +0.033201 | 0.008872 | 3.742 | 0.0002*** |

## Comparison with individual tests

| Test | Hypotheses tested | Correction needed | Result |
|---|---|---|---|
| Individual t-tests | 9 | Romano-Wolf | 0/9 significant |
| Joint F-test | 1 | None | significant (p=0.0000) |
| Difference test | 1 | None | significant (p=0.0000) |

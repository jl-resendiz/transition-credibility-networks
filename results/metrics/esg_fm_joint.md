# ESG Horse Race under Fama-MacBeth Inference

Sample restricted to firms with Refinitiv Environmental Score (14731 firm-event observations across 165 events). For each event, runs cross-sectional OLS within the ESG-covered subsample; aggregates with Newey-West (lag 4) HAC standard errors.

## Spec: ESG only

| Variable | Mean | SE (NW lag=4) | t | N events |
|---|---|---|---|---|
| intercept | +0.2075 | 0.0444 | +4.68 | 165 |
| esg_score | -0.2760 | 0.0540 | -5.11 | 165 |

## Spec: Spatial only (fuel + geo)

| Variable | Mean | SE (NW lag=4) | t | N events |
|---|---|---|---|---|
| intercept | +0.0404 | 0.0161 | +2.51 | 165 |
| w_fuel | -1.6231 | 1.9958 | -0.81 | 165 |
| w_geo | -463.2732 | 207.3774 | -2.23 | 165 |

## Spec: Both (ESG + spatial)

| Variable | Mean | SE (NW lag=4) | t | N events |
|---|---|---|---|---|
| intercept | +0.2374 | 0.0507 | +4.68 | 165 |
| esg_score | -0.2861 | 0.0569 | -5.03 | 165 |
| w_fuel | -4.8171 | 2.3135 | -2.08 | 165 |
| w_geo | -386.5638 | 168.9175 | -2.29 | 165 |

## Joint test (Both spec): H0: gamma_ESG = 0 AND gamma_fuel = 0

- Estimated mean ESG coefficient: -0.2861
- Estimated mean fuel coefficient: -4.8171
- Wald statistic: 25.925 (chi-squared with 2 df)
- Approximate p-value: 0.0000

**Reject the joint null** at the 5% level: at least one of the two coefficients is non-zero in the FM cross-section.

Both individually survive marginal significance at the 5% level. The paper claim that "both ESG and fuel-mix similarity survive a joint test" under FM inference is supported.

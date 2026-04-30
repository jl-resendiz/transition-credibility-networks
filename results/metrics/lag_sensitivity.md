# Newey-West Lag Sensitivity for Fama-MacBeth Coefficients

Time series of event-level FM coefficients: T = 117 events.

Standard lag rules of thumb suggest lags of 4 to 5 for this T. However, the event windows are 5 months wide and 63% of months host multiple overlapping windows, which can extend serial correlation in the FM time series beyond rule-of-thumb lags. We therefore report lags {4, 8, 12, 18} for transparency.

Driscoll-Kraay (1998) HAC standard errors are designed for panel data with contemporaneous cross-sectional dependence. The Fama-MacBeth time series has already collapsed the cross-section into a scalar coefficient per event, so Driscoll-Kraay applied to this series reduces algebraically to Newey-West with the same lag. Only Newey-West is reported.

## Lag sensitivity table

| Channel | Mean | SE iid | SE NW(4) | SE NW(8) | SE NW(12) | SE NW(18) | t NW(4) | t NW(8) | t NW(12) | t NW(18) |
|---|---|---|---|---|---|---|---|---|---|---|
| $\gamma_{\text{fuel}}$ | -4.7656 | 0.4705 | 0.6508 | 0.7305 | 0.7916 | 0.8743 | -7.32 | -6.52 | -6.02 | -5.45 |
| $\gamma_{\text{geo}}$ | -0.5427 | 0.2315 | 0.3090 | 0.3237 | 0.3345 | 0.3261 | -1.76 | -1.68 | -1.62 | -1.66 |
| $\gamma_{\text{reg}}$ | +2.6975 | 0.6554 | 0.9518 | 0.9475 | 0.9058 | 0.9026 | +2.83 | +2.85 | +2.98 | +2.99 |
| $\gamma_{\text{geo}} - \gamma_{\text{fuel}}$ | +4.2229 | 0.4759 | 0.7076 | 0.7886 | 0.8354 | 0.8936 | +5.97 | +5.35 | +5.05 | +4.73 |

## Interpretation

A coefficient that loses statistical significance only at long lags suggests serial correlation in the FM time series that the baseline lag may understate. A coefficient whose t-statistic is stable across lags indicates a result that does not depend on lag choice.

For the headline fuel coefficient, |t| ranges from 5.45 to 7.32 across lags [4, 8, 12, 18]. All lags reject the null at the 5% level, so the result does not depend on the lag choice within this range.

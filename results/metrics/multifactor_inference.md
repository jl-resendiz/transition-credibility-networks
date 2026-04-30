# Multi-Factor Abnormal Returns: FM Cross-Sectional Inference

Replaces single-factor (market-adjusted) abnormal returns with a 4-factor model (FF3 Mkt-RF, SMB, HML + sample-constructed utility industry excess return) estimated firm-by-firm on a 24-month pre-event window. CAR is the within-window prediction error over the headline window [-1, +3].

## Headline comparison (FM + Newey-West, lag=4)

| Channel | Single-factor | Multi-factor | Shrinkage |
|---|---|---|---|
| $\gamma_{\text{fuel}}$ | -4.7656 (0.6508) [t=-7.32] | -3.1043 (0.6898) [t=-4.50] | +34.9% |
| $\gamma_{\text{geo}}$ | -0.5427 (0.3090) [t=-1.76] | -0.0679 (0.5306) [t=-0.13] | +87.5% |
| $\gamma_{\text{reg}}$ | +2.6975 (0.9518) [t=+2.83] | +2.4749 (0.8247) [t=+3.00] | +8.3% |
| $\gamma_{\text{geo}} - \gamma_{\text{fuel}}$ | +4.2229 (0.7076) [t=+5.97] | +3.0364 (0.9176) [t=+3.31] | +28.1% |

## Specification

Pre-event regression for each (firm, event):

$$(r_{it} - rf_t) = \alpha_i + \beta_M (Mkt-RF)_t + \beta_S SMB_t + \beta_V HML_t + \beta_U (UTL_t - rf_t) + \epsilon_{it}$$

estimated on the 24-month window pre-event (minimum 12 obs). Abnormal return at month $t$ is $AR_{it} = (r_{it}-rf_t) - \hat{\alpha}_i - \hat{\beta} \cdot factors_t$.

## Notes

- Events with valid FM regressions: 117
- Total firm-event observations: 55474
- Utility factor: equal-weighted mean of all sample firm returns at each month (require >= 30 firms).
- Currency factor: omitted. The FF Mkt-RF is a US factor; for non-US firms, a USD trade-weighted index would be appropriate but is not in the repo. Limitation documented.
- Honest DID and lag sensitivity should be re-run on these multi-factor CARs in a future revision.

# 5-Factor Inference: FF3 + UMD + Utility Industry

Adds Carhart momentum (UMD) to the existing 4-factor specification:

  AR_it = r_it − rf_t − [α_i + β_M(Mkt−RF)_t + β_S SMB_t + β_H HML_t
                          + β_W UMD_t + β_U(UTL−rf)_t]

Firm-by-firm betas estimated on a 24-month pre-event window. CAR is
the within-window prediction error over [-1, +3] months.

Events with valid FM regression: 117
Total firm-event observations: 55474
Avg firms per event: 244.2
Avg within-event R²: 0.0430

## Headline coefficients (FM + NW lag 4)

| Variable | Mean β | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.0099 | 0.0060 | +1.637 | 0.1016 |
| w_geo | -0.2434 | 0.4334 | -0.562 | 0.5743 |
| w_fuel | -2.8493 | 0.8018 | -3.554 | 0.0004*** |
| w_reg | +2.6542 | 0.9400 | +2.824 | 0.0047*** |
| same_sector | +0.0085 | 0.0065 | +1.308 | 0.1907 |

**Difference (γ_geo − γ_fuel):** +2.6058, NW SE = 1.0160, t = +2.565

## Interpretation

A negative and statistically significant γ_fuel under 5-factor
adjustment indicates that the channel is NOT absorbed by:
- Market (Mkt-RF)
- Size (SMB)
- Value (HML)
- Momentum (UMD) — the new factor in this 5F spec, addressing the
  referee concern that coal-heavy peers carry persistent negative
  momentum and the channel may be a "low-momentum trap".
- Utility-industry portfolio (UTL−rf)

Comparison with the existing 4-factor result (in multifactor_inference.md):
4-factor (FF3 + Utility): γ_fuel = -3.10, t = -4.50.
If 5-factor preserves significance, the channel adds explanatory power
beyond standard risk factors AND momentum.

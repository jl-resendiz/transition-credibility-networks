# Honest DID Sensitivity (Multi-Factor CARs)

Recalibration of the placebo-CAR-window Rambachan-Roth analysis using 4-factor (FF3 + sample-constructed utility industry) abnormal returns, replacing the single-factor market-adjusted CARs used in `honest_did.md`.

## Window-by-window results

| Window | beta_fuel | SE (NW lag=4) | t-stat | N events | Role |
|---|---|---|---|---|---|
| [-12, -8] | -1.0980 | 0.5807 | -1.89 | 117 | placebo |
| [-9, -5] | +0.4919 | 0.4053 | +1.21 | 117 | placebo |
| [-6, -2] | +2.8445 | 0.5635 | +5.05 | 117 | placebo |
| [-1, +3] | -3.1043 | 0.6898 | -4.50 | 117 | **HEADLINE** |

## Recalibrated M-bar

- Headline beta_fuel (-1, +3): -3.1043
- Headline SE: 0.6898
- Headline t-stat: -4.50
- Pre-period max |beta_fuel|: 2.8445 (window [-6, -2])
- **M-bar (5% breakdown, multi-factor): 0.62**

For comparison, M-bar on single-factor CARs was 1.26 (reported in `honest_did.md`).

**Moderate** robustness under the multi-factor CARs.

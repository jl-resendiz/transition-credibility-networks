# Anomaly vs Risk: Post-Formation Decay Test

Tests whether the fuel-similarity premium accumulates with horizon (risk-premium hypothesis) or plateaus / decays (mispricing hypothesis).

## Cumulative effect by post-event horizon

| Window | $\hat\beta_{\mathrm{fuel}}$ | SE (NW lag=4) | t-stat | per-month | N events |
|---|---|---|---|---|---|
| [-1, +1] | -2.8930 | 0.6907 | -4.19 | -1.4465 | 117 |
| [-1, +3] | -4.7656 | 0.6508 | -7.32 | -1.1914 | 117 |
| [-1, +6] | -4.6961 | 0.9261 | -5.07 | -0.6709 | 117 |
| [-1, +12] | -3.1786 | 3.1945 | -1.00 | -0.2445 | 117 |
| [-1, +24] | -2.4497 | 6.1414 | -0.40 | -0.0980 | 117 |

## Interpretation

- $\hat\beta(1) = -2.8930$, $\hat\beta(3) = -4.7656$, $\hat\beta(6) = -4.6961$, $\hat\beta(12) = -3.1786$, $\hat\beta(24) = -2.4497$.

- Ratio $\hat\beta(24)/\hat\beta(3) = 0.51$. The effect decays at longer horizons, consistent with a mispricing interpretation: institutional arbitrage corrects the initial under-reaction over time.

## Caveats

- This is one of three demarcation tests. The institutional ownership split and the characteristic-matched (DGTW) benchmark, flagged by referees as natural complements, require additional data not included in the replication package; these are reserved for future versions.
- Longer horizons are subject to greater confounding from subsequent retirement events, so the H=24 result should be read as suggestive rather than definitive.

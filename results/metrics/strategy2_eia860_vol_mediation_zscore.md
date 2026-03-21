# EIA860 Volatility Mediation

- event_scope: all_matched
- transform: zscore
- N: 25384

## Volatility change ~ exposure
beta_w=+0.0005, se=0.0001, t=3.89, clusters=91

## CAR ~ exposure
beta_w=-0.0019, se=0.0108, t=-0.17, clusters=91

## CAR ~ exposure + vol_change
beta_w=-0.0028, se=0.0109, t=-0.25, clusters=91
beta_vol=+1.9178, se=0.5819, t=3.30, clusters=91

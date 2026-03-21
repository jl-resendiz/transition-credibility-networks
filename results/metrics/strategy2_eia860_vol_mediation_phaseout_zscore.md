# EIA860 Volatility Mediation

- event_scope: all_matched
- transform: zscore
- N: 129420

## Volatility change ~ exposure
beta_w=-0.0006, se=0.0002, t=-2.81, clusters=26

## CAR ~ exposure
beta_w=-0.0674, se=0.0096, t=-7.02, clusters=26

## CAR ~ exposure + vol_change
beta_w=-0.0691, se=0.0098, t=-7.04, clusters=26
beta_vol=-2.6579, se=1.8677, t=-1.42, clusters=26

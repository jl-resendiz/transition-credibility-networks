# EIA860 Volatility Mediation

- event_scope: all_matched
- transform: log1p
- N: 25384

## Volatility change ~ exposure
beta_w=+0.0610, se=0.0153, t=3.98, clusters=91

## CAR ~ exposure
beta_w=-0.2572, se=1.3998, t=-0.18, clusters=91

## CAR ~ exposure + vol_change
beta_w=-0.3741, se=1.4173, t=-0.26, clusters=91
beta_vol=+1.9178, se=0.5819, t=3.30, clusters=91

# EIA860 Volatility Mediation

- event_scope: all_matched
- transform: log1p
- N: 129420

## Volatility change ~ exposure
beta_w=-0.1040, se=0.0381, t=-2.73, clusters=26

## CAR ~ exposure
beta_w=-11.4034, se=1.5736, t=-7.25, clusters=26

## CAR ~ exposure + vol_change
beta_w=-11.6800, se=1.5991, t=-7.30, clusters=26
beta_vol=-2.6606, se=1.8670, t=-1.43, clusters=26

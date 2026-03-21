# Phase-out CAR Robustness

- events: 14 (tier_filter=1, binding_only=True)
- transform: log1p
- CAR window: [-1,20]
- N obs: 20431

## Main regression (event-clustered)
beta=-13.3306, se=7.6104, t=-1.75, clusters=14, N=20431

## Leave-one-out (event)
- beta min=-17.0942, max=-6.7233, mean=-13.3329
- t min=-2.31, max=-1.36, mean=-1.70

## Placebo timing (shifted months)
- shift -6: beta=-4.5299, se=2.6477, t=-1.71, clusters=14, N=21717
- shift +6: beta=-2.5542, se=1.4230, t=-1.79, clusters=14, N=20932
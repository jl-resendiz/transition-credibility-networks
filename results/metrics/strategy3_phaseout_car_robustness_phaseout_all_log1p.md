# Phase-out CAR Robustness

- events: 27 (tier_filter=, binding_only=False)
- transform: log1p
- CAR window: [-1,20]
- N obs: 214871

## Main regression (event-clustered)
beta=-5.7052, se=1.1279, t=-5.06, clusters=27, N=214871

## Leave-one-out (event)
- beta min=-6.1556, max=-5.1511, mean=-5.7085
- t min=-5.39, max=-4.57, mean=-4.97

## Placebo timing (shifted months)
- shift -6: beta=-1.7532, se=0.3199, t=-5.48, clusters=27, N=224974
- shift +6: beta=-1.7854, se=0.3403, t=-5.25, clusters=27, N=219137
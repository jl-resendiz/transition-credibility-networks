# Phase-out CAR Robustness

- events: 27 (tier_filter=, binding_only=False)
- transform: base
- CAR window: [-1,20]
- N obs: 214871

## Main regression (event-clustered)
beta=-4.4052, se=0.9168, t=-4.80, clusters=27, N=214871

## Leave-one-out (event)
- beta min=-4.7704, max=-3.9672, mean=-4.4081
- t min=-5.08, max=-4.35, mean=-4.72

## Placebo timing (shifted months)
- shift -6: beta=-1.3486, se=0.2555, t=-5.28, clusters=27, N=224974
- shift +6: beta=-1.3923, se=0.2788, t=-4.99, clusters=27, N=219137
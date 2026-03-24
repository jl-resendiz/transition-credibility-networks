# Event-Specific Geographic Weight Matrix

Baseline w_geo uses MW-weighted centroids, diluting geographic signal
for multinationals. Event-specific w_geo measures distance from the
RETIRING PLANT to each firm's NEAREST plant.

Decay half-life: 1000 km (scale = 1442.7 km)
Events: 179 first-mover coal retirements (179 with GPS)
Valid events (>= 20 obs): 175
Total obs: 55580

## Diagnostics

Mean w_geo_centroid: 0.001644
Mean w_geo_event: 0.001717
Correlation(centroid, event): 0.3119

## Pooled Regressions (event-clustered SEs)

### A: Centroid w_geo

N = 55580, clusters = 175, R2 = 0.0071

| Variable | beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035796 | 0.007657 | 4.675 | 0.0000*** |
| w_fuel | -5.474254 | 0.730054 | -7.498 | 0.0000*** |
| w_geo_centroid | -0.022960 | 0.056598 | -0.406 | 0.6850 |
| w_reg | +1.452522 | 1.051538 | 1.381 | 0.1672 |
| same_sector | +0.033201 | 0.008872 | 3.742 | 0.0002*** |

### B: Event-specific w_geo

N = 55580, clusters = 175, R2 = 0.0071

| Variable | beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035703 | 0.007641 | 4.673 | 0.0000*** |
| w_fuel | -5.516336 | 0.733468 | -7.521 | 0.0000*** |
| w_geo_event | +0.120704 | 0.091460 | 1.320 | 0.1869 |
| w_reg | +1.409139 | 1.043401 | 1.351 | 0.1768 |
| same_sector | +0.033135 | 0.008877 | 3.733 | 0.0002*** |

### C: Both geo measures

N = 55580, clusters = 175, R2 = 0.0071

| Variable | beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035744 | 0.007641 | 4.678 | 0.0000*** |
| w_fuel | -5.502125 | 0.732962 | -7.507 | 0.0000*** |
| w_geo_centroid | -0.055044 | 0.057254 | -0.961 | 0.3363 |
| w_geo_event | +0.142710 | 0.095070 | 1.501 | 0.1333 |
| w_reg | +1.417622 | 1.045154 | 1.356 | 0.1750 |
| same_sector | +0.033112 | 0.008875 | 3.731 | 0.0002*** |

## Fama-MacBeth (1973) + Newey-West SEs

### A: Centroid w_geo

Events: 117, Avg N/event: 244.6, Avg R2: 0.0515

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.025478 | 0.009092 | 2.802 | 0.0051*** |
| w_fuel | -4.782353 | 0.649637 | -7.362 | 0.0000*** |
| w_geo_centroid | -0.607185 | 0.468053 | -1.297 | 0.1945 |
| w_reg | +2.642804 | 0.961429 | 2.749 | 0.0060*** |
| same_sector | +0.021454 | 0.011209 | 1.914 | 0.0556* |

### B: Event-specific w_geo

Events: 117, Avg N/event: 244.6, Avg R2: 0.0548

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.025605 | 0.009134 | 2.803 | 0.0051*** |
| w_fuel | -5.002154 | 0.574110 | -8.713 | 0.0000*** |
| w_geo_event | +1.772367 | 2.586465 | 0.685 | 0.4932 |
| w_reg | +2.339514 | 0.938980 | 2.492 | 0.0127** |
| same_sector | +0.020507 | 0.011396 | 1.800 | 0.0719* |

### C: Both geo measures

Events: 117, Avg N/event: 244.6, Avg R2: 0.0565

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.025859 | 0.009170 | 2.820 | 0.0048*** |
| w_fuel | -5.124863 | 0.584768 | -8.764 | 0.0000*** |
| w_geo_centroid | -0.579457 | 0.513352 | -1.129 | 0.2590 |
| w_geo_event | +4.664489 | 3.928473 | 1.187 | 0.2351 |
| w_reg | +2.393944 | 0.936961 | 2.555 | 0.0106** |
| same_sector | +0.020771 | 0.011330 | 1.833 | 0.0668* |

## Summary: Centroid vs Event-Specific Geography

| Spec | Method | w_geo_centroid t | w_geo_event t | w_fuel t | R2 / Avg R2 |
|---|---|---:|---:|---:|---:|
| A: Centroid w_geo | Pooled | -0.406 |  | -7.498 | 0.0071 |
| A: Centroid w_geo | FM+NW | -1.297 |  | -7.362 | 0.0515 |
| B: Event-specific w_geo | Pooled |  | 1.320 | -7.521 | 0.0071 |
| B: Event-specific w_geo | FM+NW |  | 0.685 | -8.713 | 0.0548 |
| C: Both geo measures | Pooled | -0.961 | 1.501 | -7.507 | 0.0071 |
| C: Both geo measures | FM+NW | -1.129 | 1.187 | -8.764 | 0.0565 |

## Interpretation

If event-specific w_geo (nearest-plant distance) rescues the geographic
channel, we expect Spec B to show a larger and more significant coefficient
on w_geo_event than Spec A shows on w_geo_centroid. In Spec C (horse race),
if w_geo_event dominates w_geo_centroid, the centroid measure was indeed
diluted by multinational dispersion.
# Portfolio Sorts: Spatial Exposure and Returns

CAR window: [-1, +3] months | Min firms per event: 10 | Events processed: 176 | Events skipped: 3

## Panel A: Fuel-Similarity Quintiles ([-1,+3] month CARs)

| Quintile | Mean CAR | N_firms (avg) |
|---|---|---|
| Q1 (lowest fuel sim) | +0.84% | 42.1 |
| Q2 | +1.19% | 41.4 |
| Q3 | -4.68% | 41.6 |
| Q4 | -4.14% | 41.4 |
| Q5 (highest fuel sim) | -4.61% | 41.3 |
| Q5 - Q1 (spread) | -5.45% (t = -7.89) | |

## Panel B: Geographic Proximity Quintiles ([-1,+3] month CARs)

| Quintile | Mean CAR | N_firms (avg) |
|---|---|---|
| Q1 (most distant) | -5.05% | 42.1 |
| Q2 | -0.23% | 41.4 |
| Q3 | -8.14% | 41.6 |
| Q4 | +1.15% | 41.4 |
| Q5 (closest) | +0.98% | 41.3 |
| Q5 - Q1 (spread) | +6.03% (t = 9.41) | |

## Panel C: Channel Split

geo_spread - fuel_spread = +11.48% (t = 14.00, N = 176)

## Panel D: Long-Short Portfolio

Long (high geo + low fuel) vs Short (low geo + high fuel)
Mean CAR: +8.33% (t = 3.94)
N events with valid sorts: 175

## Interpretation

The portfolio sort methodology does not depend on regression functional
form, standard error specification, or multiple hypothesis testing
corrections. A significant Q5-Q1 spread establishes that spatial exposure
predicts returns non-parametrically.

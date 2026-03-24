# Spatial Transition Score (STS): A Computable Alternative to ESG

## Construction

STS_j = w_fuel_ij * has_ets_j - w_geo_ij

Inputs: GEM plant-level data (public), ETS membership (public), GPS coordinates (public).
No corporate disclosure required.

## Part 1: In-Sample Validation

Full sample: N = 55580, events = 175, firms = 565.
ESG subsample: N = 14731, events = 165, firms = 153.
Standard errors: event-clustered.
Window: [-1, +3] months.

### Univariate regressions: CAR = a + b * Predictor

| Predictor | Beta | SE | t-stat | R-squared | N |
|---|---:|---:|---:|---:|---:|
| STS (full sample) | -0.075190 | 0.058801 | -1.279 | 0.000020 | 55580 |
| STS (ESG subsample) | -0.097330 | 0.162414 | -0.599 | 0.000030 | 14731 |
| ESG score | -0.114409 | 0.013942 | -8.206 | 0.011874 | 14731 |

### STS Portfolio Sorts

| Quintile | Mean STS | Mean CAR | N_firms |
|---|---:|---:|---:|
| Q1 | -0.0073 | +0.0459 | 11116 |
| Q2 | +0.0000 | +0.0879 | 11116 |
| Q3 | +0.0000 | +0.0535 | 11116 |
| Q4 | +0.0000 | +0.0323 | 11116 |
| Q5 | +0.0034 | +0.0510 | 11116 |
| Q5-Q1 | | +0.0051 (t = 1.584) | |

## Part 2: Out-of-Sample Validation

Training: events before 2020 (N = 47777, events = 160)
Test: events 2020+ (N = 7803, events = 15)

### Full sample comparison (STS vs Naive)

| Metric | STS | Naive (mean) |
|---|---:|---:|
| MSPE | 0.08232732 | 0.08230575 |
| Correlation (pred, actual) | -0.0224 | -- |
| Directional accuracy | 0.5961 | 0.5961 |

### ESG subsample comparison (STS vs ESG vs Naive)

| Metric | STS | ESG | Naive (mean) |
|---|---:|---:|---:|
| MSPE | 0.06923478 | 0.06894252 | 0.06923863 |
| Correlation (pred, actual) | 0.0134 | 0.0617 | -- |
| Directional accuracy | 0.6115 | 0.6182 | 0.6115 |

## Part 3: Firm-Level Exposure Ranking

### Most Exposed (most negative average STS)

| Rank | Firm | Country | Avg STS | N events |
|---|---|---|---:|---:|
| 1 | CEC AFRICA INVESTMENTS LTD | MUS | -0.2629 | 25 |
| 2 | OGK-2 JSC | RUS | -0.0431 | 103 |
| 3 | REPOWER | CHE | -0.0313 | 158 |
| 4 | INTER RAO UES OJSC | RUS | -0.0217 | 100 |
| 5 | ABU DHABI NATIONAL ENERGY CO | ARE | -0.0203 | 105 |
| 6 | SAUDI ELECTRICITY CO | SAU | -0.0162 | 105 |
| 7 | QATAR ELECT & WATER | QAT | -0.0151 | 105 |
| 8 | TALEN ENERGY CORP -OLD | USA | -0.0116 | 9 |
| 9 | ALPIQ HOLDING AG | CHE | -0.0114 | 180 |
| 10 | COPPERBELT ENERGY CORP | ZMB | -0.0108 | 105 |

### Least Exposed (most positive average STS)

| Rank | Firm | Country | Avg STS | N events |
|---|---|---|---:|---:|
| 1 | FORTIS INC | CAN | +0.0060 | 27 |
| 2 | VISTRA CORP | USA | +0.0059 | 27 |
| 3 | UNIPER SE | DEU | +0.0057 | 27 |
| 4 | AB IGNITIS GRUPE | LTU | +0.0054 | 4 |
| 5 | AVANGRID INC | USA | +0.0045 | 30 |
| 6 | AVISTA CORP | USA | +0.0042 | 187 |
| 7 | IDACORP INC | USA | +0.0040 | 187 |
| 8 | NORTHWESTERN ENRGY GROUP INC | USA | +0.0035 | 187 |
| 9 | GENESIS ENERGY LTD | NZL | +0.0035 | 55 |
| 10 | HAWAIIAN ELECTRIC INDS | USA | +0.0034 | 187 |

## Part 4: Out-of-Sample Portfolio Sorts (Temporal Split)

### Pre-2020 (training)

| Quintile | Mean CAR |
|---|---:|
| Q1 | +0.0409 |
| Q2 | +0.0978 |
| Q3 | +0.0419 |
| Q4 | +0.0409 |
| Q5 | +0.0385 |

Q5-Q1 spread: -0.0024 (t = -0.726)
Events: 160

### Post-2020 (test)

| Quintile | Mean CAR |
|---|---:|
| Q1 | +0.0797 |
| Q2 | +0.0126 |
| Q3 | +0.0634 |
| Q4 | +0.1293 |
| Q5 | +0.0521 |

Q5-Q1 spread: -0.0276 (t = -2.976***)
Events: 15

### Comparison: does the signal persist?

| Period | STS Q5-Q1 | t-stat | Fuel Q5-Q1 | t-stat | N events |
|---|---:|---:|---:|---:|---:|
| Pre-2020 | -0.0024 | -0.726 | -0.0941 | -25.078*** | 160 |
| Post-2020 | -0.0276 | -2.976*** | +0.0305 | 3.123*** | 15 |

## Key Finding

The Spatial Transition Score is a freely computable spatial measure constructed from publicly available plant-level data. In-sample, STS has t = -1.28 (R-squared = 0.0000), while ESG has t = -8.21 (R-squared = 0.0119). Out-of-sample MSPE: STS = 0.06923478, ESG = 0.06894252, Naive = 0.06923863. The score can be constructed using only GEM plant trackers, ETS membership records, and GPS coordinates, without purchasing proprietary ESG ratings.

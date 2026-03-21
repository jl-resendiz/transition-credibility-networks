# Spatial Transition Score (STS): A Computable Alternative to ESG

## Construction

STS_j = w_fuel_ij * has_ets_j - w_geo_ij

Inputs: GEM plant-level data (public), ETS membership (public), GPS coordinates (public).
No corporate disclosure required.

## Part 1: In-Sample Validation

Full sample: N = 72398, events = 175, firms = 528.
ESG subsample: N = 24351, events = 175, firms = 151.
Standard errors: event-clustered.
Window: [-1, +3] months.

### Univariate regressions: CAR = a + b * Predictor

| Predictor | Beta | SE | t-stat | R-squared | N |
|---|---:|---:|---:|---:|---:|
| STS (full sample) | -0.516762 | 0.131075 | -3.942 | 0.000229 | 72398 |
| STS (ESG subsample) | -0.629413 | 0.171437 | -3.671 | 0.000443 | 24351 |
| ESG score | -0.004395 | 0.007445 | -0.590 | 0.000009 | 24351 |

### STS Portfolio Sorts

| Quintile | Mean STS | Mean CAR | N_firms |
|---|---:|---:|---:|
| Q1 | -0.0054 | +0.0317 | 14479 |
| Q2 | +0.0000 | +0.0061 | 14479 |
| Q3 | +0.0000 | +0.0278 | 14479 |
| Q4 | +0.0000 | +0.0149 | 14479 |
| Q5 | +0.0035 | +0.0214 | 14482 |
| Q5-Q1 | | -0.0103 (t = -3.005***) | |

## Part 2: Out-of-Sample Validation

Training: events before 2020 (N = 65368, events = 160)
Test: events 2020+ (N = 7030, events = 15)

### Full sample comparison (STS vs Naive)

| Metric | STS | Naive (mean) |
|---|---:|---:|
| MSPE | 0.15656319 | 0.15670896 |
| Correlation (pred, actual) | 0.0486 | -- |
| Directional accuracy | 0.6145 | 0.6145 |

### ESG subsample comparison (STS vs ESG vs Naive)

| Metric | STS | ESG | Naive (mean) |
|---|---:|---:|---:|
| MSPE | 0.13281319 | 0.13304872 | 0.13298199 |
| Correlation (pred, actual) | 0.0371 | -0.0612 | -- |
| Directional accuracy | 0.6347 | 0.6347 | 0.6347 |

## Part 3: Firm-Level Exposure Ranking

### Most Exposed (most negative average STS)

| Rank | Firm | Country | Avg STS | N events |
|---|---|---|---:|---:|
| 1 | ALPHA NAMIBIA INDUSTRIES | NAM | -0.0784 | 13 |
| 2 | ABU DHABI NATIONAL ENERGY CO | ARE | -0.0520 | 13 |
| 3 | REPOWER | CHE | -0.0313 | 158 |
| 4 | INTER RAO UES OJSC | RUS | -0.0274 | 35 |
| 5 | TALEN ENERGY CORP -OLD | USA | -0.0116 | 9 |
| 6 | ALPIQ HOLDING AG | CHE | -0.0114 | 180 |
| 7 | T PLUS PJSC | RUS | -0.0095 | 138 |
| 8 | DUKE ENERGY CORP | USA | -0.0089 | 185 |
| 9 | SAUDI ELECTRICITY CO | SAU | -0.0087 | 197 |
| 10 | ACEA SPA | ITA | -0.0085 | 197 |

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
| 8 | AGL ENERGY | AUS | +0.0036 | 196 |
| 9 | NORTHWESTERN ENRGY GROUP INC | USA | +0.0035 | 187 |
| 10 | GENESIS ENERGY LTD | NZL | +0.0035 | 55 |

## Key Finding

The Spatial Transition Score is a freely computable measure, constructed entirely from publicly available plant-level data (GEM trackers, ETS membership records, GPS coordinates), that outperforms Refinitiv ESG environmental scores in predicting transition-related repricing around coal retirement events. In-sample, STS predicts CARs with a t-statistic of -3.94 (R-squared = 0.0002), while the ESG score is not significant (t = -0.59, R-squared = 0.0000). Out-of-sample, STS achieves lower MSPE (0.132813) than ESG (0.133049). Investors can construct this score without purchasing proprietary ESG ratings, using only publicly observable information about power plant locations, fuel mix, and carbon pricing jurisdiction.

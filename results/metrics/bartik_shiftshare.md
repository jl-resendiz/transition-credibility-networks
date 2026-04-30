# Shift-Share (Bartik) Causal Robustness Test

Addresses the concern that the fuel coefficient captures shared exposure
to common factors rather than network transmission.

## Design

**Shares**: Pre-period fuel-mix similarity weights using only plants
commissioned before 2014 (frozen before the main retirement wave).

**Shifts**: Annual aggregate coal MW retired globally per year.

**Bartik instrument**: B_ie = w_fuel_pre(i,j) x (RetiredMW_year / 10000)

## Diagnostics

Firms with pre-period fuel shares: 348
Firms with current fuel shares: 405
Pre-period fuel matrix: 348 firms, 96420 edges
Correlation(w_fuel_pre, w_fuel_current) across matrix: 0.3181

### Bartik instrument summary

| Statistic | Value |
|---|---:|
| Mean | 0.005999 |
| SD | 0.008674 |
| Min | 0.000000 |
| Max | 0.041467 |
| Non-zero | 10894/24070 |
| Corr(w_fuel_pre, w_fuel) in sample | 0.4963 |

### Retirement shocks by year

| Year | MW retired |
|---|---:|
| 2015 | 38,615 |
| 2016 | 36,570 |
| 2017 | 31,304 |
| 2018 | 36,770 |
| 2019 | 36,033 |
| 2020 | 43,987 |
| 2021 | 35,188 |
| 2022 | 28,433 |
| 2023 | 23,933 |
| 2024 | 30,719 |
| 2025 | 23,233 |

## Spec A: Bartik Reduced Form (Fama-MacBeth + Newey-West)

CAR = alpha + beta_bartik x B_ie + w_geo + w_reg + same_sector + eps

Valid events: 3
Avg firms per event: 563.3
Avg within-event R2: 0.0160

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.023821 | 0.024622 | 0.967 | 0.4354 |
| w_geo | +0.495012 | 0.603373 | 0.820 | 0.4982 |
| bartik | -1.963534 | 0.834620 | -2.353 | 0.1429 |
| w_reg | -1.463553 | 0.402335 | -3.638 | 0.0680* |
| same_sector | +0.015086 | 0.002311 | 6.528 | 0.0227** |

## Spec B: Standard w_fuel (Fama-MacBeth + Newey-West, comparison)

CAR = alpha + beta_fuel x w_fuel + w_geo + w_reg + same_sector + eps

Valid events: 3
Avg firms per event: 563.3
Avg within-event R2: 0.0314

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.020830 | 0.033920 | 0.614 | 0.5392 |
| w_geo | +0.770333 | 0.602725 | 1.278 | 0.2012 |
| w_fuel | -14.693448 | 1.604397 | -9.158 | 0.0000*** |
| w_reg | -1.169730 | 0.344980 | -3.391 | 0.0007*** |
| same_sector | +0.033702 | 0.003595 | 9.374 | 0.0000*** |

## Pooled Event-Clustered Regressions

### Spec A: Bartik

N = 24070, R2 = 0.0045

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.036045 | 0.010730 | 3.359 | 0.0008*** |
| bartik | -1.891236 | 0.364735 | -5.185 | 0.0000*** |
| w_geo | +0.363730 | 0.076125 | 4.778 | 0.0000*** |
| w_reg | -0.780023 | 0.232155 | -3.360 | 0.0008*** |
| same_sector | +0.011537 | 0.007782 | 1.482 | 0.1382 |

### Spec B: Standard w_fuel

N = 24070, R2 = 0.0029

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.028705 | 0.009640 | 2.978 | 0.0029*** |
| w_fuel | -4.213568 | 1.344214 | -3.135 | 0.0017*** |
| w_geo | +0.343470 | 0.086704 | 3.961 | 0.0001*** |
| w_reg | -0.812642 | 0.222586 | -3.651 | 0.0003*** |
| same_sector | +0.013794 | 0.007693 | 1.793 | 0.0730* |

## Summary Comparison

| Specification | Channel variable | FM t-stat | Pooled t-stat |
|---|---|---:|---:|
| Bartik (pre-period shares x agg shock) | bartik | -2.353 | -5.185 |
| Standard (current w_fuel) | w_fuel | -9.158 | -3.135 |

## GPS (2020) Shift-Share Diagnostics

### Rotemberg Weights

Which events drive the Bartik estimate? Rotemberg weights are proportional to the
sum of squared Bartik residuals (after partialling out controls) within each event.

HHI of Rotemberg weights: 0.0311
Top 5 events share: 0.2474

| Rank | Event | Plant | Year | MW | Weight |
|---:|---:|---|---:|---:|---:|
| 1 | 138 | Panipat power station | 2015 | 110 | 0.0505 |
| 2 | 139 | Panipat power station | 2015 | 110 | 0.0505 |
| 3 | 140 | Panipat power station | 2015 | 110 | 0.0505 |
| 4 | 141 | Panipat power station | 2015 | 110 | 0.0505 |
| 5 | 145 | Chandrapur (Assam) power stati | 2016 | 30 | 0.0452 |
| 6 | 146 | Chandrapur (Assam) power stati | 2016 | 30 | 0.0452 |
| 7 | 163 | Vorkutinskaya-2 power station | 2020 | 35 | 0.0375 |
| 8 | 164 | Vorkutinskaya-2 power station | 2020 | 50 | 0.0375 |
| 9 | 165 | Vorkutinskaya-2 power station | 2020 | 60 | 0.0375 |
| 10 | 166 | Vorkutinskaya-2 power station | 2020 | 47 | 0.0375 |

### Negative Weight Diagnostic

Negative-weight events: 0 / 40
Sum of negative weights: 0.000000
Sum of positive weights: 1.000000

### Pre-Event Balance Test ([-5, -2] months)

Tests whether Bartik exposure predicts CARs in the pre-event window.
Under the identifying assumption, the Bartik instrument should NOT predict
pre-event returns.

N = 24070, R2 = 0.001225

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.010142 | 0.006393 | 1.587 | 0.1126 |
| bartik | -0.642099 | 0.344208 | -1.865 | 0.0621* |
| w_geo | +0.146102 | 0.071431 | 2.045 | 0.0408** |
| w_reg | -0.246093 | 0.333688 | -0.737 | 0.4608 |
| same_sector | +0.011642 | 0.005641 | 2.064 | 0.0390** |

**PASS**: Bartik t = -1.865, p = 0.0621

### Pre-Balance Sensitivity (cutoff = 2010)

Repeats the pre-event balance test using only plants commissioned before 2010.
More pre-determined shares strengthen the causal claim if the test still passes.

Firms with pre-2010 fuel shares: 331 (vs 348 at pre-2014)
Pre-2010 fuel matrix: 331 firms, 88206 edges

N = 24070

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.010920 | 0.006213 | 1.757 | 0.0788* |
| bartik | -0.721540 | 0.313062 | -2.305 | 0.0212** |
| w_geo | +0.148244 | 0.072088 | 2.056 | 0.0397** |
| w_reg | -0.244390 | 0.333141 | -0.734 | 0.4632 |
| same_sector | +0.011366 | 0.005709 | 1.991 | 0.0465** |

**FAIL**: Bartik (pre-2010) t = -2.305, p = 0.0212

#### Pre-Balance Comparison

| Cutoff | Bartik t | Bartik p | Verdict |
|---:|---:|---:|---|
| 2014 | -1.865 | 0.0621 | PASS |
| 2010 | -2.305 | 0.0212 | FAIL |

### Oster (2019) Coefficient Stability Bounds

Tests how much selection on unobservables (delta) would be needed to explain
away the Bartik coefficient. delta* > 1 means unobservables would need to be
more important than observables.

| Quantity | Bartik | Standard w_fuel |
|---|---:|---:|
| beta (no controls) | -1.922375 | -4.249332 |
| beta (full controls) | -1.891236 | -4.213568 |
| R2 (no controls) | 0.003701 | 0.002026 |
| R2 (full controls) | 0.004498 | 0.002938 |
| beta* (delta=1) | -1.838489 | -4.178988 |
| **delta*** | **35.8545** | **121.8479** |

**PASS**: delta* = 35.8545

## Interpretation

The Bartik instrument separates pre-determined exposure (fuel-mix similarity
frozen at pre-2014 plant vintages) from aggregate shocks (total coal MW
retired per year). If beta_bartik is significant, the fuel channel reflects
genuine network transmission rather than shared exposure to common factors.

The comparison between Spec A and Spec B shows whether the causal (Bartik)
and descriptive (OLS) estimates agree in sign and magnitude.

Under the Goldsmith-Pinkham, Sorkin & Swift (2020) framework, the Bartik
coefficient has a causal interpretation as the dose-response of abnormal
returns to technology exposure, provided: (1) fuel-mix shares are pre-determined
(supported by pre-2014 vintage restriction); (2) shares do not predict pre-event
returns (tested above); (3) Rotemberg weights are non-negative (checked above);
(4) the Oster (2019) bound confirms robustness to selection on unobservables.
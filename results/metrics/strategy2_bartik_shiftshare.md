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
Avg within-event R2: 0.0158

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.023943 | 0.024554 | 0.975 | 0.3295 |
| w_geo | +0.529782 | 0.518486 | 1.022 | 0.3069 |
| bartik | -1.977274 | 0.853558 | -2.317 | 0.0205** |
| w_reg | -1.432015 | 0.398365 | -3.595 | 0.0003*** |
| same_sector | +0.014932 | 0.002166 | 6.894 | 0.0000*** |

## Spec B: Standard w_fuel (Fama-MacBeth + Newey-West, comparison)

CAR = alpha + beta_fuel x w_fuel + w_geo + w_reg + same_sector + eps

Valid events: 3
Avg firms per event: 563.3
Avg within-event R2: 0.0312

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.020836 | 0.033908 | 0.614 | 0.5389 |
| w_geo | +0.736243 | 0.511382 | 1.440 | 0.1499 |
| w_fuel | -14.710990 | 1.676663 | -8.774 | 0.0000*** |
| w_reg | -1.119956 | 0.342836 | -3.267 | 0.0011*** |
| same_sector | +0.033655 | 0.003539 | 9.509 | 0.0000*** |

## Pooled Event-Clustered Regressions

### Spec A: Bartik

N = 24070, R2 = 0.0043

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.036221 | 0.010718 | 3.379 | 0.0007*** |
| bartik | -1.880015 | 0.364328 | -5.160 | 0.0000*** |
| w_geo | +0.165618 | 0.043799 | 3.781 | 0.0002*** |
| w_reg | -0.749760 | 0.210528 | -3.561 | 0.0004*** |
| same_sector | +0.011565 | 0.007782 | 1.486 | 0.1372 |

### Spec B: Standard w_fuel

N = 24070, R2 = 0.0028

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.028914 | 0.009626 | 3.004 | 0.0027*** |
| w_fuel | -4.188841 | 1.345616 | -3.113 | 0.0019*** |
| w_geo | +0.154240 | 0.049963 | 3.087 | 0.0020*** |
| w_reg | -0.783553 | 0.203255 | -3.855 | 0.0001*** |
| same_sector | +0.013807 | 0.007690 | 1.795 | 0.0726* |

## Summary Comparison

| Specification | Channel variable | FM t-stat | Pooled t-stat |
|---|---|---:|---:|
| Bartik (pre-period shares x agg shock) | bartik | -2.317 | -5.160 |
| Standard (current w_fuel) | w_fuel | -8.774 | -3.113 |

## GPS (2020) Shift-Share Diagnostics

### Rotemberg Weights

Which events drive the Bartik estimate? Rotemberg weights are proportional to the
sum of squared Bartik residuals (after partialling out controls) within each event.

HHI of Rotemberg weights: 0.0311
Top 5 events share: 0.2476

| Rank | Event | Plant | Year | MW | Weight |
|---:|---:|---|---:|---:|---:|
| 1 | 138 | Panipat power station | 2015 | 110 | 0.0506 |
| 2 | 139 | Panipat power station | 2015 | 110 | 0.0506 |
| 3 | 140 | Panipat power station | 2015 | 110 | 0.0506 |
| 4 | 141 | Panipat power station | 2015 | 110 | 0.0506 |
| 5 | 145 | Chandrapur (Assam) power stati | 2016 | 30 | 0.0453 |
| 6 | 146 | Chandrapur (Assam) power stati | 2016 | 30 | 0.0453 |
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

N = 24070, R2 = 0.001351

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.010086 | 0.006368 | 1.584 | 0.1132 |
| bartik | -0.646513 | 0.345754 | -1.870 | 0.0615* |
| w_geo | +0.172461 | 0.045667 | 3.776 | 0.0002*** |
| w_reg | -0.249264 | 0.333972 | -0.746 | 0.4554 |
| same_sector | +0.011693 | 0.005649 | 2.070 | 0.0385** |

**PASS**: Bartik t = -1.870, p = 0.0615

### Pre-Balance Sensitivity (cutoff = 2010)

Repeats the pre-event balance test using only plants commissioned before 2010.
More pre-determined shares strengthen the causal claim if the test still passes.

Firms with pre-2010 fuel shares: 331 (vs 348 at pre-2014)
Pre-2010 fuel matrix: 331 firms, 88206 edges

N = 24070

| Variable | beta | SE(cl) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.010859 | 0.006187 | 1.755 | 0.0792* |
| bartik | -0.725037 | 0.314402 | -2.306 | 0.0211** |
| w_geo | +0.173359 | 0.046216 | 3.751 | 0.0002*** |
| w_reg | -0.247398 | 0.333389 | -0.742 | 0.4580 |
| same_sector | +0.011418 | 0.005717 | 1.997 | 0.0458** |

**FAIL**: Bartik (pre-2010) t = -2.306, p = 0.0211

#### Pre-Balance Comparison

| Cutoff | Bartik t | Bartik p | Verdict |
|---:|---:|---:|---|
| 2014 | -1.870 | 0.0615 | PASS |
| 2010 | -2.306 | 0.0211 | FAIL |

### Oster (2019) Coefficient Stability Bounds

Tests how much selection on unobservables (delta) would be needed to explain
away the Bartik coefficient. delta* > 1 means unobservables would need to be
more important than observables.

| Quantity | Bartik | Standard w_fuel |
|---|---:|---:|
| beta (no controls) | -1.922375 | -4.249332 |
| beta (full controls) | -1.880015 | -4.188841 |
| R2 (no controls) | 0.003701 | 0.002026 |
| R2 (full controls) | 0.004307 | 0.002764 |
| beta* (delta=1) | -1.789627 | -4.120860 |
| **delta*** | **20.7993** | **61.6174** |

**PASS**: delta* = 20.7993

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
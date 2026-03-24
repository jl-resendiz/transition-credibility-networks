# Geographic Diversification Test

Tests whether multinational diversification attenuates the geographic
competitive channel (w_geo). If geographic effects are real at the plant
level, single-country utilities should show a stronger geo coefficient.

Events: 179 first-mover coal retirements
Window: [-1, +3] months, vwretd market-adjusted returns
Minimum obs per event: 20
Single-country definition: >=90% of MW in one country

## Firm-Level Diversification Summary

- Firms with plant data: 414
- Single-country firms: 301
- Multi-country firms: 113
- HHI (geographic concentration): mean = 0.8516, median = 1.0000

| n_countries | Firms |
|---:|---:|
| 1 | 249 |
| 2 | 44 |
| 3 | 20 |
| 4 | 13 |
| 5 | 11 |
| 6 | 11 |
| 7 | 8 |
| 8 | 4 |
| 9 | 3 |
| 10 | 16 |
| 11 | 5 |
| 12 | 11 |
| 13 | 4 |
| 14 | 3 |
| 15 | 1 |
| 16 | 2 |
| 20 | 6 |
| 22 | 1 |
| 27 | 1 |
| 31 | 1 |

## Specification 1: Full Sample Baseline

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + e

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0515

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4743 | (0.7301) | -7.50 | 0.000*** | -4.7824 | (0.6496) | -7.36 | 0.000*** |
| w_geo                  | -0.0230 | (0.0566) | -0.41 | 0.685 | -0.6072 | (0.4681) | -1.30 | 0.195 |
| w_reg                  | +1.4525 | (1.0515) | 1.38 | 0.167 | +2.6428 | (0.9614) | 2.75 | 0.006*** |
| same_sector            | +0.0332 | (0.0089) | 3.74 | 0.000*** | +0.0215 | (0.0112) | 1.91 | 0.056* |

## Specification 2: Single-Country Subsample

Same specification, restricted to firms with >=90% MW in one country.

Events with >= 20 single-country obs: 175
Total single-country obs: 49308

N = 49308, Events = 116, R2(pooled) = 0.0071, R2(FM avg) = 0.0506

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.6784 | (0.7806) | -7.27 | 0.000*** | -2.6615 | (1.1310) | -2.35 | 0.019** |
| w_geo                  | +0.0537 | (0.0926) | 0.58 | 0.562 | -1.2525 | (1.2022) | -1.04 | 0.298 |
| w_reg                  | +1.3466 | (1.0770) | 1.25 | 0.211 | +2.2767 | (0.9516) | 2.39 | 0.017** |
| same_sector            | +0.0330 | (0.0090) | 3.65 | 0.000*** | +0.0224 | (0.0117) | 1.91 | 0.056* |

### Multi-Country Subsample (for comparison)

N = 5751, Events = 36, R2(pooled) = 0.0065, R2(FM avg) = 0.1407

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -4.3420 | (1.1524) | -3.77 | 0.000*** | -7.7130 | (1.8233) | -4.23 | 0.000*** |
| w_geo                  | -0.0327 | (0.0314) | -1.04 | 0.299 | +1252.9092 | (1320.3892) | 0.95 | 0.343 |
| w_reg                  | +2.2737 | (1.4035) | 1.62 | 0.105 | +6.2270 | (3.6664) | 1.70 | 0.089* |
| same_sector            | +0.0312 | (0.0099) | 3.14 | 0.002*** | +0.0661 | (0.0132) | 4.99 | 0.000*** |

## Specification 3a: Diversification Interaction (n_countries)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x n_countries) + e

If b2 < 0 and b5 > 0: geo benefit exists but weakens with diversification

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0565

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4825 | (0.7304) | -7.51 | 0.000*** | -4.8317 | (0.6679) | -7.23 | 0.000*** |
| w_geo                  | +0.0466 | (0.0965) | 0.48 | 0.629 | +48.9767 | (49.0926) | 1.00 | 0.318 |
| w_reg                  | +1.4440 | (1.0550) | 1.37 | 0.171 | +2.6300 | (0.9747) | 2.70 | 0.007*** |
| same_sector            | +0.0332 | (0.0089) | 3.75 | 0.000*** | +0.0217 | (0.0113) | 1.93 | 0.053* |
| w_geo_x_nc             | -0.0240 | (0.0184) | -1.31 | 0.191 | -49.5760 | (49.0919) | -1.01 | 0.313 |

## Specification 3b: Diversification Interaction (log n_countries)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x log(n_countries)) + e

Log version: captures concave attenuation (first additional country matters most)

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0567

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4760 | (0.7302) | -7.50 | 0.000*** | -4.9203 | (0.6757) | -7.28 | 0.000*** |
| w_geo                  | -0.0153 | (0.0927) | -0.17 | 0.869 | -1.1298 | (0.8568) | -1.32 | 0.187 |
| w_reg                  | +1.4510 | (1.0577) | 1.37 | 0.170 | +2.7084 | (0.9998) | 2.71 | 0.007*** |
| same_sector            | +0.0332 | (0.0089) | 3.74 | 0.000*** | +0.0219 | (0.0113) | 1.94 | 0.053* |
| w_geo_x_log_nc         | -0.0106 | (0.0700) | -0.15 | 0.880 | -12.6477 | (12.5792) | -1.01 | 0.315 |

## Specification 3c: Diversification Interaction (HHI)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x HHI) + e

HHI = sum(share_k^2) where share_k = MW_in_country_k / total_MW.
HHI = 1 for single-country firms; HHI close to 0 for diversified firms.

Prediction (Lemma 2): b5 > 0 (higher concentration = stronger geo effect).
Equivalently: b_coeff on w_geo x (1 - HHI) should be negative.

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0562

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4902 | (0.7292) | -7.53 | 0.000*** | -4.7535 | (0.6453) | -7.37 | 0.000*** |
| w_geo                  | -0.2068 | (0.1129) | -1.83 | 0.067* | +195.4096 | (201.1883) | 0.97 | 0.331 |
| w_reg                  | +1.4409 | (1.0536) | 1.37 | 0.171 | +2.6553 | (0.9967) | 2.66 | 0.008*** |
| same_sector            | +0.0332 | (0.0089) | 3.75 | 0.000*** | +0.0223 | (0.0113) | 1.97 | 0.048** |
| w_geo_x_hhi            | +0.2335 | (0.1788) | 1.31 | 0.191 | -196.3137 | (201.1911) | -0.98 | 0.329 |

## Interpretation

- Full sample w_geo: beta = -0.607185, t = -1.297, p = 0.1945
- Single-country w_geo: beta = -1.252475, t = -1.042, p = 0.2975
- w_geo is NOT stronger in the single-country subsample.
- Interaction (w_geo x n_countries): beta = -49.575980, t = -1.010, p = 0.3126
- Interaction (w_geo x log(n_countries)): beta = -12.647735, t = -1.005, p = 0.3147
- Interaction (w_geo x HHI): beta = -196.313680, t = -0.976, p = 0.3292
  Pattern: negative HHI interaction -> geo effect does NOT strengthen with concentration (inconsistent with Lemma 2).

## Diagnostics

- Mean w_geo (single-country): 0.001367
- Mean w_geo (multi-country): 0.003823
- Mean CAR (single-country): 0.054167
- Mean CAR (multi-country): 0.053966
- Mean n_countries across obs: 2.07
- Mean HHI across obs: 0.9389
- HHI range: [0.0788, 1.0000]
- HHI = 1 (single-country): 45412 obs (81.7%)

# Geographic Diversification Test

Tests whether multinational diversification attenuates the geographic
competitive channel (w_geo). If geographic effects are real at the plant
level, single-country utilities should show a stronger geo coefficient.

Events: 179 first-mover-matched (175 used in pooled regression below; 117 with ≥20 firms qualify for FM)
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

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0518

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4884 | (0.7279) | -7.54 | 0.000*** | -4.7656 | (0.6508) | -7.32 | 0.000*** |
| w_geo                  | +0.0176 | (0.1007) | 0.17 | 0.861 | -0.5427 | (0.3090) | -1.76 | 0.079* |
| w_reg                  | +1.4412 | (1.0503) | 1.37 | 0.170 | +2.6975 | (0.9518) | 2.83 | 0.005*** |
| same_sector            | +0.0332 | (0.0089) | 3.74 | 0.000*** | +0.0215 | (0.0112) | 1.92 | 0.055* |

## Specification 2: Single-Country Subsample

Same specification, restricted to firms with >=90% MW in one country.

Events with >= 20 single-country obs: 175
Total single-country obs: 49308

N = 49308, Events = 116, R2(pooled) = 0.0072, R2(FM avg) = 0.0506

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.7192 | (0.7760) | -7.37 | 0.000*** | -2.6483 | (1.1303) | -2.34 | 0.019** |
| w_geo                  | +0.1772 | (0.1767) | 1.00 | 0.316 | -0.6744 | (0.6731) | -1.00 | 0.316 |
| w_reg                  | +1.3079 | (1.0727) | 1.22 | 0.223 | +2.2932 | (0.9323) | 2.46 | 0.014** |
| same_sector            | +0.0330 | (0.0091) | 3.65 | 0.000*** | +0.0224 | (0.0117) | 1.91 | 0.056* |

### Multi-Country Subsample (for comparison)

N = 5751, Events = 36, R2(pooled) = 0.0065, R2(FM avg) = 0.1455

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -4.3416 | (1.1524) | -3.77 | 0.000*** | -7.7472 | (1.8367) | -4.22 | 0.000*** |
| w_geo                  | -0.0610 | (0.0507) | -1.20 | 0.230 | +5.7886 | (31.0664) | 0.19 | 0.852 |
| w_reg                  | +2.2777 | (1.4032) | 1.62 | 0.105 | +7.0183 | (3.5702) | 1.97 | 0.049** |
| same_sector            | +0.0312 | (0.0099) | 3.14 | 0.002*** | +0.0656 | (0.0136) | 4.83 | 0.000*** |

## Specification 3a: Diversification Interaction (n_countries)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x n_countries) + e

If b2 < 0 and b5 > 0: geo benefit exists but weakens with diversification

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0577

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4987 | (0.7280) | -7.55 | 0.000*** | -4.6200 | (0.7023) | -6.58 | 0.000*** |
| w_geo                  | +0.1638 | (0.1684) | 0.97 | 0.331 | -2.5698 | (2.9383) | -0.87 | 0.382 |
| w_reg                  | +1.4232 | (1.0546) | 1.35 | 0.177 | +2.6518 | (0.9611) | 2.76 | 0.006*** |
| same_sector            | +0.0333 | (0.0089) | 3.75 | 0.000*** | +0.0218 | (0.0112) | 1.94 | 0.053* |
| w_geo_x_nc             | -0.0491 | (0.0277) | -1.77 | 0.076* | +2.2517 | (2.8905) | 0.78 | 0.436 |

## Specification 3b: Diversification Interaction (log n_countries)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x log(n_countries)) + e

Log version: captures concave attenuation (first additional country matters most)

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0573

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.4998 | (0.7276) | -7.56 | 0.000*** | -4.8273 | (0.6924) | -6.97 | 0.000*** |
| w_geo                  | +0.0871 | (0.1749) | 0.50 | 0.618 | -0.6852 | (0.4707) | -1.46 | 0.145 |
| w_reg                  | +1.4277 | (1.0588) | 1.35 | 0.178 | +2.7322 | (0.9894) | 2.76 | 0.006*** |
| same_sector            | +0.0333 | (0.0089) | 3.75 | 0.000*** | +0.0219 | (0.0113) | 1.94 | 0.053* |
| w_geo_x_log_nc         | -0.0937 | (0.1190) | -0.79 | 0.431 | +3.7663 | (4.6864) | 0.80 | 0.422 |

## Specification 3c: Diversification Interaction (HHI)

CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + b5 (w_geo x HHI) + e

HHI = sum(share_k^2) where share_k = MW_in_country_k / total_MW.
HHI = 1 for single-country firms; HHI close to 0 for diversified firms.

Prediction (Lemma 2): b5 > 0 (higher concentration = stronger geo effect).
Equivalently: b_coeff on w_geo x (1 - HHI) should be negative.

N = 55580, Events = 117, R2(pooled) = 0.0071, R2(FM avg) = 0.0564

| Variable | beta(OLS) | SE(cl) | t | p | beta(FM) | SE(NW) | t | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| w_fuel                 | -5.5168 | (0.7258) | -7.60 | 0.000*** | -4.6883 | (0.6614) | -7.09 | 0.000*** |
| w_geo                  | -0.4002 | (0.1657) | -2.42 | 0.016** | -3.2322 | (7.2924) | -0.44 | 0.658 |
| w_reg                  | +1.4137 | (1.0524) | 1.34 | 0.179 | +2.6750 | (0.9896) | 2.70 | 0.007*** |
| same_sector            | +0.0333 | (0.0089) | 3.75 | 0.000*** | +0.0222 | (0.0112) | 1.97 | 0.048** |
| w_geo_x_hhi            | +0.5351 | (0.3060) | 1.75 | 0.080* | +2.7164 | (7.3761) | 0.37 | 0.713 |

## Interpretation

- Full sample w_geo: beta = -0.542675, t = -1.756, p = 0.0791
- Single-country w_geo: beta = -0.674446, t = -1.002, p = 0.3163
- w_geo is NOT stronger in the single-country subsample.
- Interaction (w_geo x n_countries): beta = +2.251721, t = 0.779, p = 0.4360
  Pattern: negative base geo + positive interaction -> geo benefit weakens with more countries.
- Interaction (w_geo x log(n_countries)): beta = +3.766300, t = 0.804, p = 0.4216
- Interaction (w_geo x HHI): beta = +2.716445, t = 0.368, p = 0.7127
  Pattern: positive HHI interaction -> geo effect strengthens with geographic concentration (consistent with Lemma 2).

## Diagnostics

- Mean w_geo (single-country): 0.001331
- Mean w_geo (multi-country): 0.003502
- Mean CAR (single-country): 0.054167
- Mean CAR (multi-country): 0.053966
- Mean n_countries across obs: 2.07
- Mean HHI across obs: 0.9389
- HHI range: [0.0788, 1.0000]
- HHI = 1 (single-country): 45412 obs (81.7%)

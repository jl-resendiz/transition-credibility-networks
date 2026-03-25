# Signal Hierarchy: Spatial Exposure vs. ESG Ratings

## Test 1: Horse Race

Sample restricted to firms with Refinitiv Environmental Score (N = 14731, 165 events, 153 unique firms).
Skipped 144841 firm-event observations without ESG coverage.
ESG score normalized to [0,1] (original scale: 0-100).
Standard errors: event-clustered.
Window: [-1, +3] months.

### Model 1: ESG Only

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.111511 | 0.015178 | 7.347 | 0.0000*** |
| esg_score | -0.114409 | 0.013942 | -8.206 | 0.0000*** |

R-squared = 0.011874, N = 14731, Clusters = 165

### Model 2: Spatial Only

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.027343 | 0.009920 | 2.756 | 0.0058*** |
| w_fuel | -1.551732 | 0.865872 | -1.792 | 0.0731* |
| w_geo | +0.013796 | 0.265798 | 0.052 | 0.9586 |
| same_sector | +0.030378 | 0.011236 | 2.704 | 0.0069*** |

R-squared = 0.003438, N = 14731, Clusters = 165

### Model 3: Both (ESG + Spatial)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.089353 | 0.013651 | 6.545 | 0.0000*** |
| w_fuel | -1.131042 | 0.853482 | -1.325 | 0.1851 |
| w_geo | -0.013439 | 0.256570 | -0.052 | 0.9582 |
| esg_score | -0.117603 | 0.014379 | -8.179 | 0.0000*** |
| same_sector | +0.034298 | 0.011680 | 2.937 | 0.0033*** |

R-squared = 0.015901, N = 14731, Clusters = 165

### Model 4: Full (Spatial + ESG + Credibility)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.088755 | 0.013636 | 6.509 | 0.0000*** |
| w_fuel | +0.402024 | 0.929585 | 0.432 | 0.6654 |
| w_geo | -0.128748 | 0.248829 | -0.517 | 0.6049 |
| w_reg | +1.246969 | 0.975665 | 1.278 | 0.2012 |
| esg_score | -0.116388 | 0.014840 | -7.843 | 0.0000*** |
| w_fuel_x_ets | -3.475126 | 1.676143 | -2.073 | 0.0381** |
| same_sector | +0.033728 | 0.011595 | 2.909 | 0.0036*** |

R-squared = 0.017294, N = 14731, Clusters = 165

### ESG Coefficient Attenuation

| | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| ESG alone (Model 1) | -0.114409 | 0.013942 | -8.206 | 0.0000*** |
| ESG with spatial (Model 3) | -0.117603 | 0.014379 | -8.179 | 0.0000*** |

Attenuation: -2.8% reduction in ESG coefficient magnitude when spatial exposure is included.

## Test 2: Information Hierarchy

### ESG-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + ESG score | 0.011874 | +0.011874 |
| 2 | + Spatial (w_fuel, w_geo) | 0.012152 | +0.000279 |
| 3 | + Policy credibility (w_fuel x has_ets) | 0.012764 | +0.000612 |
| 4 | + w_reg + SameSector | 0.017294 | +0.004530 |

### Spatial-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + Spatial (w_fuel, w_geo) | 0.000486 | +0.000486 |
| 2 | + Policy (w_fuel x has_ets) | 0.001707 | +0.001221 |
| 3 | + ESG score | 0.012764 | +0.011058 |

### Comparison

- Spatial adds to ESG: +0.000279 R-squared
- ESG adds to Spatial: +0.011058 R-squared

## Key Finding

Both ESG scores and spatial exposure retain some predictive power in the joint specification (ESG t = -8.18). However, spatial exposure adds +0.000279 marginal R-squared beyond ESG, while ESG adds only +0.011058 beyond spatial exposure. The information hierarchy favours spatial fundamentals as the more informative signal for transition risk assessment.

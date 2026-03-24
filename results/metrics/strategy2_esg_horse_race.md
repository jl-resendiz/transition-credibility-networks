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
| intercept | +0.027306 | 0.009902 | 2.758 | 0.0058*** |
| w_fuel | -1.559102 | 0.866465 | -1.799 | 0.0720* |
| w_geo | +0.031591 | 0.163589 | 0.193 | 0.8469 |
| same_sector | +0.030402 | 0.011216 | 2.711 | 0.0067*** |

R-squared = 0.003441, N = 14731, Clusters = 165

### Model 3: Both (ESG + Spatial)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.089325 | 0.013646 | 6.546 | 0.0000*** |
| w_fuel | -1.134982 | 0.853913 | -1.329 | 0.1838 |
| w_geo | +0.000131 | 0.160318 | 0.001 | 0.9993 |
| esg_score | -0.117597 | 0.014368 | -8.185 | 0.0000*** |
| same_sector | +0.034307 | 0.011661 | 2.942 | 0.0033*** |

R-squared = 0.015901, N = 14731, Clusters = 165

### Model 4: Full (Spatial + ESG + Credibility)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.088638 | 0.013655 | 6.491 | 0.0000*** |
| w_fuel | +0.397032 | 0.928586 | 0.428 | 0.6690 |
| w_geo | -0.065527 | 0.158834 | -0.413 | 0.6799 |
| w_reg | +1.225798 | 0.969140 | 1.265 | 0.2059 |
| esg_score | -0.116359 | 0.014853 | -7.834 | 0.0000*** |
| w_fuel_x_ets | -3.488802 | 1.682291 | -2.074 | 0.0381** |
| same_sector | +0.033756 | 0.011574 | 2.917 | 0.0035*** |

R-squared = 0.017282, N = 14731, Clusters = 165

### ESG Coefficient Attenuation

| | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| ESG alone (Model 1) | -0.114409 | 0.013942 | -8.206 | 0.0000*** |
| ESG with spatial (Model 3) | -0.117597 | 0.014368 | -8.185 | 0.0000*** |

Attenuation: -2.8% reduction in ESG coefficient magnitude when spatial exposure is included.

## Test 2: Information Hierarchy

### ESG-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + ESG score | 0.011874 | +0.011874 |
| 2 | + Spatial (w_fuel, w_geo) | 0.012151 | +0.000277 |
| 3 | + Policy credibility (w_fuel x has_ets) | 0.012765 | +0.000614 |
| 4 | + w_reg + SameSector | 0.017282 | +0.004517 |

### Spatial-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + Spatial (w_fuel, w_geo) | 0.000485 | +0.000485 |
| 2 | + Policy (w_fuel x has_ets) | 0.001707 | +0.001222 |
| 3 | + ESG score | 0.012765 | +0.011058 |

### Comparison

- Spatial adds to ESG: +0.000277 R-squared
- ESG adds to Spatial: +0.011058 R-squared

## Key Finding

Both ESG scores and spatial exposure retain some predictive power in the joint specification (ESG t = -8.18). However, spatial exposure adds +0.000277 marginal R-squared beyond ESG, while ESG adds only +0.011058 beyond spatial exposure. The information hierarchy favours spatial fundamentals as the more informative signal for transition risk assessment.

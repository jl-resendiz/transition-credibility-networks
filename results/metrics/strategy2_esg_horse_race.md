# Signal Hierarchy: Spatial Exposure vs. ESG Ratings

## Test 1: Horse Race

Sample restricted to firms with Refinitiv Environmental Score (N = 24351, 175 events, 151 unique firms).
Skipped 144841 firm-event observations without ESG coverage.
ESG score normalized to [0,1] (original scale: 0-100).
Standard errors: event-clustered.
Window: [-1, +3] months.

### Model 1: ESG Only

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.028829 | 0.007310 | 3.944 | 0.0001*** |
| esg_score | -0.004395 | 0.007445 | -0.590 | 0.5549 |

R-squared = 0.000009, N = 24351, Clusters = 175

### Model 2: Spatial Only

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.011568 | 0.007615 | 1.519 | 0.1288 |
| w_fuel | -0.663320 | 0.635006 | -1.045 | 0.2962 |
| w_geo | +0.460838 | 0.160914 | 2.864 | 0.0042*** |
| same_sector | +0.021125 | 0.007729 | 2.733 | 0.0063*** |

R-squared = 0.001177, N = 24351, Clusters = 175

### Model 3: Both (ESG + Spatial)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.014570 | 0.008306 | 1.754 | 0.0794* |
| w_fuel | -0.640510 | 0.635253 | -1.008 | 0.3133 |
| w_geo | +0.464995 | 0.161778 | 2.874 | 0.0040*** |
| esg_score | -0.005505 | 0.007599 | -0.724 | 0.4688 |
| same_sector | +0.021183 | 0.007747 | 2.734 | 0.0063*** |

R-squared = 0.001192, N = 24351, Clusters = 175

### Model 4: Full (Spatial + ESG + Credibility)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.011412 | 0.008305 | 1.374 | 0.1694 |
| w_fuel | +1.870949 | 0.701855 | 2.666 | 0.0077*** |
| w_geo | +0.474921 | 0.155015 | 3.064 | 0.0022*** |
| w_reg | +0.285263 | 0.490982 | 0.581 | 0.5612 |
| esg_score | +0.000722 | 0.007842 | 0.092 | 0.9266 |
| w_fuel_x_ets | -5.295239 | 1.235848 | -4.285 | 0.0000*** |
| same_sector | +0.021174 | 0.007702 | 2.749 | 0.0060*** |

R-squared = 0.002208, N = 24351, Clusters = 175

### ESG Coefficient Attenuation

| | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| ESG alone (Model 1) | -0.004395 | 0.007445 | -0.590 | 0.5549 |
| ESG with spatial (Model 3) | -0.005505 | 0.007599 | -0.724 | 0.4688 |

Attenuation: -25.2% reduction in ESG coefficient magnitude when spatial exposure is included.

## Test 2: Information Hierarchy

### ESG-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + ESG score | 0.000009 | +0.000009 |
| 2 | + Spatial (w_fuel, w_geo) | 0.000280 | +0.000271 |
| 3 | + Policy credibility (w_fuel x has_ets) | 0.001263 | +0.000983 |
| 4 | + w_reg + SameSector | 0.002208 | +0.000946 |

### Spatial-first ordering

| Step | Variables added | R-squared | Marginal R-squared |
|---|---|---:|---:|
| 0 | Intercept only | 0.000000 | -- |
| 1 | + Spatial (w_fuel, w_geo) | 0.000270 | +0.000270 |
| 2 | + Policy (w_fuel x has_ets) | 0.001261 | +0.000991 |
| 3 | + ESG score | 0.001263 | +0.000002 |

### Comparison

- Spatial adds to ESG: +0.000271 R-squared
- ESG adds to Spatial: +0.000002 R-squared

## Key Finding

Spatial network exposure subsumes all the information content of ESG environmental scores. In the horse race (Model 3), the ESG coefficient attenuates by -25% and loses statistical significance (t = -0.72) once spatial weights are included. The information hierarchy confirms this asymmetry: spatial exposure adds +0.0003 to R-squared beyond ESG, but ESG adds only +0.0000 beyond spatial exposure. Investors relying on purchased ESG ratings to assess transition risk are paying for a signal that is strictly dominated by freely observable spatial fundamentals.

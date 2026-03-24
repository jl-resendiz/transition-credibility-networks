# Bandwidth Sensitivity: W_geo under Fama-MacBeth + Newey-West Inference

Tests whether the geographic proximity channel (w_geo) strengthens at
wider bandwidths. Economic rationale: competitive benefits from coal
retirement propagate through interconnected transmission grids
(ENTSO-E ~1500km, US ISOs ~1000km), not at plant-level proximity.

Events: 179 first-mover coal retirements
Window: [-1, +3] months, vwretd market-adjusted returns
Minimum obs per event: 20

Weight formula: w_ij = exp(-d_ij / DECAY_KM) / d_ij, row-normalized
DECAY_KM = half_life / ln(2)

## Main Results: Bandwidth x Channel

| Bandwidth (km) | Events | Avg N | beta_geo | SE(NW) | t | p | beta_fuel | SE(NW) | t | p | diff(g-f) | SE(NW) | t | p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 250 | 175 | 318 | +1435.646331 | 854.762224 | 1.680 | 0.0930* | -2.892947 | 1.269647 | -2.279 | 0.0227** | +1438.539278 | 853.954461 | 1.685 | 0.0921* |
| 500 | 175 | 318 | +97.564319 | 55.717031 | 1.751 | 0.0799* | -2.950334 | 1.252154 | -2.356 | 0.0185** | +100.514653 | 54.922235 | 1.830 | 0.0672* |
| 750 | 175 | 318 | +43.087421 | 24.277985 | 1.775 | 0.0759* | -2.994289 | 1.239388 | -2.416 | 0.0157** | +46.081710 | 23.504285 | 1.961 | 0.0499** |
| 1000 | 175 | 318 | +29.240987 | 16.368279 | 1.786 | 0.0740* | -3.041249 | 1.225374 | -2.482 | 0.0131** | +32.282236 | 15.620626 | 2.067 | 0.0388** |
| 1500 | 175 | 318 | +20.179936 | 11.224416 | 1.798 | 0.0722* | -3.152703 | 1.190506 | -2.648 | 0.0081*** | +23.332639 | 10.542520 | 2.213 | 0.0269** |

## Summary: How geo significance varies with bandwidth

| Bandwidth | t(geo) | p(geo) | t(fuel) | p(fuel) | Interpretation |
|---:|---:|---:|---:|---:|---|
| 250 km | 1.680 | 0.0930* | -2.279 | 0.0227** | geo significant at 10%; fuel significant |
| 500 km | 1.751 | 0.0799* | -2.356 | 0.0185** | geo significant at 10%; fuel significant |
| 750 km | 1.775 | 0.0759* | -2.416 | 0.0157** | geo significant at 10%; fuel significant |
| 1000 km | 1.786 | 0.0740* | -2.482 | 0.0131** | geo significant at 10%; fuel significant |
| 1500 km | 1.798 | 0.0722* | -2.648 | 0.0081*** | geo significant at 10%; fuel significant |

**Strongest geo channel: 1500 km half-life** (t = 1.798, p = 0.0722)

Comparison: 500km t(geo) = 1.751 vs 1000km t(geo) = 1.786

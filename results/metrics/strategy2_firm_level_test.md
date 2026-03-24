# Firm-Level Collapsed Specification and Corrected F-Tests

This analysis addresses two issues:
1. The pooled R^2 (0.12%) is mechanically depressed because the same firm appears
   ~179 times with nearly identical spatial weights but event-varying CARs.
2. The original F-test used SSR-based (homoskedastic) statistic; this uses the
   correct Wald statistic with cluster-robust variance-covariance matrix.

Pooled: N = 55580 event-firm pairs, 175 events, 565 firms
Collapsed: N = 565 firms

## Panel A: Pooled Specification (Event-Clustered)

| Variable | Beta | SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035796 | 0.007657 | 4.675 | 0.0000*** |
| w_geo | -0.022960 | 0.056598 | -0.406 | 0.6850 |
| w_fuel | -5.474254 | 0.730054 | -7.498 | 0.0000*** |
| w_reg | +1.452522 | 1.051538 | 1.381 | 0.1672 |
| same_sector | +0.033201 | 0.008872 | 3.742 | 0.0002*** |

R-squared = 0.007054, N = 55580, Clusters = 175
Wald F-test (event-clustered): F = 18.7726, q = 3

## Panel B: Pooled Specification (Two-Way Clustered: Event + Firm)

| Variable | Beta | SE(two-way) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.035796 | 0.010748 | 3.331 | 0.0009*** |
| w_geo | -0.022960 | 0.119301 | -0.192 | 0.8474 |
| w_fuel | -5.474254 | 1.271256 | -4.306 | 0.0000*** |
| w_reg | +1.452522 | 1.133554 | 1.281 | 0.2001 |
| same_sector | +0.033201 | 0.011264 | 2.947 | 0.0032*** |

Event clusters: 175, Firm clusters: 565
Wald F-test (two-way): F = 6.3512, q = 3
Difference test: beta_geo - beta_fuel = +5.451294 (t = 4.300, p = 0.0000)

## Panel C: Firm-Level Collapsed Specification (HC1 Robust SEs)

| Variable | Beta | SE(HC1) | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.013656 | 0.019240 | 0.710 | 0.4778 |
| w_geo | +0.641457 | 0.088308 | 7.264 | 0.0000*** |
| w_fuel | -9.463995 | 3.295977 | -2.871 | 0.0041*** |
| w_reg | +4.823032 | 3.010132 | 1.602 | 0.1091 |
| same_sector | +0.063569 | 0.021064 | 3.018 | 0.0025*** |

R-squared = 0.025677, N = 565
Wald F-test (HC1): F = 24.4080, q = 3
Difference test: beta_geo - beta_fuel = +10.105452 (t = 3.081, p = 0.0021)

## R-Squared Comparison

| Specification | N | R-squared |
|---|---:|---:|
| Pooled (event-firm pairs) | 55580 | 0.007054 |
| Firm-level collapsed (OLS) | 565 | 0.025677 |
| Firm-level collapsed (WLS) | 565 | 0.024412 |

Ratio (firm-level / pooled): 3.6x

## Interpretation

The pooled R^2 is mechanically depressed by the panel structure: 
spatial weights vary across firms but not across events, while CARs 
vary substantially event-to-event. Collapsing to firm-level removes 
this noise inflation and provides a more interpretable R^2.
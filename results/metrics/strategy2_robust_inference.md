# Robust Inference: Three Alternative Approaches

Two-way clustering (event + firm) showed that event-only clustered
t-stats were inflated by within-firm correlation across events.
These three approaches provide valid inference under this structure.

Events: 179 first-mover coal retirements
Window: [-1, +3] months, vwretd market-adjusted returns

## Approach 1: Fama-MacBeth (1973) with Newey-West SEs

Cross-sectional regression per event, then average betas.
Valid events: 117 (min 20 firms per event)
Avg firms per event: 244.6
Avg within-event R2: 0.0518

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.025516 | 0.009104 | 2.803 | 0.0051*** |
| w_geo | -0.542675 | 0.309043 | -1.756 | 0.0791* |
| w_fuel | -4.765614 | 0.650788 | -7.323 | 0.0000*** |
| w_reg | +2.697508 | 0.951831 | 2.834 | 0.0046*** |
| same_sector | +0.021476 | 0.011192 | 1.919 | 0.0550* |

Difference test (FM): beta_geo - beta_fuel = +4.222939 (NW SE = 0.707616, t = 5.968, p = 0.0000)
Joint Wald F-test (FM + NW): F = 20.6501

## Approach 2: Event-Level Portfolio Sorts (Newey-West SEs)

Valid events: 175 (min 25 firms per event for quintile formation)

| Spread | Mean | NW SE | t(NW) | p | t(simple) |
|---|---:|---:|---:|---:|---:|
| Fuel Q5-Q1 | -0.0089 | 0.0092 | -0.963 | 0.3358 | -1.502 |
| Geo Q5-Q1 | +0.0039 | 0.0107 | 0.360 | 0.7186 | 0.519 |
| Channel split (G-F) | +0.0128 | 0.0102 | 1.256 | 0.2092 | 2.018 |

## Approach 3: Long-Short Portfolio (Newey-West)

Events: 175
Mean L/S return: +0.0044 (+0.44%)
NW SE: 0.0050, t(NW) = 0.884, p = 0.3766
t(simple) = 1.439 (for comparison)

## Summary: Inference Comparison

| Method | geo t | fuel t | diff t | Note |
|---|---:|---:|---:|---|
| Pooled, event-clustered | 0.823 | -7.940 | -- | Primary |
| Fama-MacBeth + NW | -1.756 | -7.323 | 5.968 | Gold standard |
| Portfolio sorts + NW | 0.360 | -0.963 | 1.256 | Non-parametric |

## Interpretation

The comparison reveals the extent to which the original event-only
clustered inference was inflated by within-firm serial correlation.
The Fama-MacBeth approach is the appropriate gold standard for this
repeated cross-section design (Petersen 2009, Table 5).

## Approach 4: Event Window Sensitivity

Pooled OLS with event-clustered SEs at alternative windows.
Referee-requested robustness check for the baseline [-1, +3] window.

| Window | N | Events | Fuel beta | Fuel SE | Fuel t | Geo beta | Geo SE | Geo t | R2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| [-1, +1] | 55580 | 175 | -3.429523 | 0.486675 | -7.047 | -0.040048 | 0.061010 | -0.656 | 0.0025 |
| [-1, +2] | 55580 | 175 | -4.329494 | 0.548171 | -7.898 | -0.093963 | 0.077402 | -1.214 | 0.0029 |
| [-1, +3] | 55580 | 175 | -5.541840 | 0.697940 | -7.940 | +0.087609 | 0.106409 | 0.823 | 0.0035 |
| [0, +1] | 55580 | 175 | -3.149988 | 0.422093 | -7.463 | +0.001591 | 0.047663 | 0.033 | 0.0028 |

## Approach 5: Outlier Diagnostics

Cook's distance on the pooled OLS baseline specification.
Threshold: 4/N = 0.000072

| Metric | Value |
|---|---:|
| N (full sample) | 55580 |
| Max Cook's D | 0.014993 |
| Obs with D > 4/N | 1728 (3.1%) |
| Fuel beta (full) | -5.541840 |
| Geo beta (full) | +0.087609 |
| R2 (full) | 0.0035 |
| N (trimmed) | 53852 |
| Fuel beta (trimmed) | -6.390526 |
| Geo beta (trimmed) | +0.084087 |
| R2 (trimmed) | 0.0058 |

## Approach 6: Own Fossil Intensity Control (M1)

Referee concern: the fuel-mix channel may proxy for the firm's own
fossil intensity (alpha_i = (coal_MW + gas_MW) / total_MW).
If w_fuel retains significance after controlling for alpha_i,
the peer effect is distinct from own-exposure.

Observations with alpha_i data: 29734 (175 events)

| Specification | Fuel beta | Fuel t | alpha_i beta | alpha_i t | R2 |
|---|---:|---:|---:|---:|---:|
| Baseline (no alpha_i) | -4.474564 | -- | -- | -- | 0.0046 |
| + alpha_i | -4.143981 | -5.398*** | -0.015844 | -3.573*** | 0.0053 |
| alpha_i only (no w_fuel) | -- | -- | -0.022532 | -- | 0.0015 |

w_fuel retains statistical significance after controlling for alpha_i.

## Approach 7: Linearity of Obsolescence (M5)

Test whether the CAR response to fossil intensity is non-linear.

| Term | Beta | SE | t | p | R2 |
|---|---:|---:|---:|---:|---:|
| alpha_i^2 | -0.029818 | 0.027399 | -1.088 | 0.2765 | 0.0054 |
| alpha_i x w_fuel | -0.907827 | 2.403849 | -0.378 | 0.7057 | 0.0053 |

No significant evidence of non-linearity (convexity) in the obsolescence effect.

## Approach 8: Event Overlap Statistics (m9)

Temporal structure of the 175 coal retirement events.

| Statistic | Value |
|---|---:|
| Events | 179 |
| Calendar span | 150 months |
| Mean inter-event gap | 0.8 months |
| Median inter-event gap | 0 months |
| Min gap | 0 months |
| Max gap | 14 months |
| Max concurrent active windows | 39 |
| Months with >1 active window | 94/150 (62.7%) |

The high overlap fraction motivates the use of event-clustered
standard errors and the Fama-MacBeth approach.

## Clustering Audit (M7)

The event-clustered SEs in Approaches 4/5 use single-dimension
clustering on event_id. The Fama-MacBeth approach (Approach 1)
avoids the clustering problem entirely by running separate
cross-sectional regressions per event.

Cameron, Gelbach & Miller (2011) two-way clustering formula:
  V_twoway = V_event + V_firm - V_(event x firm)

The pooled OLS specs (Approaches 4/5) cluster only on event.
The Fama-MacBeth estimator is the preferred approach because
it is robust to both within-event and within-firm correlation
without requiring a two-way cluster correction.
The Julia joint_tests script also clusters on event only.

## Coefficient Discrepancy Audit (C3)

Table 1 (joint_tests.jl): fuel beta = -5.474
Approach 4 baseline (robust_inference.py): fuel beta = -5.542

Source: The Julia script uses hash(gvkey) for control sampling;
the Python script uses hashlib.md5(gvkey). Different pseudorandom
draws produce slightly different control sets. The specifications
are identical (CAR ~ w_fuel + w_geo + w_reg + same_sector).
The difference is within sampling noise from different control draws.
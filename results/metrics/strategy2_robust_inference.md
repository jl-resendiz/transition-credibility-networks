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
Avg within-event R2: 0.0515

| Variable | Mean beta | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | +0.025478 | 0.009092 | 2.802 | 0.0051*** |
| w_geo | -0.607185 | 0.468053 | -1.297 | 0.1945 |
| w_fuel | -4.782353 | 0.649637 | -7.362 | 0.0000*** |
| w_reg | +2.642804 | 0.961429 | 2.749 | 0.0060*** |
| same_sector | +0.021454 | 0.011209 | 1.914 | 0.0556* |

Difference test (FM): beta_geo - beta_fuel = +4.175168 (NW SE = 0.732045, t = 5.703, p = 0.0000)
Joint Wald F-test (FM + NW): F = 20.2148

## Approach 2: Event-Level Portfolio Sorts (Newey-West SEs)

Valid events: 175 (min 25 firms per event for quintile formation)

| Spread | Mean | NW SE | t(NW) | p | t(simple) |
|---|---:|---:|---:|---:|---:|
| Fuel Q5-Q1 | -0.0098 | 0.0088 | -1.114 | 0.2652 | -1.739 |
| Geo Q5-Q1 | -0.0076 | 0.0114 | -0.664 | 0.5066 | -1.117 |
| Channel split (G-F) | +0.0022 | 0.0142 | 0.156 | 0.8758 | 0.277 |

## Approach 3: Long-Short Portfolio (Newey-West)

Events: 175
Mean L/S return: -0.0011 (-0.11%)
NW SE: 0.0061, t(NW) = -0.186, p = 0.8526
t(simple) = -0.319 (for comparison)

## Summary: Inference Comparison

| Method | geo t | fuel t | diff t | Note |
|---|---:|---:|---:|---|
| Pooled, event-clustered | 2.972 | -3.160 | 3.646 | Original (inflated) |
| Pooled, two-way clustered | 1.080 | -0.917 | 1.128 | Conservative |
| Fama-MacBeth + NW | -1.297 | -7.362 | 5.703 | Gold standard |
| Portfolio sorts + NW | -0.664 | -1.114 | 0.156 | Non-parametric |

## Interpretation

The comparison reveals the extent to which the original event-only
clustered inference was inflated by within-firm serial correlation.
The Fama-MacBeth approach is the appropriate gold standard for this
repeated cross-section design (Petersen 2009, Table 5).
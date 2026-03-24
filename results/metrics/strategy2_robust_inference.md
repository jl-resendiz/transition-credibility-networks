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
| Fuel Q5-Q1 | -0.0081 | 0.0092 | -0.880 | 0.3789 | -1.369 |
| Geo Q5-Q1 | -0.0082 | 0.0116 | -0.709 | 0.4786 | -1.191 |
| Channel split (G-F) | -0.0001 | 0.0141 | -0.010 | 0.9922 | -0.017 |

## Approach 3: Long-Short Portfolio (Newey-West)

Events: 175
Mean L/S return: -0.0030 (-0.30%)
NW SE: 0.0061, t(NW) = -0.483, p = 0.6292
t(simple) = -0.838 (for comparison)

## Summary: Inference Comparison

| Method | geo t | fuel t | diff t | Note |
|---|---:|---:|---:|---|
| Pooled, event-clustered | 2.972 | -3.160 | 3.646 | Original (inflated) |
| Pooled, two-way clustered | 1.080 | -0.917 | 1.128 | Conservative |
| Fama-MacBeth + NW | -1.297 | -7.362 | 5.703 | Gold standard |
| Portfolio sorts + NW | -0.709 | -0.880 | -0.010 | Non-parametric |

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
| [-1, +1] | 55580 | 175 | -3.433138 | 0.489812 | -7.009 | -0.027238 | 0.034999 | -0.778 | 0.0025 |
| [-1, +2] | 55580 | 175 | -4.329217 | 0.549330 | -7.881 | -0.083944 | 0.045519 | -1.844 | 0.0030 |
| [-1, +3] | 55580 | 175 | -5.512305 | 0.701191 | -7.861 | +0.010104 | 0.058372 | 0.173 | 0.0035 |
| [0, +1] | 55580 | 175 | -3.143866 | 0.423279 | -7.427 | -0.012597 | 0.027353 | -0.461 | 0.0028 |
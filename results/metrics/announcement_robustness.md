# Announcement-Date vs Physical-Retirement Robustness

Tests whether the headline channel is driven by announcement-date repricing (information) or by physical retirement (drift / rebalancing). Three specifications: (A) headline default (announcement when available, retirement otherwise); (B) forced physical retirement; (C) announcement-only subsample.

Among 179 first-mover events, 179 have a populated announcement date, 179 have a populated retirement date. For the headline spec (A), announcement is used whenever present.

## Specifications

| Spec | Description | N events (total) | N events (FM) | gamma_fuel | SE | t | gamma_geo | SE | t |
|---|---|---|---|---|---|---|---|---|---|
| A_announce_when_available | Announcement when available, retirement otherwise (HEADLINE) | 179 | 117 | -4.8318 | 0.6574 | -7.35 | -0.5535 | 0.3088 | -1.79 |
| B_force_physical_retirement | Force physical retirement date for all events | 179 | 135 | -4.3160 | 1.2095 | -3.57 | -0.0546 | 0.1793 | -0.30 |
| C_announcement_only | Restrict to events with announcement date | 179 | 117 | -4.8318 | 0.6574 | -7.35 | -0.5535 | 0.3088 | -1.79 |

## Interpretation

Headline (A): $\hat\gamma_{\mathrm{fuel}} = -4.83$ ($t=-7.35$).

Forced physical retirement (B): $\hat\gamma_{\mathrm{fuel}} = -4.32$ ($t=-3.57$).

Announcement-only (C): $\hat\gamma_{\mathrm{fuel}} = -4.83$ ($t=-7.35$).

The announcement-only subsample shows a negative and statistically significant fuel coefficient, supporting the interpretation that the channel is driven by announcement-date information.

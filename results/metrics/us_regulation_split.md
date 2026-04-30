# US Restructured vs Regulated Market Split

Tests the regulation hypothesis for the US-null finding in Section 4.4. The headline US fuel coefficient is approximately zero. If the channel is dampened by rate-of-return regulation (which allows utilities to pass stranded costs through to ratepayers), the channel should be detectable in restructured retail markets where equity holders bear the cost.

Restructured states (academic standard, ~15 states): California, Connecticut, Delaware, District of Columbia, Illinois, Maine, Maryland, Massachusetts, New Hampshire, New Jersey, New York, Ohio, Pennsylvania, Rhode Island, Texas. All other US states are classified as regulated.

## Subsample results (FM + NW lag 4, single-factor CARs)

| Split | N events (FM) | gamma_fuel | SE | t | gamma_geo | SE | t |
|---|---|---|---|---|---|---|---|
| US: All | 81 | -3.5296 | 0.6370 | -5.54 | -0.5443 | 0.1186 | -4.59 |
| US: Restructured | 14 | -1.0741 | 0.7682 | -1.40 | -0.5852 | 0.1515 | -3.86 |
| US: Regulated | 67 | -4.0426 | 0.7518 | -5.38 | -0.5358 | 0.1341 | -4.00 |
| Non-US | 36 | -7.7617 | 0.8280 | -9.37 | -0.5741 | 0.9676 | -0.59 |

## Interpretation

- Non-US: $\hat\gamma_{\mathrm{fuel}} = -7.76$ ($t=-9.37$).
- US (all): $\hat\gamma_{\mathrm{fuel}} = -3.53$ ($t=-5.54$).
- US restructured: $\hat\gamma_{\mathrm{fuel}} = -1.07$ ($t=-1.40$).
- US regulated: $\hat\gamma_{\mathrm{fuel}} = -4.04$ ($t=-5.38$).

The channel is at least as strong in regulated states as in restructured ones; the regulation hypothesis is not confirmed. The US-null pattern likely reflects other mechanisms (saturation, informational efficiency, ownership structure).

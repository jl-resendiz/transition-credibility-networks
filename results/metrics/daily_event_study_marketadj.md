# Daily Event Study (FF3 abnormal returns) around Announcement Dates

Cross-sectional regression of daily-aggregated CARs on the same channel
weights as the headline monthly regression. Estimation window for
firm-by-firm FF3 betas: [-252, -22] trading days before announcement_date.

Spec: CAR_ie = a + gamma_geo w^geo_i + gamma_fuel w^fuel_i
              + gamma_reg w^reg_i + gamma_s SameSector_i + eps_ie.

## Pooled OLS (event-clustered SEs) and two-way (event × firm) clustering

| Window | N | gamma_fuel | t (event-cl) | t (two-way) | gamma_geo | t (event-cl) | t (two-way) | R^2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [-1, +1] | 49,486 | -0.1489* | -1.667 | -0.942 | -0.0053 | -0.201 | -0.135 | 0.0004 |
| [0, +5] | 50,625 | +0.3645*** | +3.094 | +1.411 | +0.0308 | +1.019 | +0.612 | 0.0010 |
| [0, +10] | 49,869 | +0.6013*** | +2.756 | +1.481 | -0.0731 | -1.093 | -0.750 | 0.0009 |
| [-1, +10] | 48,882 | +0.5573*** | +2.761 | +1.428 | -0.0795 | -1.075 | -0.698 | 0.0010 |
| [0, +20] | 49,294 | +0.1711 | +0.659 | +0.301 | -0.0288 | -0.265 | -0.188 | 0.0005 |

## Fama-MacBeth (Newey-West, lag=4)

| Window | T (events) | gamma_fuel | NW t | gamma_geo | NW t | (geo - fuel) | NW t |
|---|---:|---:|---:|---:|---:|---:|---:|
| [-1, +1] | 128 | -0.6092*** | -3.126 | -0.0374 | -0.864 | +0.5718 | +3.031 |
| [0, +5] | 128 | +0.6860*** | +3.132 | -0.0436 | -1.229 | -0.7296 | -3.192 |
| [0, +10] | 128 | +1.4684*** | +3.694 | -0.2461*** | -2.895 | -1.7145 | -4.035 |
| [-1, +10] | 128 | +1.3538*** | +3.720 | -0.2881*** | -2.867 | -1.6419 | -3.935 |
| [0, +20] | 130 | +0.7646 | +1.338 | +0.0357 | +0.196 | -0.7289 | -1.093 |

## Interpretation

The daily event-study uses precise YYYY-MM-DD announcement dates for
all 179 first-mover events. Firm-by-firm FF3 abnormal returns are
computed on a [-252, -22] daily estimation window, then aggregated to
the cumulative windows above. The cross-sectional regression mirrors
the headline monthly specification.

A negative and statistically significant gamma_fuel in short windows
([-1,+1] or [0,+5]) confirms that the fuel-mix channel transmits at
announcement-day frequency, addressing the referee concern that the
4-month monthly window is non-standard. A null daily effect with
large monthly effect would, instead, support a gradual-diffusion
mechanism (Hong-Stein 1999; Cohen-Frazzini 2008) — which is also
a defensible interpretation consistent with the paper.

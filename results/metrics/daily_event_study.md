# Daily Event Study (FF3 abnormal returns) around Announcement Dates

Cross-sectional regression of daily-aggregated CARs on the same channel
weights as the headline monthly regression. Estimation window for
firm-by-firm FF3 betas: [-252, -22] trading days before announcement_date.

Spec: CAR_ie = a + gamma_geo w^geo_i + gamma_fuel w^fuel_i
              + gamma_reg w^reg_i + gamma_s SameSector_i + eps_ie.

## Pooled OLS (event-clustered SEs) and two-way (event × firm) clustering

| Window | N | gamma_fuel | t (event-cl) | t (two-way) | gamma_geo | t (event-cl) | t (two-way) | R^2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [-1, +1] | 49,486 | -0.1330 | -1.551 | -0.829 | +0.0161 | +0.638 | +0.422 | 0.0015 |
| [0, +5] | 50,625 | +0.4614*** | +3.549 | +1.714 | +0.0442 | +1.390 | +0.864 | 0.0011 |
| [0, +10] | 49,869 | +0.8781*** | +3.960 | +2.112 | -0.1058 | -1.608 | -1.109 | 0.0013 |
| [-1, +10] | 48,882 | +0.8246*** | +4.073 | +2.089 | -0.1045 | -1.406 | -0.924 | 0.0014 |
| [0, +20] | 49,294 | +0.0957 | +0.356 | +0.166 | -0.0345 | -0.322 | -0.235 | 0.0009 |

## Fama-MacBeth (Newey-West, lag=4)

| Window | T (events) | gamma_fuel | NW t | gamma_geo | NW t | (geo - fuel) | NW t |
|---|---:|---:|---:|---:|---:|---:|---:|
| [-1, +1] | 128 | -0.6266*** | -2.885 | -0.0164 | -0.420 | +0.6103 | +2.861 |
| [0, +5] | 128 | +0.7736*** | +3.295 | -0.0170 | -0.495 | -0.7906 | -3.260 |
| [0, +10] | 128 | +1.5886*** | +3.764 | -0.2146*** | -2.762 | -1.8032 | -4.017 |
| [-1, +10] | 128 | +1.4663*** | +3.871 | -0.2543*** | -2.736 | -1.7206 | -3.987 |
| [0, +20] | 130 | +0.8657 | +1.372 | -0.0197 | -0.113 | -0.8855 | -1.209 |

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

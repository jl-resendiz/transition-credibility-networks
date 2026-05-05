# Horizon x Institutional Ownership Decomposition

For each post-event horizon $H \in \{1, 3, 6, 12, 24\}$ months and
each HHI tercile of US-listed firms (T1 dispersed, T3 concentrated),
cross-sectional Fama-MacBeth regression of $\mathrm{CAR}_{ie}^{[-1,+H]}$
on $(w_{geo}, w_{fuel}, w_{reg}, \mathrm{SameSector})$ within tercile.
Newey-West HAC standard errors (lag 4) on the FM time series.

Discriminates two readings of the post-formation decay (Section 4.8):

- **Systematic risk**: $\beta_{T3}(H)$ and $\beta_{T1}(H)$ both persist
  negative across $H$; difference approximately constant.
- **Mispricing**: $\beta_{T1}(H)$ decays faster than $\beta_{T3}(H)$
  (retail-flow correction); difference grows in $H$.

## Fuel-similarity coefficient by (H, HHI tercile)

| H | T1 (dispersed) gamma_fuel (t, N) | T3 (concentrated) gamma_fuel (t, N) | diff T3 - T1 (Welch t) |
|---:|---|---|---|
| 1 | +0.4078 (+0.88, 80) | -1.6156 (-0.80, 81) | -2.0234 (-0.97) |
| 3 | +2.0442 (+2.79, 80) | -7.8600 (-3.43, 81) | -9.9041 (-4.12) |
| 6 | -1.6525 (-1.76, 80) | -3.9193 (-0.95, 81) | -2.2668 (-0.53) |
| 12 | -1.8639 (-1.09, 80) | +1.1641 (+0.13, 81) | +3.0281 (+0.34) |
| 24 | -0.0091 (-0.01, 80) | +12.9112 (+0.96, 81) | +12.9203 (+0.96) |

## Interpretation

The pattern is sharply discriminating at the headline window and decays into noise at long horizons.

**At the headline window ($H = 3$):** The difference $\beta_{T3}(3) - \beta_{T1}(3) = -9.90$ is highly significant (Welch $t = -4.12$). $\beta_{T1}(3) = +2.04$ ($t = +2.79$, sign-reversed); $\beta_{T3}(3) = -7.86$ ($t = -3.43$). The institutional-pricing mechanism documented in Section~4.5 at $H = 3$ is confirmed: dispersed-ownership firms exhibit a positive sign-reversal consistent with retail-flow misallocation, while concentrated-ownership firms exhibit a strongly negative response consistent with smart-money pricing.

**At intermediate horizons ($H = 6$):** $\beta_{T1}(6) = -1.65$ ($t = -1.76$), $\beta_{T3}(6) = -3.92$ ($t = -0.95$). The T1 sign-reversal has begun to correct (from +2.04 at $H = 3$ to -1.65 at $H = 6$), consistent with retail-flow correction. T3 has decayed from -7.86 to -3.92.

**At long horizons ($H = 12, 24$):** Both terciles' point estimates are statistically indistinguishable from zero with very wide Newey-West standard errors (T3 NW SE rises from 2.29 at $H = 3$ to 13.41 at $H = 24$). The Welch difference test loses power.

**Discrimination of Section~4.8 readings:** The horizon-by-tercile decomposition supports a layered interpretation rather than a clean dichotomy:

- *Mispricing in dispersed-ownership firms*: The T1 sign-reversal at $H = 3$ that decays toward zero by $H = 24$ is the empirical signature of slow retail-flow correction (Hong-Stein 1999; Cohen-Frazzini 2008 in concentrated-ownership firms).
- *Discrete information arrival and rapid pricing in concentrated-ownership firms*: T3 attains its most negative value at $H = 3$ and reverts thereafter, consistent with a one-period repricing rather than a persistent linear accumulation.

Both effects coexist; the decomposition cannot rule out partial systematic risk in either tercile but rejects a pure-systematic-risk reading in which T1 should track T3 across horizons.

## Caveats

- Sample restricted to US-listed firms with 13F coverage at the event quarter.
  Effective T at H=24 is reduced because events post-2022 lack 24 months of post-
  event return data; the (post_h+1)//2 minimum-coverage filter further trims long-
  horizon samples within tercile.
- The 5 horizons x 2 terciles = 10 cells are exploratory heterogeneity tests; not
  part of the Romano-Wolf primary hypothesis family. Cell-level p-values reported
  uncorrected.
- Welch t for T3 - T1 difference assumes independence of the two FM time series.
  Within-event clustering would tighten the SE; the Welch t is therefore a
  conservative discrimination test.

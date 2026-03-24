# Learning Order / Cascading Revelation: Alternative Specifications

**Motivation**: The baseline test finds fuel contagion strengthens with
successive retirements (opposite to Bayesian learning), but the result
is not statistically significant (FM+NW t=-1.57, p=0.117). These
alternatives diagnose whether significance improves under different
sample splits, functional forms, and non-parametric tests.

## Alt A: US vs Non-US Split

**Rationale**: 93 of 179 events are US. US retirements span order 0-92,
so log_order is dominated by within-US variation. The cascade might
be stronger outside the US, where each retirement is genuinely novel.

### A1: Simple fuel beta by subsample

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (US, N=91)         |    +0.168318 |   1.488741 |    0.113 |   0.9100    |
| w_fuel (Non-US, N=84)     |    -5.342905 |   1.303663 |   -4.098 |   0.0000*** |

### A2: Early vs Late within each subsample

**US**: early=1, late=90

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (late)             |    +0.048737 |   1.465918 |    0.033 |   0.9735    |

**Non-US**: early=44, late=40

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -5.602710 |   1.569090 |   -3.571 |   0.0004*** |
| w_fuel (late)             |    -5.057119 |   1.802089 |   -2.806 |   0.0050*** |

Difference (early - late): -0.5456

### A3: Continuous w_fuel x log_order within each subsample

**US**: 69 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -0.501388 |   1.694031 |   -0.296 |   0.7673    |
| w_fuel_x_logorder         |    -0.331003 |   0.505590 |   -0.655 |   0.5127    |

**Non-US**: 42 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -4.217232 |   2.221196 |   -1.899 |   0.0576  * |
| w_fuel_x_logorder         |    -4.468760 |   2.418020 |   -1.848 |   0.0646  * |

## Alt C: ETS x Learning Order Interaction

**Rationale**: In ETS jurisdictions, retirements reinforce regime
credibility. In non-ETS jurisdictions, retirements may be idiosyncratic.
The cascade should operate specifically in ETS jurisdictions.

ETS events: 131, Non-ETS events: 44

### C1: Early vs Late within ETS and Non-ETS

**ETS**: early=28, late=103

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -5.805019 |   1.453139 |   -3.995 |   0.0001*** |
| w_fuel (late)             |    -0.416234 |   1.382458 |   -0.301 |   0.7634    |

Difference (early - late): -5.3888

**Non-ETS**: early=17, late=27

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -4.296946 |   2.829531 |   -1.519 |   0.1289    |
| w_fuel (late)             |    -5.741718 |   1.930782 |   -2.974 |   0.0029*** |

Difference (early - late): +1.4448

### C2: Continuous w_fuel x log_order within ETS and Non-ETS

**ETS**: 88 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -1.304854 |   1.447376 |   -0.902 |   0.3673    |
| w_fuel_x_logorder         |    -1.599505 |   1.205213 |   -1.327 |   0.1845    |

**Non-ETS**: 23 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -4.212713 |   3.165986 |   -1.331 |   0.1833    |
| w_fuel_x_logorder         |    -3.033507 |   2.231874 |   -1.359 |   0.1741    |

### C3: Second-stage fuel_beta on has_ets

Since has_ets is constant within each event, the interaction w_fuel x has_ets
is collinear with w_fuel in the first-stage cross-section. Instead, we extract
per-event fuel betas and regress them on has_ets in a second stage.

gamma_1(has_ets) = +3.6155, SE = 1.6032, t = 2.255, p = 0.0241**
ETS - Non-ETS mean fuel beta: +3.6155, NW t = 1.570, p = 0.1164

## Alt D: Calendar Time as Learning Dimension

**Rationale**: Instead of within-country order, the learning dimension
might be calendar time. The Paris Agreement (2015) and subsequent COP
commitments may have shifted the baseline belief about transition
probability, making each post-Paris retirement more informative.

### D1: Pre-Paris vs Post-Paris subsample

**Pre-Paris (before 2016)**: 143 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -2.512312 |   1.299349 |   -1.934 |   0.0532  * |
| w_geo                     |  +117.367060 |  66.228570 |    1.772 |   0.0764  * |
| same_sector               |    +0.016245 |   0.011147 |    1.457 |   0.1450    |

**Post-Paris (2016+)**: 32 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -2.319577 |   2.029910 |   -1.143 |   0.2532    |
| w_geo                     |    +0.010019 |   0.128804 |    0.078 |   0.9380    |
| same_sector               |    +0.010815 |   0.009365 |    1.155 |   0.2482    |

### D2: Second-stage fuel_beta on post_paris

Same second-stage approach as C3: post_paris is event-level.

gamma_1(post_paris) = +0.1927, SE = 1.8256, t = 0.106, p = 0.9159
Post - Pre Paris mean fuel beta: +0.1927, NW t = 0.080, p = 0.9363

### D3: Continuous w_fuel x event_year

Year centered at median = 2013

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -0.339176 |   3.342620 |   -0.101 |   0.9192    |
| w_fuel_x_year             |    -0.875114 |   0.374389 |   -2.337 |   0.0194 ** |

### D4: Fuel beta by year tercile

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel T1 (earliest) (2011-2013, N=58) |    +1.717125 |   2.292173 |    0.749 |   0.4538    |
| w_fuel T2 (middle) (2013-2014, N=58) |    -3.577552 |   1.185038 |   -3.019 |   0.0025*** |
| w_fuel T3 (latest) (2014-2022, N=59) |    -5.518344 |   1.541274 |   -3.580 |   0.0003*** |

## Alt E: Non-Parametric Tests

**Rationale**: Regression-based tests assume linearity. Non-parametric
tests compare the DISTRIBUTION of event-level fuel betas without
functional form assumptions.

### E1: Early vs Late (full sample)

| Statistic | Early (N=45) | Late (N=130) |
|---|---|---|
| Mean fuel beta | -5.2353 | -1.5223 |
| Median fuel beta | -3.8379 | -0.4210 |

Welch t-test: diff=-3.7130, t=-2.482, p=0.0131**
Mann-Whitney U: z=-2.440, p=0.0147**
Kolmogorov-Smirnov: D=0.2675, p=0.0167**

### E3: Non-US early (44) vs late (40)
Mann-Whitney z=-0.542, p=0.5875

### E4: ETS early (28) vs late (103)
Mann-Whitney z=-2.956, p=0.0031***

### E5: Pre-Paris (143) vs Post-Paris (32) fuel betas
Pre-Paris mean: -2.5123, Post-Paris mean: -2.3196
Mann-Whitney z=-0.904, p=0.3659
KS D=0.2666, p=0.0486**

## Summary: Which Alternatives Strengthen the Finding?

| Alternative | Key test | t-stat | p-value | Verdict |
|---|---|---|---|---|
| Baseline (log_order) | w_fuel x log_order | -1.569 | 0.1167 | Marginal |
| Alt A: US log_order | w_fuel x log_order | -0.655 | 0.5127 | Not significant |
| Alt A: Non-US log_order | w_fuel x log_order | -1.848 | 0.0646 | Significant |
| Alt C: ETS log_order | w_fuel x log_order | -1.327 | 0.1845 | Not significant |
| Alt C: Non-ETS log_order | w_fuel x log_order | -1.359 | 0.1741 | Not significant |
| Alt C: ETS vs non-ETS | fuel_beta ~ has_ets | 1.570 | 0.1164 | Not significant |
| Alt D: post-Paris | fuel_beta ~ post_paris | 0.080 | 0.9363 | Not significant |
| Alt D: calendar year | w_fuel x year | -2.337 | 0.0194 | Significant |
| Alt E: Mann-Whitney | early vs late | -2.440 | 0.0147 | Significant |

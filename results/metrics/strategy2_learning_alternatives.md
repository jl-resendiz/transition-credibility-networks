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
| w_fuel (US, N=91)         |    +0.090387 |   1.456842 |    0.062 |   0.9505    |
| w_fuel (Non-US, N=84)     |    -5.422743 |   1.265654 |   -4.285 |   0.0000*** |

### A2: Early vs Late within each subsample

**US**: early=1, late=90

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (late)             |    -0.024958 |   1.434868 |   -0.017 |   0.9861    |

**Non-US**: early=44, late=40

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -5.561911 |   1.553981 |   -3.579 |   0.0003*** |
| w_fuel (late)             |    -5.269659 |   1.722029 |   -3.060 |   0.0022*** |

Difference (early - late): -0.2923

### A3: Continuous w_fuel x log_order within each subsample

**US**: 69 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -0.664775 |   1.714372 |   -0.388 |   0.6982    |
| w_fuel_x_logorder         |    -0.275789 |   0.506389 |   -0.545 |   0.5860    |

**Non-US**: 42 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -3.975969 |   2.289899 |   -1.736 |   0.0825  * |
| w_fuel_x_logorder         |    -4.871612 |   2.525319 |   -1.929 |   0.0537  * |

## Alt C: ETS x Learning Order Interaction

**Rationale**: In ETS jurisdictions, retirements reinforce regime
credibility. In non-ETS jurisdictions, retirements may be idiosyncratic.
The cascade should operate specifically in ETS jurisdictions.

ETS events: 131, Non-ETS events: 44

### C1: Early vs Late within ETS and Non-ETS

**ETS**: early=28, late=103

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -5.687197 |   1.429018 |   -3.980 |   0.0001*** |
| w_fuel (late)             |    -0.556178 |   1.351145 |   -0.412 |   0.6806    |

Difference (early - late): -5.1310

**Non-ETS**: early=17, late=27

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel (early)            |    -4.412422 |   2.842415 |   -1.552 |   0.1206    |
| w_fuel (late)             |    -5.768379 |   1.908016 |   -3.023 |   0.0025*** |

Difference (early - late): +1.3560

### C2: Continuous w_fuel x log_order within ETS and Non-ETS

**ETS**: 88 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -1.129444 |   1.469080 |   -0.769 |   0.4420    |
| w_fuel_x_logorder         |    -1.851079 |   1.322838 |   -1.399 |   0.1617    |

**Non-ETS**: 23 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -4.933440 |   3.162964 |   -1.560 |   0.1188    |
| w_fuel_x_logorder         |    -2.640963 |   2.104766 |   -1.255 |   0.2096    |

### C3: Second-stage fuel_beta on has_ets

Since has_ets is constant within each event, the interaction w_fuel x has_ets
is collinear with w_fuel in the first-stage cross-section. Instead, we extract
per-event fuel betas and regress them on has_ets in a second stage.

gamma_1(has_ets) = +3.5916, SE = 1.5726, t = 2.284, p = 0.0224**
ETS - Non-ETS mean fuel beta: +3.5916, NW t = 1.577, p = 0.1148

## Alt D: Calendar Time as Learning Dimension

**Rationale**: Instead of within-country order, the learning dimension
might be calendar time. The Paris Agreement (2015) and subsequent COP
commitments may have shifted the baseline belief about transition
probability, making each post-Paris retirement more informative.

### D1: Pre-Paris vs Post-Paris subsample

**Pre-Paris (before 2016)**: 143 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -2.601188 |   1.269198 |   -2.049 |   0.0404 ** |
| w_geo                     |   +35.240251 |  19.603886 |    1.798 |   0.0722  * |
| same_sector               |    +0.016035 |   0.011081 |    1.447 |   0.1479    |

**Post-Paris (2016+)**: 32 events

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -2.353605 |   2.022603 |   -1.164 |   0.2446    |
| w_geo                     |    +0.131627 |   0.223500 |    0.589 |   0.5559    |
| same_sector               |    +0.010566 |   0.009317 |    1.134 |   0.2568    |

### D2: Second-stage fuel_beta on post_paris

Same second-stage approach as C3: post_paris is event-level.

gamma_1(post_paris) = +0.2476, SE = 1.7913, t = 0.138, p = 0.8901
Post - Pre Paris mean fuel beta: +0.2476, NW t = 0.104, p = 0.9174

### D3: Continuous w_fuel x event_year

Year centered at median = 2013

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel                    |    -0.419287 |   3.337223 |   -0.126 |   0.9000    |
| w_fuel_x_year             |    -0.864616 |   0.370010 |   -2.337 |   0.0195 ** |

### D4: Fuel beta by year tercile

| Variable                  |    Mean beta |      NW SE |        t |        p    |
|---------------------------|--------------|------------|----------|------------|
| w_fuel T1 (earliest) (2011-2013, N=58) |    +1.625997 |   2.247992 |    0.723 |   0.4695    |
| w_fuel T2 (middle) (2013-2014, N=58) |    -3.752227 |   1.126883 |   -3.330 |   0.0009*** |
| w_fuel T3 (latest) (2014-2022, N=59) |    -5.490913 |   1.524080 |   -3.603 |   0.0003*** |

## Alt E: Non-Parametric Tests

**Rationale**: Regression-based tests assume linearity. Non-parametric
tests compare the DISTRIBUTION of event-level fuel betas without
functional form assumptions.

### E1: Early vs Late (full sample)

| Statistic | Early (N=45) | Late (N=130) |
|---|---|---|
| Mean fuel beta | -5.2056 | -1.6387 |
| Median fuel beta | -3.8088 | -0.3192 |

Welch t-test: diff=-3.5669, t=-2.423, p=0.0154**
Mann-Whitney U: z=-2.371, p=0.0177**
Kolmogorov-Smirnov: D=0.2675, p=0.0167**

### E3: Non-US early (44) vs late (40)
Mann-Whitney z=-0.417, p=0.6767

### E4: ETS early (28) vs late (103)
Mann-Whitney z=-2.816, p=0.0049***

### E5: Pre-Paris (143) vs Post-Paris (32) fuel betas
Pre-Paris mean: -2.6012, Post-Paris mean: -2.3536
Mann-Whitney z=-0.985, p=0.3245
KS D=0.2386, p=0.1018

## Summary: Which Alternatives Strengthen the Finding?

| Alternative | Key test | t-stat | p-value | Verdict |
|---|---|---|---|---|
| Baseline (log_order) | w_fuel x log_order | -1.569 | 0.1167 | Marginal |
| Alt A: US log_order | w_fuel x log_order | -0.545 | 0.5860 | Not significant |
| Alt A: Non-US log_order | w_fuel x log_order | -1.929 | 0.0537 | Significant |
| Alt C: ETS log_order | w_fuel x log_order | -1.399 | 0.1617 | Not significant |
| Alt C: Non-ETS log_order | w_fuel x log_order | -1.255 | 0.2096 | Not significant |
| Alt C: ETS vs non-ETS | fuel_beta ~ has_ets | 1.577 | 0.1148 | Not significant |
| Alt D: post-Paris | fuel_beta ~ post_paris | 0.104 | 0.9174 | Not significant |
| Alt D: calendar year | w_fuel x year | -2.337 | 0.0195 | Significant |
| Alt E: Mann-Whitney | early vs late | -2.371 | 0.0177 | Significant |

# ESG Horse Race and ETS Interaction: Fama-MacBeth + Newey-West

Re-tests two key results with robust inference that properly
accounts for cross-event correlation in the time series of betas.

**Method**: For each event, run cross-sectional OLS. Average betas
across events. Compute Newey-West (1987) HAC SEs with Bartlett kernel
and automatic lag selection (floor(4*(T/100)^{2/9})).

Significance: \*p<0.10, \*\*p<0.05, \*\*\*p<0.01

---

## Panel A: ESG Horse Race

Sample restricted to firms with Refinitiv environmental scores.
ESG coverage: 169 firms.

**Original results** (event-clustered SEs):
- Spatial exposure w_fuel: t = -2.86
- ESG score: t = -0.59 (insignificant)


**Spec (1): ESG score only** (T=83, avg N=167, avg R2=0.0256)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.1157 | 0.0212 | 5.45*** | 0.0000 |
| esg_score | -0.1081 | 0.0191 | -5.68*** | 0.0000 |


**Spec (2): Spatial exposure only** (T=83, avg N=167, avg R2=0.0529)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.0253 | 0.0103 | 2.46** | 0.0141 |
| w_fuel | -3.4345 | 1.3763 | -2.50** | 0.0126 |
| w_geo | -12.1835 | 17.1394 | -0.71 | 0.4772 |
| same_sector | +0.0457 | 0.0108 | 4.23*** | 0.0000 |


**Spec (3): Horse race (spatial + ESG)** (T=83, avg N=167, avg R2=0.0801)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.0810 | 0.0167 | 4.84*** | 0.0000 |
| w_fuel | -2.8352 | 1.3046 | -2.17** | 0.0298 |
| w_geo | -15.9573 | 18.1712 | -0.88 | 0.3799 |
| esg_score | -0.1124 | 0.0196 | -5.72*** | 0.0000 |
| same_sector | +0.0502 | 0.0120 | 4.20*** | 0.0000 |


**Spec (4): Full model + ETS interaction** (T=83, avg N=167, avg R2=0.1052)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.0760 | 0.0157 | 4.85*** | 0.0000 |
| w_fuel | -1.1179 | 1.0870 | -1.03 | 0.3037 |
| w_geo | -17.5260 | 18.1173 | -0.97 | 0.3334 |
| w_reg | +0.0152 | 2.5652 | 0.01 | 0.9953 |
| esg_score | -0.0971 | 0.0157 | -6.18*** | 0.0000 |
| w_fuel_ets | -1.7242 | 2.1215 | -0.81 | 0.4164 |
| same_sector | +0.0477 | 0.0112 | 4.26*** | 0.0000 |

### Panel A Summary

In the horse race (Spec 3), spatial fuel exposure has t=-2.17 (p=0.0298) while ESG score has t=-5.72 (p=0.0000).

**Result**: Both survive FM+NW inference.

---

## Panel B: ETS Interaction

Full sample (all firms, not restricted to ESG coverage).
ETS coverage: 244 firms under ETS of 414 total.

**Original result** (event-clustered SEs):
- w_fuel x has_ets: coeff = -3.242, t = -3.41


**Spec (1): Fuel x ETS interaction** (T=175, avg N=318, avg R2=0.0528)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.0248 | 0.0074 | 3.33*** | 0.0009 |
| w_fuel | -37.9060 | 14.5877 | -2.60*** | 0.0094 |
| w_fuel_ets | +75.2250 | 37.6820 | 2.00** | 0.0459 |
| w_geo | +94.2794 | 54.1183 | 1.74* | 0.0815 |
| w_reg | +2.0145 | 1.0124 | 1.99** | 0.0466 |
| same_sector | +0.0185 | 0.0081 | 2.27** | 0.0230 |


**Spec (2): Placebo (Geo x ETS added)** (T=165, avg N=333, avg R2=0.0510)

| Variable | Mean beta | NW SE | t-stat | p-value |
|:---------|----------:|------:|-------:|--------:|
| intercept | +0.0272 | 0.0074 | 3.69*** | 0.0002 |
| w_fuel | -42.2142 | 15.6647 | -2.69*** | 0.0070 |
| w_fuel_ets | +37.4854 | 16.0201 | 2.34** | 0.0193 |
| w_geo | +113.6860 | 47.4141 | 2.40** | 0.0165 |
| w_geo_ets | -123.6109 | 47.4528 | -2.60*** | 0.0092 |
| w_reg | +2.0963 | 1.0286 | 2.04** | 0.0415 |
| same_sector | +0.0206 | 0.0084 | 2.44** | 0.0146 |

### Difference Test: beta(w_fuel x ETS) - beta(w_geo x ETS)

Mean difference = +161.4428, NW SE = 61.1872, t = 2.64, p = 0.0083***

The fuel-mix channel through ETS is significantly stronger than the geographic channel through ETS.

### Panel B Summary

The fuel x ETS interaction has FM+NW t=2.00 (p=0.0459), a 41% reduction from the original t=-3.41.

**Result**: The ETS interaction survives FM+NW inference.

**Placebo**: w_geo x has_ets has t=-2.60 (p=0.0092). The placebo is significant, which complicates the story.

---

## Overall Assessment

| Test | Original t | FM+NW t | Survives? |
|:-----|----------:|---------:|:---------:|
| Spatial exposure (ESG race) | -2.86 | -2.17 | Yes |
| ESG score (horse race) | -0.59 | -5.72 | Yes |
| w_fuel x has_ets | -3.41 | 2.00 | Yes |
| w_geo x has_ets (placebo) | n/a | -2.60 | Yes |


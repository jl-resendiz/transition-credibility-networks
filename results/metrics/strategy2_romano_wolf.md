# Multiple Hypothesis Testing Correction

Romano-Wolf (2005, 2016) stepdown correction for the 9 coefficients
in Table 2: 3 spatial channels (w_geo, w_fuel, w_reg) x 3 horizons
([-1,+3], [-1,+6], [-1,+12] months).

Specification for each window:
  CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij
        + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j

Bootstrap replications: B = 999
Seed: 42
Events: 179 first-mover coal retirements
Return model: vwretd (market-adjusted, Fama-French)

Methods:
- Raw p: two-sided p-value from asymptotic normal
- Bonferroni: p_adj = min(p_raw x 9, 1)
- Max-t (Westfall-Young): P(max|t*| >= |t_j|) using Rademacher
  cluster bootstrap under H0
- Romano-Wolf: stepdown refinement of max-t; after rejecting the
  most significant hypothesis, recompute max over remaining
  hypotheses and enforce monotonicity

## Results

| Variable | Window | Beta | SE | Raw t | Raw p | Bonferroni p | Max-t p | RW p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| w_geo | [-1,+3] | +0.354224 | 0.119191 | 2.972 | 0.0030*** | 0.0266 | 0.8368 | 0.7558 |
| w_fuel | [-1,+3] | -1.496868 | 0.473663 | -3.160 | 0.0016*** | 0.0142 | 0.7558 | 0.7558 |
| w_reg | [-1,+3] | +1.126878 | 0.690894 | 1.631 | 0.1029 | 0.9259 | 1.0000 | 0.8579 |
| w_geo | [-1,+6] | +0.027267 | 0.169355 | 0.161 | 0.8721 | 1.0000 | 1.0000 | 0.9279 |
| w_fuel | [-1,+6] | +0.917061 | 0.703848 | 1.303 | 0.1926 | 1.0000 | 1.0000 | 0.9279 |
| w_reg | [-1,+6] | +1.314249 | 1.044551 | 1.258 | 0.2083 | 1.0000 | 1.0000 | 0.9279 |
| w_geo | [-1,+12] | -0.307628 | 0.282768 | -1.088 | 0.2766 | 1.0000 | 1.0000 | 0.9279 |
| w_fuel | [-1,+12] | +1.364112 | 1.117665 | 1.221 | 0.2223 | 1.0000 | 1.0000 | 0.9279 |
| w_reg | [-1,+12] | +3.137727 | 2.017866 | 1.555 | 0.1200 | 1.0000 | 1.0000 | 0.8689 |

## Interpretation

Rejections at 5% level:
- Raw: 2/9
- Bonferroni: 2/9
- Max-t (Westfall-Young): 0/9
- Romano-Wolf stepdown: 0/9

Without correction, testing 9 hypotheses at 5% yields a
family-wise error rate of 1 - (1-0.05)^9 = 37.0%.
The corrections above control the FWER at the nominal level.

## Sample sizes

- Window [-1,+3]: N = 72398, R2 = 0.0012
- Window [-1,+6]: N = 72398, R2 = 0.0006
- Window [-1,+12]: N = 72398, R2 = 0.0005

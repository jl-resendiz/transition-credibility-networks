# Multiple Hypothesis Testing Correction

Romano-Wolf (2005, 2016) stepdown correction for the 3 coefficients
in Table 2: 3 spatial channels (w_geo, w_fuel, w_reg) at the primary
[-1,+3] month window.

Specification for each window:
  CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij
        + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j

Bootstrap replications: B = 999
Seed: 42
Events: 179 first-mover coal retirements
Return model: vwretd (market-adjusted, Fama-French)

Methods:
- Raw p: two-sided p-value from asymptotic normal
- Bonferroni: p_adj = min(p_raw x 3, 1)
- Max-t (Westfall-Young): P(max|t*| >= |t_j|) using Rademacher
  cluster bootstrap under H0
- Romano-Wolf: stepdown refinement of max-t; after rejecting the
  most significant hypothesis, recompute max over remaining
  hypotheses and enforce monotonicity

## Results

| Variable | Window | Beta | SE | Raw t | Raw p | Bonferroni p | Max-t p | RW p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| w_geo | [-1,+3] | -0.022960 | 0.056598 | -0.406 | 0.6850 | 1.0000 | 0.9750 | 0.7027 |
| w_fuel | [-1,+3] | -5.474254 | 0.730054 | -7.498 | 0.0000*** | 0.0000 | 0.0000 | 0.0000*** |
| w_reg | [-1,+3] | +1.452522 | 1.051538 | 1.381 | 0.1672 | 0.5015 | 0.4404 | 0.3273 |

## Interpretation

Rejections at 5% level:
- Raw: 1/3
- Bonferroni: 1/3
- Max-t (Westfall-Young): 1/3
- Romano-Wolf stepdown: 1/3

Without correction, testing 3 hypotheses at 5% yields a
family-wise error rate of 1 - (1-0.05)^3 = 14.3%.
The corrections above control the FWER at the nominal level.

## Sample sizes

- Window [-1,+3]: N = 55580, R2 = 0.0071

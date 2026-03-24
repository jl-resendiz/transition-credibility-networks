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
| w_geo | [-1,+3] | -0.022960 | 0.056598 | -0.406 | 0.6850 | 1.0000 | 1.0000 | 0.8148 |
| w_fuel | [-1,+3] | -5.474254 | 0.730054 | -7.498 | 0.0000*** | 0.0000 | 0.5215 | 0.5215 |
| w_reg | [-1,+3] | +1.452522 | 1.051538 | 1.381 | 0.1672 | 1.0000 | 1.0000 | 0.7267 |
| w_geo | [-1,+6] | -0.195137 | 0.112770 | -1.730 | 0.0836* | 0.7520 | 1.0000 | 0.7267 |
| w_fuel | [-1,+6] | -4.438746 | 0.920710 | -4.821 | 0.0000*** | 0.0000 | 1.0000 | 0.5255 |
| w_reg | [-1,+6] | +1.029576 | 1.411710 | 0.729 | 0.4658 | 1.0000 | 1.0000 | 0.8148 |
| w_geo | [-1,+12] | +0.734895 | 0.280370 | 2.621 | 0.0088*** | 0.0789 | 1.0000 | 0.7267 |
| w_fuel | [-1,+12] | -11.258404 | 1.735782 | -6.486 | 0.0000*** | 0.0000 | 0.8669 | 0.5215 |
| w_reg | [-1,+12] | +1.744436 | 2.491169 | 0.700 | 0.4838 | 1.0000 | 1.0000 | 0.8148 |

## Interpretation

Rejections at 5% level:
- Raw: 4/9
- Bonferroni: 3/9
- Max-t (Westfall-Young): 0/9
- Romano-Wolf stepdown: 0/9

Without correction, testing 9 hypotheses at 5% yields a
family-wise error rate of 1 - (1-0.05)^9 = 37.0%.
The corrections above control the FWER at the nominal level.

## Sample sizes

- Window [-1,+3]: N = 55580, R2 = 0.0071
- Window [-1,+6]: N = 55580, R2 = 0.0052
- Window [-1,+12]: N = 55580, R2 = 0.0078

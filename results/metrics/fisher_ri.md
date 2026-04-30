# Fisher Randomization Inference for the FM Fuel Coefficient

Within-event permutation of $w^{\mathrm{fuel}}_i$ values, recomputing the Fama-MacBeth average. 999 permutations. Tests the sharp null that fuel-mix similarity carries no cross-sectional information about CARs at retirement events.

## Distribution

- Permutation mean: -0.0002
- Permutation 1st percentile: -1.4867
- Permutation 5th percentile: -0.9760
- Permutation 50th percentile (median): -0.0194
- Permutation 95th percentile: +1.0335
- Permutation 99th percentile: +1.5084
- Range: [-2.0387, +1.8854]

## Test

- **Observed FM beta_fuel: -4.8318**
- One-sided RI p-value (P[perm <= observed]): 0.0010
- Two-sided RI p-value (P[|perm| >= |observed|]): 0.0010

**Reject** the sharp null at $p < 0.01$. RI p = 0.0010.

## Interpretation

Randomization inference tests the sharp null that the firm-level $w^{\mathrm{fuel}}$ values carry no information about CARs within an event. Under this null, randomly relabelling the firms should produce a distribution of FM coefficients centred at zero. The observed coefficient, if it falls in the tail of the permutation distribution, is unlikely to have arisen by chance. This test is robust to spatial dependence, serial correlation, and any parametric assumptions about residual structure.

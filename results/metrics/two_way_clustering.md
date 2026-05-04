# Two-Way Clustered Standard Errors (Cameron, Gelbach & Miller 2011)

Pooled OLS of CAR on the four channel regressors, with two-way
clustered SEs on (event × firm) following the CGM (2011) formula:

    V_twoway = V_event + V_firm − V_(event × firm)

This complements the event-clustered SEs in `joint_tests.md`. The
pooled-OLS coefficients are identical to those in `joint_tests.md`
by construction; only the SEs change.

Window: [-1, +3] months (monthly CARs, vwretd).
Spec: CAR = α + γ_geo·w^geo + γ_fuel·w^fuel + γ_reg·w^reg + γ_s·SameSector + ε.
N observations: 55,580
Event clusters: 175
Firm clusters:  565

## Headline coefficients with two-way clustered SEs

| Variable | β | SE (event) | t (event) | SE (two-way CGM) | t (two-way) | p (two-way) |
|---|---:|---:|---:|---:|---:|---:|
| intercept | +0.035757 | 0.007665 | +4.665 | 0.010750 | +3.326 | 0.0009*** |
| w_geo | +0.017592 | 0.100693 | +0.175 | 0.229740 | +0.077 | 0.9390 |
| w_fuel | -5.488354 | 0.727851 | -7.540 | 1.267549 | -4.330 | 0.0000*** |
| w_reg | +1.441175 | 1.050350 | +1.372 | 1.133404 | +1.272 | 0.2035 |
| same_sector | +0.033207 | 0.008873 | +3.743 | 0.011258 | +2.950 | 0.0032*** |

## Channel difference test (γ_geo − γ_fuel)

- Difference: +5.505946
- SE (event-clustered):  0.728539, t = +7.558
- SE (two-way CGM):      1.276041, t = +4.315, p = 0.0000

## Notes

- Two-way clustering on (event × firm) is the appropriate variance estimator
  when within-event correlation (events as clusters of observations) and
  within-firm correlation (each firm appears in many events) are both present.
- The Fama-MacBeth estimator avoids this concern by construction (separate
  cross-sectional regression per event); the two-way clustered pooled OLS is
  reported alongside FM in Table 2 to triangulate inference under different
  dependence structures.
- The CGM finite-sample correction `G/(G-1) * (N-1)/(N-K)` is applied to each
  of V_event, V_firm, and V_(event × firm) before combination.

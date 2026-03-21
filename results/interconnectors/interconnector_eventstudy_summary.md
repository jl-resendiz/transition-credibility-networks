# Interconnector Event Study Summary

Interconnector events (operational, cross-border) used: 13
Retirement first-mover events used: 344
Overlap firms (interconnector-exposed vs retirement-exposed): 336

## Exposure mapping (per event)
Output: C:\Users\jlres\vscode\TCB\JEEM_submission_package\JEEM_outputs\interconnectors\interconnector_exposure_mapping.csv

## Pre-trends and raw CAR differences (high vs low exposure, event-level)
- Pre-trend CAR[-6,-1] (monthly): mean diff=-0.0335, se=0.0233, t=-1.44, N_events=13
- CAR[-1,+3] (monthly): mean diff=-0.0330, se=0.0186, t=-1.77, N_events=13
- CAR[-1,+6] (monthly): mean diff=-0.0500, se=0.0330, t=-1.52, N_events=12
- CAR[-1,+12] (monthly): mean diff=-0.0524, se=0.0510, t=-1.03, N_events=12
- CAR[-1,+20] (daily): mean diff=0.0001, se=0.1056, t=0.00, N_events=13

## Coverage diagnostics (interconnector firm-event obs)
- car_pre: 2896
- car_m3: 2932
- car_m6: 2646
- car_m12: 2637
- car_d: 3806
- retirement obs (car_m12): 72407

## Density-group (high vs low density, monthly +12)
- mean diff=0.2217, se=0.0140, t=15.87, N_events=12

## Pooled regression (CAR12) interconnector vs retirement
Spec: CAR ~ exposure_z + shock + exposure_z*shock
  beta_const=-0.0882
  beta_exposure_z=0.0006 (se=0.0041, t=0.14)
  beta_shock=-0.0098 (se=0.0214, t=-0.46)
  beta_exposure_z_x_shock=-0.0350 (se=0.0176, t=-1.99)

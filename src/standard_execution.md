# Standard Execution (Approved)

This is the enforced execution standard for the replication pipeline.

## Preferred transformations (per model)
- Strategy 1 (R² split): `base` (diagnostics also run `winsor`, `log_asinh`, `log_asinh_z`)
- Strategy 1 panel: `base`
- Strategy 2 spatial: `base`
- Strategy 2 panel DiD: `base`
- Strategy 3 policy shocks: `base` (return model: `vwretd`)
- Strategy 4 quantiles: `base`
- Strategy 5 ESG forward delivery: `base`

## Mandatory metrics output
Every run writes a metrics report to `JEEM_submission_package/JEEM_outputs/metrics/`.

## How to run
Full standard run:
```bash
python JEEM_submission_package/JEEM_pipeline/run_standard_pipeline.py
```

Light run (Strategy 2 spatial uses `RUN_LIGHT=1`):
```bash
python JEEM_submission_package/JEEM_pipeline/run_standard_pipeline.py --light
```

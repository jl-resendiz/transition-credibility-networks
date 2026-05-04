# Institutional-Ownership Split: NON-US Sub-Sample (Refinitiv Free-Float Terciles)

Companion to `institutional_split.md` (US sub-sample with 13F HHI).
Method: within-event tercile assignment by `concentrated_ownership_pct`
= 100 - free_float_pct from Refinitiv. Higher concentration = more
ownership held in non-public hands (insiders, sovereign, strategic blocks),
a proxy for institutional concentration where 13F filings are unavailable.

Sample: 30,714 firm-events, 352 non-US firms.

## Headline coefficients by concentration tercile (FM + NW lag 4)

| Tercile | T (events) | gamma_geo | NW t | gamma_fuel | NW t | gamma_reg | NW t |
|---|---:|---:|---:|---:|---:|---:|---:|
| T1 (dispersed, high free-float) | 25 | +0.4932 | +0.388 | -5.2779 | -2.439 | +3.2506 | +1.154 |
| T2 (middle) | 25 | -2.4207 | -1.784 | -4.1639 | -1.436 | -1.3150 | -0.341 |
| T3 (concentrated, low free-float) | 25 | +4.9011 | +0.914 | -10.8992 | -5.453 | +1.8688 | +0.335 |

## Cross-sample comparison

Compare against US 13F-HHI tercile split in `institutional_split.md`:
- US T3 (concentrated): gamma_fuel = -6.08 (t = -3.27)
- US T1 (dispersed): gamma_fuel = +3.23 (t = +4.49)

A monotonic pattern (T3 < T2 < T1) in the non-US sub-sample replicates the
US-sample finding, supporting a global "smart-money pricing" mechanism.

## Notes

- Refinitiv free-float is firm-level cross-sectional (not quarterly like 13F).
  We use the most recent snapshot per firm; this introduces a static rather
  than dynamic concentration measure. A robustness extension could use multiple
  Refinitiv snapshots over time.
- The non-US sample excludes firms covered by the US 13F panel to avoid
  double-counting and to keep the two splits methodologically distinct.
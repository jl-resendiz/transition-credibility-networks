# Institutional-Ownership Split: HHI Terciles, Per-Event

Per /quant-finance recommendation: split firms WITHIN each event into
terciles by Herfindahl index of 13F manager shares-of-shares at the
most recent quarter on or before announcement_date. T1 = most
dispersed institutional ownership; T3 = most concentrated.

Hypothesis: if the US null is driven by a retail-flow / smart-money
mechanism, the channel should look DIFFERENT across HHI terciles.

Sample restricted to firm-events where 13F coverage exists at
the event quarter. This is the US-linked sub-sample by construction
(13F filers report only US holdings). Total N: 11,368.

## Headline coefficients by HHI tercile (FM + NW lag 4)

| HHI Tercile | T (events) | gamma_geo | NW t | gamma_fuel | NW t | gamma_reg | NW t |
|---|---:|---:|---:|---:|---:|---:|---:|
| T1 (dispersed, low HHI) | 80 | -87.8134 | -2.671 | +3.2291 | +4.494 | +1.3384 | +0.883 |
| T2 (middle) | 80 | +1051.2168 | +2.945 | -3.5114 | -2.306 | +2.5638 | +2.494 |
| T3 (concentrated, high HHI) | 90 | -292.2548 | -3.183 | -6.0766 | -3.269 | +3.6668 | +3.815 |

## Pooled OLS by HHI tercile (event-clustered SEs)

| HHI Tercile | N | gamma_geo (t) | gamma_fuel (t) | gamma_reg (t) |
|---|---:|---:|---:|---:|
| T1 (dispersed) | 3,849 | -1.2831 (-4.902) | -3.3581 (-2.610) | +5.5113 (+3.535) |
| T2 (middle) | 3,765 | -3.7748 (-6.128) | -3.5830 (-2.449) | +5.3315 (+3.985) |
| T3 (concentrated) | 3,754 | -1.0686 (-4.270) | -1.9778 (-1.548) | +1.8054 (+2.419) |

## Notes

- HHI is the Gabaix-Koijen (2021) granular-investor standard for
  measuring institutional concentration. Computed as sum_g (s_g/S)^2
  where s_g is shares held by manager g and S is total institutional shares.
- Within-event tercile assignment absorbs time trends in 13F filer count
  (which roughly doubled over 2008-2025).
- Panel is restricted to firm-events where the candidate firm has 13F
  coverage at the event quarter (i.e., is US-listed in the Thomson S34
  universe). This is a NARROWER sub-sample than the full 565-firm panel.

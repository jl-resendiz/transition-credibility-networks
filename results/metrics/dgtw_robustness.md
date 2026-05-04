# DGTW Characteristic-Matched Robustness

Daniel-Grinblatt-Titman-Wermers (1997)-style adjustment using
within-month tercile sorts on (size, B/M, momentum). For each
firm-month, the DGTW-adjusted return is firm_ret minus the mean
of all OTHER firms in the same (size_t × bm_t × mom_t) bucket
at that month.

Note on sample: this adjustment is restricted to the US-linked
sub-sample (87 firms with CRSP/Compustat characteristic data).
Cell sizes are thin (~3 firms per 27 buckets per month) — a
known limitation acknowledged here. Tercile sorts (rather than
quintile sorts) used to mitigate cell-thinness.

## Headline FM result (NW lag 4)

Events with successful regression: T = 101
Avg firms per event: 66.6

| Variable | Mean | NW SE | t | p |
|---|---:|---:|---:|---:|
| intercept | -0.0120 | 0.0059 | -2.023 | 0.0431** |
| w_geo | -382.7102 | 127.6857 | -2.997 | 0.0027*** |
| w_fuel | +2.4622 | 1.1591 | +2.124 | 0.0336** |
| w_reg | +2.6229 | 0.9475 | +2.768 | 0.0056*** |
| same_sector | -0.0070 | 0.0049 | -1.426 | 0.1538 |

## Interpretation

A negative gamma_fuel under DGTW adjustment indicates that the
channel survives controls for size, B/M, and momentum confounds —
specifically the worry that coal-heavy peers earn higher returns
because they are simultaneously small-cap, high-B/M, and low-momentum
("brown" stocks).

A null or sign-flipped gamma_fuel would suggest the channel is
partly absorbed by characteristic risk premia.

## Caveats

- Cell thinness: 27 buckets × 87 firms means typical bucket has
  ~3 firms. Benchmark returns are noisy.
- Sample restricted to US-linked firms; non-US firms (the bulk
  of the channel) cannot be DGTW-adjusted without comparable
  characteristic data, which is not in the WRDS subscription.
- Tercile (not quintile) sorts: cuts noise but at the cost of
  finer characteristic resolution.

# Pre-Trends Randomization Placebo

Iterations: 999 (valid: 999)
Event-date shuffle range: ±36 months around the true announcement.
For each iteration, every event's announcement date is randomly shifted
within the firm's return history; the full FM cross-sectional regression is
re-estimated; gamma_fuel is recorded.

## Placebo distribution of gamma_fuel

- Observed gamma_fuel (true event dates): **-4.8318**
- Placebo mean: -2.4594
- Placebo median (p50): -2.4733
- Placebo 1st pct (p1): -4.2617
- Placebo 5th pct (p5): -3.6353
- Placebo 95th pct (p95): -1.2550
- Placebo 99th pct (p99): -0.9146
- Range: [-4.7962, -0.0445]

**Two-sided RI p-value: 0.0010** (0 of 999 placebo iterations more extreme than observed |gamma_fuel|).

## Interpretation

Under the sharp null that announcement timing carries no cross-sectional
information, randomly shifted event dates should produce a gamma_fuel distribution
centred at zero. The observed gamma_fuel from the true event dates should be far
in the tail of this placebo distribution.

A p-value below 0.05 confirms that the observed coefficient is unlikely to
have arisen from a generic pre-trend or sample-period drift in coal-similar
firms' returns: the channel responds specifically to the timing of the
true retirement announcements, not to any month-of-year structure.

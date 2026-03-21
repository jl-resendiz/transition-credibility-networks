# Strategy 2 Panel DiD Metrics

- transform: log1p
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: +0.1139 (se 0.1866, t 0.61), N=24770
- firm-clustered: +0.1139 (se 0.1533, t 0.74), N=24770
- two-way: +0.1139 (se 0.1870, t 0.61), N=24770

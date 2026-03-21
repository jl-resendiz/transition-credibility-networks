# Strategy 2 Panel DiD Metrics

- transform: log1p
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: +0.0150 (se 0.1310, t 0.11), N=24712
- firm-clustered: +0.0150 (se 0.1609, t 0.09), N=24712
- two-way: +0.0150 (se 0.1573, t 0.10), N=24712

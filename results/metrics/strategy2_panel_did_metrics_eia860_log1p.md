# Strategy 2 Panel DiD Metrics

- transform: log1p
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: +0.0153 (se 0.1314, t 0.12), N=24712
- firm-clustered: +0.0153 (se 0.1607, t 0.10), N=24712
- two-way: +0.0153 (se 0.1576, t 0.10), N=24712

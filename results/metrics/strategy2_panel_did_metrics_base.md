# Strategy 2 Panel DiD Metrics

- transform: base
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: +0.0165 (se 0.1251, t 0.13), N=24712
- firm-clustered: +0.0165 (se 0.1529, t 0.11), N=24712
- two-way: +0.0165 (se 0.1488, t 0.11), N=24712

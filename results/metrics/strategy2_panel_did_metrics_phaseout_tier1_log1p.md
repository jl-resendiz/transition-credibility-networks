# Strategy 2 Panel DiD Metrics

- transform: log1p
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: +0.0543 (se 0.0861, t 0.63), N=30034
- firm-clustered: +0.0543 (se 0.1385, t 0.39), N=30034
- two-way: +0.0543 (se 0.1013, t 0.54), N=30034

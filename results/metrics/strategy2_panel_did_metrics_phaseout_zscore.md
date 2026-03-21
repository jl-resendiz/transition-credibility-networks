# Strategy 2 Panel DiD Metrics

- transform: zscore
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: -0.0003 (se 0.0006, t -0.56), N=48656
- firm-clustered: -0.0003 (se 0.0006, t -0.53), N=48656
- two-way: -0.0003 (se 0.0005, t -0.59), N=48656

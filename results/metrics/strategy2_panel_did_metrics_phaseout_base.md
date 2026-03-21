# Strategy 2 Panel DiD Metrics

- transform: base
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: -0.0587 (se 0.1490, t -0.39), N=48656
- firm-clustered: -0.0587 (se 0.1578, t -0.37), N=48656
- two-way: -0.0587 (se 0.1648, t -0.36), N=48656

# Strategy 2 Panel DiD Metrics

- transform: log1p
- event_scope: all_matched
- exact_only: False
- foreign_only: False
- tau: [-6,12]
- post: [0,12]
- overlap: nearest

## exp_post coefficient
- event-clustered: -0.2346 (se 0.0816, t -2.87), N=48656
- firm-clustered: -0.2346 (se 0.0424, t -5.54), N=48656
- two-way: -0.2346 (se 0.0567, t -4.14), N=48656

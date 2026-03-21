# Wild Cluster Bootstrap: Coal Phase-Out DiD

## Specification

AR_{j,t} = firm_FE + month_FE + beta * (coal_share_j * Post_t) + eps

- Exposure: coal_share from firm_alpha_panel (treated firms only; controls = 0)
- Post: months [0, +12] relative to event
- Event window: tau in [-6, 12]
- Overlap rule: nearest
- Tier filter: 1 (binding only: True)

## Observed estimates

- beta(exp_post): -0.027594
- se(cluster): 0.012678
- t-stat: -2.176
- N: 3411
- Clusters (G): 14

## Wild cluster bootstrap

- B: 999
- Weight distribution: Webb 6-point
- Seed: 42
- Bootstrap p-value: 0.0930
- |t*| >= |t_obs|: 92 / 999

### Bootstrap t-stat distribution

- 2.5th percentile: -2.5571
- 5.0th percentile: -2.1346
- 50.0th percentile: +0.0134
- 95.0th percentile: +2.0939
- 97.5th percentile: +2.4975

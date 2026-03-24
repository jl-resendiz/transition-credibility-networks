# Wild Cluster Bootstrap: Coal Phase-Out DiD

## Specification

AR_{j,t} = firm_FE + month_FE + beta * (coal_share_j * Post_t) + eps

- Exposure: coal_share from firm_alpha_panel (treated firms only; controls = 0)
- Post: months [0, +12] relative to event
- Event window: tau in [-6, 12]
- Overlap rule: nearest
- Tier filter: 1 (binding only: True)

## Observed estimates

- beta(exp_post): -0.018143
- se(cluster): 0.008651
- t-stat: -2.097
- N: 3568
- Clusters (G): 14

## Wild cluster bootstrap

- B: 999
- Weight distribution: Webb 6-point
- Seed: 42
- Bootstrap p-value: 0.0910
- |t*| >= |t_obs|: 90 / 999

### Bootstrap t-stat distribution

- 2.5th percentile: -2.1759
- 5.0th percentile: -2.0591
- 50.0th percentile: -0.0024
- 95.0th percentile: +2.0828
- 97.5th percentile: +2.2146

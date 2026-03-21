# Conley Standard Errors: Channel Decomposition

N = 37109 event-firm pairs

| Variable | Coeff | Event-clustered SE (t) | Conley 250km SE (t) | Conley 500km SE (t) | Conley 1000km SE (t) |
|---|---|---|---|---|---|
| w_geo | 0.132 | 0.118 (1.12) | 0.275 (0.48) | 0.266 (0.50) | 0.279 (0.47) |
| w_fuel | -4.578 | 0.833 (-5.50) | 1.192 (-3.84) | 1.114 (-4.11) | 1.100 (-4.16) |
| w_reg | -0.954 | 0.930 (-1.02) | 0.775 (-1.23) | 0.972 (-0.98) | 1.527 (-0.62) |
| same_sector | 0.024 | 0.011 (2.19) | 0.010 (2.33) | 0.012 (2.11) | 0.015 (1.62) |

Notes: Conley (1999) SEs use Bartlett kernel with spherical (Haversine) distances.
Firm coordinates are capacity-weighted plant centroids from GEM/GeoAsset.
Dependent variable: cumulative abnormal return over [-1, +3] months.
Market adjustment: vwretd from Fama-French factors.

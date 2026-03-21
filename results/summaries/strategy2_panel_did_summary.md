# Strategy 2 Panel DiD Summary

**Script:** `strategy2_panel_did.py`  
**Run date:** 2026-02-20  

## Transformations
- Monthly abnormal returns: `AR_it = R_it - vwretd_t` (Fama-French market return).
- Event-time window: default `tau ? [-6, +12]` unless specified.
- Post indicator: default `Post = 1{tau ? [0,12]}`.
- Exposure: `exp_post = w_ij × Post` using geographic W.
- Two-way FE: firm + calendar-month (two-way demeaning).
- Controls: non-connected firms sampled at 5× neighbor count.

## Results (beta on `exp_post`)

### Baseline (nearest overlaps)
| Clustering | Beta | SE | t | N |
|---|---:|---:|---:|---:|
| Event | +0.0228 | 0.0233 | 0.98 | 18,515 |
| Firm | +0.0228 | 0.0212 | 1.08 | 18,515 |
| Two-way (event × firm) | +0.0228 | 0.0182 | 1.25 | 18,515 |

### Overlap rule = drop (tau [-6,+12])
| Clustering | Beta | SE | t | N |
|---|---:|---:|---:|---:|
| Event | +0.1224 | 0.1251 | 0.98 | 4,725 |
| Firm | +0.1224 | 0.1402 | 0.87 | 4,725 |
| Two-way (event × firm) | +0.1224 | 0.1436 | 0.85 | 4,725 |

### Exact-only (nearest overlaps, tau [-6,+12])
| Clustering | Beta | SE | t | N |
|---|---:|---:|---:|---:|
| Event | +0.0051 | 0.0178 | 0.29 | 17,527 |
| Firm | +0.0051 | 0.0157 | 0.33 | 17,527 |
| Two-way (event × firm) | +0.0051 | 0.0120 | 0.43 | 17,527 |

### Short window (nearest overlaps, tau [-3,+6], post [0,6])
| Clustering | Beta | SE | t | N |
|---|---:|---:|---:|---:|
| Event | +0.0089 | 0.0318 | 0.28 | 11,178 |
| Firm | +0.0089 | 0.0229 | 0.39 | 11,178 |
| Two-way (event × firm) | +0.0089 | 0.0223 | 0.40 | 11,178 |

### Exact-only + short window (nearest overlaps, tau [-3,+6], post [0,6])
| Clustering | Beta | SE | t | N |
|---|---:|---:|---:|---:|
| Event | +0.0077 | 0.0219 | 0.35 | 10,805 |
| Firm | +0.0077 | 0.0152 | 0.50 | 10,805 |
| Two-way (event × firm) | +0.0077 | 0.0122 | 0.63 | 10,805 |

Notes: Event scope = first-mover retirements. Announcement dates used where available.

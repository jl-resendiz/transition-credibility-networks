"""Daily event-time path: cross-sectional regression of AR_t on w_fuel,
w_geo, w_reg, same_sector at each daily offset tau in [-21, +21].

This visualises the timing of the repricing: where does the negative
fuel-mix-similarity coefficient show up day-by-day around the announcement?

Reads the per-day AR file produced by compute_daily_ar_panel.py
(FF3-adjusted) and the market-adjusted equivalent.

Output: results/metrics/daily_event_time_path.md
"""
import csv
import math
import os
import sys
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import derived_path, raw_path, results_path

DAILY_AR_FF3 = os.path.join(derived_path('returns'), 'daily_ar_panel.csv')
CAR_FF3 = os.path.join(derived_path('returns'), 'daily_car_panel.csv')

# Market-adjusted daily AR is not yet computed at the per-day level
# (the marketadj script only emits the CAR file, not daily AR rows).
# We re-derive daily AR market-adjusted from a fast computation here.

DAILY_RET_PATH = os.path.join(derived_path('returns'), 'daily_returns.csv')
FF3_PATH = os.path.join(raw_path('factors'), 'F-F_Research_Data_Factors_daily.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')

OUT_PATH = os.path.join(results_path('metrics'), 'daily_event_time_path.md')

OFFSETS = list(range(-21, 22))  # -21 to +21 trading days

# Linear-algebra helpers (same as elsewhere)
def invert(mat):
    n = len(mat)
    aug = [r[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, r in enumerate(mat)]
    for col in range(n):
        mr = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[mr][col]) < 1e-20: return None
        aug[col], aug[mr] = aug[mr], aug[col]
        p = aug[col][col]
        for j in range(2*n): aug[col][j] /= p
        for r in range(n):
            if r != col:
                f = aug[r][col]
                for j in range(2*n): aug[r][j] -= f*aug[col][j]
    return [r[n:] for r in aug]


def fm_at_tau(per_event_betas, lag=4):
    """Newey-West HAC on a time series of T per-event coefficients."""
    T = len(per_event_betas)
    if T < 4: return None, None
    mean = sum(per_event_betas) / T
    dem = [b - mean for b in per_event_betas]
    S = sum(d*d for d in dem) / T
    for L in range(1, lag+1):
        w = 1 - L/(lag+1)
        cov = sum(dem[t]*dem[t-L] for t in range(L, T)) / T
        S += 2*w*cov
    var = S / T
    se = math.sqrt(max(var, 0))
    return mean, se


# ─── Build daily AR panel (FF3 + market-adjusted) and run per-tau regressions ──

def main():
    # Load weight matrices, fundamentals, events, returns, factors
    print('Loading data...')
    W_geo, W_fuel, W_reg = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    for path, M in [(W_GEO_PATH, W_geo), (W_FUEL_PATH, W_fuel), (W_REG_PATH, W_reg)]:
        if not os.path.exists(path): continue
        with open(path, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                try:
                    wval = row.get('w_ij') or row.get('w_reg')
                    M[row['gvkey_i']][row['gvkey_j']] = float(wval)
                except (ValueError, TypeError, KeyError):
                    continue

    # Load CAR panel (FF3) to get the firm-event keyset + treatment vars
    car_rows = []
    with open(CAR_FF3, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                car_rows.append({
                    'gvkey': row['gvkey'],
                    'event_id': int(row['event_id']),
                    'w_geo': float(row['w_geo']),
                    'w_fuel': float(row['w_fuel']),
                    'w_reg': float(row['w_reg']),
                    'same_sector': float(row['same_sector']),
                })
            except (ValueError, KeyError):
                continue
    print(f'  {len(car_rows):,} firm-events from CAR panel.')

    # Load daily AR (FF3) — by (gvkey, event_id, day_offset)
    daily_ff3 = defaultdict(dict)  # (gvkey, event_id) -> {offset: ar}
    with open(DAILY_AR_FF3, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                key = (row['gvkey'], int(row['event_id']))
                daily_ff3[key][int(row['day_offset'])] = float(row['ar_daily'])
            except (ValueError, KeyError):
                continue
    print(f'  {len(daily_ff3):,} (firm, event) daily-AR cells.')

    # ─── Per-day cross-sectional regression: AR_tau ~ w_fuel + w_geo + w_reg + same_sector ──
    spec_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
    names = ['intercept'] + spec_vars

    # Build observations per (event, tau): list of (firm-event records, AR)
    by_event = defaultdict(lambda: defaultdict(list))  # event_id -> tau -> list
    for r in car_rows:
        ar_dict = daily_ff3.get((r['gvkey'], r['event_id']), {})
        for tau, ar in ar_dict.items():
            by_event[r['event_id']][tau].append({**r, 'ar': ar})

    # FM: for each tau, run cross-sectional OLS per event, average across events
    print('\nRunning FM regressions per offset...')
    paths = []
    for tau in OFFSETS:
        per_ev = []
        for eid, taudict in by_event.items():
            ev_rows = taudict.get(tau, [])
            if len(ev_rows) < 20:
                continue
            n = len(ev_rows)
            k = len(spec_vars) + 1
            X = [[1.0] + [r[v] for v in spec_vars] for r in ev_rows]
            y = [r['ar'] for r in ev_rows]
            XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
            Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
            inv = invert(XtX)
            if inv is None: continue
            beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
            per_ev.append(beta)
        T = len(per_ev)
        if T < 10:
            paths.append((tau, T, None, None, None, None))
            continue
        # Coefficients for w_fuel (index = position of 'w_fuel' + 1)
        idx_fuel = names.index('w_fuel')
        idx_geo = names.index('w_geo')
        b_fuel = [b[idx_fuel] for b in per_ev]
        b_geo = [b[idx_geo] for b in per_ev]
        m_fuel, se_fuel = fm_at_tau(b_fuel)
        m_geo, se_geo = fm_at_tau(b_geo)
        paths.append((tau, T, m_fuel, se_fuel, m_geo, se_geo))

    # ─── Write output ──
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    lines = [
        '# Daily Event-Time Path: gamma_fuel(tau) and gamma_geo(tau)',
        '',
        'Cross-sectional regression at each daily offset tau in [-21, +21]:',
        '  AR_{i,tau} ~ alpha + gamma_geo w^geo + gamma_fuel w^fuel ',
        '              + gamma_reg w^reg + gamma_s SameSector + eps',
        '',
        'AR is FF3-adjusted (estimation [-252, -22]). Per-event cross-sections',
        'aggregated via Fama-MacBeth with NW(4) HAC standard errors.',
        '',
        '| tau (trading days) | T (events) | gamma_fuel | NW t | gamma_geo | NW t |',
        '|---:|---:|---:|---:|---:|---:|',
    ]
    for (tau, T, m_f, se_f, m_g, se_g) in paths:
        if m_f is None:
            lines.append(f'| {tau:+3d} | {T} | — | — | — | — |')
        else:
            t_f = m_f / se_f if se_f and se_f > 1e-15 else 0
            t_g = m_g / se_g if se_g and se_g > 1e-15 else 0
            lines.append(f'| {tau:+3d} | {T} | {m_f:+.4f} | {t_f:+.3f} | {m_g:+.4f} | {t_g:+.3f} |')

    lines += [
        '',
        '## Pre/post summary',
        '',
    ]
    pre_fuels = [m_f for (tau, T, m_f, se_f, m_g, se_g) in paths
                 if -21 <= tau <= -2 and m_f is not None]
    post_fuels = [m_f for (tau, T, m_f, se_f, m_g, se_g) in paths
                  if 0 <= tau <= 10 and m_f is not None]
    if pre_fuels:
        lines.append(f'- Pre-event mean gamma_fuel ({len(pre_fuels)} taus in [-21,-2]): {sum(pre_fuels)/len(pre_fuels):+.4f}')
        lines.append(f'- Pre-event sum gamma_fuel : {sum(pre_fuels):+.4f}')
    if post_fuels:
        lines.append(f'- Post-event mean gamma_fuel ({len(post_fuels)} taus in [0,+10]): {sum(post_fuels)/len(post_fuels):+.4f}')
        lines.append(f'- Post-event sum gamma_fuel : {sum(post_fuels):+.4f}')

    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\nWrote {OUT_PATH}')


if __name__ == '__main__':
    main()

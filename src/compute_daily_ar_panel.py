"""Compute daily abnormal returns for the FULL event-firm regression panel
(neighbours + controls, all 565 firms × all 175 events) around each event's
announcement date.

Mirrors the panel-construction logic in `joint_tests.py` and
`two_way_clustering.py` to ensure the daily-AR sample matches the monthly-CAR
sample one-to-one (modulo daily-availability dropouts).

Estimation window: [-252, -22] trading days BEFORE announcement_date
(approximately 1 year of daily returns, ending 22 days before the event).
Event window: [-5, +20] trading days.

Factor model: Fama-French 3 (Mkt-RF, SMB, HML) using the daily FF3 file
already downloaded to data/raw/factors/F-F_Research_Data_Factors_daily.csv.
For non-US firms we acknowledge the FF3 factors are US-priced (a documented
limitation that already applies to the monthly multi-factor analysis).

Output:
  data/derived/returns/daily_car_panel.csv
    cols: gvkey, event_id, w_geo, w_fuel, w_reg, same_sector,
          car_m1_p1, car_0_p5, car_0_p10, car_m1_p10, car_0_p20

Also includes the daily AR series:
  data/derived/returns/daily_ar_panel.csv
    cols: gvkey, event_id, day_offset, date, ar_daily

The output preserves the (event_id, gvkey) key structure used by
joint_tests so that downstream daily-event-study scripts can
plug in identical cross-sectional regressions.
"""
import csv
import hashlib
import os
import random
import sys
from collections import defaultdict
from datetime import datetime, timedelta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import raw_path, derived_path

DAILY_RET_PATH = os.path.join(derived_path('returns'), 'daily_returns.csv')
FF3_PATH = os.path.join(raw_path('factors'), 'F-F_Research_Data_Factors_daily.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')
FUND_PATH = os.path.join(derived_path('fundamentals'), 'firm_fundamentals.csv')

OUT_DAILY = os.path.join(derived_path('returns'), 'daily_ar_panel.csv')
OUT_CAR = os.path.join(derived_path('returns'), 'daily_car_panel.csv')

EST_PRE_DAYS = 252
EST_GAP_DAYS = 22
EST_MIN_DAYS = 100  # minimum estimation observations
EVENT_PRE = 22       # extended to capture pre-event drift / anticipation
EVENT_POST = 22


# ─── Linear-algebra helpers (stdlib OLS) ──────────────────────────────

def ols(y, X):
    n = len(y)
    if n < 12:
        return None
    k = len(X[0])
    XtX = [[0.0] * k for _ in range(k)]
    Xty = [0.0] * k
    for i in range(n):
        for a in range(k):
            Xty[a] += X[i][a] * y[i]
            for b in range(a, k):
                XtX[a][b] += X[i][a] * X[i][b]
    for a in range(k):
        for b in range(a):
            XtX[a][b] = XtX[b][a]
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(k)] for i, row in enumerate(XtX)]
    for col in range(k):
        max_r = max(range(col, k), key=lambda r: abs(aug[r][col]))
        if abs(aug[max_r][col]) < 1e-14:
            return None
        aug[col], aug[max_r] = aug[max_r], aug[col]
        piv = aug[col][col]
        for j in range(2 * k):
            aug[col][j] /= piv
        for r in range(k):
            if r != col:
                f = aug[r][col]
                for j in range(2 * k):
                    aug[r][j] -= f * aug[col][j]
    inv = [r[k:] for r in aug]
    return [sum(inv[a][b] * Xty[b] for b in range(k)) for a in range(k)]


# ─── Date helpers ─────────────────────────────────────────────────────

def parse_date(s):
    return datetime.strptime(s[:10], '%Y-%m-%d').date()


# ─── Data loaders ─────────────────────────────────────────────────────

def load_ff3():
    factors = {}
    with open(FF3_PATH, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit():
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5:
                continue
            try:
                d = datetime.strptime(parts[0], '%Y%m%d').date()
                mktrf = float(parts[1]) / 100.0
                smb = float(parts[2]) / 100.0
                hml = float(parts[3]) / 100.0
                rf = float(parts[4]) / 100.0
            except ValueError:
                continue
            factors[d] = (mktrf, smb, hml, rf)
    return factors


def load_daily_returns():
    """{gvkey: {date: ret}}"""
    by_gvkey = defaultdict(dict)
    with open(DAILY_RET_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                gvkey = str(row['gvkey']).strip().split('.')[0].zfill(6)
                d = parse_date(row['datadate'])
                ret = float(row['ret_daily'])
            except (ValueError, KeyError):
                continue
            by_gvkey[gvkey][d] = ret
    return by_gvkey


def load_weight_matrix(path):
    """{gvkey_i: {gvkey_j: w}}"""
    M = defaultdict(dict)
    if not os.path.exists(path):
        return M
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                wval = row.get('w_ij') or row.get('w_reg')
                M[row['gvkey_i']][row['gvkey_j']] = float(wval)
            except (ValueError, TypeError, KeyError):
                continue
    return M


def load_fundamentals_sic():
    """{gvkey: sic4}"""
    sic = {}
    with open(FUND_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gk = row['gvkey']
            s = row.get('sic')
            if s and gk not in sic:
                sic[gk] = s[:4]
    return sic


def load_events():
    """First-mover events with announcement dates (precise YYYY-MM-DD)."""
    events = []
    with open(EVENTS_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if row.get('is_first_mover') != 'True':
                continue
            ann = row.get('announcement_date', '').strip()
            if len(ann) < 10:
                continue
            try:
                ann_date = parse_date(ann)
            except ValueError:
                continue
            matched = row.get('matched_gvkeys', '').strip()
            if not matched:
                continue
            events.append({
                'event_id': i,
                'announcement_date': ann_date,
                'gvkeys': [g.strip().zfill(6) for g in matched.split(';')],
                'plant': row.get('plant_name', ''),
            })
    return events


# ─── Main ─────────────────────────────────────────────────────────────

def main():
    print('Loading FF3 daily factors...')
    factors = load_ff3()
    print(f'  {len(factors):,} dates.')

    print('Loading daily returns (Compustat panel)...')
    rets = load_daily_returns()
    print(f'  {len(rets)} firms with daily returns.')

    print('Loading weight matrices...')
    W_geo = load_weight_matrix(W_GEO_PATH)
    W_fuel = load_weight_matrix(W_FUEL_PATH)
    W_reg = load_weight_matrix(W_REG_PATH)
    print(f'  W_geo: {len(W_geo)} firms, W_fuel: {len(W_fuel)}, W_reg: {len(W_reg)}')

    print('Loading firm SIC codes...')
    firm_sic = load_fundamentals_sic()
    print(f'  {len(firm_sic)} firms with SIC.')

    print('Loading first-mover events...')
    events = load_events()
    print(f'  {len(events)} first-mover events with announcement_date.')

    # Build candidate panels event-by-event (mirroring joint_tests.py)
    print('\nBuilding candidate panels...')
    universe_gvkeys = list(set(rets.keys()) | set(firm_sic.keys()))

    car_rows = []
    daily_rows = []
    failed_events = 0
    failed_firms_no_data = 0
    failed_firms_short_window = 0
    n_processed = 0

    # Pre-build ordered date arrays for speed: sorted FF3 trading days
    ff3_dates = sorted(factors.keys())

    for ev in events:
        t0 = ev['announcement_date']
        event_gvkeys = set(ev['gvkeys'])

        # First-mover SIC
        fm_sic4 = None
        for gk in event_gvkeys:
            fm_sic4 = firm_sic.get(gk)
            if fm_sic4:
                break

        # Build the list of candidate firms following the same logic
        # as joint_tests.py / two_way_clustering.py
        candidates_for_event = []
        for fm_gk in event_gvkeys:
            if fm_gk not in W_geo:
                continue
            neighbors = W_geo[fm_gk]
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in universe_gvkeys
                             if gk not in event_gvkeys and gk not in neighbors]
            stable_seed = int(hashlib.md5(
                str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
            random.seed(stable_seed)
            n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
            ctrl_sample = (random.sample(non_connected, n_ctrl)
                           if len(non_connected) > n_ctrl else non_connected)
            for gk in list(neighbor_gks) + ctrl_sample:
                w_geo = neighbors.get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
                j_sic = firm_sic.get(gk)
                same_sector = 1.0 if (fm_sic4 and j_sic and fm_sic4 == j_sic) else 0.0
                candidates_for_event.append({
                    'gvkey': gk, 'w_geo': w_geo, 'w_fuel': w_fuel, 'w_reg': w_reg,
                    'same_sector': same_sector,
                })

        # Find event date in trading-day index (use ff3_dates as canonical)
        idx_t0 = next((i for i, d in enumerate(ff3_dates) if d >= t0), None)
        if idx_t0 is None or idx_t0 < EST_PRE_DAYS:
            failed_events += 1
            continue

        est_start_idx = max(0, idx_t0 - EST_PRE_DAYS - EST_GAP_DAYS)
        est_end_idx = idx_t0 - EST_GAP_DAYS
        est_dates = ff3_dates[est_start_idx:est_end_idx]
        evt_dates = ff3_dates[max(0, idx_t0 - EVENT_PRE): idx_t0 + EVENT_POST + 1]

        for cand in candidates_for_event:
            gk = cand['gvkey']
            firm_rets = rets.get(gk)
            if not firm_rets:
                failed_firms_no_data += 1
                continue

            # Build estimation y, X
            y, X = [], []
            for d in est_dates:
                if d in firm_rets and d in factors:
                    mktrf, smb, hml, rf = factors[d]
                    y.append(firm_rets[d] - rf)
                    X.append([1.0, mktrf, smb, hml])
            if len(y) < EST_MIN_DAYS:
                failed_firms_short_window += 1
                continue
            beta = ols(y, X)
            if beta is None:
                continue
            a, bM, bS, bH = beta

            # Compute event-window ARs
            cum = {'m1_p1': 0.0, '0_p5': 0.0, '0_p10': 0.0, 'm1_p10': 0.0, '0_p20': 0.0}
            cnt = {k: 0 for k in cum}
            for d in evt_dates:
                if d not in firm_rets or d not in factors:
                    continue
                offset = ff3_dates.index(d) - idx_t0  # cheap because dates list is small
                mktrf, smb, hml, rf = factors[d]
                ar = (firm_rets[d] - rf) - (a + bM * mktrf + bS * smb + bH * hml)
                daily_rows.append({
                    'gvkey': gk, 'event_id': ev['event_id'],
                    'day_offset': offset, 'date': d.strftime('%Y-%m-%d'),
                    'ar_daily': ar,
                })
                if -1 <= offset <= 1:  cum['m1_p1']  += ar; cnt['m1_p1']  += 1
                if  0 <= offset <= 5:  cum['0_p5']   += ar; cnt['0_p5']   += 1
                if  0 <= offset <= 10: cum['0_p10']  += ar; cnt['0_p10']  += 1
                if -1 <= offset <= 10: cum['m1_p10'] += ar; cnt['m1_p10'] += 1
                if  0 <= offset <= 20: cum['0_p20']  += ar; cnt['0_p20']  += 1

            row = {
                'gvkey': gk, 'event_id': ev['event_id'],
                'w_geo': cand['w_geo'], 'w_fuel': cand['w_fuel'],
                'w_reg': cand['w_reg'], 'same_sector': cand['same_sector'],
                'car_m1_p1':  cum['m1_p1']  if cnt['m1_p1']  >= 2 else None,
                'car_0_p5':   cum['0_p5']   if cnt['0_p5']   >= 4 else None,
                'car_0_p10':  cum['0_p10']  if cnt['0_p10']  >= 8 else None,
                'car_m1_p10': cum['m1_p10'] if cnt['m1_p10'] >= 9 else None,
                'car_0_p20':  cum['0_p20']  if cnt['0_p20']  >= 16 else None,
            }
            car_rows.append(row)
            n_processed += 1
            if n_processed % 5000 == 0:
                print(f'  {n_processed:,} firm-events processed; '
                      f'last event {ev["event_id"]} date {t0}')

    print('\n=== Summary ===')
    print(f'Events: {len(events)}')
    print(f'Failed events (no pre-window): {failed_events}')
    print(f'Firm-events with no daily data : {failed_firms_no_data:,}')
    print(f'Firm-events with short window  : {failed_firms_short_window:,}')
    print(f'Firm-events processed (CAR rows): {len(car_rows):,}')
    print(f'Daily-AR rows                   : {len(daily_rows):,}')

    os.makedirs(os.path.dirname(OUT_CAR), exist_ok=True)
    with open(OUT_CAR, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(car_rows[0].keys()))
        w.writeheader(); w.writerows(car_rows)
    print(f'Wrote {OUT_CAR}')

    with open(OUT_DAILY, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gvkey','event_id','day_offset','date','ar_daily'])
        w.writeheader(); w.writerows(daily_rows)
    print(f'Wrote {OUT_DAILY}')


if __name__ == '__main__':
    main()

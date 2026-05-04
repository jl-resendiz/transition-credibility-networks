"""Compute daily abnormal returns using MARKET-ADJUSTED (vwretd) returns,
mirroring the monthly headline specification (no firm-specific betas).

This is the H3 diagnostic recommended by /econometrics: the FF3 daily AR
spec correlates the factor loadings with the w_fuel treatment, biasing
the cross-sectional regression. Market-adjusted returns avoid this by
using only the value-weighted market return as a benchmark.

Spec:
    AR_{i,t} = r_{i,t} - vwretd_t          (no firm-specific beta)
    pre_mean_AR_i = mean( AR over [-252, -22] daily )
    CAR_{i, [tau1, tau2]} = sum_{t in window} (AR_{i,t} - pre_mean_AR_i)

vwretd is taken from FF3 daily as Mkt-RF + RF, matching the monthly
construction in compute_returns.py.

Output: data/derived/returns/daily_car_marketadj_panel.csv
"""
import csv
import hashlib
import os
import random
import sys
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import raw_path, derived_path

DAILY_RET_PATH = os.path.join(derived_path('returns'), 'daily_returns.csv')
FF3_PATH = os.path.join(raw_path('factors'), 'F-F_Research_Data_Factors_daily.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')
FUND_PATH = os.path.join(derived_path('fundamentals'), 'firm_fundamentals.csv')

OUT_CAR = os.path.join(derived_path('returns'), 'daily_car_marketadj_panel.csv')

EST_PRE_DAYS = 252
EST_GAP_DAYS = 22
EST_MIN_DAYS = 100
EVENT_PRE = 5
EVENT_POST = 20


def parse_date(s):
    return datetime.strptime(s[:10], '%Y-%m-%d').date()


def load_market_index_from_ff3():
    """vwretd_t = MktRF_t + RF_t  (matches monthly construction)."""
    vw = {}
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
                rf = float(parts[4]) / 100.0
            except ValueError:
                continue
            vw[d] = mktrf + rf
    return vw


def load_daily_returns():
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
            })
    return events


def main():
    print('Loading vwretd (FF3 Mkt-RF + RF, daily)...')
    vw = load_market_index_from_ff3()
    print(f'  {len(vw):,} dates.')

    print('Loading daily returns (Compustat panel)...')
    rets = load_daily_returns()
    print(f'  {len(rets)} firms.')

    print('Loading weight matrices...')
    W_geo = load_weight_matrix(W_GEO_PATH)
    W_fuel = load_weight_matrix(W_FUEL_PATH)
    W_reg = load_weight_matrix(W_REG_PATH)

    print('Loading firm SIC codes...')
    firm_sic = load_fundamentals_sic()

    print('Loading events...')
    events = load_events()
    print(f'  {len(events)} first-mover events with announcement_date.')

    universe_gvkeys = list(set(rets.keys()) | set(firm_sic.keys()))
    vw_dates = sorted(vw.keys())

    car_rows = []
    n_proc = 0
    skipped_no_data = 0
    skipped_short_window = 0

    for ev in events:
        t0 = ev['announcement_date']
        event_gvkeys = set(ev['gvkeys'])
        fm_sic4 = next((firm_sic.get(gk) for gk in event_gvkeys if firm_sic.get(gk)), None)

        # Build candidates
        candidates = []
        for fm_gk in event_gvkeys:
            if fm_gk not in W_geo:
                continue
            neighbors = W_geo[fm_gk]
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in universe_gvkeys
                             if gk not in event_gvkeys and gk not in neighbors]
            seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
            random.seed(seed)
            n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
            ctrl_sample = (random.sample(non_connected, n_ctrl)
                           if len(non_connected) > n_ctrl else non_connected)
            for gk in list(neighbor_gks) + ctrl_sample:
                w_geo = neighbors.get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
                j_sic = firm_sic.get(gk)
                same_sector = 1.0 if (fm_sic4 and j_sic and fm_sic4 == j_sic) else 0.0
                candidates.append({
                    'gvkey': gk, 'w_geo': w_geo, 'w_fuel': w_fuel,
                    'w_reg': w_reg, 'same_sector': same_sector,
                })

        idx_t0 = next((i for i, d in enumerate(vw_dates) if d >= t0), None)
        if idx_t0 is None or idx_t0 < EST_PRE_DAYS:
            continue

        est_start_idx = max(0, idx_t0 - EST_PRE_DAYS - EST_GAP_DAYS)
        est_end_idx = idx_t0 - EST_GAP_DAYS
        est_dates = vw_dates[est_start_idx:est_end_idx]
        evt_dates = vw_dates[max(0, idx_t0 - EVENT_PRE): idx_t0 + EVENT_POST + 1]

        for cand in candidates:
            gk = cand['gvkey']
            firm_rets = rets.get(gk)
            if not firm_rets:
                skipped_no_data += 1
                continue

            # Pre-event mean of (r - vwretd)
            pre_ars = []
            for d in est_dates:
                if d in firm_rets and d in vw:
                    pre_ars.append(firm_rets[d] - vw[d])
            if len(pre_ars) < EST_MIN_DAYS:
                skipped_short_window += 1
                continue
            pre_mean = sum(pre_ars) / len(pre_ars)

            # Event-window CAR (subtract pre_mean per day)
            cum = {'m1_p1': 0.0, '0_p5': 0.0, '0_p10': 0.0, 'm1_p10': 0.0, '0_p20': 0.0}
            cnt = {k: 0 for k in cum}
            for d in evt_dates:
                if d not in firm_rets or d not in vw:
                    continue
                offset = vw_dates.index(d) - idx_t0
                ar = (firm_rets[d] - vw[d]) - pre_mean
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
            n_proc += 1
            if n_proc % 5000 == 0:
                print(f'  {n_proc:,} firm-events; last event {ev["event_id"]} {t0}')

    print(f'\n=== Summary ===')
    print(f'CAR rows: {len(car_rows):,}')
    print(f'Skipped (no daily data): {skipped_no_data:,}')
    print(f'Skipped (short window) : {skipped_short_window:,}')

    os.makedirs(os.path.dirname(OUT_CAR), exist_ok=True)
    with open(OUT_CAR, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(car_rows[0].keys()))
        w.writeheader(); w.writerows(car_rows)
    print(f'Wrote {OUT_CAR}')


if __name__ == '__main__':
    main()

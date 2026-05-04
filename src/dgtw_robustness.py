"""DGTW-style characteristic-matched robustness for the US-linked sub-sample.

Adapted Daniel-Grinblatt-Titman-Wermers (1997) implementation:
- For each (firm, month-end), classify into within-month tercile of:
    * size (log market cap)
    * book-to-market (log BM)
    * 12-month momentum (months [t-12, t-2])
- Build benchmark monthly return per (size_t, bm_t, mom_t) bucket
  as the equal-weighted mean of all firms in the bucket EXCLUDING
  the focal firm.
- DGTW-adjusted return = firm_ret - benchmark_ret.
- Compute DGTW-adjusted CAR over the [-1, +3] monthly headline window.
- Run cross-sectional regression on DGTW-adjusted CARs.

If the channel survives DGTW adjustment, it is not absorbed by
size/value/momentum confounds. With 87 US firms in 27 buckets, cells
will be thin (~3 firms each) — a known limitation we acknowledge.

Inputs:
  - data/derived/dgtw/dgtw_chars_us.csv  (firm × month × terciles)
  - data/derived/returns/monthly_returns.csv (firm × month × ret)
  - data/derived/events/coal_retirement_events.csv

Output:
  - results/metrics/dgtw_robustness.md
"""
import csv
import hashlib
import math
import os
import random
import sys
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import raw_path, derived_path, results_path

CHARS_PATH = os.path.join(derived_path('dgtw'), 'dgtw_chars_us.csv')
RET_PATH = os.path.join(derived_path('returns'), 'monthly_returns.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')
FUND_PATH = os.path.join(derived_path('fundamentals'), 'firm_fundamentals.csv')

OUT_PATH = os.path.join(results_path('metrics'), 'dgtw_robustness.md')

POST_MONTHS = 3
PRE_MONTHS = 24


def parse_date(s):
    return datetime.strptime(s[:10], '%Y-%m-%d').date()


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


def fm_with_nw(per_event_betas, lag=4):
    T = len(per_event_betas)
    if T < 4: return None
    k = len(per_event_betas[0])
    means = [sum(b[a] for b in per_event_betas) / T for a in range(k)]
    nw_se = []
    for a in range(k):
        x = [per_event_betas[t][a] - means[a] for t in range(T)]
        S = sum(xt*xt for xt in x) / T
        for L in range(1, lag+1):
            w = 1 - L/(lag+1)
            cov = sum(x[t]*x[t-L] for t in range(L, T)) / T
            S += 2*w*cov
        nw_se.append(math.sqrt(max(S/T, 0)))
    return {'mean': means, 'nw_se': nw_se, 'T': T}


def main():
    # Load characteristic deciles + monthly returns
    print('Loading DGTW characteristics...')
    chars = defaultdict(dict)  # gvkey -> month_str -> (size_t, bm_t, mom_t)
    with open(CHARS_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                gk = str(row['gvkey']).split('.')[0].zfill(6)
                me = row['month_end']
                month_str = me[:7]
                size_t = int(row['size_t']) if row['size_t'] not in ('', 'None') else None
                bm_t = int(row['bm_t']) if row['bm_t'] not in ('', 'None') else None
                mom_t = int(row['mom_t']) if row['mom_t'] not in ('', 'None') else None
                if size_t and bm_t and mom_t:
                    chars[gk][month_str] = (size_t, bm_t, mom_t)
            except (ValueError, KeyError):
                continue
    print(f'  {len(chars)} firms with DGTW chars; '
          f'{sum(len(v) for v in chars.values())} firm-months total.')

    print('\nLoading monthly returns...')
    rets = defaultdict(dict)
    with open(RET_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                gk = str(row['gvkey']).split('.')[0].zfill(6)
                ym = row['datadate'][:7]
                ret = float(row['ret_monthly'])
            except (ValueError, KeyError):
                continue
            rets[gk][ym] = ret

    # Build benchmarks: for each month and (size_t, bm_t, mom_t), mean return
    # of all firms in that bucket. Note: we exclude the focal firm at usage time.
    print('\nBuilding DGTW benchmarks (per month × bucket)...')
    bucket_returns = defaultdict(lambda: defaultdict(list))  # ym -> bucket -> [(gvkey, ret)]
    for gk, mmap in chars.items():
        for ym, bucket in mmap.items():
            if ym in rets.get(gk, {}):
                bucket_returns[ym][bucket].append((gk, rets[gk][ym]))

    # Compute DGTW-adjusted returns: AR = firm_ret - mean(bucket excluding self)
    print('\nComputing DGTW-adjusted monthly returns...')
    dgtw_ar = defaultdict(dict)  # gvkey -> ym -> dgtw_adjusted_return
    for gk, mmap in chars.items():
        for ym, bucket in mmap.items():
            if ym not in rets.get(gk, {}):
                continue
            firm_ret = rets[gk][ym]
            bucket_obs = bucket_returns[ym][bucket]
            others = [r for (g, r) in bucket_obs if g != gk]
            if len(others) < 2:  # need at least 2 other firms in bucket
                continue
            bench = sum(others) / len(others)
            dgtw_ar[gk][ym] = firm_ret - bench

    n_dgtw = sum(len(v) for v in dgtw_ar.values())
    print(f'  {n_dgtw:,} DGTW-adjusted firm-months computed.')

    # ─── Build event-firm panel using same logic as joint_tests.py ──
    # Restricted to firms that HAVE DGTW-adjusted returns (i.e., US-linked
    # with characteristic data).
    print('\nLoading network matrices, events, fundamentals...')
    W_geo = defaultdict(dict)
    W_fuel = defaultdict(dict)
    W_reg = defaultdict(dict)
    for path, M in [(W_GEO_PATH, W_geo), (W_FUEL_PATH, W_fuel), (W_REG_PATH, W_reg)]:
        if not os.path.exists(path): continue
        with open(path, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                try:
                    wval = row.get('w_ij') or row.get('w_reg')
                    M[row['gvkey_i']][row['gvkey_j']] = float(wval)
                except (ValueError, TypeError, KeyError):
                    continue

    firm_sic = {}
    with open(FUND_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']; s = row.get('sic')
            if s and gk not in firm_sic: firm_sic[gk] = s[:4]

    events = []
    with open(EVENTS_PATH, 'r', encoding='utf-8') as f:
        for i, row in enumerate(csv.DictReader(f)):
            if row.get('is_first_mover') != 'True': continue
            ann = row.get('announcement_date', '').strip()
            if len(ann) < 10: continue
            try:
                ann_date = parse_date(ann)
            except ValueError:
                continue
            matched = row.get('matched_gvkeys', '').strip()
            if not matched: continue
            events.append({
                'event_id': i,
                'ann_date': ann_date,
                'gvkeys': [g.strip().zfill(6) for g in matched.split(';')],
            })
    print(f'  {len(events)} first-mover events.')

    universe = list(set(rets.keys()) | set(firm_sic.keys()))

    # Compute DGTW-adjusted CAR for [-1, +3] window
    def dgtw_car(gk, event_month, post=POST_MONTHS, pre=PRE_MONTHS):
        if gk not in dgtw_ar: return None
        all_months = sorted(dgtw_ar[gk].keys())
        idx = next((i for i, m in enumerate(all_months) if m >= event_month), None)
        if idx is None: return None
        # pre-event mean of dgtw_ar
        pre_idx = max(0, idx - pre)
        pre_ar = [dgtw_ar[gk][m] for m in all_months[pre_idx:idx]]
        if len(pre_ar) < 6: return None
        pre_mean = sum(pre_ar) / len(pre_ar)
        car = 0.0
        cnt = 0
        for offset in range(-1, post+1):
            i2 = idx + offset
            if 0 <= i2 < len(all_months):
                m = all_months[i2]
                if m in dgtw_ar[gk]:
                    car += dgtw_ar[gk][m] - pre_mean
                    cnt += 1
        return car if cnt >= 3 else None

    print('\nBuilding DGTW-adjusted CAR panel...')
    panel = []
    for ev in events:
        ann_month = ev['ann_date'].strftime('%Y-%m')
        event_gvkeys = set(ev['gvkeys'])
        fm_sic4 = next((firm_sic.get(gk) for gk in event_gvkeys if firm_sic.get(gk)), None)

        for fm_gk in event_gvkeys:
            if fm_gk not in W_geo: continue
            neighbors = W_geo[fm_gk]
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in universe
                             if gk not in event_gvkeys and gk not in neighbors]
            seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
            random.seed(seed)
            n_ctrl = min(len(non_connected), max(5*len(neighbor_gks), 20))
            ctrl_sample = (random.sample(non_connected, n_ctrl)
                           if len(non_connected) > n_ctrl else non_connected)
            for gk in list(neighbor_gks) + ctrl_sample:
                car = dgtw_car(gk, ann_month)
                if car is None: continue
                panel.append({
                    'event_id': ev['event_id'], 'gvkey': gk, 'car': car,
                    'w_geo': neighbors.get(gk, 0.0),
                    'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                    'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                    'same_sector': 1.0 if (fm_sic4 and firm_sic.get(gk) == fm_sic4) else 0.0,
                })
    print(f'  Panel size (DGTW-eligible firm-events): {len(panel):,}')
    print(f'  Unique gvkeys: {len(set(r["gvkey"] for r in panel))}')
    print(f'  Unique events: {len(set(r["event_id"] for r in panel))}')

    # ─── FM regression on DGTW-adjusted CARs ──
    spec_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

    by_event = defaultdict(list)
    for r in panel:
        by_event[r['event_id']].append(r)

    per_event_betas = []
    event_ns = []
    for eid in sorted(by_event.keys()):
        evrows = by_event[eid]
        if len(evrows) < 8: continue
        n = len(evrows); k = len(spec_vars) + 1
        X = [[1.0] + [r[v] for v in spec_vars] for r in evrows]
        y = [r['car'] for r in evrows]
        XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
        inv = invert(XtX)
        if inv is None: continue
        beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
        per_event_betas.append(beta)
        event_ns.append(n)

    print(f'\n{len(per_event_betas)} events with successful FM regression')
    fm = fm_with_nw(per_event_betas, lag=4)
    if fm is None:
        sys.exit('FM aggregation failed.')

    print('\n=== DGTW-adjusted FM result ===')
    names = ['intercept'] + spec_vars
    for i, nm in enumerate(names):
        m = fm['mean'][i]; se = fm['nw_se'][i]
        t = m/se if se > 1e-15 else 0
        print(f'  {nm:12s}: mean = {m:+.4f}, SE_NW = {se:.4f}, t = {t:+.3f}')

    # ─── Write output ──
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    lines = [
        '# DGTW Characteristic-Matched Robustness',
        '',
        'Daniel-Grinblatt-Titman-Wermers (1997)-style adjustment using',
        'within-month tercile sorts on (size, B/M, momentum). For each',
        'firm-month, the DGTW-adjusted return is firm_ret minus the mean',
        'of all OTHER firms in the same (size_t × bm_t × mom_t) bucket',
        'at that month.',
        '',
        'Note on sample: this adjustment is restricted to the US-linked',
        'sub-sample (87 firms with CRSP/Compustat characteristic data).',
        'Cell sizes are thin (~3 firms per 27 buckets per month) — a',
        'known limitation acknowledged here. Tercile sorts (rather than',
        'quintile sorts) used to mitigate cell-thinness.',
        '',
        '## Headline FM result (NW lag 4)',
        '',
        f'Events with successful regression: T = {fm["T"]}',
        f'Avg firms per event: {sum(event_ns)/len(event_ns):.1f}',
        '',
        '| Variable | Mean | NW SE | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for i, nm in enumerate(names):
        m = fm['mean'][i]; se = fm['nw_se'][i]
        t = m/se if se > 1e-15 else 0
        p = 2*(1 - 0.5*(1 + math.erf(abs(t)/math.sqrt(2))))
        stars = '***' if p<0.01 else '**' if p<0.05 else '*' if p<0.10 else ''
        lines.append(f'| {nm} | {m:+.4f} | {se:.4f} | {t:+.3f} | {p:.4f}{stars} |')

    lines += [
        '',
        '## Interpretation',
        '',
        'A negative gamma_fuel under DGTW adjustment indicates that the',
        'channel survives controls for size, B/M, and momentum confounds —',
        'specifically the worry that coal-heavy peers earn higher returns',
        'because they are simultaneously small-cap, high-B/M, and low-momentum',
        '("brown" stocks).',
        '',
        'A null or sign-flipped gamma_fuel would suggest the channel is',
        'partly absorbed by characteristic risk premia.',
        '',
        '## Caveats',
        '',
        '- Cell thinness: 27 buckets × 87 firms means typical bucket has',
        '  ~3 firms. Benchmark returns are noisy.',
        '- Sample restricted to US-linked firms; non-US firms (the bulk',
        '  of the channel) cannot be DGTW-adjusted without comparable',
        '  characteristic data, which is not in the WRDS subscription.',
        '- Tercile (not quintile) sorts: cuts noise but at the cost of',
        '  finer characteristic resolution.',
        '',
    ]

    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\nWrote {OUT_PATH}')


if __name__ == '__main__':
    main()

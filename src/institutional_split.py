"""Institutional-ownership split test: does the fuel-mix channel exist
within HIGH vs LOW institutional-ownership concentration tiers?

Per /quant-finance recommendation: split on Herfindahl (HHI) of 13F manager
shares-of-shares at the EVENT QUARTER (most recent 13F snapshot prior to
announcement; no look-ahead). Terciles within each event quarter (so the
split absorbs time trends in 13F filer counts).

For each US-linked firm, look up its inst-ownership HHI at the snapshot
quarter just before the event date, then sort firms within each event
into HHI terciles. Run cross-sectional regression separately on T3
(concentrated) and T1 (dispersed). Report difference of fuel-mix
coefficients.

Inputs:
  - data/derived/institutional/institutional_ownership.csv
  - data/raw/wrds/ccm_link_us.csv (gvkey ↔ permno mapping for US firms)
  - data/derived/events/coal_retirement_events.csv (announcement_date)
  - existing event-firm panel: rebuild from scratch in this script using same
    logic as joint_tests.py to maintain consistency

Output:
  results/metrics/institutional_split.md
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

INST_PATH = os.path.join(derived_path('institutional'), 'institutional_ownership.csv')
DAILY_RET_PATH = os.path.join(derived_path('returns'), 'monthly_returns.csv')  # monthly for the headline
FF3_PATH = os.path.join(raw_path('factors'), 'F-F_Research_Data_Factors.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')
FUND_PATH = os.path.join(derived_path('fundamentals'), 'firm_fundamentals.csv')

OUT_PATH = os.path.join(results_path('metrics'), 'institutional_split.md')

POST_MONTHS = 3
PRE_MONTHS = 24


def parse_date(s):
    return datetime.strptime(s[:10], '%Y-%m-%d').date()


# ─── OLS / linear-algebra helpers ─────────────────────────────────────

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
    """Per-event betas (list of vectors) → FM mean + NW(lag) HAC SEs."""
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
        var = S / T
        nw_se.append(math.sqrt(max(var, 0)))
    return {'mean': means, 'nw_se': nw_se, 'T': T}


# ─── Build CAR panel mirroring joint_tests.py ──────────────────────────

def load_market_index_monthly():
    """vwretd_t = MktRF_t + RF_t monthly."""
    vw = {}
    with open(FF3_PATH, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit():
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5: continue
            try:
                d = parts[0]
                if not d.isdigit() or len(d) != 6: continue
                mktrf = float(parts[1]) / 100.0
                rf = float(parts[4]) / 100.0
            except ValueError:
                continue
            vw[f'{d[:4]}-{d[4:6]}'] = mktrf + rf
    return vw


def load_monthly_returns():
    by_gvkey = defaultdict(dict)
    with open(DAILY_RET_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                gk = str(row['gvkey']).split('.')[0].zfill(6)
                ym = row['datadate'][:7]
                ret = float(row['ret_monthly'])
            except (ValueError, KeyError):
                continue
            by_gvkey[gk][ym] = ret
    return by_gvkey


def load_weight_matrix(path):
    M = defaultdict(dict)
    if not os.path.exists(path): return M
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                wval = row.get('w_ij') or row.get('w_reg')
                M[row['gvkey_i']][row['gvkey_j']] = float(wval)
            except (ValueError, TypeError, KeyError):
                continue
    return M


def load_firm_sic():
    sic = {}
    with open(FUND_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            s = row.get('sic')
            if s and gk not in sic: sic[gk] = s[:4]
    return sic


def load_events():
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
                'announcement_date': ann_date,
                'gvkeys': [g.strip().zfill(6) for g in matched.split(';')],
            })
    return events


def load_inst_ownership():
    """{gvkey: sorted [(fdate, hhi, n_managers, inst_pct, top5)]}"""
    by_gvkey = defaultdict(list)
    with open(INST_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                gk = str(row['gvkey']).split('.')[0].zfill(6)
                d = parse_date(row['fdate'])
                hhi = float(row['hhi_managers'])
                nm = int(row['n_managers'])
                ipc = float(row['inst_pct_capped'])
                t5 = float(row['top5_concentration'])
            except (ValueError, KeyError):
                continue
            by_gvkey[gk].append((d, hhi, nm, ipc, t5))
    for gk in by_gvkey:
        by_gvkey[gk].sort()
    return by_gvkey


def get_inst_at(gvkey, target_date, inst):
    """Most recent inst quarter on or before target_date."""
    history = inst.get(gvkey, [])
    eligible = [t for t in history if t[0] <= target_date]
    return eligible[-1] if eligible else None


def compute_monthly_car(gvkey, event_month, returns, vw, post=POST_MONTHS, pre=PRE_MONTHS):
    if gvkey not in returns: return None
    months = sorted(returns[gvkey].keys())
    idx = next((i for i, m in enumerate(months) if m >= event_month), None)
    if idx is None: return None
    pre_idx = max(0, idx - pre)
    pre_rets = [returns[gvkey][m] for m in months[pre_idx:idx]]
    if len(pre_rets) < 12: return None
    ar_list = []
    for m in months[pre_idx:idx]:
        if m in returns[gvkey] and m in vw:
            ar_list.append(returns[gvkey][m] - vw[m])
    pre_mean = sum(ar_list) / len(ar_list) if ar_list else 0.0
    car = 0.0
    for offset in range(-1, post + 1):
        i2 = idx + offset
        if 0 <= i2 < len(months):
            m = months[i2]
            if m in returns[gvkey] and m in vw:
                car += (returns[gvkey][m] - vw[m]) - pre_mean
    return car


# ─── Main ──────────────────────────────────────────────────────────────

def main():
    print('Loading data...')
    vw = load_market_index_monthly()
    rets = load_monthly_returns()
    W_geo = load_weight_matrix(W_GEO_PATH)
    W_fuel = load_weight_matrix(W_FUEL_PATH)
    W_reg = load_weight_matrix(W_REG_PATH)
    firm_sic = load_firm_sic()
    events = load_events()
    inst = load_inst_ownership()
    print(f'  Events: {len(events)}, Inst-coverage gvkeys: {len(inst)}')

    universe_gvkeys = list(set(rets.keys()) | set(firm_sic.keys()))

    # Build the panel mirroring joint_tests.py
    print('Building panel...')
    panel = []
    for ev in events:
        t0 = ev['announcement_date']
        event_gvkeys = set(ev['gvkeys'])
        fm_sic4 = next((firm_sic.get(gk) for gk in event_gvkeys if firm_sic.get(gk)), None)
        ann_year = t0.year
        ann_month = f'{ann_year}-{t0.month:02d}'

        for fm_gk in event_gvkeys:
            if fm_gk not in W_geo: continue
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
                car = compute_monthly_car(gk, ann_month, rets, vw)
                if car is None: continue
                inst_rec = get_inst_at(gk, t0, inst)
                if inst_rec is None: continue
                _, hhi, nm, ipc, t5 = inst_rec
                panel.append({
                    'event_id': ev['event_id'],
                    'gvkey': gk,
                    'car': car,
                    'w_geo': neighbors.get(gk, 0.0),
                    'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                    'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                    'same_sector': 1.0 if (fm_sic4 and firm_sic.get(gk) == fm_sic4) else 0.0,
                    'hhi': hhi, 'n_mgr': nm, 'inst_pct': ipc, 'top5': t5,
                })
    print(f'  Panel rows with inst-coverage: {len(panel):,}')
    print(f'  Unique gvkeys: {len(set(r["gvkey"] for r in panel))}')
    print(f'  Unique events: {len(set(r["event_id"] for r in panel))}')

    # Within-event tercile assignment for HHI
    print('\nAssigning HHI terciles within event...')
    by_event = defaultdict(list)
    for r in panel:
        by_event[r['event_id']].append(r)

    panel_with_tercile = []
    for eid, rows in by_event.items():
        # Sort by HHI ascending within event
        rows_sorted = sorted(rows, key=lambda r: r['hhi'])
        n = len(rows_sorted)
        if n < 6:  # need at least 6 firms for terciles to be informative
            continue
        for i, r in enumerate(rows_sorted):
            if i < n / 3:
                r['hhi_tercile'] = 'T1_dispersed'
            elif i < 2 * n / 3:
                r['hhi_tercile'] = 'T2_mid'
            else:
                r['hhi_tercile'] = 'T3_concentrated'
            panel_with_tercile.append(r)

    print(f'  Rows after tercile assignment: {len(panel_with_tercile):,}')

    # ─── Run cross-sectional FM regression separately on T1 and T3 ──
    spec_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
    names = ['intercept'] + spec_vars

    def per_event_ols(rows, min_n=10):
        by_e = defaultdict(list)
        for r in rows:
            by_e[r['event_id']].append(r)
        per = []
        for eid, evrows in by_e.items():
            if len(evrows) < min_n: continue
            n = len(evrows)
            k = len(spec_vars) + 1
            X = [[1.0] + [r[v] for v in spec_vars] for r in evrows]
            y = [r['car'] for r in evrows]
            XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
            Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
            inv = invert(XtX)
            if inv is None: continue
            beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
            per.append(beta)
        return per

    print('\nRunning FM regressions per HHI tercile...')
    results_by_tercile = {}
    for ter in ('T1_dispersed', 'T2_mid', 'T3_concentrated'):
        rows = [r for r in panel_with_tercile if r['hhi_tercile'] == ter]
        per = per_event_ols(rows, min_n=8)  # lower threshold within tercile
        fm = fm_with_nw(per, lag=4) if per else None
        if fm:
            print(f'  {ter}: T = {fm["T"]} events, '
                  f'gamma_fuel = {fm["mean"][2]:+.4f} (NW SE = {fm["nw_se"][2]:.4f}, '
                  f't = {fm["mean"][2]/fm["nw_se"][2]:.3f})' if fm["nw_se"][2] > 1e-15 else '  no SE')
        results_by_tercile[ter] = fm

    # Pooled OLS within tercile (event-clustered)
    def pooled_ols_evcl(rows):
        if len(rows) < 50: return None
        n = len(rows); k = len(spec_vars) + 1
        X = [[1.0] + [r[v] for v in spec_vars] for r in rows]
        y = [r['car'] for r in rows]
        XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
        inv = invert(XtX)
        if inv is None: return None
        beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
        yhat = [sum(X[i][a]*beta[a] for a in range(k)) for i in range(n)]
        resid = [y[i] - yhat[i] for i in range(n)]
        # Event-clustered SE
        cluster_keys = [r['event_id'] for r in rows]
        cmap = defaultdict(list)
        for i, c in enumerate(cluster_keys): cmap[c].append(i)
        S = [[0.0]*k for _ in range(k)]
        for idxs in cmap.values():
            xu = [0.0]*k
            for i in idxs:
                ri = resid[i]
                for a in range(k): xu[a] += X[i][a]*ri
            for a in range(k):
                for b in range(a, k):
                    v = xu[a]*xu[b]; S[a][b] += v
                    if a != b: S[b][a] += v
        # V = inv * S * inv
        from itertools import product
        V = [[0.0]*k for _ in range(k)]
        for a, b in product(range(k), range(k)):
            for x in range(k):
                for z in range(k):
                    V[a][b] += inv[a][x]*S[x][z]*inv[z][b]
        G = len(cmap)
        if G > 1:
            sc = (G/(G-1)) * ((n-1)/(n-k))
            for a in range(k):
                for b in range(k): V[a][b] *= sc
        se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0 for a in range(k)]
        t = [beta[a]/se[a] if se[a] > 1e-15 else 0 for a in range(k)]
        return {'beta': beta, 'se': se, 't': t, 'n': n, 'G': G}

    pool_by_tercile = {}
    for ter in ('T1_dispersed', 'T2_mid', 'T3_concentrated'):
        rows = [r for r in panel_with_tercile if r['hhi_tercile'] == ter]
        pool_by_tercile[ter] = pooled_ols_evcl(rows)

    # ─── Write output ──
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    lines = [
        '# Institutional-Ownership Split: HHI Terciles, Per-Event',
        '',
        'Per /quant-finance recommendation: split firms WITHIN each event into',
        'terciles by Herfindahl index of 13F manager shares-of-shares at the',
        'most recent quarter on or before announcement_date. T1 = most',
        'dispersed institutional ownership; T3 = most concentrated.',
        '',
        'Hypothesis: if the US null is driven by a retail-flow / smart-money',
        'mechanism, the channel should look DIFFERENT across HHI terciles.',
        '',
        f'Sample restricted to firm-events where 13F coverage exists at',
        f'the event quarter. This is the US-linked sub-sample by construction',
        f'(13F filers report only US holdings). Total N: {len(panel):,}.',
        '',
        '## Headline coefficients by HHI tercile (FM + NW lag 4)',
        '',
        '| HHI Tercile | T (events) | gamma_geo | NW t | gamma_fuel | NW t | gamma_reg | NW t |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for ter, label in [('T1_dispersed', 'T1 (dispersed, low HHI)'),
                       ('T2_mid', 'T2 (middle)'),
                       ('T3_concentrated', 'T3 (concentrated, high HHI)')]:
        fm = results_by_tercile.get(ter)
        if fm is None:
            lines.append(f'| {label} | — | — | — | — | — | — | — |')
            continue
        T = fm['T']; m = fm['mean']; se = fm['nw_se']
        def t_(idx):
            return f'{m[idx]/se[idx]:+.3f}' if se[idx] > 1e-15 else '—'
        lines.append(
            f'| {label} | {T} | {m[1]:+.4f} | {t_(1)} | '
            f'{m[2]:+.4f} | {t_(2)} | {m[3]:+.4f} | {t_(3)} |'
        )

    lines += [
        '',
        '## Pooled OLS by HHI tercile (event-clustered SEs)',
        '',
        '| HHI Tercile | N | gamma_geo (t) | gamma_fuel (t) | gamma_reg (t) |',
        '|---|---:|---:|---:|---:|',
    ]
    for ter, label in [('T1_dispersed', 'T1 (dispersed)'),
                       ('T2_mid', 'T2 (middle)'),
                       ('T3_concentrated', 'T3 (concentrated)')]:
        p = pool_by_tercile.get(ter)
        if p is None:
            lines.append(f'| {label} | — | — | — | — |')
            continue
        b, t = p['beta'], p['t']
        lines.append(
            f'| {label} | {p["n"]:,} | {b[1]:+.4f} ({t[1]:+.3f}) | '
            f'{b[2]:+.4f} ({t[2]:+.3f}) | {b[3]:+.4f} ({t[3]:+.3f}) |'
        )

    lines += [
        '',
        '## Notes',
        '',
        '- HHI is the Gabaix-Koijen (2021) granular-investor standard for',
        '  measuring institutional concentration. Computed as sum_g (s_g/S)^2',
        '  where s_g is shares held by manager g and S is total institutional shares.',
        '- Within-event tercile assignment absorbs time trends in 13F filer count',
        '  (which roughly doubled over 2008-2025).',
        '- Panel is restricted to firm-events where the candidate firm has 13F',
        '  coverage at the event quarter (i.e., is US-listed in the Thomson S34',
        '  universe). This is a NARROWER sub-sample than the full 565-firm panel.',
        '',
    ]

    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\nWrote {OUT_PATH}')


if __name__ == '__main__':
    main()

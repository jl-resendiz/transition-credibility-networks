"""Institutional-ownership split for the NON-US sub-sample, using
Refinitiv free-float concentration in lieu of 13F manager Herfindahl.

PRE-STAGED for Phase 6 of EXECUTION_PLAN.md. Run AFTER:
  1. pull_refinitiv_extra.py  (produces refinitiv_extra.csv)
  2. build_nonus_institutional_panel.py  (produces nonus_concentration.csv)

This is the non-US analogue of `institutional_split.py`, which used the 13F
manager-Herfindahl on the US sub-sample. Method:

For each NON-US event, sort the candidate firms in the panel by their
`concentrated_ownership_pct` (= 100 - free_float_pct from Refinitiv).
Within-event tercile assignment. Run cross-sectional FM regression on T1
(most dispersed = high free-float) vs T3 (most concentrated = low free-float).

Output: results/metrics/institutional_split_nonus.md

Decision logic (Phase 6.4):
  - If T3 (concentrated) shows MORE NEGATIVE gamma_fuel than T1 (dispersed),
    SAME monotonic pattern as the US HHI-tercile split → GLOBAL MECHANISM
    finding. Update manuscript per Scenario A in PHASE4_REVIEW.md / EXECUTION_PLAN.md.
  - If both terciles show similar negative gamma_fuel → channel uniform across
    ownership concentration in non-US. Scenario B in PHASE4_REVIEW.md.
  - If pattern reverses (T1 more negative than T3) → unexpected, spawn
    /quant-finance agent for interpretation.
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
from _paths import derived_path, raw_path, results_path

CONC_PATH = os.path.join(derived_path('institutional'), 'nonus_concentration.csv')
DAILY_RET_PATH = os.path.join(derived_path('returns'), 'monthly_returns.csv')
FF3_PATH = os.path.join(raw_path('factors'), 'F-F_Research_Data_Factors.csv')
EVENTS_PATH = os.path.join(derived_path('events'), 'coal_retirement_events.csv')
W_GEO_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_geo.csv')
W_FUEL_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_fuel.csv')
W_REG_PATH = os.path.join(derived_path('networks'), 'weight_matrix_W_regulatory.csv')
FUND_PATH = os.path.join(derived_path('fundamentals'), 'firm_fundamentals.csv')

OUT_PATH = os.path.join(results_path('metrics'), 'institutional_split_nonus.md')

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
    if not os.path.exists(CONC_PATH):
        print(f'SKIP: {CONC_PATH} not found.')
        print('  Run build_nonus_institutional_panel.py first (requires Refinitiv pull).')
        print('  This is the optional Phase 6 non-US institutional extension.')
        sys.exit(0)

    # ─── Load ownership concentration data (non-US only) ──
    print('Loading non-US ownership concentration...')
    concentration = {}
    with open(CONC_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if row.get('us_listed') == '1':
                continue  # Skip US-listed firms (covered by institutional_split.py)
            try:
                gk = str(row['gvkey']).split('.')[0].zfill(6)
                conc = float(row['concentrated_ownership_pct'])
            except (ValueError, KeyError):
                continue
            concentration[gk] = conc
    print(f'  Non-US firms with concentration data: {len(concentration)}')

    # ─── Load returns + factors + matrices + events (mirror of institutional_split.py) ──
    monthly_ret = defaultdict(dict)
    with open(DAILY_RET_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = str(row['gvkey']).split('.')[0].zfill(6)
            ym = row['datadate'][:7]
            try:
                monthly_ret[gk][ym] = float(row['ret_monthly'])
            except ValueError:
                pass

    vw = {}
    with open(FF3_PATH, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit(): continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5: continue
            d = parts[0]
            if not d.isdigit() or len(d) != 6: continue
            try:
                mktrf = float(parts[1])/100; rf = float(parts[4])/100
            except ValueError:
                continue
            vw[f'{d[:4]}-{d[4:6]}'] = mktrf + rf

    W_geo = defaultdict(dict)
    with open(W_GEO_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
    W_fuel = defaultdict(dict)
    with open(W_FUEL_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
    W_reg = defaultdict(dict)
    if os.path.exists(W_REG_PATH):
        with open(W_REG_PATH, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                wval = row.get('w_ij') or row.get('w_reg')
                try:
                    W_reg[row['gvkey_i']][row['gvkey_j']] = float(wval)
                except (ValueError, TypeError):
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
            ret = row.get('event_date', '').strip()
            ed = ann if ann else ret
            if not ed or len(ed) < 7: continue
            events.append({
                'event_id': i,
                'event_month': ed[:7],
                'gvkeys': [g.strip().zfill(6) for g in row['matched_gvkeys'].split(';')],
            })
    print(f'  {len(events)} events')

    universe = list(set(monthly_ret.keys()) | set(firm_sic.keys()))


    def compute_monthly_car(gvkey, event_month):
        if gvkey not in monthly_ret: return None
        months = sorted(monthly_ret[gvkey].keys())
        idx = next((i for i, m in enumerate(months) if m >= event_month), None)
        if idx is None: return None
        pre = max(0, idx - PRE_MONTHS)
        pre_ar = []
        for m in months[pre:idx]:
            if m in monthly_ret[gvkey] and m in vw:
                pre_ar.append(monthly_ret[gvkey][m] - vw[m])
        if len(pre_ar) < 12: return None
        pre_mean = sum(pre_ar)/len(pre_ar)
        car = 0.0; cnt = 0
        for off in range(-1, POST_MONTHS+1):
            i2 = idx + off
            if 0 <= i2 < len(months):
                m = months[i2]
                if m in monthly_ret[gvkey] and m in vw:
                    car += (monthly_ret[gvkey][m] - vw[m]) - pre_mean; cnt += 1
        return car if cnt >= 3 else None


    # ─── Build panel restricted to non-US firms with concentration data ──
    print('\nBuilding panel...')
    panel = []
    for ev in events:
        event_gvkeys = set(ev['gvkeys'])
        em = ev['event_month']
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
            ctrl = (random.sample(non_connected, n_ctrl)
                    if len(non_connected) > n_ctrl else non_connected)
            for gk in list(neighbor_gks) + ctrl:
                if gk not in concentration:
                    continue  # Skip firms without ownership data (i.e., US-listed or unmatched)
                car = compute_monthly_car(gk, em)
                if car is None: continue
                panel.append({
                    'event_id': ev['event_id'], 'gvkey': gk, 'car': car,
                    'w_geo': neighbors.get(gk, 0.0),
                    'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                    'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                    'same_sector': 1.0 if (fm_sic4 and firm_sic.get(gk) == fm_sic4) else 0.0,
                    'concentration': concentration[gk],
                })

    print(f'  Panel size (non-US, concentration-eligible): {len(panel):,}')
    print(f'  Unique gvkeys: {len(set(r["gvkey"] for r in panel))}')
    print(f'  Unique events: {len(set(r["event_id"] for r in panel))}')

    if len(panel) < 100:
        sys.exit('Panel too small for meaningful split. Check coverage.')

    # ─── Within-event tercile assignment ──
    by_event = defaultdict(list)
    for r in panel:
        by_event[r['event_id']].append(r)

    panel_with_tercile = []
    for eid, rows in by_event.items():
        rows_sorted = sorted(rows, key=lambda r: r['concentration'])
        n = len(rows_sorted)
        if n < 6: continue
        for i, r in enumerate(rows_sorted):
            if i < n / 3:
                r['tercile'] = 'T1_dispersed'
            elif i < 2 * n / 3:
                r['tercile'] = 'T2_mid'
            else:
                r['tercile'] = 'T3_concentrated'
            panel_with_tercile.append(r)

    print(f'  Rows with tercile assignment: {len(panel_with_tercile):,}')


    # ─── Run FM per tercile ──
    SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

    def per_event_ols(rows, min_n=8):
        by_e = defaultdict(list)
        for r in rows:
            by_e[r['event_id']].append(r)
        per = []
        for eid, evrows in by_e.items():
            if len(evrows) < min_n: continue
            n = len(evrows); k = len(SPEC_VARS) + 1
            X = [[1.0] + [r[v] for v in SPEC_VARS] for r in evrows]
            y = [r['car'] for r in evrows]
            XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
            Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
            inv = invert(XtX)
            if inv is None: continue
            beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
            per.append(beta)
        return per

    print('\nRunning FM per tercile...')
    results = {}
    for ter in ('T1_dispersed', 'T2_mid', 'T3_concentrated'):
        rows = [r for r in panel_with_tercile if r['tercile'] == ter]
        per = per_event_ols(rows, min_n=8)
        fm = fm_with_nw(per, lag=4) if per else None
        if fm:
            t = fm['mean'][2] / fm['nw_se'][2] if fm['nw_se'][2] > 1e-15 else 0
            print(f'  {ter}: T = {fm["T"]} events, '
                  f'gamma_fuel = {fm["mean"][2]:+.4f} (NW SE = {fm["nw_se"][2]:.4f}, t = {t:+.3f})')
        results[ter] = fm


    # ─── Output ──
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    lines = [
        '# Institutional-Ownership Split: NON-US Sub-Sample (Refinitiv Free-Float Terciles)',
        '',
        'Companion to `institutional_split.md` (US sub-sample with 13F HHI).',
        'Method: within-event tercile assignment by `concentrated_ownership_pct`',
        '= 100 - free_float_pct from Refinitiv. Higher concentration = more',
        'ownership held in non-public hands (insiders, sovereign, strategic blocks),',
        'a proxy for institutional concentration where 13F filings are unavailable.',
        '',
        f'Sample: {len(panel):,} firm-events, {len(set(r["gvkey"] for r in panel))} non-US firms.',
        '',
        '## Headline coefficients by concentration tercile (FM + NW lag 4)',
        '',
        '| Tercile | T (events) | gamma_geo | NW t | gamma_fuel | NW t | gamma_reg | NW t |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for ter, label in [('T1_dispersed', 'T1 (dispersed, high free-float)'),
                       ('T2_mid', 'T2 (middle)'),
                       ('T3_concentrated', 'T3 (concentrated, low free-float)')]:
        fm = results.get(ter)
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
        '## Cross-sample comparison',
        '',
        'Compare against US 13F-HHI tercile split in `institutional_split.md`:',
        '- US T3 (concentrated): gamma_fuel = -6.08 (t = -3.27)',
        '- US T1 (dispersed): gamma_fuel = +3.23 (t = +4.49)',
        '',
        'A monotonic pattern (T3 < T2 < T1) in the non-US sub-sample replicates the',
        'US-sample finding, supporting a global "smart-money pricing" mechanism.',
        '',
        '## Notes',
        '',
        '- Refinitiv free-float is firm-level cross-sectional (not quarterly like 13F).',
        '  We use the most recent snapshot per firm; this introduces a static rather',
        '  than dynamic concentration measure. A robustness extension could use multiple',
        '  Refinitiv snapshots over time.',
        '- The non-US sample excludes firms covered by the US 13F panel to avoid',
        '  double-counting and to keep the two splits methodologically distinct.',
    ]

    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\nWrote {OUT_PATH}')


if __name__ == '__main__':
    main()

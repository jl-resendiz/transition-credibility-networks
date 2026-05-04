"""Pre-trends randomization placebo: shuffle event dates within firm,
recompute gamma_fuel, verify the headline coefficient is in the tail of
the placebo distribution.

Mirrors the panel-construction logic of joint_tests.py / two_way_clustering.py
exactly, but with the event date randomized for each iteration. If the
coefficient distribution under randomized event dates is centred at zero
and the observed coefficient is far in the negative tail, the channel
is not a spurious pre-trend artifact.

Method: For each placebo iteration b:
1. Randomize each event's announcement date by adding a random offset
   in [-36, +36] months within the same firm's return history.
2. Re-run the cross-sectional FM regression.
3. Record gamma_fuel.
After 999 iterations, compare observed gamma_fuel = -4.77 to the placebo
distribution.

Output: results/metrics/pretrends_placebo.md
"""
import csv
import hashlib
import math
import os
import random
import sys
from collections import defaultdict
from datetime import datetime, timedelta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import derived_path, raw_path, results_path

POST_MONTHS = 3
PRE_MONTHS = 24
N_PLACEBO = 999
SEED = 42
SHUFFLE_RANGE_MONTHS = 36  # ±3 years


def _print(msg=''):
    print(msg); sys.stdout.flush()


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


# ─── Load returns + factors + matrices + events (mirrors joint_tests.py) ──

_print('Loading data...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = str(row['gvkey']).split('.')[0].zfill(6)
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass

vw = {}
with open(raw_path('factors', 'F-F_Research_Data_Factors.csv'), 'r', encoding='utf-8') as f:
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
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
W_fuel = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_fuel.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
W_reg = defaultdict(dict)
reg_path = derived_path('networks', 'weight_matrix_W_regulatory.csv')
if os.path.exists(reg_path):
    with open(reg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            wval = row.get('w_ij') or row.get('w_reg')
            try:
                W_reg[row['gvkey_i']][row['gvkey_j']] = float(wval)
            except (ValueError, TypeError):
                continue

firm_sic = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']; s = row.get('sic')
        if s and gk not in firm_sic: firm_sic[gk] = s[:4]

events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
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
_print(f'  {len(events)} first-mover events')

universe = list(set(monthly_ret.keys()) | set(firm_sic.keys()))


# ─── CAR helper (vwretd-adjusted, pre-event demeaning) ──

def compute_monthly_car(gvkey, event_month, post=POST_MONTHS):
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
    for off in range(-1, post+1):
        i2 = idx + off
        if 0 <= i2 < len(months):
            m = months[i2]
            if m in monthly_ret[gvkey] and m in vw:
                car += (monthly_ret[gvkey][m] - vw[m]) - pre_mean; cnt += 1
    return car if cnt >= 3 else None


# ─── Build per-event candidate panels ONCE (panel structure is fixed) ──

def build_candidate_panels():
    """Returns dict: event_id → list of candidate firm dicts (without CAR)."""
    panels = {}
    for ev in events:
        event_gvkeys = set(ev['gvkeys'])
        fm_sic4 = next((firm_sic.get(gk) for gk in event_gvkeys if firm_sic.get(gk)), None)
        cands = []
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
                cands.append({
                    'gvkey': gk,
                    'w_geo': neighbors.get(gk, 0.0),
                    'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                    'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                    'same_sector': 1.0 if (fm_sic4 and firm_sic.get(gk) == fm_sic4) else 0.0,
                })
        panels[ev['event_id']] = cands
    return panels


_print('Building candidate panels (cached)...')
panels = build_candidate_panels()
_print(f'  Panels for {len(panels)} events')


# ─── Compute one FM gamma_fuel given a (possibly shuffled) event_month per event ──

SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

def fm_gamma_fuel(event_months_by_id, min_firms=20):
    """For each event, run cross-section regression at the given event month;
    return the FM mean gamma_fuel."""
    per_event_betas = []
    for eid, ev_month in event_months_by_id.items():
        cands = panels.get(eid, [])
        if len(cands) < min_firms: continue
        # Compute CARs for each candidate at this (possibly placebo) event month
        rows = []
        for c in cands:
            car = compute_monthly_car(c['gvkey'], ev_month)
            if car is not None:
                rows.append({**c, 'car': car})
        if len(rows) < min_firms: continue

        n = len(rows); k = len(SPEC_VARS) + 1
        X = [[1.0] + [r[v] for v in SPEC_VARS] for r in rows]
        y = [r['car'] for r in rows]
        XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
        inv = invert(XtX)
        if inv is None: continue
        beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
        per_event_betas.append(beta[2])  # gamma_fuel index = 2 (after intercept, w_geo, w_fuel)

    if len(per_event_betas) < 5: return None
    return sum(per_event_betas) / len(per_event_betas), len(per_event_betas)


# ─── Observed coefficient (no shuffle) ──

_print('\nComputing observed gamma_fuel (no shuffle)...')
true_event_months = {ev['event_id']: ev['event_month'] for ev in events}
observed, T_obs = fm_gamma_fuel(true_event_months)
_print(f'  Observed gamma_fuel = {observed:+.4f} (T = {T_obs} events)')


# ─── Placebo: shuffle event months within feasible range, repeat ──

_print(f'\nRunning {N_PLACEBO} placebo iterations (event-date randomization)...')
random.seed(SEED)

def shuffle_event_month(orig_ym, rng):
    """Add a random ±SHUFFLE_RANGE_MONTHS offset to the original month."""
    y, m = int(orig_ym[:4]), int(orig_ym[5:7])
    total = y * 12 + (m - 1) + rng.randint(-SHUFFLE_RANGE_MONTHS, SHUFFLE_RANGE_MONTHS)
    new_y = total // 12
    new_m = (total % 12) + 1
    if new_m == 13:
        new_y += 1; new_m = 1
    return f'{new_y:04d}-{new_m:02d}'

placebo_results = []
for b in range(N_PLACEBO):
    rng = random.Random(SEED * 10 + b)
    placebo_event_months = {
        eid: shuffle_event_month(orig_em, rng)
        for eid, orig_em in true_event_months.items()
    }
    res = fm_gamma_fuel(placebo_event_months)
    if res:
        placebo_results.append(res[0])
    if (b+1) % 100 == 0:
        _print(f'  iteration {b+1}/{N_PLACEBO}: '
               f'last gamma= {res[0]:+.4f}' if res else f'  iteration {b+1} skipped')

_print(f'\nValid placebo iterations: {len(placebo_results)}/{N_PLACEBO}')

placebo_results.sort()
n = len(placebo_results)
p1 = placebo_results[max(0, n//100)]
p5 = placebo_results[max(0, n//20)]
p50 = placebo_results[n//2]
p95 = placebo_results[min(n-1, 19*n//20)]
p99 = placebo_results[min(n-1, 99*n//100)]
mn = placebo_results[0]
mx = placebo_results[-1]

# Two-sided RI p-value: P(|placebo| >= |observed|)
n_extreme = sum(1 for x in placebo_results if abs(x) >= abs(observed))
p_two = (n_extreme + 1) / (n + 1)


# ─── Output ──

OUT_PATH = os.path.join(results_path('metrics'), 'pretrends_placebo.md')
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
lines = [
    '# Pre-Trends Randomization Placebo',
    '',
    f'Iterations: {N_PLACEBO} (valid: {n})',
    f'Event-date shuffle range: ±{SHUFFLE_RANGE_MONTHS} months around the true announcement.',
    'For each iteration, every event\'s announcement date is randomly shifted',
    'within the firm\'s return history; the full FM cross-sectional regression is',
    're-estimated; gamma_fuel is recorded.',
    '',
    '## Placebo distribution of gamma_fuel',
    '',
    f'- Observed gamma_fuel (true event dates): **{observed:+.4f}**',
    f'- Placebo mean: {sum(placebo_results)/n:+.4f}',
    f'- Placebo median (p50): {p50:+.4f}',
    f'- Placebo 1st pct (p1): {p1:+.4f}',
    f'- Placebo 5th pct (p5): {p5:+.4f}',
    f'- Placebo 95th pct (p95): {p95:+.4f}',
    f'- Placebo 99th pct (p99): {p99:+.4f}',
    f'- Range: [{mn:+.4f}, {mx:+.4f}]',
    '',
    f'**Two-sided RI p-value: {p_two:.4f}** ({n_extreme} of {n} placebo iterations more extreme than observed |gamma_fuel|).',
    '',
    '## Interpretation',
    '',
    'Under the sharp null that announcement timing carries no cross-sectional',
    'information, randomly shifted event dates should produce a gamma_fuel distribution',
    'centred at zero. The observed gamma_fuel from the true event dates should be far',
    'in the tail of this placebo distribution.',
    '',
    'A p-value below 0.05 confirms that the observed coefficient is unlikely to',
    'have arisen from a generic pre-trend or sample-period drift in coal-similar',
    'firms\' returns: the channel responds specifically to the timing of the',
    'true retirement announcements, not to any month-of-year structure.',
    '',
]

with open(OUT_PATH, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
_print(f'\nWrote {OUT_PATH}')

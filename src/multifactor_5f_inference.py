"""5-factor (FF3 + UMD + Utility) extension of multifactor_inference.py.

Adds Carhart momentum (UMD) to the existing 4-factor model:
  (r_it - rf_t) = alpha_i + b_M (Mkt-RF)_t + b_S SMB_t + b_V HML_t
                  + b_W UMD_t + b_U (UTL_excess)_t + epsilon_it

Direct response to /quant-finance referee anticipation: a finance referee
will ask whether the channel survives momentum control, since coal-heavy
firms tend to have negative momentum (low past returns → high beta_W).
If gamma_fuel survives 5-factor adjustment, the channel is not absorbed
by momentum either.

Inputs:  data/derived/returns/monthly_returns.csv
         data/raw/factors/F-F_Research_Data_Factors.csv
         data/raw/factors/F-F_Momentum_Factor.csv  (NEW)
         data/derived/networks/weight_matrix_W_*.csv
         data/derived/events/coal_retirement_events.csv

Output:  results/metrics/multifactor_5f_inference.md
"""
import csv
import os
import sys
import math
import random
import hashlib
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import derived_path, raw_path, results_path


def _print(msg=''):
    print(msg); sys.stdout.flush()


# ─── OLS helpers (reused) ─────────────────────────────────────────────

def invert_matrix(mat):
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


def ols_simple(y, X):
    n = len(y)
    if n < len(X[0]) + 1: return None
    k = len(X[0])
    XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
    inv = invert_matrix(XtX)
    if inv is None: return None
    beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
    return {'beta': beta}


def newey_west_se(series, lag=4):
    n = len(series)
    if n < 2: return float('nan')
    mean = sum(series) / n
    dev = [x - mean for x in series]
    g0 = sum(d*d for d in dev) / n
    var_nw = g0
    for L in range(1, min(lag, n-1) + 1):
        w = 1 - L/(lag+1)
        c = sum(dev[t]*dev[t-L] for t in range(L, n)) / n
        var_nw += 2*w*c
    if var_nw <= 0: return float('nan')
    return math.sqrt(var_nw / n)


# ─── Load data ────────────────────────────────────────────────────────

_print('Loading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = str(row['gvkey']).split('.')[0].zfill(6)
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass


def load_ff3_monthly(path):
    f3 = {}
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The '): continue
            if line.startswith(','): continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5: continue
            d = parts[0]
            if not d.isdigit() or len(d) != 6: continue
            try:
                mktrf = float(parts[1]) / 100.0
                smb = float(parts[2]) / 100.0
                hml = float(parts[3]) / 100.0
                rf = float(parts[4]) / 100.0
            except ValueError:
                continue
            f3[f'{d[:4]}-{d[4:6]}'] = {
                'mkt_rf': mktrf, 'smb': smb, 'hml': hml, 'rf': rf,
            }
    return f3


def load_umd_monthly(path):
    """Ken French UMD monthly file: yyyymm, mom (in percent)."""
    umd = {}
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 2: continue
            d = parts[0]
            if not d.isdigit() or len(d) != 6: continue
            try:
                mom = float(parts[1]) / 100.0
            except ValueError:
                continue
            umd[f'{d[:4]}-{d[4:6]}'] = mom
    return umd


_print('Loading FF3 monthly factors...')
ff3 = load_ff3_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
_print(f'  FF3 months: {len(ff3)}')

_print('Loading UMD monthly factor...')
umd = load_umd_monthly(raw_path('factors', 'F-F_Momentum_Factor.csv'))
_print(f'  UMD months: {len(umd)}')

# Construct utility industry factor (same as multifactor_inference.py)
month_returns = defaultdict(list)
for gk, dates in monthly_ret.items():
    for m, r in dates.items():
        month_returns[m].append(r)
util_factor = {m: sum(rs)/len(rs) for m, rs in month_returns.items() if len(rs) >= 30}

# Combine all into factors_panel: 5 factors per month
factors_panel = {}
for m, ff in ff3.items():
    if m in util_factor and m in umd:
        factors_panel[m] = {
            **ff,
            'umd': umd[m],
            'utl_excess': util_factor[m] - ff['rf'],
        }
_print(f'  Months with all 5 factors: {len(factors_panel)}')


# ─── Networks, fundamentals, events ──────────────────────────────────

W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
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

fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']; fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    return f['sic'][:4] if (f and f.get('sic')) else None


_print('Loading events...')
all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'): continue
        if row.get('is_first_mover') != 'True': continue
        ann = row.get('announcement_date', '').strip()
        ret = row.get('event_date', '').strip()
        ed = ann if ann else ret
        if not ed or len(ed) < 7: continue
        all_events.append({
            'plant': row['plant_name'], 'event_month': ed[:7],
            'gvkeys': row['matched_gvkeys'].split(';'),
        })
_print(f'  First-mover events: {len(all_events)}')


# ─── 5-factor abnormal-return CAR ─────────────────────────────────────

PRE_MONTHS = 24
PRE_MIN = 12
POST_MONTHS = 3


def compute_car_5f(gvkey, event_month):
    if gvkey not in monthly_ret: return None
    months = sorted(monthly_ret[gvkey].keys())
    idx = next((i for i, m in enumerate(months) if m >= event_month), None)
    if idx is None: return None

    y_pre, X_pre = [], []
    for i in range(max(0, idx - PRE_MONTHS), idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in factors_panel:
            f = factors_panel[m]
            y_pre.append(monthly_ret[gvkey][m] - f['rf'])
            X_pre.append([1.0, f['mkt_rf'], f['smb'], f['hml'], f['umd'], f['utl_excess']])
    if len(y_pre) < PRE_MIN:
        return None
    res = ols_simple(y_pre, X_pre)
    if res is None: return None
    a, bM, bS, bH, bW, bU = res['beta']

    car = 0.0; n = 0
    for off in range(-1, POST_MONTHS + 1):
        i2 = idx + off
        if 0 <= i2 < len(months):
            m = months[i2]
            if m in monthly_ret[gvkey] and m in factors_panel:
                f = factors_panel[m]
                expected = (a + bM*f['mkt_rf'] + bS*f['smb'] + bH*f['hml']
                            + bW*f['umd'] + bU*f['utl_excess'])
                actual = monthly_ret[gvkey][m] - f['rf']
                car += actual - expected; n += 1
    if n < 3: return None
    return car


# ─── Build event-firm panel and run FM ────────────────────────────────

_print('\nBuilding 5-factor event-firm panel...')
MIN_OBS_PER_EVENT = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

event_datasets = {}
for event_id, ev in enumerate(all_events):
    event_gvkeys = set(ev['gvkeys'])
    em = ev['event_month']
    fm_sic4 = next((get_sic4(gk) for gk in event_gvkeys if get_sic4(gk)), None)

    obs = []
    for fm_gk in event_gvkeys:
        if fm_gk not in W_geo: continue
        neighbors = W_geo[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals
                         if gk not in event_gvkeys and gk not in neighbors]
        seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
        random.seed(seed)
        n_ctrl = min(len(non_connected), max(5*len(neighbor_gks), 20))
        ctrl = (random.sample(non_connected, n_ctrl)
                if len(non_connected) > n_ctrl else non_connected)
        for gk in list(neighbor_gks) + ctrl:
            j_sic4 = get_sic4(gk)
            same_sec = 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0
            car = compute_car_5f(gk, em)
            if car is None: continue
            obs.append({
                'car': car,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': same_sec, 'gvkey': gk,
            })
    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'  Valid events: {n_valid}, total obs: {total_obs}')


# ─── FM cross-sectional regressions ──────────────────────────────────

_print('\nRunning FM regressions...')

event_betas = defaultdict(list)
event_r2s, event_ns = [], []
for eid in sorted(event_datasets.keys()):
    obs = event_datasets[eid]
    n = len(obs)
    k = len(SPEC_VARS) + 1
    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in SPEC_VARS] for o in obs]
    XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
    inv = invert_matrix(XtX)
    if inv is None: continue
    beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
    yhat = [sum(X[i][a]*beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - yhat[i] for i in range(n)]
    ss_tot = sum((yi - sum(y)/n)**2 for yi in y)
    ss_res = sum(r**2 for r in resid)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    event_r2s.append(r2)
    event_ns.append(n)
    names = ['intercept'] + SPEC_VARS
    for i, nm in enumerate(names):
        event_betas[nm].append(beta[i])

T = len(event_r2s)
_print(f'  Events with FM beta: {T}')

names = ['intercept'] + SPEC_VARS
fm_results = {}
for nm in names:
    series = event_betas[nm]
    mean = sum(series)/len(series)
    se = newey_west_se(series, lag=4)
    t = mean/se if se > 1e-15 else 0
    p = 2*(1 - 0.5*(1 + math.erf(abs(t)/math.sqrt(2))))
    fm_results[nm] = {'mean': mean, 'se': se, 't': t, 'p': p}
    _print(f'  {nm:12s}: mean = {mean:+.4f}, SE_NW = {se:.4f}, t = {t:+.3f}, p = {p:.4f}')

# Difference test (geo - fuel)
diff_series = [g - f for g, f in zip(event_betas['w_geo'], event_betas['w_fuel'])]
diff_mean = sum(diff_series) / len(diff_series)
diff_se = newey_west_se(diff_series, lag=4)
diff_t = diff_mean / diff_se if diff_se > 1e-15 else 0
_print(f'\nDifference (geo - fuel): {diff_mean:+.4f}, SE_NW = {diff_se:.4f}, t = {diff_t:+.3f}')

# ─── Write output ──
OUT_PATH = os.path.join(results_path('metrics'), 'multifactor_5f_inference.md')
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

lines = [
    '# 5-Factor Inference: FF3 + UMD + Utility Industry',
    '',
    'Adds Carhart momentum (UMD) to the existing 4-factor specification:',
    '',
    '  AR_it = r_it − rf_t − [α_i + β_M(Mkt−RF)_t + β_S SMB_t + β_H HML_t',
    '                          + β_W UMD_t + β_U(UTL−rf)_t]',
    '',
    'Firm-by-firm betas estimated on a 24-month pre-event window. CAR is',
    'the within-window prediction error over [-1, +3] months.',
    '',
    f'Events with valid FM regression: {T}',
    f'Total firm-event observations: {total_obs}',
    f'Avg firms per event: {sum(event_ns)/T:.1f}',
    f'Avg within-event R²: {sum(event_r2s)/T:.4f}',
    '',
    '## Headline coefficients (FM + NW lag 4)',
    '',
    '| Variable | Mean β | NW SE | t | p |',
    '|---|---:|---:|---:|---:|',
]
for nm in names:
    r = fm_results[nm]
    stars = '***' if r['p']<0.01 else '**' if r['p']<0.05 else '*' if r['p']<0.10 else ''
    lines.append(f'| {nm} | {r["mean"]:+.4f} | {r["se"]:.4f} | {r["t"]:+.3f} | {r["p"]:.4f}{stars} |')

lines += [
    '',
    f'**Difference (γ_geo − γ_fuel):** {diff_mean:+.4f}, NW SE = {diff_se:.4f}, t = {diff_t:+.3f}',
    '',
    '## Interpretation',
    '',
    'A negative and statistically significant γ_fuel under 5-factor',
    'adjustment indicates that the channel is NOT absorbed by:',
    '- Market (Mkt-RF)',
    '- Size (SMB)',
    '- Value (HML)',
    '- Momentum (UMD) — the new factor in this 5F spec, addressing the',
    '  referee concern that coal-heavy peers carry persistent negative',
    '  momentum and the channel may be a "low-momentum trap".',
    '- Utility-industry portfolio (UTL−rf)',
    '',
    'Comparison with the existing 4-factor result (in multifactor_inference.md):',
    '4-factor (FF3 + Utility): γ_fuel = -3.10, t = -4.50.',
    'If 5-factor preserves significance, the channel adds explanatory power',
    'beyond standard risk factors AND momentum.',
    '',
]

with open(OUT_PATH, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
_print(f'\nWrote {OUT_PATH}')

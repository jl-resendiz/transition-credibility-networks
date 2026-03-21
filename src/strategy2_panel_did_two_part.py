"""Panel DiD two-part exposure model.

Estimate: AR_{j,t} = alpha_j + gamma_t + b1*(D_ij * Post) + b2*(log1p(w_ij) * Post) + eps
where D_ij = 1{w_ij > 0}.

Outputs metrics in results/metrics.
"""
import csv
import math
import os
import random
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path

def _env_int(name, default):
    v = os.getenv(name)
    if v is None or v == '':
        return default
    try:
        return int(v)
    except ValueError:
        return default


def _env_bool(name, default):
    v = os.getenv(name)
    if v is None or v == '':
        return default
    return v.strip().lower() in ('1', 'true', 'yes', 'y')


EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'all_matched')
EVENTS_PATH = os.getenv('EVENTS_PATH', '')
EXACT_ONLY = _env_bool('EXACT_ONLY', False)
TAU_START = _env_int('TAU_START', -6)
TAU_END = _env_int('TAU_END', 12)
POST_START = _env_int('POST_START', 0)
POST_END = _env_int('POST_END', 12)
CONTROL_MULT = _env_int('CONTROL_MULT', 5)
OVERLAP_RULE = os.getenv('OVERLAP_RULE', 'nearest')
FOREIGN_ONLY = _env_bool('FOREIGN_ONLY', False)
WRITE_METRICS = _env_bool('WRITE_METRICS', True)
METRICS_SUFFIX = os.getenv('METRICS_SUFFIX', '').strip()

TIER_FILTER = os.getenv('TIER_FILTER', '').strip()
BINDING_ONLY = _env_bool('BINDING_ONLY', False)
EXCLUDE_YEARS = os.getenv('EXCLUDE_YEARS', '').strip()


def add_months(ym, delta):
    y, m = ym.split('-')
    y = int(y)
    m = int(m)
    m += delta
    while m > 12:
        y += 1
        m -= 12
    while m < 1:
        y -= 1
        m += 12
    return f"{y:04d}-{m:02d}"


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The ') or line.startswith(','):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5:
                continue
            date = parts[0]
            if not date.isdigit() or len(date) != 6:
                continue
            try:
                mktrf_val = float(parts[1])
                rf_val = float(parts[4])
            except ValueError:
                continue
            vw = (mktrf_val + rf_val) / 100.0
            vwretd[f'{date[:4]}-{date[4:6]}'] = vw
    return vwretd


def invert_matrix(mat):
    n = len(mat)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(mat)]
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[max_row][col]) < 1e-20:
            return None
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pivot
        for row in range(n):
            if row != col:
                factor = aug[row][col]
                for j in range(2 * n):
                    aug[row][j] -= factor * aug[col][j]
    return [row[n:] for row in aug]


def mat_mul(a, b):
    rows = len(a)
    cols = len(b[0])
    mid = len(b)
    out = [[0.0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = a[i][k]
            if aik == 0:
                continue
            for j in range(cols):
                out[i][j] += aik * b[k][j]
    return out


def _cluster_cov(X, resid, clusters):
    k = len(X[0])
    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for _, idxs in clusters.items():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]
    return S, len(clusters)


def ols_cluster(data, y_var, x_vars, cluster_var):
    n = len(data)
    k = len(x_vars) + 1
    X = [[1.0] + [row[v] for v in x_vars] for row in data]
    y = [row[y_var] for row in data]
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    resid = [y[i] - sum(beta[a] * X[i][a] for a in range(k)) for i in range(n)]

    clusters = {}
    for i, d in enumerate(data):
        cid = d.get(cluster_var)
        clusters.setdefault(cid, []).append(i)
    S, G = _cluster_cov(X, resid, clusters)
    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
    se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]
    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'n': n,
    }


def two_way_demean(data, y_var, x_vars, firm_key='gvkey', time_key='ym'):
    vars_all = [y_var] + x_vars
    overall = {v: 0.0 for v in vars_all}
    firm_sums = defaultdict(lambda: defaultdict(float))
    firm_counts = defaultdict(int)
    time_sums = defaultdict(lambda: defaultdict(float))
    time_counts = defaultdict(int)

    for row in data:
        f = row[firm_key]
        t = row[time_key]
        firm_counts[f] += 1
        time_counts[t] += 1
        for v in vars_all:
            val = row[v]
            overall[v] += val
            firm_sums[f][v] += val
            time_sums[t][v] += val

    n = len(data)
    if n == 0:
        return []
    for v in vars_all:
        overall[v] /= n

    firm_means = {f: {v: firm_sums[f][v] / firm_counts[f] for v in vars_all}
                  for f in firm_counts}
    time_means = {t: {v: time_sums[t][v] / time_counts[t] for v in vars_all}
                  for t in time_counts}

    out = []
    for row in data:
        f = row[firm_key]
        t = row[time_key]
        d = {firm_key: f, time_key: t}
        for v in vars_all:
            d[v] = row[v] - firm_means[f][v] - time_means[t][v] + overall[v]
        d['event_id'] = row.get('event_id')
        out.append(d)
    return out


# ----------------------------- Load data -----------------------------

monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass

vwretd = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
if not vwretd:
    raise RuntimeError('Missing F-F monthly factors for vwretd.')

W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# ----------------------------- Build panel -----------------------------

all_events = []
events_path = EVENTS_PATH if EVENTS_PATH else derived_path('events', 'coal_retirement_events.csv')
with open(events_path, 'r', encoding='utf-8') as f:
    tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()]) if TIER_FILTER else None
    exclude_years = set([y.strip() for y in EXCLUDE_YEARS.split(',') if y.strip()]) if EXCLUDE_YEARS else set()
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        ann_src = row.get('announcement_source', '').strip()
        if EXACT_ONLY:
            if not ann_date:
                continue
            if ann_src and 'proxy' in ann_src.lower():
                continue
        effective_date = ann_date if ann_date else ret_date
        if effective_date and len(effective_date) >= 7:
            event_month = effective_date[:7]
        else:
            event_month = f'{row["ret_year"]}-07'
        event_year = int(event_month[:4]) if event_month[:4].isdigit() else None

        if tiers is not None and row.get('exogeneity_tier', ''):
            if row.get('exogeneity_tier', '').strip() not in tiers:
                continue
        if BINDING_ONLY and row.get('binding', ''):
            if row.get('binding', '').strip().lower() != 'yes':
                continue
        if exclude_years and event_month[:4] in exclude_years:
            continue

        all_events.append({
            'event_month': event_month,
            'event_year': event_year,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

events = [e for e in all_events if e['is_first_mover']] if EVENT_SCOPE == 'first_mover' else list(all_events)

obs = []
for event_id, event in enumerate(events):
    event_month = event['event_month']
    if not event_month:
        continue
    event_gvkeys = set(event['gvkeys'])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        fm_fic = None
        frow_fm = fundamentals.get(fm_gk)
        if frow_fm and frow_fm.get('fic'):
            fm_fic = frow_fm.get('fic')
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
        stable_seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
        random.seed(stable_seed)
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = random.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
        candidate_firms = list(neighbor_gks) + ctrl_sample

        for gk in candidate_firms:
            w_ij = neighbors.get(gk, 0.0)
            if FOREIGN_ONLY and fm_fic:
                fic = fundamentals.get(gk, {}).get('fic')
                if fic == fm_fic:
                    w_ij = 0.0
            if gk not in monthly_ret:
                continue
            for tau in range(TAU_START, TAU_END + 1):
                ym = add_months(event_month, tau)
                if ym not in monthly_ret[gk]:
                    continue
                ar = monthly_ret[gk][ym] - vwretd.get(ym, 0.0)
                post = 1.0 if (tau >= POST_START and tau <= POST_END) else 0.0
                d = 1.0 if w_ij > 0 else 0.0
                wpos = math.log1p(w_ij) if w_ij > 0 else 0.0
                obs.append({
                    'gvkey': gk,
                    'ym': ym,
                    'ar': ar,
                    'event_id': event_id,
                    'exp_post_d': d * post,
                    'exp_post_w': wpos * post,
                })

# Overlap handling: keep nearest event by |tau|
if OVERLAP_RULE == 'nearest':
    grouped = defaultdict(list)
    for row in obs:
        key = (row['gvkey'], row['ym'])
        grouped[key].append(row)
    obs = []
    for _, rows in grouped.items():
        # keep one (all rows in same month); pick max exp_post_d then exp_post_w
        rows_sorted = sorted(rows, key=lambda r: (r['exp_post_d'], r['exp_post_w']), reverse=True)
        obs.append(rows_sorted[0])

# Two-way demean (firm + month)
obs_dm = two_way_demean(obs, 'ar', ['exp_post_d', 'exp_post_w'], firm_key='gvkey', time_key='ym')

res_event = ols_cluster(obs_dm, 'ar', ['exp_post_d', 'exp_post_w'], 'event_id')
res_firm = ols_cluster(obs_dm, 'ar', ['exp_post_d', 'exp_post_w'], 'gvkey')

# Two-way cluster: build pair clusters (event, gvkey)
clusters = defaultdict(list)
for i, row in enumerate(obs_dm):
    clusters[(row['event_id'], row['gvkey'])].append(i)

def ols_cluster_custom(data, y_var, x_vars, clusters_map):
    n = len(data)
    k = len(x_vars) + 1
    X = [[1.0] + [row[v] for v in x_vars] for row in data]
    y = [row[y_var] for row in data]
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    resid = [y[i] - sum(beta[a] * X[i][a] for a in range(k)) for i in range(n)]
    S, G = _cluster_cov(X, resid, clusters_map)
    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
    se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]
    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'n': n,
    }

res_tw = ols_cluster_custom(obs_dm, 'ar', ['exp_post_d', 'exp_post_w'], clusters)

def fmt(res, label):
    if not res:
        return f'{label}: failed'
    b1 = res['beta'].get('exp_post_d', 0.0)
    b2 = res['beta'].get('exp_post_w', 0.0)
    se1 = res['se'].get('exp_post_d', 0.0)
    se2 = res['se'].get('exp_post_w', 0.0)
    t1 = res['t'].get('exp_post_d', 0.0)
    t2 = res['t'].get('exp_post_w', 0.0)
    return f'{label}: d={b1:+.4f} (se {se1:.4f}, t {t1:.2f}); w={b2:+.4f} (se {se2:.4f}, t {t2:.2f}); N={res["n"]}'

print(fmt(res_event, 'Event-clustered'))
print(fmt(res_firm, 'Firm-clustered'))
print(fmt(res_tw, 'Two-way'))

if WRITE_METRICS:
    suffix = f'_{METRICS_SUFFIX}' if METRICS_SUFFIX else ''
    out_path = os.path.join('JEEM_submission_package', 'JEEM_outputs', 'metrics',
                            f'strategy2_panel_did_two_part{suffix}.md')
    lines = [
        '# Strategy 2 Panel DiD (Two-Part Exposure)',
        '',
        f'- event_scope: {EVENT_SCOPE}',
        f'- tier_filter: {TIER_FILTER}',
        f'- binding_only: {BINDING_ONLY}',
        f'- exclude_years: {EXCLUDE_YEARS}',
        f'- N: {res_event["n"] if res_event else 0}',
        '',
        '## exp_post_d and exp_post_w',
        fmt(res_event, 'Event-clustered'),
        fmt(res_firm, 'Firm-clustered'),
        fmt(res_tw, 'Two-way'),
        '',
    ]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'Wrote: {out_path}')

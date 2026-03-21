"""Strategy 2 portfolio sorts: Fama-French methodology for spatial exposure.

Non-parametric test of whether spatial exposure predicts returns around
first-mover coal retirement events.  Sorts firms into quintiles by
fuel-similarity (W_fuel) and geographic proximity (W_geo), then computes
equal-weighted mean CARs per quintile.

Panel A: Fuel-similarity quintile sorts
Panel B: Geographic-proximity quintile sorts
Panel C: Channel split (geo spread minus fuel spread)
Panel D: Long-short portfolio (high geo + low fuel vs low geo + high fuel)
"""
import csv
import math
import os
from collections import defaultdict

from _paths import raw_path, derived_path, results_path


# ── Helpers ─────────────────────────────────────────────────────────

def add_months(ym, delta):
    """Shift a YYYY-MM string by delta months."""
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
    """Load Fama-French monthly factors; return (mktrf, rf, vwretd) dicts."""
    if not os.path.exists(path):
        return None
    mktrf = {}
    rf = {}
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The '):
                continue
            if line.startswith(','):
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
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f'{date[:4]}-{date[4:6]}'
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


def is_exact_source(src):
    if not src:
        return False
    s = src.lower()
    if 'proxy' in s or 'approx' in s or 'mid' in s or 'month' in s:
        return False
    return True


def assign_quintiles(values):
    """Given a list of (key, value) pairs, return dict key -> quintile (1..5).

    Q1 = lowest values, Q5 = highest values.
    """
    if not values:
        return {}
    sorted_vals = sorted(values, key=lambda x: x[1])
    n = len(sorted_vals)
    assignment = {}
    for rank, (key, _) in enumerate(sorted_vals):
        q = min(int(rank * 5 / n) + 1, 5)
        assignment[key] = q
    return assignment


# ── Load data ───────────────────────────────────────────────────────

print('Loading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass
print(f'  Monthly returns: {len(monthly_ret)} firms')

print('Loading Fama-French monthly factors...')
ff_monthly_path = raw_path('factors', 'F-F_Research_Data_Factors.csv')
ff_factors = load_ff_factors_monthly(ff_monthly_path)
if not ff_factors:
    raise RuntimeError('Missing Fama-French monthly factors for vwretd.')
mktrf_monthly, rf_monthly, market_ret = ff_factors
print(f'  Market months (F-F vwretd): {len(market_ret)}')

print('Loading geographic weight matrix (W_geo)...')
W_geo = defaultdict(dict)
geo_path = derived_path('networks', 'weight_matrix_W_geo.csv')
with open(geo_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
print(f'  Firms in W_geo: {len(W_geo)}')

print('Loading fuel-similarity weight matrix (W_fuel)...')
W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
    print(f'  Firms in W_fuel: {len(W_fuel)}')
else:
    raise RuntimeError('Missing W_fuel. Run build_fuel_matrix.py first.')

print('Loading retirement events...')
all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        elif row.get('ret_year'):
            event_year = int(row['ret_year'])
        if effective_date and len(effective_date) >= 7:
            event_month = effective_date[:7]
        elif event_year:
            event_month = f'{event_year}-07'
        else:
            continue
        all_events.append({
            'plant': row.get('plant_name', ''),
            'year': event_year,
            'event_month': event_month,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

events = [e for e in all_events if e['is_first_mover']]
print(f'  First-mover events: {len(events)} (of {len(all_events)} total)')


# ── Compute CARs ───────────────────────────────────────────────────

# CAR window: [-1, +3] months (5-month window)
CAR_PRE = 1   # months before event
CAR_POST = 3  # months after event
MIN_FIRMS = 10  # minimum firms per event for meaningful quintile sorts

print(f'\nComputing monthly CARs [{-CAR_PRE}, +{CAR_POST}]...')


def compute_monthly_car(gvkey, event_month):
    """Compute CAR[-1, +3] as sum of abnormal returns (vwretd model)."""
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    # Find event_idx: first month >= event_month
    event_idx = None
    for i, m in enumerate(months):
        if m >= event_month:
            event_idx = i
            break
    if event_idx is None:
        return None
    # Require some pre-event data for validity
    if event_idx < 6:
        return None
    car = 0.0
    valid_months = 0
    for offset in range(-CAR_PRE, CAR_POST + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in monthly_ret[gvkey] and m in market_ret:
                ar = monthly_ret[gvkey][m] - market_ret[m]
                car += ar
                valid_months += 1
    if valid_months < (CAR_PRE + CAR_POST + 1) // 2:
        return None
    return car


# ── Event-level portfolio sorts ────────────────────────────────────

# Accumulators for event-level results
fuel_quintile_cars = {q: [] for q in range(1, 6)}  # q -> list of event-level mean CARs
fuel_quintile_counts = {q: [] for q in range(1, 6)}
geo_quintile_cars = {q: [] for q in range(1, 6)}
geo_quintile_counts = {q: [] for q in range(1, 6)}
fuel_spreads = []   # Q5-Q1 spread per event (fuel)
geo_spreads = []    # Q5-Q1 spread per event (geo)
channel_diffs = []  # geo_spread - fuel_spread per event
long_short_cars = []  # long-short portfolio CAR per event
events_processed = 0
events_skipped = 0

print(f'Processing {len(events)} first-mover events...')

for event_id, event in enumerate(events):
    event_month = event['event_month']
    event_gvkeys = set(event['gvkeys'])

    # Collect all neighbor firms and their exposures + CARs
    firm_data = {}  # gvkey -> {car, w_fuel, w_geo}

    for fm_gk in event_gvkeys:
        # Get all firms connected to this first-mover via either matrix
        all_neighbors = set()
        if fm_gk in W_geo:
            all_neighbors.update(W_geo[fm_gk].keys())
        if fm_gk in W_fuel:
            all_neighbors.update(W_fuel[fm_gk].keys())
        # Exclude event firms themselves
        all_neighbors -= event_gvkeys

        for gk in all_neighbors:
            if gk in firm_data:
                # If firm appears for multiple first-movers in same event,
                # keep the maximum exposure (strongest link)
                w_fuel_val = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                w_geo_val = W_geo.get(fm_gk, {}).get(gk, 0.0)
                firm_data[gk]['w_fuel'] = max(firm_data[gk]['w_fuel'], w_fuel_val)
                firm_data[gk]['w_geo'] = max(firm_data[gk]['w_geo'], w_geo_val)
                continue
            car = compute_monthly_car(gk, event_month)
            if car is None:
                continue
            w_fuel_val = W_fuel.get(fm_gk, {}).get(gk, 0.0)
            w_geo_val = W_geo.get(fm_gk, {}).get(gk, 0.0)
            firm_data[gk] = {
                'car': car,
                'w_fuel': w_fuel_val,
                'w_geo': w_geo_val,
            }

    n_firms = len(firm_data)
    if n_firms < MIN_FIRMS:
        events_skipped += 1
        continue

    events_processed += 1
    gvkeys_list = list(firm_data.keys())

    # --- Panel A: Fuel-similarity quintile sorts ---
    fuel_pairs = [(gk, firm_data[gk]['w_fuel']) for gk in gvkeys_list]
    fuel_q = assign_quintiles(fuel_pairs)

    quintile_car_fuel = {q: [] for q in range(1, 6)}
    for gk in gvkeys_list:
        q = fuel_q[gk]
        quintile_car_fuel[q].append(firm_data[gk]['car'])

    fuel_means = {}
    for q in range(1, 6):
        vals = quintile_car_fuel[q]
        if vals:
            fuel_means[q] = sum(vals) / len(vals)
            fuel_quintile_cars[q].append(fuel_means[q])
            fuel_quintile_counts[q].append(len(vals))
        else:
            fuel_means[q] = None

    if fuel_means.get(5) is not None and fuel_means.get(1) is not None:
        fuel_spread = fuel_means[5] - fuel_means[1]
        fuel_spreads.append(fuel_spread)
    else:
        fuel_spread = None

    # --- Panel B: Geographic-proximity quintile sorts ---
    geo_pairs = [(gk, firm_data[gk]['w_geo']) for gk in gvkeys_list]
    geo_q = assign_quintiles(geo_pairs)

    quintile_car_geo = {q: [] for q in range(1, 6)}
    for gk in gvkeys_list:
        q = geo_q[gk]
        quintile_car_geo[q].append(firm_data[gk]['car'])

    geo_means = {}
    for q in range(1, 6):
        vals = quintile_car_geo[q]
        if vals:
            geo_means[q] = sum(vals) / len(vals)
            geo_quintile_cars[q].append(geo_means[q])
            geo_quintile_counts[q].append(len(vals))
        else:
            geo_means[q] = None

    if geo_means.get(5) is not None and geo_means.get(1) is not None:
        geo_spread = geo_means[5] - geo_means[1]
        geo_spreads.append(geo_spread)
    else:
        geo_spread = None

    # --- Panel C: Channel split ---
    if geo_spread is not None and fuel_spread is not None:
        channel_diffs.append(geo_spread - fuel_spread)

    # --- Panel D: Long-short portfolio ---
    # Long: top quintile geo AND bottom quintile fuel
    # Short: bottom quintile geo AND top quintile fuel
    long_cars = []
    short_cars = []
    for gk in gvkeys_list:
        g_q = geo_q[gk]
        f_q = fuel_q[gk]
        if g_q == 5 and f_q == 1:
            long_cars.append(firm_data[gk]['car'])
        elif g_q == 1 and f_q == 5:
            short_cars.append(firm_data[gk]['car'])

    if long_cars and short_cars:
        long_mean = sum(long_cars) / len(long_cars)
        short_mean = sum(short_cars) / len(short_cars)
        long_short_cars.append(long_mean - short_mean)

    if events_processed % 10 == 0:
        print(f'  Processed {events_processed} events ({n_firms} firms in latest)...')

print(f'\nEvents processed: {events_processed}, skipped (< {MIN_FIRMS} firms): {events_skipped}')


# ── Compute statistics ─────────────────────────────────────────────

def ts_mean_tstat(vals):
    """Time-series mean and t-statistic for a list of event-level values."""
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    n = len(vals)
    if n == 0:
        return None, None, 0
    mean = sum(vals) / n
    if n < 2:
        return mean, None, n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1)
    sd = math.sqrt(var) if var > 0 else 0.0
    t = mean / (sd / math.sqrt(n)) if sd > 0 else 0.0
    return mean, t, n


def fmt_pct(val):
    """Format a return value as percentage string."""
    if val is None:
        return 'N/A'
    return f'{val * 100:+.2f}%'


def fmt_t(t):
    """Format a t-statistic."""
    if t is None:
        return 'N/A'
    return f'{t:.2f}'


# ── Panel A: Fuel-similarity results ──────────────────────────────

print('\n=== Panel A: Fuel-Similarity Quintiles ===')
panel_a_rows = []
for q in range(1, 6):
    mean_car, _, n_events = ts_mean_tstat(fuel_quintile_cars[q])
    avg_firms = sum(fuel_quintile_counts[q]) / len(fuel_quintile_counts[q]) if fuel_quintile_counts[q] else 0
    label = f'Q{q}'
    if q == 1:
        label += ' (lowest fuel sim)'
    elif q == 5:
        label += ' (highest fuel sim)'
    print(f'  {label}: mean CAR = {fmt_pct(mean_car)}, avg N_firms = {avg_firms:.1f}')
    panel_a_rows.append((label, mean_car, avg_firms))

fuel_spread_mean, fuel_spread_t, fuel_spread_n = ts_mean_tstat(fuel_spreads)
print(f'  Q5-Q1 spread: {fmt_pct(fuel_spread_mean)} (t = {fmt_t(fuel_spread_t)}, N = {fuel_spread_n})')

# ── Panel B: Geographic-proximity results ─────────────────────────

print('\n=== Panel B: Geographic Proximity Quintiles ===')
panel_b_rows = []
for q in range(1, 6):
    mean_car, _, n_events = ts_mean_tstat(geo_quintile_cars[q])
    avg_firms = sum(geo_quintile_counts[q]) / len(geo_quintile_counts[q]) if geo_quintile_counts[q] else 0
    label = f'Q{q}'
    if q == 1:
        label += ' (most distant)'
    elif q == 5:
        label += ' (closest)'
    print(f'  {label}: mean CAR = {fmt_pct(mean_car)}, avg N_firms = {avg_firms:.1f}')
    panel_b_rows.append((label, mean_car, avg_firms))

geo_spread_mean, geo_spread_t, geo_spread_n = ts_mean_tstat(geo_spreads)
print(f'  Q5-Q1 spread: {fmt_pct(geo_spread_mean)} (t = {fmt_t(geo_spread_t)}, N = {geo_spread_n})')

# ── Panel C: Channel split ───────────────────────────────────────

print('\n=== Panel C: Channel Split ===')
chan_mean, chan_t, chan_n = ts_mean_tstat(channel_diffs)
print(f'  geo_spread - fuel_spread = {fmt_pct(chan_mean)} (t = {fmt_t(chan_t)}, N = {chan_n})')

# ── Panel D: Long-short portfolio ────────────────────────────────

print('\n=== Panel D: Long-Short Portfolio ===')
ls_mean, ls_t, ls_n = ts_mean_tstat(long_short_cars)
print(f'  Long (high geo + low fuel) - Short (low geo + high fuel)')
print(f'  Mean CAR: {fmt_pct(ls_mean)} (t = {fmt_t(ls_t)})')
print(f'  Events with valid sorts: {ls_n}')


# ── Write output ────────────────────────────────────────────────────

out_path = results_path('metrics', 'strategy2_portfolio_sorts.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = []
lines.append('# Portfolio Sorts: Spatial Exposure and Returns')
lines.append('')
lines.append(f'CAR window: [{-CAR_PRE}, +{CAR_POST}] months | '
             f'Min firms per event: {MIN_FIRMS} | '
             f'Events processed: {events_processed} | '
             f'Events skipped: {events_skipped}')
lines.append('')

# Panel A
lines.append('## Panel A: Fuel-Similarity Quintiles ([-1,+3] month CARs)')
lines.append('')
lines.append('| Quintile | Mean CAR | N_firms (avg) |')
lines.append('|---|---|---|')
for label, mean_car, avg_firms in panel_a_rows:
    lines.append(f'| {label} | {fmt_pct(mean_car)} | {avg_firms:.1f} |')
lines.append(f'| Q5 - Q1 (spread) | {fmt_pct(fuel_spread_mean)} (t = {fmt_t(fuel_spread_t)}) | |')
lines.append('')

# Panel B
lines.append('## Panel B: Geographic Proximity Quintiles ([-1,+3] month CARs)')
lines.append('')
lines.append('| Quintile | Mean CAR | N_firms (avg) |')
lines.append('|---|---|---|')
for label, mean_car, avg_firms in panel_b_rows:
    lines.append(f'| {label} | {fmt_pct(mean_car)} | {avg_firms:.1f} |')
lines.append(f'| Q5 - Q1 (spread) | {fmt_pct(geo_spread_mean)} (t = {fmt_t(geo_spread_t)}) | |')
lines.append('')

# Panel C
lines.append('## Panel C: Channel Split')
lines.append('')
lines.append(f'geo_spread - fuel_spread = {fmt_pct(chan_mean)} (t = {fmt_t(chan_t)}, N = {chan_n})')
lines.append('')

# Panel D
lines.append('## Panel D: Long-Short Portfolio')
lines.append('')
lines.append('Long (high geo + low fuel) vs Short (low geo + high fuel)')
lines.append(f'Mean CAR: {fmt_pct(ls_mean)} (t = {fmt_t(ls_t)})')
lines.append(f'N events with valid sorts: {ls_n}')
lines.append('')

# Interpretation
lines.append('## Interpretation')
lines.append('')
lines.append('The portfolio sort methodology does not depend on regression functional')
lines.append('form, standard error specification, or multiple hypothesis testing')
lines.append('corrections. A significant Q5-Q1 spread establishes that spatial exposure')
lines.append('predicts returns non-parametrically.')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')

print(f'\nResults written to {out_path}')

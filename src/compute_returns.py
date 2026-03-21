"""Compute daily and monthly equity returns.

Data hierarchy (CRSP preferred for US, Datastream for non-US if available):
  1. CRSP/Compustat Merged daily + monthly (US-listed firms, gold standard)
  2. Datastream monthly total returns (non-US, if provided)
  3. Compustat Global Security daily + monthly (non-US fallback)

For firms present in both, CRSP takes priority.

CRSP daily:   total return = (prccd_t * trfd_t) / (prccd_{t-1} * trfd_{t-1}) - 1
CRSP monthly: total return = trt1m / 100  (already provided by CRSP)
Compustat daily:  same formula via trfd
Compustat monthly: price return = prccm_t / prccm_{t-1} - 1
Datastream monthly (if available): use provided total return series

Filters: price >= $1, daily volume >= 1000, returns capped at +/- 50%.
"""
import csv, os
from collections import defaultdict

from _paths import raw_path, derived_path

# CRSP/Compustat Merged files (may be split across batches)
CRSP_DAILY_FILES = [
    'crsp_ccm_daily_part1.csv',   # batch 1
    'crsp_ccm_daily_part2.csv',   # batch 2
]
CRSP_MONTHLY_FILES = [
    'crsp_ccm_monthly_part1.csv',   # batch 1
    'crsp_ccm_monthly_part2.csv',   # batch 2
]
# Compustat Global (non-US fallback)
GLOBAL_DAILY = raw_path('compustat', 'compustat_global_security_daily.csv')
GLOBAL_MONTHLY = raw_path('compustat', 'compustat_global_security_monthly.csv')

# Optional Datastream monthly total returns (non-US) with gvkey mapping
DATASTREAM_MONTHLY_FILES = [
    raw_path('datastream', 'datastream_monthly_returns.csv'),
]


# ====================================================================
# DAILY RETURNS
# ====================================================================
print('=== DAILY RETURNS ===')

# Step 1: Load CRSP daily prices (GVKEY column is uppercase)
crsp_daily = defaultdict(list)  # gvkey -> [(date, price, trfd, volume)]
crsp_gvkeys = set()

for fname in CRSP_DAILY_FILES:
    fpath = raw_path('crsp_compustat', fname)
    if not os.path.exists(fpath):
        continue
    print(f'  Loading CRSP daily: {fname}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['GVKEY']
            date = row['datadate']
            prccd = row.get('prccd', '')
            trfd = row.get('trfd', '')
            cshtrd = row.get('cshtrd', '')
            if prccd and trfd and date:
                try:
                    p = float(prccd)
                    t = float(trfd)
                    v = float(cshtrd) if cshtrd else 0.0
                    crsp_daily[gk].append((date, p, t, v))
                    crsp_gvkeys.add(gk)
                except ValueError:
                    pass

print(f'  CRSP daily firms: {len(crsp_daily)}')

# Step 2: Load Compustat Global daily (skip gvkeys already in CRSP)
global_daily = defaultdict(list)
fpath = GLOBAL_DAILY
if os.path.exists(fpath):
    print(f'  Loading Compustat Global daily: {GLOBAL_DAILY}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            if gk in crsp_gvkeys:
                continue  # CRSP takes priority
            date = row['datadate']
            prccd = row.get('prccd', '')
            trfd = row.get('trfd', '')
            cshtrd = row.get('cshtrd', '')
            if prccd and trfd and date:
                try:
                    p = float(prccd)
                    t = float(trfd)
                    v = float(cshtrd) if cshtrd else 0.0
                    global_daily[gk].append((date, p, t, v))
                except ValueError:
                    pass
    print(f'  Compustat Global daily firms (non-CRSP): {len(global_daily)}')

# Step 3: Compute daily returns from both sources (same formula)
daily_returns = []
all_daily = {**crsp_daily, **global_daily}
crsp_obs = 0
global_obs = 0

for gk in all_daily:
    data = sorted(all_daily[gk], key=lambda x: x[0])
    is_crsp = gk in crsp_gvkeys
    for i in range(1, len(data)):
        date_t, p_t, trf_t, v_t = data[i]
        date_tm1, p_tm1, trf_tm1, v_tm1 = data[i - 1]
        if p_tm1 <= 0 or trf_tm1 <= 0 or p_t <= 0 or trf_t <= 0:
            continue
        if p_t < 1 or p_tm1 < 1:
            continue
        if v_t < 1000 or v_tm1 < 1000:
            continue
        ri_t = (p_t * trf_t) / (p_tm1 * trf_tm1) - 1
        if ri_t > 0.5:
            ri_t = 0.5
        elif ri_t < -0.5:
            ri_t = -0.5
        daily_returns.append({
            'gvkey': gk,
            'datadate': date_t,
            'prccd': f'{p_t:.4f}',
            'ret_daily': f'{ri_t:.6f}',
        })
        if is_crsp:
            crsp_obs += 1
        else:
            global_obs += 1

outpath = derived_path('returns', 'daily_returns.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['gvkey', 'datadate', 'prccd', 'ret_daily'])
    w.writeheader()
    w.writerows(daily_returns)

n_firms_d = len(set(r['gvkey'] for r in daily_returns))
print(f'\nDaily returns: {len(daily_returns):,} obs, {n_firms_d} firms')
print(f'  CRSP source:     {crsp_obs:,} obs ({len(crsp_daily)} firms)')
print(f'  Global source:   {global_obs:,} obs ({len(global_daily)} firms)')
print(f'Saved {outpath}')

# Distribution
rets = sorted(float(r['ret_daily']) for r in daily_returns)
n = len(rets)
if n > 0:
    print(f'\nDaily return distribution:')
    print(f'  Min={rets[0]:.4f}  P5={rets[5*n//100]:.4f}  '
          f'Median={rets[n//2]:.4f}  P95={rets[95*n//100]:.4f}  Max={rets[-1]:.4f}')


# ====================================================================
# MONTHLY RETURNS
# ====================================================================
print('\n=== MONTHLY RETURNS ===')

# Step 1: Load CRSP monthly (uses trt1m = CRSP total return, already computed)
crsp_monthly_ret = []  # list of dicts
crsp_monthly_gvkeys = set()

for fname in CRSP_MONTHLY_FILES:
    fpath = raw_path('crsp_compustat', fname)
    if not os.path.exists(fpath):
        continue
    print(f'  Loading CRSP monthly: {fname}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['GVKEY']
            date = row['datadate']
            prccm = row.get('prccm', '')
            trt1m = row.get('trt1m', '')
            if not date or not prccm:
                continue
            try:
                p = float(prccm)
            except ValueError:
                continue
            if p < 1:
                continue
            # Use trt1m (CRSP total return in %) if available
            ret = None
            if trt1m and trt1m.strip():
                try:
                    ret = float(trt1m) / 100.0  # CRSP returns are in percent
                except ValueError:
                    pass
            if ret is not None:
                if ret > 0.5:
                    ret = 0.5
                elif ret < -0.5:
                    ret = -0.5
                crsp_monthly_ret.append({
                    'gvkey': gk,
                    'datadate': date,
                    'price': f'{p:.4f}',
                    'ret_monthly': f'{ret:.6f}',
                })
                crsp_monthly_gvkeys.add(gk)

print(f'  CRSP monthly: {len(crsp_monthly_ret):,} obs, {len(crsp_monthly_gvkeys)} firms')

# Step 2: Load optional Datastream monthly total returns (non-US)
def _find_col(row, candidates):
    for c in candidates:
        if c in row and row[c] != '':
            return c
    return None


def load_datastream_monthly(files, crsp_gvkeys):
    ds_ret = []
    ds_gvkeys = set()
    for fname in files:
        fpath = fname
        if not os.path.exists(fpath):
            continue
        print(f'  Loading Datastream monthly: {fname}...')
        with open(fpath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Require gvkey mapping
                gk = row.get('gvkey') or row.get('GVKEY') or row.get('gvkey_i')
                if not gk or gk in crsp_gvkeys:
                    continue
                date = row.get('datadate') or row.get('date') or row.get('DATE')
                if not date:
                    continue
                # Determine return column
                ret_col = _find_col(row, ['ret_monthly', 'ret', 'return', 'trt1m', 'retm'])
                if not ret_col:
                    continue
                try:
                    ret = float(row[ret_col])
                except ValueError:
                    continue
                # If in percent, scale
                if abs(ret) > 2.0:
                    ret = ret / 100.0
                # Format date to YYYY-MM-DD for compatibility
                date_str = str(date).strip()
                if date_str.isdigit() and len(date_str) == 6:
                    date_fmt = f'{date_str[:4]}-{date_str[4:6]}-01'
                elif date_str.isdigit() and len(date_str) == 8:
                    date_fmt = f'{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}'
                else:
                    # Assume already YYYY-MM-DD or YYYY-MM; standardize to YYYY-MM-01 if needed
                    date_fmt = date_str if len(date_str) >= 10 else f'{date_str}-01'
                if ret > 0.5:
                    ret = 0.5
                elif ret < -0.5:
                    ret = -0.5
                ds_ret.append({
                    'gvkey': gk,
                    'datadate': date_fmt,
                    'price': '',
                    'ret_monthly': f'{ret:.6f}',
                })
                ds_gvkeys.add(gk)
    return ds_ret, ds_gvkeys


datastream_monthly_ret, datastream_monthly_gvkeys = load_datastream_monthly(
    DATASTREAM_MONTHLY_FILES, crsp_monthly_gvkeys
)
if datastream_monthly_ret:
    print(f'  Datastream monthly: {len(datastream_monthly_ret):,} obs, {len(datastream_monthly_gvkeys)} firms')
else:
    print('  Datastream monthly: none found (skipping)')

# Step 3: Load Compustat Global monthly (price-based returns, skip CRSP and Datastream firms)
global_monthly_prices = defaultdict(list)
fpath = GLOBAL_MONTHLY
if os.path.exists(fpath):
    print(f'  Loading Compustat Global monthly: {GLOBAL_MONTHLY}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            if gk in crsp_monthly_gvkeys or gk in datastream_monthly_gvkeys:
                continue  # CRSP takes priority
            date = row['datadate']
            prccm = row.get('prccm', '')
            if prccm and date:
                try:
                    p = float(prccm)
                    if p >= 1:
                        global_monthly_prices[gk].append((date, p))
                except ValueError:
                    pass
    print(f'  Compustat Global monthly firms (non-CRSP): {len(global_monthly_prices)}')

# Compute price returns for Global firms
global_monthly_ret = []
for gk in global_monthly_prices:
    prices = sorted(global_monthly_prices[gk], key=lambda x: x[0])
    for i in range(1, len(prices)):
        date_t, p_t = prices[i]
        date_tm1, p_tm1 = prices[i - 1]
        if p_tm1 >= 1 and p_t >= 1:
            ret = (p_t / p_tm1) - 1
            if ret > 0.5:
                ret = 0.5
            elif ret < -0.5:
                ret = -0.5
            global_monthly_ret.append({
                'gvkey': gk,
                'datadate': date_t,
                'price': f'{p_t:.4f}',
                'ret_monthly': f'{ret:.6f}',
            })

# Step 4: Merge and write
monthly_returns = crsp_monthly_ret + datastream_monthly_ret + global_monthly_ret

outpath = derived_path('returns', 'monthly_returns.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['gvkey', 'datadate', 'price', 'ret_monthly'])
    w.writeheader()
    w.writerows(monthly_returns)

n_firms_m = len(set(r['gvkey'] for r in monthly_returns))
n_crsp_m = len(crsp_monthly_gvkeys)
n_ds_m = len(datastream_monthly_gvkeys)
n_global_m = len(global_monthly_prices)
print(f'\nMonthly returns: {len(monthly_returns):,} obs, {n_firms_m} firms')
print(f'  CRSP source (trt1m):    {len(crsp_monthly_ret):,} obs ({n_crsp_m} firms)')
print(f'  Datastream source:      {len(datastream_monthly_ret):,} obs ({n_ds_m} firms)')
print(f'  Global source (price):  {len(global_monthly_ret):,} obs ({n_global_m} firms)')
print(f'Saved {outpath}')

# Distribution
rets_m = sorted(float(r['ret_monthly']) for r in monthly_returns)
nm = len(rets_m)
if nm > 0:
    print(f'\nMonthly return distribution:')
    print(f'  Min={rets_m[0]:.4f}  P5={rets_m[5*nm//100]:.4f}  '
          f'Median={rets_m[nm//2]:.4f}  P95={rets_m[95*nm//100]:.4f}  Max={rets_m[-1]:.4f}')

# Summary
print(f'\n=== SUMMARY ===')
print(f'Daily:   {len(daily_returns):,} obs, {n_firms_d} firms '
      f'(CRSP: {len(crsp_daily)}, Global: {len(global_daily)})')
print(f'Monthly: {len(monthly_returns):,} obs, {n_firms_m} firms '
      f'(CRSP: {n_crsp_m}, Datastream: {n_ds_m}, Global: {n_global_m})')

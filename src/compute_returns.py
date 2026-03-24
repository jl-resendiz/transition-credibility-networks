"""Compute daily and monthly equity returns.

Data sources:
  Daily:  CRSP (US) + Compustat Global Security (non-US, deduplicated)
  Monthly: CRSP total returns (US) + Eikon TR.TotalReturn (non-US)

Monthly returns use total return series (including dividends) where available:
  - CRSP trt1m (US, gold standard)
  - Eikon TR.TotalReturn via pull_eikon_returns.py (non-US, preferred)
  - Compustat Global Security price returns (non-US fallback, deduplicated)
    Used only for ~100 firms not covered by Eikon. Share-class duplicates
    removed per Ince & Porter (2006, JFE).

Daily returns use CRSP + Compustat Global Security (deduplicated) because
Eikon does not provide daily total return series via the batch API.

Filters per Ince & Porter (2006) and Griffin, Kelly & Nardari (2010):
  - Price >= $1 (penny stock filter)
  - Daily volume >= 1000 shares
  - Monthly returns capped at +/- 100% (daily at +/- 50%)
  - Eikon: end-of-month only, illiquid firms removed (>50% months at cap)
  - Compustat Global daily: deduplicated by gvkey-date (highest volume)
"""
import csv, os
from collections import defaultdict

from _paths import raw_path, derived_path

# Return caps per Ince & Porter (2006): +/-50% daily, +/-100% monthly
DAILY_CAP = 0.50
MONTHLY_CAP = 1.00

# CRSP/Compustat Merged files
CRSP_DAILY_FILES = [
    'crsp_ccm_daily_part1.csv',
    'crsp_ccm_daily_part2.csv',
]
CRSP_MONTHLY_FILES = [
    'crsp_ccm_monthly_part1.csv',
    'crsp_ccm_monthly_part2.csv',
]
# Compustat Global daily only (for non-US daily returns)
GLOBAL_DAILY = raw_path('compustat', 'compustat_global_security_daily.csv')

# Eikon monthly total returns (non-US)
EIKON_MONTHLY_FILES = [
    raw_path('eikon', 'eikon_monthly_returns.csv'),
]


# ── Helper: deduplicate Compustat Global share classes (daily only) ──
def deduplicate_global(rows_by_gvkey):
    """For each gvkey, deduplicate by date keeping the highest-volume observation."""
    deduped = {}
    for gk, rows in rows_by_gvkey.items():
        by_date = defaultdict(list)
        for row in rows:
            by_date[row[0]].append(row)
        clean = []
        for date, obs_list in sorted(by_date.items()):
            if len(obs_list) == 1:
                clean.append(obs_list[0])
            else:
                best = max(obs_list, key=lambda x: (x[3], x[1]))
                clean.append(best)
        deduped[gk] = clean
    return deduped


# ── Helper: clean Eikon returns ──
def clean_eikon(ds_ret):
    """End-of-month only + remove illiquid firms (>50% months at cap)."""
    by_gk = defaultdict(list)
    for r in ds_ret:
        by_gk[r['gvkey']].append(r)

    cleaned = []
    removed_firms = 0
    for gk, obs in by_gk.items():
        by_month = defaultdict(list)
        for r in obs:
            ym = r['datadate'][:7]
            if len(ym) == 6:
                ym = ym[:4] + '-' + ym[4:6]
            by_month[ym].append(r)

        eom = []
        for ym, month_obs in sorted(by_month.items()):
            best = max(month_obs, key=lambda x: x['datadate'])
            eom.append(best)

        n_total = len(eom)
        if n_total < 6:
            removed_firms += 1
            continue
        n_at_cap = sum(1 for r in eom if abs(float(r['ret_monthly'])) >= MONTHLY_CAP * 0.99)
        if n_at_cap / n_total > 0.50:
            removed_firms += 1
            continue

        cleaned.extend(eom)

    if removed_firms > 0:
        print(f'    Removed {removed_firms} illiquid/short-history firms')
    return cleaned


# ====================================================================
# DAILY RETURNS (CRSP + Compustat Global deduplicated)
# ====================================================================
print('=== DAILY RETURNS ===')

crsp_daily = defaultdict(list)
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

global_daily_raw = defaultdict(list)
fpath = GLOBAL_DAILY
if os.path.exists(fpath):
    print(f'  Loading Compustat Global daily: {GLOBAL_DAILY}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            if gk in crsp_gvkeys:
                continue
            date = row['datadate']
            prccd = row.get('prccd', '')
            trfd = row.get('trfd', '')
            cshtrd = row.get('cshtrd', '')
            if prccd and trfd and date:
                try:
                    p = float(prccd)
                    t = float(trfd)
                    v = float(cshtrd) if cshtrd else 0.0
                    global_daily_raw[gk].append((date, p, t, v))
                except ValueError:
                    pass

    n_before = sum(len(v) for v in global_daily_raw.values())
    global_daily = deduplicate_global(global_daily_raw)
    n_after = sum(len(v) for v in global_daily.values())
    print(f'  Compustat Global daily: {len(global_daily)} firms '
          f'(deduped {n_before:,} -> {n_after:,} obs)')
else:
    global_daily = {}

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
        ri_t = max(-DAILY_CAP, min(DAILY_CAP, ri_t))
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
print(f'  CRSP:   {crsp_obs:,} obs ({len(crsp_daily)} firms)')
print(f'  Global: {global_obs:,} obs ({len(global_daily)} firms)')
print(f'Saved {outpath}')

rets = sorted(float(r['ret_daily']) for r in daily_returns)
n = len(rets)
if n > 0:
    print(f'\nDaily return distribution:')
    print(f'  Min={rets[0]:.4f}  P5={rets[5*n//100]:.4f}  '
          f'Median={rets[n//2]:.4f}  P95={rets[95*n//100]:.4f}  Max={rets[-1]:.4f}')


# ====================================================================
# MONTHLY RETURNS (CRSP + Eikon total returns ONLY)
# ====================================================================
print('\n=== MONTHLY RETURNS ===')

# Step 1: CRSP monthly total returns (US)
crsp_monthly_ret = []
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
            ret = None
            if trt1m and trt1m.strip():
                try:
                    ret = float(trt1m) / 100.0
                except ValueError:
                    pass
            if ret is not None:
                ret = max(-MONTHLY_CAP, min(MONTHLY_CAP, ret))
                crsp_monthly_ret.append({
                    'gvkey': gk,
                    'datadate': date,
                    'price': f'{p:.4f}',
                    'ret_monthly': f'{ret:.6f}',
                })
                crsp_monthly_gvkeys.add(gk)

print(f'  CRSP monthly: {len(crsp_monthly_ret):,} obs, {len(crsp_monthly_gvkeys)} firms')

# Step 2: Eikon monthly total returns (non-US)
def _find_col(row, candidates):
    for c in candidates:
        if c in row and row[c] != '':
            return c
    return None


def load_eikon_monthly(files, crsp_gvkeys):
    ds_ret = []
    ds_gvkeys = set()
    for fpath in files:
        if not os.path.exists(fpath):
            continue
        print(f'  Loading Eikon monthly: {fpath}...')
        with open(fpath, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                gk = row.get('gvkey') or row.get('GVKEY')
                if not gk or gk in crsp_gvkeys:
                    continue
                date = row.get('datadate') or row.get('date')
                if not date:
                    continue
                ret_col = _find_col(row, ['ret_monthly', 'ret', 'return', 'trt1m', 'retm'])
                if not ret_col:
                    continue
                try:
                    ret = float(row[ret_col])
                except ValueError:
                    continue
                date_str = str(date).strip()
                if date_str.isdigit() and len(date_str) == 6:
                    date_fmt = f'{date_str[:4]}-{date_str[4:6]}-01'
                elif date_str.isdigit() and len(date_str) == 8:
                    date_fmt = f'{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}'
                else:
                    date_fmt = date_str if len(date_str) >= 10 else f'{date_str}-01'
                ret = max(-MONTHLY_CAP, min(MONTHLY_CAP, ret))
                ds_ret.append({
                    'gvkey': gk,
                    'datadate': date_fmt,
                    'price': '',
                    'ret_monthly': f'{ret:.6f}',
                })
                ds_gvkeys.add(gk)
    return ds_ret, ds_gvkeys


eikon_monthly_raw, eikon_monthly_gvkeys = load_eikon_monthly(
    EIKON_MONTHLY_FILES, crsp_monthly_gvkeys
)
if eikon_monthly_raw:
    print(f'  Eikon raw: {len(eikon_monthly_raw):,} obs, {len(eikon_monthly_gvkeys)} firms')
    eikon_monthly_ret = clean_eikon(eikon_monthly_raw)
    eikon_monthly_gvkeys = set(r['gvkey'] for r in eikon_monthly_ret)
    print(f'  Eikon clean: {len(eikon_monthly_ret):,} obs, {len(eikon_monthly_gvkeys)} firms')
else:
    eikon_monthly_ret = []
    print('  Eikon monthly: none found')

# Step 3: Compustat Global monthly DEDUPLICATED (fallback for firms without CRSP or Eikon)
# Price returns only — used for ~100 firms not covered by Eikon
GLOBAL_MONTHLY = raw_path('compustat', 'compustat_global_security_monthly.csv')
global_monthly_raw = defaultdict(list)
if os.path.exists(GLOBAL_MONTHLY):
    print(f'  Loading Compustat Global monthly (fallback): {GLOBAL_MONTHLY}...')
    with open(GLOBAL_MONTHLY, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            if gk in crsp_monthly_gvkeys or gk in eikon_monthly_gvkeys:
                continue  # CRSP and Eikon take priority
            date = row['datadate']
            prccm = row.get('prccm', '')
            cshom = row.get('cshom', '')
            if prccm and date:
                try:
                    p = float(prccm)
                    vol = float(cshom) if cshom else 0.0
                    if p >= 1:
                        global_monthly_raw[gk].append((date, p, vol))
                except ValueError:
                    pass

    # Deduplicate by gvkey-date (keep highest volume/price per Ince & Porter 2006)
    global_monthly_prices = {}
    n_raw = sum(len(v) for v in global_monthly_raw.values())
    for gk, rows in global_monthly_raw.items():
        by_date = defaultdict(list)
        for row in rows:
            by_date[row[0]].append(row)
        clean = []
        for date, obs_list in sorted(by_date.items()):
            if len(obs_list) == 1:
                clean.append(obs_list[0])
            else:
                best = max(obs_list, key=lambda x: (x[2], x[1]))
                clean.append(best)
        global_monthly_prices[gk] = clean
    n_clean = sum(len(v) for v in global_monthly_prices.values())
    print(f'  Compustat Global monthly (fallback): {len(global_monthly_prices)} firms '
          f'(deduped {n_raw:,} -> {n_clean:,} obs)')
else:
    global_monthly_prices = {}

global_monthly_ret = []
for gk in global_monthly_prices:
    prices = sorted(global_monthly_prices[gk], key=lambda x: x[0])
    for i in range(1, len(prices)):
        date_t, p_t, _ = prices[i]
        date_tm1, p_tm1, _ = prices[i - 1]
        if p_tm1 >= 1 and p_t >= 1:
            ret = (p_t / p_tm1) - 1
            ret = max(-MONTHLY_CAP, min(MONTHLY_CAP, ret))
            global_monthly_ret.append({
                'gvkey': gk,
                'datadate': date_t,
                'price': f'{p_t:.4f}',
                'ret_monthly': f'{ret:.6f}',
            })

# Step 4: Merge and write
monthly_returns = crsp_monthly_ret + eikon_monthly_ret + global_monthly_ret

outpath = derived_path('returns', 'monthly_returns.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['gvkey', 'datadate', 'price', 'ret_monthly'])
    w.writeheader()
    w.writerows(monthly_returns)

n_firms_m = len(set(r['gvkey'] for r in monthly_returns))
n_crsp_m = len(crsp_monthly_gvkeys)
n_eikon_m = len(eikon_monthly_gvkeys)
n_global_m = len(global_monthly_prices)
print(f'\nMonthly returns: {len(monthly_returns):,} obs, {n_firms_m} firms')
print(f'  CRSP (total return):           {len(crsp_monthly_ret):,} obs ({n_crsp_m} firms)')
print(f'  Eikon (total return):          {len(eikon_monthly_ret):,} obs ({n_eikon_m} firms)')
print(f'  Compustat Global (fallback):   {len(global_monthly_ret):,} obs ({n_global_m} firms)')
print(f'Saved {outpath}')

rets_m = sorted(float(r['ret_monthly']) for r in monthly_returns)
nm = len(rets_m)
if nm > 0:
    print(f'\nMonthly return distribution:')
    print(f'  Min={rets_m[0]:.4f}  P5={rets_m[5*nm//100]:.4f}  '
          f'Median={rets_m[nm//2]:.4f}  P95={rets_m[95*nm//100]:.4f}  Max={rets_m[-1]:.4f}')

print(f'\n=== SUMMARY ===')
print(f'Daily:   {len(daily_returns):,} obs, {n_firms_d} firms (CRSP: {len(crsp_daily)}, Global: {len(global_daily)})')
print(f'Monthly: {len(monthly_returns):,} obs, {n_firms_m} firms '
      f'(CRSP: {n_crsp_m}, Eikon: {n_eikon_m}, Global fallback: {n_global_m})')

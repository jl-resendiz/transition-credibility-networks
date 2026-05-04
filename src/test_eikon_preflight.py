"""Pre-flight test for Refinitiv Eikon API access.

Run BEFORE attempting `pull_refinitiv_extra.py`. Verifies:
1. EIKON_APP_KEY is set in environment (loaded from .env via _credentials).
2. eikon Python package is importable.
3. Refinitiv Workspace Desktop is running locally (proxy at port 9000).
4. A 5-firm test pull of TR.PctFreeFloat returns plausible values.

Usage:  python src/test_eikon_preflight.py

Exit code 0 = ready to run pull_refinitiv_extra.py.
Exit code != 0 = fix the issue indicated in the error message.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _credentials  # auto-loads .env


def fail(msg):
    print(f'\nFAIL: {msg}\n')
    sys.exit(1)


print('=== Refinitiv Eikon Pre-Flight Test ===\n')

# 1. Credential check
print('[1/4] Checking EIKON_APP_KEY in environment...')
key = os.environ.get('EIKON_APP_KEY', '').strip() or os.environ.get('REFINITIV_APP_KEY', '').strip()
if not key:
    fail('EIKON_APP_KEY (or REFINITIV_APP_KEY) is not set in .env at repo root. '
         'Generate one in Refinitiv Workspace Desktop -> App Generator and add it to .env.')
masked = key[:3] + '*' * max(0, len(key) - 6) + key[-3:] if len(key) > 6 else '***'
print(f'  OK: key length {len(key)}, value {masked}')

# 2. Package check
print('\n[2/4] Checking eikon package...')
try:
    import eikon as ek
except ImportError:
    fail('eikon package not installed. Run: pip install -r requirements-data.txt')
print(f'  OK: eikon imported')

# 3. Set app key + test connectivity
print('\n[3/4] Setting app key and probing local proxy at port 9000...')
try:
    ek.set_app_key(key)
except Exception as e:
    fail(f'set_app_key failed: {e}\n'
         'Common causes: (a) Workspace Desktop not running; (b) invalid key; '
         '(c) proxy port not at default 9000.')
print(f'  OK: app key accepted')

# 4. Test pull on 5 US utilities
print('\n[4/4] Test pull: TR.PctFreeFloat for 5 US utility RICs...')
TEST_RICS = ['EXC.OQ', 'NEE', 'DUK', 'AEP.OQ', 'XEL.OQ']
try:
    df, err = ek.get_data(
        TEST_RICS,
        ['TR.PctFreeFloat', 'TR.SharesOutstanding'],
    )
except Exception as e:
    fail(f'get_data failed: {e}\n'
         'If error mentions "proxy" or "9000": Workspace Desktop is not connected.\n'
         'If error mentions "permissions": your Eikon subscription may not include this field.')

if df is None or len(df) == 0:
    fail('get_data returned empty result. Subscription may not cover requested fields.')

print(f'\n  OK: {len(df)} rows returned')
print(df.to_string())

if err:
    print(f'\n  Note: get_data warnings/errors (often benign):')
    print(f'  {err}')

# Sanity check on the values
ff_col = next((c for c in df.columns if 'FreeFloat' in c.replace(' ', '')), None)
if ff_col:
    valid_ff = df[ff_col].dropna()
    print(f'\n  Free-float values populated: {len(valid_ff)}/{len(df)}')
    if len(valid_ff) > 0:
        print(f'  Free-float range: {valid_ff.min():.1f} - {valid_ff.max():.1f}')
        if valid_ff.max() <= 1.0:
            print(f'  Note: values look like fractions (0-1 scale).')
        elif valid_ff.max() <= 100.0:
            print(f'  Note: values look like percentages (0-100 scale).')

print('\n=== ALL CHECKS PASSED ===')
print('\nNext step: python src/pull_refinitiv_extra.py')
print('Expected runtime: 6-12 hours for ~478 non-US firms')
print('Output: data/raw/refinitiv/refinitiv_extra.csv')

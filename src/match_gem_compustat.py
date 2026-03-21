"""Match GEM parent entities to Compustat Global gvkeys via name normalization + token matching."""
import csv, re, os
from collections import defaultdict
from _paths import raw_path, derived_path

def normalize(name):
    """Strip corporate suffixes and normalize for matching."""
    name = name.lower().strip()
    replacements = {
        '\xe9': 'e', '\xe8': 'e', '\xf1': 'n', '\xfc': 'u', '\xf6': 'o',
        '\xe4': 'a', '\xe1': 'a', '\xe0': 'a', '\xed': 'i', '\xf3': 'o',
        '\xfa': 'u', '\xe7': 'c', '\xdf': 'ss',
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    suffixes = [
        r'\bpjsc\b', r'\bjsc\b', r'\bpsc\b', r'\bojsc\b',
        r'\bcorp\b', r'\bcorporation\b', r'\binc\b', r'\bltd\b', r'\blimited\b',
        r'\bplc\b', r'\bllc\b', r'\bnv\b', r'\bbv\b',
        r'\bsa\b', r'\bspa\b', r'\bse\b', r'\bag\b', r'\baktiengesellschaft\b',
        r'\bco\b', r'\bcompany\b', r'\bgroup\b',
        r'\bholdings?\b', r'\bsoc\b', r'\bepe\b',
        r'\bpt\b', r'\btbk\b', r'\bpersero\b',
        r'\bpublic\b', r'\bues\b',
    ]
    for s in suffixes:
        name = re.sub(s, '', name)
    name = re.sub(r'[^\w\s]', '', name)
    name = re.sub(r'\s+', ' ', name).strip()
    return name

def tokenize(norm_name):
    """Split normalized name into token set."""
    return set(norm_name.split())

# Load Compustat company names (most recent year per gvkey) from both Global and NA
compustat = {}
for src_file in ['compustat_global_utilities.csv', 'compustat_na_utilities.csv']:
    fpath = raw_path('compustat', src_file)
    if not os.path.exists(fpath):
        continue
    with open(fpath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gvkey = row['gvkey']
            fyear = row['fyear']
            if gvkey not in compustat or fyear > compustat[gvkey]['fyear']:
                compustat[gvkey] = {
                    'conm': row['conm'],
                    'conml': row.get('conml', row['conm']),
                    'fic': row['fic'],
                    'sic': row['sic'],
                    'fyear': fyear,
                    'isin': row.get('isin', ''),
                    'source': src_file,
                }

# Build normalized name -> gvkey lookup + token index
norm_to_gvkey = defaultdict(list)
token_index = {}  # gvkey -> (norm_name, token_set)
for gvkey, info in compustat.items():
    for name_field in ['conm', 'conml']:
        n = normalize(info[name_field])
        if n:
            norm_to_gvkey[n].append(gvkey)
            token_index[gvkey] = (n, tokenize(n))

print(f'Compustat firms: {len(compustat)}')
print(f'Normalized name entries: {len(norm_to_gvkey)}')

# Load GEM parents
gem_parents = []
with open(derived_path('mappings', 'gem_parents_parsed.csv'), 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gem_parents.append(row)

print(f'GEM parent entities: {len(gem_parents)}')

def find_match(gem_norm):
    """Try exact match, then token-subset match."""
    # 1. Exact normalized match
    if gem_norm in norm_to_gvkey:
        return norm_to_gvkey[gem_norm], 'exact_norm'

    gem_tokens = tokenize(gem_norm)
    if len(gem_tokens) < 2:
        return [], None

    # 2. Token subset: all GEM tokens appear in Compustat name (or vice versa)
    candidates = []
    for gvkey, (comp_norm, comp_tokens) in token_index.items():
        if len(comp_tokens) < 2:
            continue
        # GEM tokens are subset of Compustat tokens
        if gem_tokens <= comp_tokens:
            candidates.append((gvkey, 'gem_subset'))
        # Compustat tokens are subset of GEM tokens
        elif comp_tokens <= gem_tokens:
            candidates.append((gvkey, 'comp_subset'))
        # High overlap (Jaccard > 0.6) with at least 2 shared tokens
        else:
            shared = gem_tokens & comp_tokens
            if len(shared) >= 2:
                jaccard = len(shared) / len(gem_tokens | comp_tokens)
                if jaccard >= 0.5:
                    candidates.append((gvkey, f'jaccard_{jaccard:.2f}'))

    if candidates:
        return [c[0] for c in candidates], candidates[0][1]
    return [], None

# Match
matches = []
unmatched_large = []

for gp in gem_parents:
    name = gp['parent_name']
    total_mw = float(gp['total_mw'])
    norm = normalize(name)

    matched_gvkeys, match_type = find_match(norm)

    if matched_gvkeys:
        matched_gvkeys = list(set(matched_gvkeys))
        for gvkey in matched_gvkeys:
            info = compustat[gvkey]
            matches.append({
                'gem_parent': name,
                'gem_norm': norm,
                'gvkey': gvkey,
                'conm': info['conm'],
                'fic': info['fic'],
                'isin': info['isin'],
                'match_type': match_type,
                'coal_mw': gp['coal_mw'],
                'gas_mw': gp['gas_mw'],
                'solar_mw': gp['solar_mw'],
                'wind_mw': gp['wind_mw'],
                'fossil_mw': gp['fossil_mw'],
                'total_mw': gp['total_mw'],
                'alpha': gp['alpha'],
            })
    elif total_mw >= 500:
        unmatched_large.append((name, norm, total_mw, gp['fossil_mw'], gp['alpha']))

# Summary
matched_parents = set(m['gem_parent'] for m in matches)
matched_gvkeys_set = set(m['gvkey'] for m in matches)
matched_mw = sum(float(m['total_mw']) for m in matches)
total_mw_all = sum(float(gp['total_mw']) for gp in gem_parents)

# Deduplicate MW (avoid double-counting same GEM parent matched to multiple gvkeys)
mw_by_parent = {}
for m in matches:
    if m['gem_parent'] not in mw_by_parent:
        mw_by_parent[m['gem_parent']] = float(m['total_mw'])
dedup_mw = sum(mw_by_parent.values())

print(f'\n=== MATCHING RESULTS ===')
print(f'Matched GEM parents: {len(matched_parents)} / {len(gem_parents)}')
print(f'Matched Compustat gvkeys: {len(matched_gvkeys_set)} / {len(compustat)}')
print(f'Matched MW (dedup): {dedup_mw:,.0f} / {total_mw_all:,.0f} ({100*dedup_mw/total_mw_all:.1f}%)')

# Match type breakdown
from collections import Counter
type_counts = Counter(m['match_type'] for m in matches)
print(f'\nMatch types: {dict(type_counts)}')

# Save matches
outpath = derived_path('mappings', 'gem_compustat_matches.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['gem_parent', 'gem_norm', 'gvkey', 'conm', 'fic', 'isin', 'match_type',
                                       'coal_mw', 'gas_mw', 'solar_mw', 'wind_mw', 'fossil_mw', 'total_mw', 'alpha'])
    w.writeheader()
    w.writerows(matches)
print(f'Saved matches to {outpath}')

# Show top unmatched
print(f'\nTop 30 UNMATCHED (>500 MW):')
unmatched_large.sort(key=lambda x: x[2], reverse=True)
fmt = '{:<50} {:>10} {:>10} {:>6}'
print(fmt.format('GEM Parent', 'Total MW', 'Fossil MW', 'Alpha'))
for name, norm, tmw, fmw, alpha in unmatched_large[:30]:
    print(fmt.format(name[:50], f'{tmw:,.0f}', fmw, alpha))

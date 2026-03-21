"""Summarize Strategy 2 spatial regression logs into CSV/Markdown tables."""
import csv
import os
import re
import sys

from _paths import results_path

DEFAULT_LOG = results_path('logs', 'strategy2_spatial_regression_run.log')
OUT_CSV = results_path('summaries', 'strategy2_spatial_regression_summary.csv')
OUT_MD = results_path('summaries', 'strategy2_spatial_regression_summary.md')


def read_lines(path):
    for enc in ('utf-8', 'utf-16'):
        try:
            with open(path, 'r', encoding=enc) as f:
                return f.read().splitlines()
        except UnicodeDecodeError:
            continue
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        return f.read().splitlines()


def parse_log(lines):
    rows = []
    scope = None
    ret_model = None
    block = None
    reg = None
    stats = {'n': None, 'clusters': None, 'r2': None, 'adj_r2': None}

    def parse_stats(line):
        out = {'n': None, 'clusters': None, 'r2': None, 'adj_r2': None}
        m = re.search(r'N=(\d+)', line)
        if m:
            out['n'] = int(m.group(1))
        m = re.search(r'clusters=([0-9, ]+)', line)
        if m:
            out['clusters'] = m.group(1).strip()
        m = re.search(r'R2=([0-9\.\-]+)', line)
        if m:
            out['r2'] = float(m.group(1))
        m = re.search(r'Adj R2=([0-9\.\-]+)', line)
        if m:
            out['adj_r2'] = float(m.group(1))
        return out

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith('EVENT SCOPE:'):
            scope = line.split(':', 1)[1].strip()
            continue
        if line.startswith('#  RETURN MODEL:'):
            ret_model = line.split(':', 1)[1].strip()
            continue
        block_headers = (
            'STRATEGY 2 SPATIAL REGRESSION',
            'CONTINUOUS CONDSCORE',
            'RANK-BASED CONDSCORE',
            'KERNEL CONDSCORE',
            'INTERACTION Z-SCORED',
            'EXTENDED SPEC',
            'ALPHA TRAJECTORY',
            'CHANNEL DECOMPOSITION',
        )
        if any(line.startswith(h) for h in block_headers):
            block = line
            continue

        # Regression label
        if line.startswith('DAILY') or line.startswith('MONTHLY'):
            # Skip the obs count lines
            if 'regression obs' in line.lower():
                continue
            reg = line
            continue

        if line.startswith('N='):
            stats = parse_stats(line)
            continue

        if line.startswith('Variable') or line.startswith('-----'):
            continue

        # Variable rows (format: var beta se t sig)
        parts = line.split()
        if len(parts) >= 4 and re.match(r'^[A-Za-z_]', parts[0]):
            var = parts[0]
            try:
                beta = float(parts[1])
                se = float(parts[2])
                t = float(parts[3])
            except ValueError:
                continue
            sig = parts[4] if len(parts) > 4 else ''
            rows.append({
                'scope': scope,
                'return_model': ret_model,
                'block': block,
                'regression': reg,
                'variable': var,
                'beta': beta,
                'se': se,
                't': t,
                'sig': sig,
                'n': stats.get('n'),
                'clusters': stats.get('clusters'),
                'r2': stats.get('r2'),
                'adj_r2': stats.get('adj_r2'),
            })

    return rows


def write_csv(rows, path):
    if not rows:
        return
    fields = ['scope', 'return_model', 'block', 'regression', 'variable',
              'beta', 'se', 't', 'sig', 'n', 'clusters', 'r2', 'adj_r2']
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_md(rows, path):
    if not rows:
        return
    header = '| Scope | Return | Block | Regression | Variable | Beta | t | N | R2 |\n'
    sep = '|---|---|---|---|---|---|---|---|---|\n'
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r['scope']} | {r['return_model']} | {r['block']} | {r['regression']} | "
            f"{r['variable']} | {r['beta']:.6f} | {r['t']:.2f} | {r['n']} | "
            f"{r['r2'] if r['r2'] is not None else ''} |"
        )
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


def main():
    log_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_LOG
    if not os.path.exists(log_path):
        raise SystemExit(f'Log not found: {log_path}')
    lines = read_lines(log_path)
    rows = parse_log(lines)
    write_csv(rows, OUT_CSV)
    write_md(rows, OUT_MD)
    print(f'Wrote {len(rows)} rows to {OUT_CSV}')
    print(f'Wrote {len(rows)} rows to {OUT_MD}')


if __name__ == '__main__':
    main()

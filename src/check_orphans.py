"""Check that all strategy*.py and build*.py in src/ are in run_all.py DAG.

Allowlist: _paths.py, run_all.py, romano_wolf_bootstrap.jl, pull_*.py, 
update_*.py, check_orphans.py, strategy2_referee_compute.py (called by run_all
via referee_tables dependency), strategy2_joint_tests.py (Julia fallback).

Usage: python src/check_orphans.py
Exit code 0 = clean, 1 = orphans found.
"""
import os, sys, re

SRC = os.path.dirname(os.path.abspath(__file__))
ALLOWLIST = {
    '_paths.py', 'run_all.py', 'romano_wolf_bootstrap.jl',
    'check_orphans.py',
}
ALLOW_PREFIXES = ('pull_', 'update_')

# Parse DAG from run_all.py
dag = set()
with open(os.path.join(SRC, 'run_all.py')) as f:
    for m in re.finditer(r"'([^']+\.(?:py|jl))'", f.read()):
        dag.add(m.group(1))

# List all scripts on disk
on_disk = set()
for fname in os.listdir(SRC):
    if fname.endswith(('.py', '.jl')):
        on_disk.add(fname)

# Find orphans
orphans = []
for f in sorted(on_disk - dag - ALLOWLIST):
    if any(f.startswith(p) for p in ALLOW_PREFIXES):
        continue
    orphans.append(f)

if orphans:
    print(f'ORPHAN SCRIPTS (not in DAG, not in allowlist):')
    for o in orphans:
        print(f'  {o}')
    print(f'\nEither add to run_all.py, add to ALLOWLIST in check_orphans.py, or delete.')
    sys.exit(1)
else:
    print('All scripts accounted for.')
    sys.exit(0)

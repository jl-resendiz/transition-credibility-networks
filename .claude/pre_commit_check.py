"""Hook: run check_orphans.py before git commit commands."""
import json, sys, subprocess, os

data = json.load(sys.stdin)
cmd = data.get('tool_input', {}).get('command', '')

if 'git commit' not in cmd:
    sys.exit(0)

src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src')
result = subprocess.run([sys.executable, os.path.join(src, 'check_orphans.py')])
if result.returncode != 0:
    print(json.dumps({
        "decision": "block",
        "reason": "Orphan scripts found in src/. Run python src/check_orphans.py for details."
    }))
sys.exit(result.returncode)

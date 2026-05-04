"""Load credentials from `.env` at the repo root into `os.environ`.

The pipeline core is stdlib-only by design, so this loader avoids
`python-dotenv` and parses `.env` itself. Lines starting with `#` are
comments; blank lines are ignored; existing environment variables are
NEVER overwritten (so command-line overrides take precedence).

Pull scripts (e.g. `pull_eikon_returns.py`, `pull_wrds_13f.py`) call
`load_env()` at startup and then read credentials with `os.getenv()`.

Usage
-----
    from _credentials import load_env, require
    load_env()                            # populate os.environ
    user = require("WRDS_USERNAME")       # raises if missing
    pwd  = require("WRDS_PASSWORD")
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
ENV_PATH = os.path.join(ROOT, '.env')


def load_env(path=ENV_PATH):
    """Read KEY=VALUE pairs from `.env` and set them in `os.environ`.

    No-op if the file does not exist (so the pipeline can run without
    credentials when re-running on previously-pulled derived data).
    """
    if not os.path.exists(path):
        return False
    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = value
    return True


def require(key, hint=None):
    """Return env var `key` or exit with an actionable message."""
    if key not in os.environ:
        load_env()
    val = os.environ.get(key, '').strip()
    if not val:
        msg = f'ERROR: environment variable {key} is not set.\n'
        msg += f'Set it in .env at the repo root (see .env.example).'
        if hint:
            msg += f'\n{hint}'
        print(msg, file=sys.stderr)
        raise SystemExit(2)
    return val


# Auto-load on import so existing pull scripts that simply call
# `os.getenv("EIKON_APP_KEY")` work without modification once they
# `import _credentials`.
load_env()

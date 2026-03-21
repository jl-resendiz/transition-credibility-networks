"""Standard execution pipeline for econometric diagnostics and models.

This script enforces a single, explicit standard:
  - Preferred transformations per model
  - Mandatory metrics output for every model
  - Deterministic execution order

Usage:
  python run_standard_pipeline.py            # full run
  python run_standard_pipeline.py --light    # light run for Strategy 2 spatial
"""
import os
import sys
import subprocess


ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, 'scripts')


def run(cmd, env=None):
    print(f'\n>>> {cmd}')
    res = subprocess.run(cmd, shell=True, env=env)
    if res.returncode != 0:
        raise SystemExit(res.returncode)


def main():
    light = '--light' in sys.argv

    base_env = os.environ.copy()
    base_env['WRITE_METRICS'] = '1'

    # Strategy 1: R² split diagnostics (multiple transforms)
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy1_r2_transform_experiments.py")}"', env=base_env)

    # Strategy 1 panel: two-way FE + density split
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy1_panel_regression.py")}"', env=base_env)

    # Strategy 2 spatial (preferred transform: base)
    env_s2 = base_env.copy()
    env_s2['TRANSFORM_SET'] = 'base'
    if light:
        env_s2['RUN_LIGHT'] = '1'
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy2_spatial_regression.py")}"', env=env_s2)

    # Strategy 2 panel DiD (preferred transform: base)
    env_s2p = base_env.copy()
    env_s2p['TRANSFORM_SET'] = 'base'
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy2_panel_did.py")}"', env=env_s2p)

    # Strategy 3 policy shocks (preferred transform: base, return model vwretd)
    env_s3 = base_env.copy()
    env_s3['TRANSFORM_SET'] = 'base'
    env_s3['RET_MODEL'] = 'vwretd'
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy3_policy_shocks.py")}"', env=env_s3)

    # Strategy 4 quantile regression
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy4_quantile_regression.py")}"', env=base_env)

    # Strategy 5 ESG forward delivery
    run(f'"{sys.executable}" "{os.path.join(SCRIPTS, "strategy5_esg_forward_delivery.py")}"', env=base_env)

    print('\nStandard pipeline complete.')


if __name__ == '__main__':
    main()

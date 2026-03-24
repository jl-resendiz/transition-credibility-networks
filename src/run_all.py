"""Pipeline orchestrator: build derived data then run all analysis scripts.

Usage:
  python src/run_all.py              # full pipeline (build + analysis)
  python src/run_all.py --analysis   # analysis only (skip build steps)

Build scripts run in dependency order:
  parse_gem -> match_gem_compustat -> build_* / compute_returns -> build_events

Analysis scripts are independent and run in any order.
"""
import os
import sys
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(ROOT, 'src')


def run(script_name):
    path = os.path.join(SRC, script_name)
    if not os.path.exists(path):
        print(f'  SKIP (not found): {script_name}')
        return
    print(f'\n{"=" * 60}')
    print(f'>>> {script_name}')
    print('=' * 60)
    result = subprocess.run(
        [sys.executable, path],
        cwd=SRC,
    )
    if result.returncode != 0:
        print(f'\n*** FAILED: {script_name} (exit code {result.returncode}) ***')
        raise SystemExit(result.returncode)


def main():
    analysis_only = '--analysis' in sys.argv

    if not analysis_only:
        # ── Stage 1: Parse raw data ──
        run('parse_gem.py')

        # ── Stage 2: Match GEM to Compustat ──
        run('match_gem_compustat.py')

        # ── Stage 3: Build derived datasets ──
        run('build_fundamentals.py')
        run('build_ets_matrix.py')
        run('build_weight_matrix.py')
        run('build_fuel_matrix.py')
        run('build_time_varying_alpha.py')
        run('compute_returns.py')

        # ── Stage 4: Build events (depends on Stage 3) ──
        run('build_retirement_events.py')
        run('build_coal_phaseout_events.py')
        run('build_eia860_announcement_events.py')

        # ── Stage 5: Summary statistics ──
        run('summary_statistics.py')

    # ── Stage 6: Analysis scripts ──
    analysis_scripts = [
        # T1/T8: Fuel channel dominance and difference tests
        'strategy2_robust_inference.py',
        'strategy2_difference_test_summary.py',
        'strategy2_bandwidth_fmb.py',
        'strategy2_joint_tests.py',
        'strategy2_firm_level_test.py',

        # T2: Geography null result
        'strategy2_event_specific_geo.py',
        'strategy2_geo_diversification.py',

        # T3/T4: ETS and carbon pricing
        'strategy2_credibility_interaction.py',
        'strategy2_esg_ets_fmb.py',

        # T5: ESG complementarity
        'strategy2_esg_horse_race.py',

        # T6: Cascading revelation
        'strategy2_learning_alternatives.py',

        # T7: Spatial Transition Score
        'strategy2_spatial_score.py',

        # T9: Phase-out wild bootstrap
        'strategy3_phaseout_wild_bootstrap.py',

        # T10: Multiple testing correction
        'strategy2_romano_wolf.py',

        # Shift-share (Bartik) causal diagnostics
        'strategy2_bartik_shiftshare.py',

        # Appendix robustness tables
        'strategy2_referee_tables.py',
    ]

    for s in analysis_scripts:
        run(s)

    print(f'\n{"=" * 60}')
    print('Pipeline complete.')
    print('=' * 60)


if __name__ == '__main__':
    main()

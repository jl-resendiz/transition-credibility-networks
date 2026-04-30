"""Pipeline orchestrator: build derived data then run all analysis scripts.

Usage:
  python src/run_all.py              # full pipeline (build + analysis)
  python src/run_all.py --analysis   # analysis only (skip build steps)

Build scripts run in dependency order:
  parse_gem -> match_gem_compustat -> build_* / compute_returns -> build_events

Analysis scripts are independent and run in any order.

Julia (.jl) scripts are used for computationally intensive bootstrap
procedures where they provide >10x speedup over Python. If Julia is
not installed, the pipeline falls back to equivalent Python scripts.

R (.R) scripts are used for spatial-econometric standard errors that
have no stdlib-Python equivalent (Conley SEs via fixest::conley()). If
Rscript is not on PATH, those steps are skipped with a warning.
"""
import os
import sys
import shutil
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(ROOT, 'src')

# Locate Julia binary (optional; pipeline falls back to Python if absent)
JULIA = shutil.which('julia')
if not JULIA:
    _julia_winapp = os.path.expanduser(
        r'~\AppData\Local\Microsoft\WindowsApps\julia.exe')
    if os.path.exists(_julia_winapp):
        JULIA = _julia_winapp

# Locate Rscript binary (optional; R-only steps skipped if absent)
RSCRIPT = shutil.which('Rscript')

# Fallback map: Julia script -> Python equivalent
_FALLBACK = {
    'joint_tests.jl': 'joint_tests.py',
}


def run(script_name):
    path = os.path.join(SRC, script_name)
    if not os.path.exists(path):
        print(f'  SKIP (not found): {script_name}')
        return
    print(f'\n{"=" * 60}', flush=True)
    print(f'>>> {script_name}', flush=True)
    print('=' * 60, flush=True)
    if script_name.endswith('.jl'):
        if JULIA:
            cmd = [JULIA, path]
        elif script_name in _FALLBACK:
            fb = _FALLBACK[script_name]
            print(f'  Julia not found, falling back to {fb}', flush=True)
            cmd = [sys.executable, os.path.join(SRC, fb)]
        else:
            print(f'  SKIP (Julia not found, no fallback): {script_name}')
            return
    elif script_name.endswith('.R'):
        if RSCRIPT:
            cmd = [RSCRIPT, path]
        else:
            print(f'  SKIP (Rscript not found): {script_name}')
            return
    else:
        cmd = [sys.executable, path]
    result = subprocess.run(cmd, cwd=SRC)
    if result.returncode != 0:
        print(f'\n*** FAILED: {script_name} (exit code {result.returncode}) ***')
        raise SystemExit(result.returncode)


def main():
    analysis_only = '--analysis' in sys.argv

    if not analysis_only:
        # ── Stage 0: Convert GEM xlsx to CSV (one-time, fast reads downstream) ──
        run('build_gem_csv.py')

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

    # ── Stage 6: Analysis scripts (9 scripts, consolidated) ──
    #
    # Each script addresses a specific identification threat:
    #   robust_inference       → Main results table (channel decomposition)
    #   joint_tests            → Joint significance + fuel ≠ geo difference
    #   esg_horse_race         → ESG complementarity (coverage argument)
    #   bartik_shiftshare      → Causal identification (GPS 2020 diagnostics)
    #   romano_wolf            → Multiple testing correction (3 hypotheses, primary window)
    #   geo_diversification    → Aggregation lemma test (theory validation)
    #   learning_alternatives  → Calendar-time learning (market efficiency)
    #   referee_tables         → LaTeX output + appendix robustness tables
    #
    analysis_scripts = [
        # Main results
        'robust_inference.py',
        'joint_tests.jl',        # ~42s Julia vs ~268s Python
        'esg_horse_race.py',

        # Identification and robustness
        'bartik_shiftshare.py',
        'romano_wolf.py',
        'geo_diversification.py',
        'learning_alternatives.py',

        # Phase 2: pre-trend transparency (event-time path + Honest DID)
        'event_time_path.py',
        'honest_did.py',

        # Phase 2: NW lag sensitivity (depends on event_level_betas.csv)
        'lag_sensitivity.py',

        # Phase 2: multi-factor abnormal returns (FF3 + sample utility industry)
        'multifactor_inference.py',
        'honest_did_mf.py',

        # Phase 2: bridge interaction (Section 2.5 augmented spec)
        'bridge_interaction.py',

        # Phase 2: announcement vs physical retirement robustness
        'announcement_robustness.py',

        # Phase 3: heterogeneity tests (US split, ESG joint, country robustness)
        'us_regulation_split.py',
        'esg_fm_joint.py',
        'country_robustness.py',

        # Post-units-fix: identification + anomaly-vs-risk demarcation
        'fisher_ri.py',
        'anomaly_vs_risk.py',

        # Conley spatial standard errors (skipped if Rscript not on PATH)
        'robustness_conley_se.R',

        # Output generation (compute exports JSON; tables reads it)
        'referee_compute.py',
        'referee_tables.py',
    ]

    for s in analysis_scripts:
        run(s)

    print(f'\n{"=" * 60}')
    print('Pipeline complete.')
    print('=' * 60)


if __name__ == '__main__':
    main()

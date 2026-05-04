"""Generate reviewer-requested robustness tables for Strategy 2.

Reads pre-computed regression results from results/json/referee_tables.json
(produced by referee_compute.py) and formats LaTeX tables.

Tasks:
1) Correlation matrix among w_geo, w_fuel, w_reg (3x3).
2) Channel decomposition with firm controls (Size, lambda, rho) for 3-month CARs.
3) Strong placebo: shuffle exposure networks (permute gvkeys) for 3-month CARs.
4) Specification progression (w_geo only, w_fuel only, both, full + controls).
5) Bandwidth sensitivity: w_geo with half-life 250km and 1000km.

Outputs LaTeX tables under results/tables/.
"""
import json
import os

from _paths import results_path


def stars(t):
    at = abs(t)
    if at >= 2.58:
        return '***'
    if at >= 1.96:
        return '**'
    if at >= 1.65:
        return '*'
    return ''


def format_coef(beta, se, t):
    return f"{beta:.3f}{stars(t)}"


def latex_var(name):
    mapping = {
        'w_geo': r'$w^{\mathrm{geo}}$',
        'w_fuel': r'$w^{\mathrm{fuel}}$',
        'w_reg': r'$w^{\mathrm{reg}}$',
        'same_sector': 'Same sector',
        'log_assets': 'Log assets',
        'lambda': r'Leverage ($\lambda$)',
        'rho': r'Return spread ($\rho$)',
    }
    return mapping.get(name, name.replace('_', r'\_'))


def latex_table_corr(corrs, out_path):
    lines = []
    lines.append('\\begin{table}[!htbp]')
    lines.append('\\centering')
    lines.append('\\caption{Correlation Matrix of Spatial Weights}')
    lines.append('\\label{tab:weight_corr}')
    lines.append('\\begin{tabular}{lccc}')
    lines.append('\\toprule')
    lines.append(' & $w^{\\mathrm{geo}}$ & $w^{\\mathrm{fuel}}$ & $w^{\\mathrm{reg}}$ \\\\')
    lines.append('\\midrule')
    for i, row_name in enumerate(['$w^{\\mathrm{geo}}$', '$w^{\\mathrm{fuel}}$', '$w^{\\mathrm{reg}}$']):
        row = [row_name]
        for j in range(3):
            val = corrs[i][j]
            row.append(f"{val:.3f}" if val is not None else '')
        lines.append(' & '.join(row) + ' \\\\')
    lines.append('\\bottomrule')
    lines.append('\\end{tabular}')
    lines.append('\\end{table}')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


def latex_table_spec(cols, out_path, caption, label, col_labels=None, extra_rows=None):
    """Write a specification-style LaTeX table.

    Each element of *cols* is a dict with:
        'n'    : observation count
        'vars' : list of variable names included in this column
        'coefficients' : {var: {'beta', 'se', 't'}}
    """
    lines = []
    lines.append('\\begin{table}[!htbp]')
    lines.append('\\centering')
    lines.append(f'\\caption{{{caption}}}')
    lines.append(f'\\label{{{label}}}')
    lines.append('\\begin{tabular}{l' + 'c' * len(cols) + '}')
    lines.append('\\toprule')
    if col_labels and len(col_labels) == len(cols):
        header = [''] + col_labels
    else:
        header = [''] + [f'({i})' for i in range(1, len(cols) + 1)]
    lines.append(' & '.join(header) + ' \\\\')
    lines.append('\\midrule')
    var_order = []
    for col in cols:
        for v in col['vars']:
            if v not in var_order:
                var_order.append(v)
    for v in var_order:
        row = [latex_var(v)]
        row_se = ['']
        for col in cols:
            if v in col['vars'] and v in col['coefficients']:
                c = col['coefficients'][v]
                beta = c['beta']
                se = c['se']
                t = c['t']
                row.append(format_coef(beta, se, t))
                row_se.append(f"({se:.3f})")
            else:
                row.append('')
                row_se.append('')
        lines.append(' & '.join(row) + ' \\\\')
        lines.append(' & '.join(row_se) + ' \\\\')
    lines.append('\\midrule')
    lines.append('N & ' + ' & '.join(str(col['n']) for col in cols) + ' \\\\')
    if extra_rows:
        for name, vals in extra_rows:
            lines.append(f'{name} & ' + ' & '.join(vals) + ' \\\\')
    lines.append('\\bottomrule')
    lines.append('\\end{tabular}')
    lines.append('\\end{table}')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


def main():
    json_path = results_path('json', 'referee_tables.json')
    if not os.path.exists(json_path):
        raise RuntimeError(
            f'Pre-computed results not found: {json_path}\n'
            'Run referee_compute.py first.'
        )
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # ── Table 1: Correlation matrix ──
    latex_table_corr(data['correlations'],
                     results_path('tables', 'table_weight_correlations.tex'))

    # ── Table 2: Channel decomposition with controls ──
    cc = data['channel_controls']
    latex_table_spec(
        [cc],
        results_path('tables', 'table_channel_controls.tex'),
        'Channel Decomposition with Firm Controls (3-Month CAR)',
        'tab:channel_controls',
        col_labels=['Controls'],
    )

    # ── Table 3: VIF ──
    vifs = data['vifs']
    vif_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets', 'lambda', 'rho']
    vif_lines = []
    vif_lines.append('\\begin{table}[!htbp]')
    vif_lines.append('\\centering')
    vif_lines.append('\\caption{Variance Inflation Factors (3-Month CAR Controls)}')
    vif_lines.append('\\label{tab:vif_controls}')
    vif_lines.append('\\begin{tabular}{lc}')
    vif_lines.append('\\toprule')
    vif_lines.append('Variable & VIF \\\\')
    vif_lines.append('\\midrule')
    for v in vif_vars:
        val = vifs.get(v)
        vif_lines.append(f'{latex_var(v)} & {val:.2f} \\\\')
    vif_lines.append('\\bottomrule')
    vif_lines.append('\\end{tabular}')
    vif_lines.append('\\end{table}')
    with open(results_path('tables', 'table_vif_controls.tex'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(vif_lines))

    # ── Table 4: Placebo shuffle ──
    pl = data['placebo']
    placebo_lines = []
    placebo_lines.append('\\begin{table}[!htbp]')
    placebo_lines.append('\\centering')
    placebo_lines.append('\\caption{Strong Placebo: Shuffled Exposure Networks (3-Month CAR)}')
    placebo_lines.append('\\label{tab:placebo_shuffle}')
    placebo_lines.append('\\begin{tabular}{lcc}')
    placebo_lines.append('\\toprule')
    placebo_lines.append(' & Mean coefficient & SD across permutations \\\\')
    placebo_lines.append('\\midrule')
    placebo_lines.append(f'$w^{{\\mathrm{{geo}}}}$ & {pl["w_geo"]["mean"]:.3f} & ({pl["w_geo"]["sd"]:.3f}) \\\\')
    placebo_lines.append(f'$w^{{\\mathrm{{fuel}}}}$ & {pl["w_fuel"]["mean"]:.3f} & ({pl["w_fuel"]["sd"]:.3f}) \\\\')
    placebo_lines.append(f'$w^{{\\mathrm{{reg}}}}$ & {pl["w_reg"]["mean"]:.3f} & ({pl["w_reg"]["sd"]:.3f}) \\\\')
    placebo_lines.append('\\bottomrule')
    placebo_lines.append('\\end{tabular}')
    placebo_lines.append('\\end{table}')
    out_placebo = results_path('tables', 'table_placebo_shuffle.tex')
    with open(out_placebo, 'w', encoding='utf-8') as f:
        f.write('\n'.join(placebo_lines))

    # ── Table 5: Specification progression ──
    sp = data['spec_progression']
    cols = [
        sp['geo_only'],
        sp['fuel_only'],
        sp['geo_fuel'],
        sp['full_controls'],
    ]
    latex_table_spec(
        cols,
        results_path('tables', 'table_spec_progression.tex'),
        'Specification Progression: 3-Month CARs',
        'tab:spec_progression',
        col_labels=['Geo only', 'Fuel only', 'Geo + Fuel', 'Full + controls'],
    )

    # ── Table 6: Controls sensitivity ──
    cs = data['controls_sensitivity']
    cols_ctrl = [
        cs['size_only'],
        cs['leverage_only'],
        cs['firm_fe'],
    ]
    latex_table_spec(
        cols_ctrl,
        results_path('tables', 'table_channel_controls_sensitivity.tex'),
        'Channel Decomposition: Controls Sensitivity (3-Month CARs)',
        'tab:channel_controls_sensitivity',
        col_labels=['Size only', 'Leverage only', 'Firm FE'],
        extra_rows=[('Firm FE', ['No', 'No', 'Yes'])],
    )

    # ── Table 7: Bandwidth sensitivity (5 half-lives) ──
    bw = data['bandwidth_sensitivity']
    bandwidths = [250, 500, 750, 1000, 1500]
    cols_bw = [bw[f'half_life_{h}km'] for h in bandwidths]
    latex_table_spec(
        cols_bw,
        results_path('tables', 'table_bandwidth_sensitivity.tex'),
        'Bandwidth Sensitivity (3-Month CARs)',
        'tab:bandwidth_sensitivity',
        col_labels=[f'Half-life {h} km' for h in bandwidths],
    )

    print('Reviewer tables written to results/tables/')


if __name__ == '__main__':
    main()

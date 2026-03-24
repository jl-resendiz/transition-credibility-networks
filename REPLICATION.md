# Replication Instructions

## Requirements

- **Python** 3.8+ (stdlib only; the only non-stdlib dependency is `openpyxl`)
- **R** 4.3+ with packages: `ggplot2`, `showtext`, `sysfonts`, `sf`, `rnaturalearth`, `rnaturalearthdata`
- **LaTeX** (MiKTeX or TeX Live) with `pdflatex` and `bibtex`

## Proprietary Data

The following data sources require institutional access and are not included in the repository:

| Source | Files | Access |
|---|---|---|
| Compustat Global/NA | `data/raw/compustat/` | Wharton Research Data Services (WRDS) |
| LSEG Eikon | `data/raw/eikon/` | Refinitiv Eikon terminal or API |
| CRSP | Returns merged via Compustat link | WRDS |

Public data sources (included):
- Global Energy Monitor plant tracker
- World Bank Carbon Pricing Dashboard
- Fama-French factor data

## Reproduce from Scratch

```bash
# Step 1: Full pipeline (build derived data + all analysis)
python src/run_all.py

# Step 2: Generate figures (requires R)
Rscript results/figures/generate_fig1.R
Rscript results/figures/generate_fig2.R
Rscript results/figures/generate_fig3.R

# Step 3: Copy figures to manuscript
cp results/figures/fig*.pdf manuscript/figures/

# Step 4: Compile manuscript (from manuscript/ directory)
cd manuscript
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Pipeline Structure

Build scripts (Stage 1-5) run in dependency order:
```
parse_gem → match_gem_compustat → build_fundamentals, build_fuel_matrix,
build_ets_matrix, build_weight_matrix → build_time_varying_alpha →
compute_returns → build_retirement_events, build_coal_phaseout_events,
build_eia860_announcement_events → summary_statistics
```

Analysis scripts (Stage 6) are independent and run in any order. See `src/run_all.py` for the complete list.

## Runtime

Full pipeline: approximately 15-20 minutes on a modern laptop.
Figure generation: approximately 1 minute.

## Output

- `results/metrics/`: Markdown summaries of all analysis results
- `results/tables/`: LaTeX tables for the appendix
- `results/summaries/`: CSV data exports
- `results/figures/`: PDF figures and R generation scripts
- `manuscript/main.tex`: The paper (references `results/tables/` and `manuscript/figures/`)

## Contact

Jose Luis Resendiz, Smith School of Enterprise and the Environment, University of Oxford.

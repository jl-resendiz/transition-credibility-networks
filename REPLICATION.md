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

# Step 1b: Conley spatial standard errors (requires R + fixest)
Rscript src/robustness_conley_se.R

# Step 2: Generate figures (requires R + ggplot2 + sf)
Rscript results/figures/generate_fig1.R
Rscript results/figures/generate_fig2.R
Rscript results/figures/generate_fig3.R
Rscript results/figures/generate_fig4.R

# Step 3: Copy figures to manuscript
cp results/figures/fig*.pdf manuscript/figures/

# Step 4: Compile manuscript (from manuscript/ directory)
cd manuscript
pdflatex when_coal_retires.tex && bibtex when_coal_retires && pdflatex when_coal_retires.tex && pdflatex when_coal_retires.tex
```

**Note on R dependencies:** The Conley spatial standard errors (Section 4.6) are computed in R using the `fixest` package, which provides built-in Conley SE support via `vcov = conley()`. This is the only analysis step that requires R. All other analysis runs in stdlib Python. R packages needed: `fixest`, `ggplot2`, `showtext`, `sysfonts`, `sf`, `rnaturalearth`, `rnaturalearthdata`.

## Pipeline Structure

Build scripts (Stage 1-5) run in dependency order:
```
parse_gem → match_gem_compustat → build_fundamentals, build_fuel_matrix,
build_ets_matrix, build_weight_matrix → build_time_varying_alpha →
compute_returns → build_retirement_events, build_coal_phaseout_events,
build_eia860_announcement_events → summary_statistics
```

Analysis scripts (Stage 6) are independent and run in any order. See `src/run_all.py` for the complete list. Recent additions:

- `event_time_path.py` — event-time path of the fuel coefficient $\hat\beta_\text{fuel}(\tau)$ for $\tau\in[-12,+6]$, exported to `results/summaries/event_time_betas.csv` (consumed by `generate_fig4.R`).
- `honest_did.py` — Rambachan-Roth (2023) sensitivity bound on parallel-trends violations using placebo CAR windows. Reports the breakdown $\bar M$ in `results/metrics/honest_did.md`.
- `lag_sensitivity.py` — Newey-West HAC lag sensitivity at lags {4, 8, 12, 18} on the FM time series of event-level coefficients.
- `multifactor_inference.py` — multi-factor abnormal returns (FF3 + sample-constructed utility industry excess return) with firm-specific betas estimated on a 24-month pre-event window. Reports the comparison vs single-factor headline in `results/metrics/multifactor_inference.md`.
- `honest_did_mf.py` — recalibration of Honest DID on the multi-factor CARs.
- `bridge_interaction.py` — empirical test of the augmented specification in Section 2.5 of the manuscript ($w^{\text{fuel}} \times (\alpha_i^{\text{pre-2014}} - \bar\alpha)$ interaction).
- `announcement_robustness.py` — comparison of announcement-default vs forced-physical event timing.
- `us_regulation_split.py` — US restructured-vs-regulated state split for the regulation hypothesis.
- `esg_fm_joint.py` — Fama-MacBeth joint test of ESG and fuel-mix similarity coefficients in the ESG-overlap subsample.
- `country_robustness.py` — non-US country event counts, leave-one-country-out diagnostic, developed-ex-US vs emerging-market split.

## Runtime

Full pipeline: approximately 15-20 minutes on a modern laptop.
Figure generation: approximately 1 minute.

## Output

- `results/metrics/`: Markdown summaries of all analysis results
- `results/tables/`: LaTeX tables for the appendix
- `results/summaries/`: CSV data exports
- `results/figures/`: PDF figures and R generation scripts
- `manuscript/when_coal_retires.tex`: The paper (references `results/tables/` and `manuscript/figures/`)

## Contact

Jose Luis Resendiz, Smith School of Enterprise and the Environment, University of Oxford.

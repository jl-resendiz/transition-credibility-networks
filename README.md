# When Coal Retires: The Propagation of Stranding Risk

**Jose Luis Resendiz**
Smith School of Enterprise and the Environment, University of Oxford

---

## Overview

This repository contains the replication package for the paper. Using plant-level data for 703 listed power utilities across 80 countries, I show that fuel-mix similarity is the dominant channel through which coal plant retirements transmit to stock prices. Geographic proximity, despite its theoretical appeal through local competition, does not transmit to equity valuations. The fuel channel is robust across all inference methods (Fama-MacBeth, event-clustered, two-way clustered, Conley spatial) and survives all multiple testing corrections (Bonferroni, max-t, Romano-Wolf). ESG scores dominate fuel-mix similarity where both exist but cover only 153 of 703 firms; fuel-mix similarity extends transition risk measurement to the remainder.

---

## Repository Structure

```
transition-credibility-networks/
├── manuscript/          Manuscript source (LaTeX + figures)
├── data/
│   ├── raw/             External input data (GEM, Compustat, CRSP, Eikon)
│   └── derived/         Intermediate datasets produced by the pipeline
├── src/                 All scripts (build + analysis)
├── results/
│   ├── metrics/         Regression output (one MD per analysis)
│   ├── tables/          LaTeX table fragments (imported by manuscript)
│   └── summaries/       Descriptive statistics
└── literature/          Literature review notes
```

---

## Data Sources

| Source | Contents | Access |
|---|---|---|
| Global Energy Monitor (GEM) | Plant-level capacity, fuel type, GPS coordinates | Public |
| Compustat Global/NA | Firm financials | Licensed |
| CRSP/Compustat Merged | US equity total returns (daily + monthly) | Licensed |
| LSEG Eikon (TR.TotalReturn) | Non-US equity total returns (monthly) | Licensed |
| LSEG Eikon (ESG) | ESG scores, governance, emissions | Licensed |
| World Bank Carbon Pricing Dashboard | Carbon prices by country-year | Public |
| EIA Form 860 | US coal retirement announcements | Public |

Monthly equity returns use total return series only: CRSP `trt1m` (US) and Eikon `TR.TotalReturn` (non-US). Compustat Global Security price returns are not used for monthly data because they exclude dividends (Ince & Porter 2006).

---

## Reproduction

**Requirements:** Python 3.8+ (standard library only). The only non-stdlib dependency is `openpyxl`.

```bash
python src/run_all.py              # full pipeline (build + analysis)
python src/run_all.py --analysis   # analysis only (skip build)
```

### Pipeline

```
Stage 1: parse_gem.py                    (GEM xlsx -> parsed CSV)
Stage 2: match_gem_compustat.py          (GEM + Compustat -> matched firms)
Stage 3: build_fundamentals.py           (financials panel)
         build_weight_matrix.py          (geographic weight matrix)
         build_fuel_matrix.py            (fuel-similarity weight matrix)
         build_ets_matrix.py             (regulatory weight matrix)
         build_time_varying_alpha.py     (fossil intensity panel)
         compute_returns.py              (CRSP + Eikon -> returns)
Stage 4: build_retirement_events.py      (coal retirement events)
         build_coal_phaseout_events.py   (phase-out policy events)
         build_eia860_announcement_events.py
Stage 5: summary_statistics.py
Stage 6: 9 analysis scripts (see src/run_all.py for the DAG)
```

---

## Citation

> Resendiz, J. L. (2026). When Coal Retires: The Propagation of Stranding Risk. Working paper, University of Oxford.

---

## Contact

Jose Luis Resendiz — Smith School of Enterprise and the Environment, University of Oxford

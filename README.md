# Pricing Transition Credibility on Spatial Networks

**Jose Luis Resendiz**  
Smith School of Enterprise and the Environment, University of Oxford  

---

## Overview

This repository contains the replication package for the paper *"Pricing Transition Credibility on Spatial Networks."* The paper develops a networked asset-pricing framework to show that climate transition risk is fundamentally a spatial network property, not an isolated firm-level attribute. Using a global sample of 389 geolocated power utilities, we document a channel split in how plant retirements transmit through networks: negative contagion for technologically similar peers (stranding risk) and positive competitive benefits for geographically proximate neighbors. We further reveal a credibility gap wherein voluntary retirements transmit as ambiguous signals (volatility without directional repricing), while binding phase-out laws force directional repricing proportional to fossil exposure.

---

## Repository Structure

```
transition-credibility-networks/
├── manuscript/          Manuscript source files (LaTeX + figures)
├── data/
│   ├── raw/             External input data (Compustat, GEM, Refinitiv, etc.)
│   └── derived/         Intermediate datasets produced by the pipeline
├── src/                 All analysis scripts
├── results/             Regression output tables, metrics, and summaries
└── literature/          Thematic literature review notes
```

---

## Data Sources

| Source | Contents | Access |
|---|---|---|
| Global Energy Monitor (GEM) | Plant-level capacity, fuel type, coordinates | Public |
| Compustat Global/NA | Firm financials, equity returns | Licensed |
| CRSP / Datastream | Daily and monthly equity returns | Licensed |
| Refinitiv | ESG scores, governance, Scope 1+2 emissions | Licensed |
| World Bank Carbon Pricing Dashboard | Carbon prices by country-year | Public |
| EIA Form 860 | US coal retirement announcements | Public |

Raw licensed data are not redistributed. See `data/raw/` for source documentation.

---

## Reproduction

**Requirements:** Python 3.8+ (standard library only — no pandas, numpy, or scipy).

### Step 1: Build all derived datasets and run all analyses

```bash
python src/run_pipeline.py
```

This runs all build and strategy scripts in dependency order. Outputs are written to `results/`.

### Step 2: Compile the manuscript

```bash
cd manuscript
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Requires a LaTeX distribution (MiKTeX or TeX Live).

---

## Citation

If you use this code or data, please cite:

> Resendiz, J. L. (2026). Pricing Transition Credibility on Spatial Networks. Working paper, University of Oxford.

---

## Contact

Jose Luis Resendiz — Smith School of Enterprise and the Environment, University of Oxford

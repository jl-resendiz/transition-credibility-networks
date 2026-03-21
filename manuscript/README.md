# Manuscript

This folder contains the LaTeX source for the paper.

## Files

| File | Description |
|---|---|
| `main.tex` | Main manuscript source |
| `references.bib` | BibTeX bibliography |
| `figures/` | Publication figures (PDF and PNG) |

## Figures

| File | Description |
|---|---|
| `fig1_global_sample` | Global sample of 389 power utilities |
| `fig2_network_structure` | Spatial network visualization |
| `fig2b_regional_networks` | Regional network decomposition |
| `fig3_retirement_events` | Coal retirement event time paths |
| `event_time_fuel_similarity` | Event-time response by fuel-similarity channel |

## Compilation

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

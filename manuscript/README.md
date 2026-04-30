# Manuscript

This folder contains the LaTeX source for the paper.

## Files

| File | Description |
|---|---|
| `when_coal_retires.tex` | Main manuscript source |
| `references.bib` | BibTeX bibliography |
| `figures/` | Publication figures (PDF) |

## Figures

| File | Description |
|---|---|
| `fig1_fuel_vs_geo.pdf` | Distribution of event-level coefficients (fuel vs geo) |
| `fig2_calendar_time.pdf` | Fuel coefficient by event-year tercile |
| `fig3_world_map.pdf` | Global map of retirements and sample firms |
| `fig4_event_time_path.pdf` | Event-time path of fuel-similarity coefficient |

## Compilation

```bash
pdflatex when_coal_retires.tex
bibtex when_coal_retires
pdflatex when_coal_retires.tex
pdflatex when_coal_retires.tex
```

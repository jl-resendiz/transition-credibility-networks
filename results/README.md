# Results

Outputs produced by `src/run_all.py`. Do not edit manually.

## Structure

| Folder | Contents |
|---|---|
| `metrics/` | Pipeline output (one `.md` per analysis script) |
| `tables/` | LaTeX table source files for the manuscript appendix |
| `json/` | Machine-readable results for `strategy2_referee_tables.py` |
| `summaries/` | Event-level betas and other CSV exports |
| `figures/` | R scripts and PDF figures |

## Regenerating

```bash
python src/run_all.py --analysis   # ~60 seconds
```

# Data

## Structure

```
data/
├── raw/          External inputs (not modified after download)
└── derived/      Intermediate datasets produced by src/ scripts
```

## Raw Data (`data/raw/`)

| Folder | Source | Contents |
|---|---|---|
| `compustat/` | Compustat Global/NA | Firm financials; daily and monthly equity returns |
| `crsp_compustat/` | CRSP-Compustat Merged | CCM daily and monthly returns (two-part splits) |
| `datastream/` | Refinitiv Datastream | Non-US equity returns |
| `events/` | GEM / EIA-860 | Coal retirement events with announcement dates |
| `factors/` | Ken French Data Library | Fama-French factors (daily and monthly) |
| `gem/` | Global Energy Monitor | Plant trackers: coal, solar, wind, oil/gas (Jan-Feb 2026) |
| `policy/` | World Bank Carbon Pricing Dashboard | Carbon prices by country-year; coal phase-out shock dates |
| `refinitiv/` | Refinitiv ESG | ESG scores, governance extras, Scope 1+2 emissions panel |
| `interconnectors/` | ENTSO-E / manual | Cross-border electricity interconnector data |

Licensed data (Compustat, CRSP, Refinitiv) are not redistributed. Contact the respective providers.

## Derived Data (`data/derived/`)

Derived datasets are produced by the build scripts in `src/` and should not be edited manually.

| Folder | Contents |
|---|---|
| `fundamentals/` | Fossil intensity panel (`firm_alpha_panel.csv`), firm financials |
| `mappings/` | GEM-Compustat matches, GvKey-RIC mapping, utility lists |
| `networks/` | Spatial weight matrices (fuel similarity, geographic, regulatory, ETS) |
| `events/` | Processed coal retirements, EIA-860 announcements, phase-out shocks |
| `returns/` | Market-adjusted daily and monthly returns |
| `trajectories/` | Time-varying fossil intensity trajectories |
| `validation/` | Alpha measurement cross-checks |

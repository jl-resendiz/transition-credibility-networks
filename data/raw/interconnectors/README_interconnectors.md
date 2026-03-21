# Cross-Border Electricity Interconnector Dataset

## Purpose
Event set for expanding Strategy 2 (spatial transmission) beyond coal retirements. Interconnector commissionings change spatial market structure exogenously to individual utility decisions, providing cleaner identification of spatial repricing.

## File
- `cross_border_interconnectors.csv` — project-level dataset of cross-border and inter-regional electricity interconnectors commissioned 2006-2025, plus planned/under-construction projects.

## Fields
| Field | Description |
|---|---|
| project_name | Common name of the interconnector |
| country_a, country_b | Countries (or regions) connected |
| region | Geographic region (Europe-North Sea, ASEAN, South Asia, etc.) |
| capacity_mw | Rated capacity in MW |
| technology | HVDC, HVAC, HVDC-B2B (back-to-back) |
| length_km | Cable/line length in km (some entries have voltage instead; needs cleaning) |
| commissioning_year | Year of commercial operation |
| commissioning_date | Exact date where available (YYYY-MM-DD) |
| status | Operational, Under construction, Planned |
| announcement_or_fid_date | Announcement or final investment decision date (to be filled) |
| source | Primary source(s) |

## Coverage
- **Operational:** 26 projects (2006-2025)
- **Under construction:** 4 projects
- **Planned:** 4 projects
- **Regions:** Europe (20), South Asia (6), ASEAN (2), East Asia (2), Oceania (1)

## Key Sources

### Structured Datasets (downloadable)
1. **Global Transmission Database (Zenodo)**
   - DOI: 10.5281/zenodo.10594108
   - URL: https://zenodo.org/records/10594108
   - Format: 4 CSV files (existing/planned, national/regional)
   - Coverage: Global, existing + planned capacities
   - Fields: MW capacity, voltage, distance; **includes expected commissioning years for planned projects**
   - By: TransitionZero, Simon Fraser University, CCG, Dartmouth

2. **ENTSO-E TYNDP 2024 Project Collection**
   - URL: https://tyndp2024.entsoe.eu/projects-map
   - Data portal: https://tyndp-data.netlify.app/maps-data/
   - Excel download (TYNDP 2020): https://www.entsoe.eu/Documents/TYNDP%20documents/TYNDP2020/201102_TYNDP2020_Portfolio_updated.xlsx
   - Coverage: 178 transmission projects + 33 storage projects (Europe)
   - Fields: Project sheets with commissioning year, CBA results

3. **ENTSO-E TYNDP 2022 Project Sheets**
   - URL: https://tyndp2022-project-platform.azurewebsites.net/projectsheets/transmission
   - Individual project sheets with commissioning targets

4. **Ember Europe Electricity Interconnection Data**
   - URL: https://ember-energy.org/data/europe-electricity-interconnection-data/
   - Data tool: https://ember-energy.org/data/europe-electricity-interconnection-data-tool/
   - Methodology: https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/europe_interconnection_data_tool/Ember%20Europe%20Electricity%20Interconnection%20Data%20Tool%20-%20Methodology.pdf

5. **ASEAN Power Grid Interconnections Project Profiles**
   - Publisher: ASEAN Centre for Energy (ACE)
   - URL: https://aseanenergy.org/publications/asean-power-grid-interconnections-project-profiles/

6. **India Ministry of Power — Interconnection with Neighbouring Countries**
   - URL: https://powermin.gov.in/en/content/interconnection-neighbouring-countries

### Reference Sources (project-level verification)
- Wikipedia: List of HVDC projects — https://en.wikipedia.org/wiki/List_of_HVDC_projects
- 4C Offshore Subsea Interconnectors Database — https://www.4coffshore.com/transmission/interconnectors.aspx (premium for full data)
- Wood Mackenzie Europe Power Cross-Border Interconnector Tracker (commercial)
- EU PCI/PMI Interactive Map — https://energy.ec.europa.eu/topics/infrastructure/projects-common-interest-and-projects-mutual-interest/key-cross-border-infrastructure-projects_en
- ACER PCI Selection — https://www.acer.europa.eu/electricity/infrastructure/projects-common-interest/pci-selection

## Data Quality Notes

1. **Exact commissioning dates** are available for approximately half the projects. For the remainder, only the commissioning year is known. Announcement/FID dates need to be filled from news archives and regulatory filings.

2. **Announcement dates are critical** for event studies. The commissioning date is when the interconnector begins operation, but markets may price the event at announcement, regulatory approval, or financial close. You will need to collect these separately from:
   - EU Official Journal (for PCI designation dates)
   - National regulatory authority decisions
   - Company press releases (for FID dates)
   - ENTSO-E project sheets (for regulatory milestone dates)

3. **The `length_km` field** has some entries showing voltage (e.g., "400kV") instead of length for HVAC lines. These need cleaning.

4. **Capacity figures** may differ between sources (rated vs. derated vs. commercial). The dataset uses rated capacity where available.

5. **Internal interconnectors** (SAPEI within Italy, Basslink within Australia, Hokkaido-Honshu within Japan, Attica-Crete within Greece) are included because they change inter-regional market structure even though they do not cross sovereign borders.

## How to Use for the Paper

### As Positive/Neutral Shocks (complement to coal retirement negative shocks)
- When a new interconnector is commissioned, previously isolated utilities become exposed to cross-border competition or gain access to cheaper imports.
- The spatial transmission prediction: utilities on both sides experience repricing depending on their fuel mix relative to the new competitive landscape.
- Sign depends on relative cost positions (observable from existing fuel mix data).

### Identification Advantage
- Interconnector approvals are regulatory decisions by grid operators/energy regulators.
- Timing driven by permitting, environmental review, and construction — largely exogenous to individual utility stock prices.
- Much cleaner than coal retirements (which are endogenous to market conditions).

### Next Steps
1. Download the Global Transmission Database from Zenodo to fill gaps in capacity and planned projects.
2. Download ENTSO-E TYNDP 2020/2022 Excel files for European project-level data with commissioning years.
3. Collect announcement/FID dates from press releases and regulatory filings for each operational project.
4. Match interconnector endpoints to utility service territories using the existing spatial weight matrix infrastructure.
5. Run event studies around commissioning (and announcement) dates using the same CAR methodology as Strategy 2.

# Coal Phase-Out Regulatory Shocks Dataset

## Purpose
Event set for a clean causal identification strategy: legally exogenous coal phase-out shocks with exact dates, suitable for event studies and panel DiD designs. Complements the existing EIA-860 plant-level retirement events with regulation-level shocks whose timing is determined by legislative/regulatory processes rather than firm decisions.

## File
- `coal_phaseout_shocks.csv` — 44 regulatory shock events across 4 tiers of exogeneity.

## Fields
| Field | Description |
|---|---|
| shock_id | Unique identifier |
| jurisdiction | Country, state, or institution |
| shock_name | Official name of law/regulation |
| event_date | Date markets first learned the final outcome (vote, signing, announcement) |
| legal_date | Date of entry into force or Official Journal publication (may differ from event_date) |
| instrument_type | Law, regulation, directive, decision, communique, etc. |
| binding | Whether the instrument carries legal force (yes/no) |
| scope | What entities are covered |
| coal_mandate | Specific coal-related requirement |
| affected_countries | ISO3 codes or descriptive scope |
| exogeneity_tier | 1 (cleanest) to 4 (weakest) — see below |
| source | Primary source reference |

## Exogeneity Tiers

### Tier 1 — Binding national/state coal laws with exact event dates (strongest)
Parliamentary votes or governor signatures on laws that directly mandate coal plant closures or prohibit coal generation. Timing determined by legislative process, not firm economics. **Recommended for the primary event study.**

Includes: Germany (Kohleausstiegsgesetz), Netherlands (coal ban), France (energy-climate law), Finland (coal ban), Greece (climate law), 10 US state laws (CA, WA, OR, NM, CO, NY, VA, IL, NC, MN), Indonesia (PR 112/2022), Vietnam (PDP8), Philippines (DOE moratorium), US EPA MATS.

### Tier 2 — EU directives, federal regulations, institutional policies
EU-level instruments and federal US regulations with identifiable adoption dates. Slightly weaker because: (a) trilogue negotiations may leak information before formal adoption, (b) compliance deadlines are years away, creating noise in short-window event studies.

Includes: IED, BAT conclusions, ETS Phase IV reform, Denmark Climate Act, Japan energy plans, Bangladesh cancellations, ADB energy policy.

### Tier 3 — Market mechanism reforms and indirect instruments
Instruments that affect coal economics indirectly through carbon pricing, taxonomy exclusion, or financing restrictions rather than direct mandates. Weaker identification because the link to specific plant closures is mediated by market prices.

Includes: EU ETS MSR, EU Taxonomy, CBAM, Denmark Energy Agreement, South Korea 10th Plan, World Bank coal policy, FERC MOPR.

### Tier 4 — Voluntary commitments and political signals
Non-binding communiques, coalition launches, and COP decisions. Useful as controls or placebo tests, not as primary event dates.

Includes: PPCA launch, G7 Carbis Bay, G20 Rome, COP26 Glasgow.

## Recommended Event Study Design

### Cleanest specification (N small, identification strong)
Use **Tier 1 binding national laws only**, focusing on the 2019 cluster:
- Finland: 2019-03-06
- New Mexico: 2019-03-22
- Colorado: 2019-05-30
- Washington CETA: 2019-05-07
- New York CLCPA: 2019-07-18
- France: 2019-09-26
- Netherlands: 2019-12-10

Plus:
- Oregon: 2016-03-11
- Germany: 2020-07-03
- Virginia: 2020-04-12
- Illinois: 2021-09-15
- Greece: 2022-05-26
- Indonesia: 2022-09-13
- Vietnam: 2023-05-15
- Philippines: 2020-10-27

This gives N ≈ 15 clean shocks with exact dates, binding force, and legislative timing exogenous to individual utility stock prices.

### Pre-trend concern
The 2019 US cluster (5 state laws in 7 months) creates potential cross-contamination. Consider:
1. Using only the first-mover (Finland 2019-03-06) and last-mover (Netherlands 2019-12-10) with exclusion windows between.
2. Treating the 2019 cluster as a single shock window (2019-03 to 2019-12).
3. Running leave-one-out diagnostics on each event.

### Which utilities to study
For each shock, the treatment group is utilities with coal exposure in the affected jurisdiction:
- **National laws** (DE, NL, FR, FI, GR, ID, VN, PH): utilities operating coal plants in that country.
- **US state laws**: utilities with coal plants in that state, or serving load in that state (for coal-by-wire bans like CA, OR).
- **EU regulations**: all EU utilities with coal capacity.

Control group: comparable utilities without coal exposure in the jurisdiction, or utilities in unaffected jurisdictions.

## Key Sources
1. Climate Change Laws of the World (Grantham Research Institute, LSE) — https://climate-laws.org/
2. EUR-Lex (EU Official Journal) — https://eur-lex.europa.eu/
3. US Federal Register — https://www.federalregister.gov/
4. State legislature websites (linked in individual entries)
5. IEA Policy Database — https://www.iea.org/policies/

## Data Quality Notes
1. **Event date vs. legal date**: For event studies, use `event_date` (when markets learned the outcome). For compliance analysis, use `legal_date`. They differ for EU instruments where publication follows adoption by weeks.
2. **US Virginia VCEA**: Signing date 2020-04-12 falls on a Sunday. Use Monday 2020-04-13 as day 0.
3. **Clean Power Plan**: Stayed by SCOTUS 2016-02-09, never took effect. Both announcement (2015-08-03) and stay (2016-02-09) are relevant events in opposite directions.
4. **FERC MOPR**: Directionally pro-coal (disadvantaged subsidized renewables). The December 2021 reversal is the anti-coal date.
5. **Bangladesh**: Partially reversed in 2025 (Matarbari revival), complicating the irreversibility assumption.
6. **2019 cluster**: Five US state coal laws signed March-December 2019. Treat carefully for cross-contamination.

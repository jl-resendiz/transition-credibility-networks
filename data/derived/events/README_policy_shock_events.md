# policy_shock_events.csv

## Contents

This file contains 20 national-level coal phase-out policy announcements
covering the period 2015 to 2022. Each row records the event date, country
ISO3 code, a short announcement description, the stated target year, the
type of policy instrument, and the precision of the date.

## Columns

| Column | Description |
|--------|-------------|
| `event_id` | Short identifier (country_YYYY_MM_DD format) |
| `event_date` | Date of the policy announcement (YYYY-MM-DD, blank if only month known) |
| `event_month` | Year-month string (YYYY-MM) |
| `country_iso3` | ISO 3166-1 alpha-3 country code(s), semicolon-separated for multi-country events |
| `announcement` | Brief description of the policy action |
| `target_year` | Year by which coal phase-out or net-zero is targeted |
| `type` | Instrument type (e.g. Legislation, Ministerial speech, International coalition) |
| `precision` | Date precision (Day-exact, Day-approx, Month-level) |
| `notes` | Additional context |

## Status: Legacy / Superseded

This file is a legacy precursor to `coal_phaseout_shocks_events.csv`, which
supersedes it in all current analysis scripts. The superseding file has a
richer structure, including:

- EU-level regulatory events (directives, taxonomy, ETS phases)
- US state-level legislation
- Binding vs. non-binding classification
- Exogeneity tier scoring (1 = binding national law, 2 = directive/plan, 3 = indirect)
- Matched gvkey lists and first-mover flags

No current analysis script in `src/` imports or reads `policy_shock_events.csv`.
All policy event study regressions use `coal_phaseout_shocks_events.csv` instead.

## Recommendation

Retain for reference and audit trail. The file documents the initial set of
country-level phase-out events that seeded the more comprehensive event
catalogue. It may be safely removed once the provenance of
`coal_phaseout_shocks_events.csv` is fully documented.

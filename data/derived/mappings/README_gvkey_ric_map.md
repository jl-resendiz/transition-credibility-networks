# gvkey_ric_map.csv

## Contents

This file maps Compustat Global/North America gvkey identifiers to Refinitiv
(formerly Thomson Reuters) RIC identifiers for the full sample of power
utilities used in the transition credibility network analysis. Each row
represents one firm.

## Columns

| Column | Description |
|--------|-------------|
| `gvkey` | Compustat permanent company identifier (6-digit, zero-padded) |
| `isin`  | International Securities Identification Number used as the primary join key during construction |
| `ric`   | Refinitiv Instrument Code, used to query Eikon/Datastream for ESG scores, governance data, and emissions panels |

## Construction

The mapping was constructed manually by cross-referencing Compustat Global
and North America company records with the Refinitiv Eikon symbology service.
The matching procedure was:

1. Extract ISINs from Compustat Security Daily tables for all gvkeys
   identified as electric utilities (SIC 4911-4941).
2. Look up each ISIN in the Refinitiv Eikon symbol-search API to retrieve
   the corresponding primary-listed RIC.
3. For cases where the ISIN lookup returned no result or an ambiguous match,
   verify by company name matching against the Refinitiv company directory.
4. Where a company has changed its primary listing (indicated by a suffix
   such as `^D24`, `^F23`, etc.), the RIC reflects the listing that was
   active during the sample period.

No automated script produces this file. It is a static, researcher-curated
input that must exist before running any of the Refinitiv pull scripts.

## Consumers

The following scripts read this file to retrieve Refinitiv data:

- `src/pull_refinitiv_esg.py`
- `src/pull_refinitiv_panel.py`
- `src/pull_refinitiv_extra.py`

## Last updated

The file contains 532 gvkey-RIC pairs. It was last updated during the
initial sample construction phase and should be reviewed if the firm sample
changes or if Refinitiv retires or reassigns any listed RICs.

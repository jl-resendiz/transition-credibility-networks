# Summary Statistics

## Sample Overview

- Total Compustat utility firms: 873
- Analysis sample (complete theta): 428 firms
- Countries (all): 81
- Countries (analysis sample): 69
- Total installed capacity (analysis sample): 2,983,226 MW = 3.0 TW

### Event-study sample (chain to manuscript abstract)

| Stage | Sample | N firms |
|---|---|---:|
| 1 | Firms with monthly returns coverage (manuscript: "703 listed utilities") | 703 |
| 2 | Returns AND at least one network layer (W_geo) | 350 |
| 3 | Returns AND all three network layers (manuscript: "565 have complete spatial weight data") | 183 |
| 4 | Returns AND ESG coverage (manuscript: "153 also have ESG environmental scores") | 169 |

Note: the "analysis sample (complete theta)" reported above is the firm universe with complete *fundamentals* (alpha, lambda, rho, kappa) used by Panels A and B. The event-study regressions in Section 4 use the monthly-returns sample (Stage 1) intersected with the network layers (Stage 3), which is the sample size cited in the abstract.

## Panel A: Analysis Sample (N = 428 firms, latest fiscal year)

| Variable                  |      N |       Mean |         SD |        Min |     Median |        Max |
|---|---|---|---|---|---|---|
| Leverage (lambda)         |    428 |      0.388 |      0.255 |      0.000 |      0.394 |      3.704 |
| Return spread (rho)       |    428 |      0.069 |      0.138 |     -1.754 |      0.077 |      0.305 |
| Cash flow adequacy (kappa) |    428 |     17.179 |    126.694 |   -259.250 |      3.512 |   1939.027 |
| Legacy intensity (alpha)  |    428 |      0.366 |      0.424 |      0.000 |      0.085 |      1.000 |

## Panel B: GEM-Matched Subsample (N = 449 firms)

| Variable                  |      N |       Mean |         SD |        Min |     Median |        Max |
|---|---|---|---|---|---|---|
| Legacy Intensity (alpha)  |    449 |      0.367 |      0.425 |      0.000 |      0.074 |      1.000 |
| Leverage (lambda)         |    443 |      1.940 |     32.766 |      0.000 |      0.390 |    690.000 |
| Operating ROA (rho)       |    439 |     -0.334 |      8.452 |   -177.000 |      0.076 |      0.305 |
| Interest Coverage (kappa) |    433 |     10.252 |    189.468 |  -2928.000 |      3.505 |   1939.027 |
| Obligation Rigidity (delta) |    436 |      0.766 |      0.231 |      0.000 |      0.846 |      1.000 |
| Network Density           |    414 |      0.085 |      0.221 |      0.001 |      0.043 |      3.091 |
| Total Assets ($M)         |    445 | 1346752.402 | 12778312.682 |      0.000 |  20431.093 | 246807795.000 |

## Panel C: Geographic Distribution

Countries in analysis sample: 69

| Country | Firms |
|---|---|
| USA | 96 |
| CHN | 45 |
| IND | 21 |
| BRA | 19 |
| TUR | 15 |
| CHL | 12 |
| JPN | 12 |
| RUS | 12 |
| ESP | 11 |
| CAN | 11 |
| PHL | 9 |
| ITA | 9 |
| THA | 9 |
| AUS | 7 |
| DEU | 7 |
| CYM | 7 |
| BMU | 6 |
| VNM | 6 |
| TWN | 6 |
| GBR | 5 |
| FRA | 5 |
| PAK | 5 |
| NZL | 5 |
| ISR | 5 |
| SGP | 4 |
| ARG | 4 |
| CHE | 4 |
| MYS | 4 |
| KOR | 3 |
| PRT | 3 |
| PER | 3 |
| COL | 3 |
| GRC | 3 |
| ARE | 3 |
| POL | 3 |
| NOR | 3 |
| MUS | 3 |
| AUT | 2 |
| NLD | 2 |
| LTU | 2 |
| ZAF | 2 |
| SAU | 2 |
| SWE | 2 |
| LKA | 2 |
| IDN | 2 |
| HKG | 1 |
| JAM | 1 |
| JOR | 1 |
| CZE | 1 |
| JEY | 1 |
| FIN | 1 |
| QAT | 1 |
| BEL | 1 |
| OMN | 1 |
| KEN | 1 |
| BGD | 1 |
| ZMB | 1 |
| IMN | 1 |
| PSE | 1 |
| UGA | 1 |
| ROU | 1 |
| HUN | 1 |
| BRB | 1 |
| VGB | 1 |
| NAM | 1 |
| MEX | 1 |
| IRL | 1 |
| NGA | 1 |
| CYP | 1 |

## Panel D: Coal Retirement Events

- Total retirements: 1844
- Matched to Compustat: 833
- First-mover events: 344
- First-mover + matched: 179
- Total retired MW: 364,784
- Matched retired MW: 209,526
- Countries with first-mover events: 35
- US first-mover events: 135

First-mover events by country:

| Country | Events |
|---|---|
| United States | 135 |
| China | 50 |
| Russia | 21 |
| India | 15 |
| South Africa | 15 |
| Australia | 11 |
| France | 11 |
| Romania | 8 |
| Germany | 6 |
| Ukraine | 6 |
| Brazil | 5 |
| Chile | 5 |
| Bulgaria | 4 |
| Poland | 4 |
| Uzbekistan | 4 |
| Japan | 4 |
| Kazakhstan | 4 |
| Netherlands | 3 |
| United Kingdom | 3 |
| Greece | 3 |
| Thailand | 3 |
| Panama | 3 |
| Italy | 2 |
| Philippines | 2 |
| Slovakia | 2 |
| Spain | 2 |
| Finland | 2 |
| Honduras | 2 |
| Canada | 2 |
| United Arab Emirates | 2 |
| Austria | 1 |
| South Korea | 1 |
| Guadeloupe | 1 |
| Indonesia | 1 |
| Israel | 1 |

## Panel E: EIA-860 Announcement Events

- Total EIA-860 events: 92

## Panel F: Coal Phase-out Events

- Total phase-out events: 27


# Country-Level Robustness for the Non-US Coefficient

## (A) Non-US event counts by country

| Country | Events |
|---|---|
| Russia | 19 |
| China | 12 |
| South Africa | 9 |
| India | 7 |
| France | 5 |
| Germany | 5 |
| Chile | 4 |
| Netherlands | 3 |
| Poland | 3 |
| Greece | 3 |
| Italy | 2 |
| Philippines | 2 |
| Slovakia | 2 |
| Spain | 2 |
| Canada | 2 |
| United Arab Emirates | 2 |
| Austria | 1 |
| Brazil | 1 |
| South Korea | 1 |
| Australia | 1 |

Total non-US events: 86

## (B) Leave-one-country-out

| Dropped | N events removed | N FM events | gamma_fuel | SE | t |
|---|---|---|---|---|---|
| (none) | 0 | 36 | -7.7617 | 0.8280 | -9.37 |
| Russia | 19 | 36 | -7.7617 | 0.8280 | -9.37 |
| China | 12 | 26 | -7.8917 | 1.0973 | -7.19 |
| South Africa | 9 | 36 | -7.7617 | 0.8280 | -9.37 |
| India | 7 | 36 | -7.7617 | 0.8280 | -9.37 |
| France | 5 | 31 | -7.3677 | 0.8652 | -8.52 |
| Germany | 5 | 33 | -7.7643 | 0.8385 | -9.26 |
| Chile | 4 | 36 | -7.7617 | 0.8280 | -9.37 |
| Netherlands | 3 | 33 | -7.1012 | 1.1051 | -6.43 |
| Poland | 3 | 33 | -8.3435 | 0.9727 | -8.58 |
| Greece | 3 | 33 | -8.5311 | 0.8431 | -10.12 |

## (C) Developed-ex-US vs Emerging vs Frontier

| Split | N events | gamma_fuel | SE | t | gamma_geo | SE | t |
|---|---|---|---|---|---|---|---|
| Developed (ex-US) | 17 | -9.6312 | 1.2422 | -7.75 | -2.2955 | 1.4978 | -1.53 |
| Emerging | 17 | -5.6241 | 1.3679 | -4.11 | +1.4806 | 0.7748 | +1.91 |
| Frontier / Other | 2 | -10.0403 | nan | +nan | -3.4070 | nan | +nan |

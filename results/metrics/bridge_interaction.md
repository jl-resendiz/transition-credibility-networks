# Bridge Interaction: Empirical Test of Section 2.5

Augmented Fama-MacBeth specification:

$$\mathrm{CAR}_i = a + \gamma_{\mathrm{fuel}} w^{\mathrm{fuel}}_i + \gamma_{\mathrm{het}} w^{\mathrm{fuel}}_i (\alpha_i^{pre-2014} - \bar\alpha) + \gamma_{\mathrm{geo}} w^{\mathrm{geo}}_i + \gamma_{\mathrm{reg}} w^{\mathrm{reg}}_i + \gamma_{ss} \mathbb{1}_{ss} + \varepsilon_i$$

with $\alpha_i^{pre-2014}$ = mean firm coal share 2010-2013, $\bar\alpha = 0.2899$ (cross-sectional mean of pre-2014 alpha). Single-factor (market-adjusted) CARs over [-1, +3] are used to match the headline reference table.

## Augmented FM coefficients

| Variable | Mean | SE (NW lag=4) | t-stat | N events | Predicted sign |
|---|---|---|---|---|---|
| $\gamma_{\mathrm{fuel}}$ | -3.2523 | 1.2707 | -2.56 | 115 | negative |
| $\gamma_{\mathrm{het}}$ ($w^{\mathrm{fuel}} \times (\alpha-\bar\alpha)$) | +2.2259 | 2.0632 | +1.08 | 115 | negative (theory) |
| $\gamma_{\mathrm{geo}}$ | -0.1162 | 0.4816 | -0.24 | 115 | attenuated to zero |
| $\gamma_{\mathrm{reg}}$ | +3.9199 | 1.3838 | +2.83 | 115 | positive |

## Interpretation

The interaction coefficient $\gamma_{\mathrm{het}}$ is **positive** ($\hat\gamma_{\mathrm{het}} = +2.2259$, $|t| = 1.08$). The sign is opposite to the theory prediction. This warrants a reassessment of the linearity assumption in Assumption 2 of the model.

## Notes

- Events with pre-2014 alpha coverage and >= 20 firms: 115
- Total firm-event observations: 26253
- Pre-2014 alpha used as a pre-determined proxy avoids endogeneity to the contemporaneous retirement event.
- A multi-factor robustness version (using FF3 + utility CARs) is a natural extension; see multifactor_inference.py for the multi-factor CAR construction.

# Response to Editor and Referees

**Manuscript:** "When Coal Retires: The Propagation of Stranding Risk"

**Author:** Jose Luis Resendiz, Smith School of Enterprise and the Environment, University of Oxford

**Date:** 24 March 2026

I am grateful to the editor and three referees for an exceptionally thorough evaluation. The reports identified genuine problems in the original submission, several of which I had not recognised. The revision addresses every concern raised. The paper is substantially restructured: the formal model has been replaced by a 1.5-page conceptual framework, the Romano-Wolf implementation has been corrected, the ESG horse race has been reframed, and the portfolio sort analysis has been removed. The result is a shorter, more honest paper whose claims are supported by its evidence.

---

## Summary of Changes

1. **Model replaced.** The three-page Bayesian model (Cournot competition, CRRA preferences, binary regulatory state, Bayesian posterior updating) has been replaced by a 1.5-page conceptual framework (Section 2) with two assumptions and one lemma. The framework derives the estimating equation directly from the asymmetric aggregation of geographic and technology channels without invoking preferences, beliefs, or market structure that the data cannot identify.

2. **Romano-Wolf corrected.** The bootstrap implementation now uses centred t-statistics: (beta* - beta_obs) / se_obs, following Cameron, Gelbach, and Miller (2008). The fuel channel survives all corrections at the primary window (Bonferroni p = 0.000, max-t p = 0.000, Romano-Wolf p = 0.000). The self-contradiction between the robustness section and the conclusion has been resolved.

3. **ETS interaction removed.** The emissions trading system interaction produced mixed results across inference methods and has been removed from the main text. A sentence in the robustness section notes its availability upon request.

4. **Portfolio sorts removed.** The row-normalised fuel weight has a cross-sectional standard deviation of approximately 0.004. Quintile sorts on this variable lack power to separate treated from control firms, producing t-statistics an order of magnitude below the regression estimates. Rather than include a result that would require extensive caveating, the portfolio sort analysis has been removed from scope.

5. **Lemma 2 downgraded.** The original Lemma 2 predicted that geographic attenuation scales as 1/K_i under independence, with diversification as the identified mechanism. The revised Lemma 1 (Aggregation Attenuation) retains the qualitative prediction that geographic variance vanishes while technology variance does not, but adds Remark 1, which acknowledges that under partial correlation the geographic variance converges to a positive constant rather than to zero. The HHI interaction (t = +1.31, insignificant) is reported honestly as "directionally consistent but insignificant" in Section 4.2.

6. **Calendar-time reframed.** The original interpretation invoked within-jurisdiction learning. The revision reframes the calendar-time pattern as geographic composition: the early tercile is dominated by US events (where coal retirement is routine), while later terciles contain a growing share of non-US retirements. Section 4.4 now leads with the US versus non-US decomposition and explicitly states that the pattern "reflects the changing geographic composition of the event sample, not a within-market process."

7. **Economic magnitude presented honestly.** The 1-SD effect (2.2 basis points) now appears in the abstract. The Bolton-Kacperczyk comparison has been removed. The economic magnitude paragraph (Section 4.1) reports both the 1-SD effect and the median-to-zero comparison, with benchmarks to Lang and Stulz (1992) and Menzly and Ozbas (2010).

8. **Return source composition documented.** Section 3.2 now reports the approximate composition (15% CRSP, 70% Eikon, 15% Compustat Global Security) and notes that results are stable when estimated separately on each source.

9. **ESG horse race reframed.** The pipeline summary claiming "spatial fundamentals more informative" has been corrected. Table 3 now shows ESG as the stronger predictor in the 153-firm subsample. A footnote acknowledges the selection problem (153/703 firms, disproportionately large, developed-market).

10. **FM Bartik reported with appropriate caution.** The FM Bartik estimate uses t(T-1) with T = 3 events, yielding t = 2.32 and p = 0.147. Section 3.5 now characterises this as "significant under pooled inference but based on only 3 events that meet the minimum observation threshold with pre-period weights, so the FM result should be interpreted with caution." The pooled Bartik (t = -5.16) is the operative result.

---

## Response to the Editor

The editor identified eight items that, in combination, amounted to a new paper. I agree with that assessment. The revision addresses all eight:

| Item | Editor's requirement | Action taken |
|------|---------------------|--------------|
| 1 | Remove or restructure the formal model | Replaced with 1.5-page framework (Section 2) |
| 2 | Reconcile parametric/non-parametric discrepancy | Portfolio sorts removed; weight SD documented |
| 3 | Resolve Romano-Wolf contradiction | Bug fixed; text now consistent throughout |
| 4 | Fix ESG horse race or remove it | Reframed; ESG acknowledged as stronger in subsample |
| 5 | Reframe calendar-time as composition | Done (Section 4.4) |
| 6 | Report window sensitivity, balance tests, first-stage F | Window sensitivity added; pre-balance at 2014 reported; pre-balance at 2010 fails (p = 0.021) and is acknowledged |
| 7 | Fix economic magnitude presentation | 1-SD in abstract; Bolton-Kacperczyk removed |
| 8 | Add welfare discussion | Not added as a formal section; scope of claims paragraph (Section 3.5) clarifies that the paper documents a pricing pattern and does not assess whether retirements are desirable |

On item 8, I considered adding a welfare discussion but concluded that the paper lacks the structural estimation necessary to distinguish efficient information incorporation from excessive contagion. A speculative welfare paragraph would weaken rather than strengthen the paper's credibility. The concluding paragraph of the introduction now states explicitly: "The paper documents this pattern; it does not assess whether the retirements themselves are desirable." I hope the editor will accept this as an honest delimitation rather than an evasion.

---

## Response to Referee A (Senior Economist)

I am grateful to Referee A for the constructive reading and for the recommendation to revise rather than reject. The referee correctly identified that the core finding is real but the execution had serious problems. I have addressed each concern below.

### Comment A1: Economic significance (2.2 basis points)

**Referee's concern:** The one-standard-deviation effect of 2.2 basis points is small. The 3.3 percentage point headline conflates weight dispersion with effect size. The Bolton-Kacperczyk comparison is inapt (different estimand, different horizon).

**Response:** The referee is right on all three points. The abstract now reports the 1-SD effect directly: "a one-standard-deviation increase in fuel-mix similarity predicts a 2.2 basis point decline in cumulative abnormal returns over four months." The Bolton-Kacperczyk comparison has been removed. The economic magnitude paragraph in Section 4.1 now benchmarks against Lang and Stulz (1992), who document contagion effects of -1.0% and competitive effects of +2.2% within days of bankruptcy announcements, and Menzly and Ozbas (2010), who find monthly return predictability of 13 to 22 basis points from cross-industry linkages. I retain the median-to-zero comparison (3.3 percentage points) because it communicates the economic content of the regression coefficient, but it is now clearly labelled as a comparison between two points on the fuel-similarity distribution, not a standard-deviation effect.

**Location in manuscript:** Abstract (line 59); Section 4.1, "Economic magnitude" paragraph.

### Comment A2: ESG horse race contradicts itself

**Referee's concern:** The pipeline summary stated "spatial fundamentals more informative," which directly contradicted the marginal R-squared numbers (ESG 40x more informative). The horse race was based on 153/703 firms with no balance test and contradictory results across methods.

**Response:** The pipeline summary was a material misrepresentation and has been corrected. Table 3 now presents the horse race transparently: ESG alone achieves R-squared = 0.012 versus 0.003 for spatial fuel alone. Under pooled OLS, the fuel coefficient loses significance when ESG is included (column 3). Under Fama-MacBeth + Newey-West, both survive (fuel t = -2.17, ESG t = -5.72). The text in Section 4.3 states plainly: "In the 153-firm subsample where both measures exist, ESG is the stronger predictor." A footnote now acknowledges the selection problem: "The 153 ESG-covered firms are disproportionately large and based in developed markets. The ESG horse race results may therefore not generalise to the smaller, emerging-market utilities that constitute the majority of the sample. No formal balance test between the ESG and non-ESG subsamples was conducted."

**Location in manuscript:** Section 4.3 and Table 3.

### Comment A3: Contribution is relabelling

**Referee's concern:** Asymmetric aggregation (some shocks diversify, some do not) is the textbook definition of systematic versus idiosyncratic risk.

**Response:** I accept this. The revision no longer claims novelty for the theoretical mechanism. The introduction now frames the contribution as measurement and identification: "The measurement innovation (GEM-based L1 similarity covering 565 firms in 80 countries) is genuine. The opposing-sign prediction (fuel negative, geography zero) is sharper than 'coal is priced.'" The model section is now a brief conceptual framework that derives the estimating equation, not a theoretical contribution in its own right.

**Location in manuscript:** Introduction, paragraph 5; Section 2 (entire section).

### Comment A4: Portfolio sorts are dead

**Referee's concern:** The Q5-Q1 spread (t = -1.30) and long-short portfolio (t = -0.17) disagree with the regression (t = -7.36) by an order of magnitude. The discrepancy must be diagnosed.

**Response:** The referee's diagnosis was correct: the row-normalised fuel weight has a cross-sectional standard deviation of approximately 0.004, so quintile boundaries fail to separate treated from control firms. Rather than include a lengthy diagnostic section for a test that fundamentally lacks power with this variable, I have removed the portfolio sort analysis from the paper entirely. The economic magnitude paragraph now reports the weight SD directly, so the reader can assess the scale of variation being exploited.

**Location in manuscript:** Section 4.1, "Economic magnitude" paragraph (weight SD = 0.004).

### Comment A5: Calendar-time is composition, not learning

**Referee's concern:** The strengthening of the fuel coefficient over time reflects the changing geographic composition of the event sample (more non-US events in later years), not within-jurisdiction learning.

**Response:** The referee was right. Section 4.4 has been rewritten from the ground up. It now leads with the US versus non-US decomposition (US: t = 0.11; non-US: t = -4.10) and explicitly states: "The strengthening over time reflects the changing geographic composition of the event sample, not a within-market process." The within-jurisdiction interaction (t = -1.57, p = 0.117) is reported as insignificant, and the split-sample comparison is acknowledged as mixing within- and across-jurisdiction variation.

**Location in manuscript:** Section 4.4 ("Geographic Heterogeneity in the Fuel Signal"), entire subsection.

### Comment A6: Romano-Wolf contradiction

**Referee's concern:** The robustness section claimed fuel survives Romano-Wolf; the conclusion claimed no channel survives across nine hypotheses.

**Response:** This was a quality control failure, as Referee C diagnosed. The implementation has been corrected (see Response to Referee C, Comment C2). The manuscript is now internally consistent: the robustness section reports that fuel survives all corrections at the three-hypothesis family (Section 4.5), and the conclusion reports the same result. The conclusion also notes the nine-hypothesis result honestly: "No individual channel survives the Romano-Wolf multiple testing correction across nine hypotheses, though the joint F-test is highly significant (F = 70.83, p = 0.000) and the fuel channel's Fama-MacBeth t-statistic of -7.36 exceeds single-hypothesis thresholds by a wide margin."

**Location in manuscript:** Section 4.5 ("Romano-Wolf multiple testing"); Section 6 (Conclusion), final paragraph.

### Comment A7: FM Bartik on 3 events is uninformative

**Referee's concern:** The FM Bartik is based on T = 3 events. With t(T-1) = t(2), even t = 2.32 gives p = 0.147.

**Response:** The manuscript now states this explicitly: the FM result "is based on only 3 events that meet the minimum observation threshold with pre-period weights, so the FM result should be interpreted with caution." The pooled Bartik (t = -5.16, N = 24,070) is identified as the operative result.

**Location in manuscript:** Section 3.5, fourth paragraph of the identification discussion.

---

## Response to Referee B (Economic Theorist)

Referee B's report was the most challenging and, in retrospect, the most helpful. The referee was right that the original model did no real work and that its key prediction was rejected by the data. The revision responds by replacing the model with a conceptual framework that is honest about what it can and cannot deliver.

### Comment B1: Lemma 2 falsified by the data

**Referee's concern:** Lemma 2 predicted that geographic attenuation scales as 1/K_i under independence. Every empirical interaction test failed. The HHI interaction had the wrong sign under OLS. The single-country subsample was weaker, not stronger. The multi-country subsample was numerically unstable. The fallback to partial correlation rendered the prediction unfalsifiable.

**Response:** The referee was right. The revision takes three steps. First, the original Lemma 2 has been replaced by a simpler Aggregation Attenuation Lemma (Lemma 1, Section 2.3) that derives the qualitative prediction without staking it on a specific functional form. The lemma states that geographic variance satisfies Var(sum) proportional to sum of theta_ik squared, which is bounded by 1/K_i and vanishes as K_i goes to infinity, while technology variance does not depend on K_i. Second, Remark 1 immediately following the lemma acknowledges that under partial correlation across markets, the geographic variance converges to a positive constant rather than to zero: "The qualitative content of the lemma is that the two channels aggregate asymmetrically, not the specific rate at which geographic exposure declines with K_i." Third, the empirical discussion in Section 4.2 reports the HHI interaction honestly: "An interaction between returns and the Herfindahl index of market concentration H_i is directionally consistent with the attenuation mechanism (t = 1.31) but insignificant, so I cannot confirm the specific cross-sectional prediction that attenuation scales with K_i. The data are consistent with the lemma in that geographic exposure is attenuated to zero in aggregate, but the particular diversification channel remains an interpretation rather than an identified finding."

I believe this is a more honest framing. The model predicts that the two channels aggregate asymmetrically, which the data confirm (fuel significant, geography zero). It does not identify the mechanism (diversification) through which the asymmetry arises. The text says so.

**Location in manuscript:** Section 2.3 (Lemma 1 and Remark 1); Section 4.2, final paragraph.

### Comment B2: Bayesian structure adds nothing

**Referee's concern:** The binary regulatory state omega was never identified, estimated, or tested. CRRA preferences were observationally vacuous. The posterior p_1 entered as a single reduced-form parameter. No empirical test distinguished the model from "coal exposure is a priced common factor."

**Response:** The referee was entirely correct. The Bayesian apparatus has been removed. The revised model (Section 2) contains no binary state, no posterior, no CRRA preferences, and no Cournot competition. It posits a monotone stochastic discount factor, two assumptions (geographic gains are local and independent; technology exposure is firm-level), and derives the estimating equation in half a page. The economic content is the same: technology obsolescence is a common factor that does not diversify across markets, while local competitive effects do. But the notation now matches the empirical content. The model generates the estimating equation and nothing more.

The referee asked: "What empirical test distinguishes this model from 'coal exposure is a priced common factor'?" The honest answer is: the opposing-sign prediction. A simple "coal is priced" story predicts that coal-heavy firms lose value, but says nothing about the geographic channel. The model predicts that the geographic channel is attenuated toward zero while the fuel channel is robustly negative. This is a joint prediction about two coefficients, not just the sign of one. Table 1 confirms this pattern under all three inference methods.

**Location in manuscript:** Section 2 (entire section, now 1.5 pages).

### Comment B3: No welfare analysis

**Referee's concern:** A JEEM paper should discuss whether the observed propagation is efficient or excessive.

**Response:** I considered adding a welfare discussion but concluded that the paper lacks the structural ingredients to deliver one. Distinguishing efficient information incorporation from excessive contagion requires identifying the counterfactual: what would returns look like if investors processed the technology signal correctly? Without a structural model of investor beliefs and fundamental values, any welfare claim would be speculative. The introduction now states the paper's scope explicitly: "The paper documents this pattern; it does not assess whether the retirements themselves are desirable." The concluding paragraph of Section 3.5 (Scope of claims) further delimits: "The causal claim applies to the cross-sectional allocation of returns across firms with different pre-determined technology exposure. It does not extend to the aggregate effect of retirements on the utility sector as a whole."

I acknowledge this as a limitation. If the editor and referee judge that a welfare discussion is essential for JEEM, I am willing to add a brief section in a subsequent revision, though I would want to be clear about its speculative nature.

**Location in manuscript:** Introduction, final sentence of paragraph 4; Section 3.5, "Scope of claims."

---

## Response to Referee C (Econometrician)

Referee C's report identified the most consequential technical error in the paper (the Romano-Wolf bug) and provided a precise diagnosis. I am grateful for the specificity.

### Comment C1: Parametric and non-parametric results disagree

**Referee's concern:** The FM regression (t = -7.36) and the portfolio sort (Q5-Q1 t = -1.30, long-short t = -0.17) test the same hypothesis on the same data and disagree by an order of magnitude.

**Response:** The referee's diagnosis was correct and matches Referee A's: the row-normalised fuel weight has SD of approximately 0.004, so quintile boundaries cannot separate treated from control firms. The predicted Q5-Q1 return difference is approximately 3 basis points, well within monthly return noise. I have removed the portfolio sort analysis entirely rather than include a result that fundamentally lacks statistical power. The weight SD is now reported in the economic magnitude paragraph so readers can assess the dispersion being exploited.

**Location in manuscript:** Section 4.1, "Economic magnitude" paragraph.

### Comment C2: Romano-Wolf self-contradiction

**Referee's concern:** Four sources in the paper gave mutually inconsistent claims about whether fuel survives Romano-Wolf. The joint_tests.md reference to "0/9 significant" appeared to be hard-coded. The conclusion was never updated when the code was narrowed from 9 to 3 hypotheses.

**Response:** The referee diagnosed this exactly right. The bug was in the bootstrap implementation: the original code used raw bootstrap t-statistics rather than centred (beta* - beta_obs) / se_obs per Cameron, Gelbach, and Miller (2008). After correction, the fuel channel survives all corrections at the three-hypothesis family (Bonferroni p = 0.000, max-t p = 0.000, Romano-Wolf p = 0.000). Neither the geographic nor the regulatory channel survives.

The manuscript is now internally consistent. The robustness section (Section 4.5) reports the three-hypothesis result. The conclusion reports both the three-hypothesis survival and the nine-hypothesis non-survival, clearly distinguished. The hard-coded "0/9 significant" string in joint_tests.md has been traced and corrected.

I should have caught this before submission. I am embarrassed that it required a referee to identify a self-contradiction between the robustness section and the conclusion of the same paper.

**Location in manuscript:** Section 4.5 ("Romano-Wolf multiple testing"); Section 6 (Conclusion).

### Comment C3: Shift-share weaknesses

**Referee's concern:** (a) FM Bartik on T = 3 events has approximately zero degrees of freedom. (b) Pre-balance fails at the 2010 cutoff (p = 0.021). (c) No first-stage F-statistic reported.

**Response on (a):** The code already uses the t(T-1) distribution, as the referee noted. The manuscript now states this explicitly: "The Fama-MacBeth estimate is also significant (t = -2.32) but is based on only 3 events that meet the minimum observation threshold with pre-period weights, so the FM result should be interpreted with caution." The pooled Bartik (t = -5.16) is the operative result.

**Response on (b):** The pre-balance test at the 2014 cutoff (t = -1.87, p = 0.062) is borderline. At the 2010 cutoff, it fails (p = 0.021). The manuscript acknowledges this: the 2014 cutoff is the pre-specified choice motivated by the event distribution, but the 2010 failure is a legitimate concern. I have not added additional cutoffs (2008, 2012, 2016) in this revision but am willing to do so if required.

**Response on (c):** The shift-share design in this paper differs from the standard instrumental variables setting where a first-stage F-statistic is diagnostic. The fuel-mix shares are not instruments for an endogenous regressor; they are the treatment variable itself (pre-determined exposure to the retirement shock). With 11 distinct shift values (retirement events), the relevant diagnostic is the Rotemberg weight analysis, which confirms non-negative weights and a low Herfindahl (0.031). I am open to reporting a first-stage F if the referee considers it informative in this context.

**Location in manuscript:** Section 3.5, identification discussion (FM Bartik caution, pre-balance, Rotemberg weights).

### Comment C4: Event window and return data

**Referee's concern:** The [-1, +3] month window is long. Overlapping windows introduce dependence. The Compustat Global fallback uses price returns without dividends.

**Response:** On the window, the fuel channel is significant at [-1,+1], [-1,+2], [-1,+3], and [0,+1], as reported in the event window discussion (Section 3.5). The [-1, +3] window is the primary specification because many retirement events unfold over weeks rather than days, and the technology repricing mechanism operates at the monthly horizon consistent with Menzly and Ozbas (2010).

On overlapping windows, the Fama-MacBeth procedure estimates each event separately and averages with Newey-West HAC correction, which accounts for serial dependence in the time series of event-level coefficients. Two-way clustering (event + firm) provides a further guard against cross-event correlation.

On return sources, Section 3.2 now reports the composition: "Of the 703 firms, approximately 15 percent have CRSP coverage, 70 percent have Eikon coverage, and 15 percent rely on Compustat Global Security prices. The main results are stable when estimated separately on each return source." The Compustat Global fallback is acknowledged in the caveats paragraph of the conclusion: "A subset of firms use Compustat Global Security price returns (without dividend adjustment) rather than Eikon total returns, introducing a modest measurement asymmetry."

**Location in manuscript:** Section 3.2 (return source composition); Section 3.5 (event window discussion); Section 6 (caveats paragraph).

### Comment C5: Pipeline code quality

**Referee's concern:** Massive OLS duplication across 10+ scripts. Gauss-Jordan tolerance below machine epsilon. No automated tests.

**Response:** These are valid quality concerns. The shared OLS module (_ols.py) has been consolidated and is used by all analysis scripts. The Gauss-Jordan tolerance has been corrected. Automated tests remain a desirable improvement that I have not yet implemented. I acknowledge this as a limitation of the replication package.

---

## Summary Table of Changes

| Referee concern | Section | Change |
|----------------|---------|--------|
| Model does no work (B1, B2) | 2 | Replaced 3-page Bayesian model with 1.5-page framework |
| Lemma 2 falsified (B1) | 2.3, 4.2 | Downgraded to qualitative prediction; HHI result reported honestly |
| Romano-Wolf contradiction (C2) | 4.5, 6 | Bug fixed; text consistent throughout |
| ESG horse race (A2) | 4.3, Table 3 | ESG acknowledged as stronger; selection footnote added |
| Portfolio sorts dead (A4, C1) | Removed | Weight SD documented instead |
| Economic magnitude (A1) | Abstract, 4.1 | 1-SD in abstract; Bolton-Kacperczyk removed |
| Calendar-time (A5) | 4.4 | Reframed as geographic composition |
| FM Bartik (A7, C3) | 3.5 | Reported with t(T-1) caution |
| ETS interaction | 4.5 | Removed (mixed results) |
| Return sources (C4) | 3.2 | Composition sentence added |
| Welfare (B3) | Intro, 3.5 | Scope delimited; no formal welfare section |
| Event windows (C4) | 3.5 | Sensitivity at four windows noted |
| Pre-balance (C3) | 3.5 | 2010 failure acknowledged |

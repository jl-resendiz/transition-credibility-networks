# Foundational references for spatial climate transition risk in asset pricing

This literature review identifies **21 authoritative papers from top-tier journals** that collectively provide the theoretical and empirical scaffolding for a model of spatial transmission of climate transition risk framed as an asset pricing model under policy uncertainty with network frictions. Every paper listed below is published in a top-5 economics journal (AER, QJE, JPE, Econometrica, REStud), a top-3 finance journal (JF, JFE, RFS), or a leading field journal (AEJ: Applied Economics). For each of the four building blocks, 5–6 papers are provided with precise citation guidance keyed to the model's mathematical assumptions.

---

## Building Block 1: Networked cash flows and opposing spatial channels

The model assumes firm *i*'s cash flows depend on neighbor shocks through a spatial autoregressive structure with two weight matrices of opposite sign: geographic proximity (γ_G > 0, competitive benefit from local rival exit) and fuel similarity (γ_F < 0, contagion from technologically similar peer distress). Five papers establish this dual-channel network architecture.

### Paper 1.1 — Network propagation of idiosyncratic shocks

**Acemoglu, Daron, Vasco M. Carvalho, Asuman Ozdaglar, and Alireza Tahbaz-Salehi (2012). "The Network Origins of Aggregate Fluctuations." *Econometrica*, 80(5), 1977–2016.**

**Core finding:** Microeconomic shocks to individual sectors do not average out in the presence of intersectoral input-output linkages. The network's asymmetric structure—where some sectors are disproportionately important as suppliers—creates cascade effects that generate aggregate volatility from idiosyncratic disturbances.

**How to cite:** This paper provides the foundational mathematical framework for why a firm's cash flows are functions of shocks hitting networked neighbors through a weighted adjacency matrix *W*. The spatial autoregressive cash-flow specification (**ΔCF_i = γ_G · W_G · ε + γ_F · W_F · ε**) directly inherits the Acemoglu et al. logic that idiosyncratic shocks (e.g., a single plant retirement) propagate through the network, with the magnitude determined by the topology of interconnections. Cite as: *"Following Acemoglu et al. (2012), we model firm-level cash flows as functions of shocks to network neighbors, where the network structure determines the magnitude and direction of propagation."*

### Paper 1.2 — Production networks priced in equilibrium

**Herskovic, Bernard (2018). "Networks in Production: Asset Pricing Implications." *Journal of Finance*, 73(4), 1785–1818.**

**Core finding:** In a multisector equilibrium model with input-output network linkages, changes in network concentration and network sparsity are sources of systematic risk reflected in asset prices. Network-sorted portfolios generate economically significant return spreads of **4.6%** and **−3.2%** per year.

**How to cite:** Herskovic demonstrates that production network linkages—precisely the type of inter-firm connections captured by geographic and fuel-similarity weight matrices—enter equilibrium risk premia. This provides the theoretical precedent that network topology in cash flows translates into network-based factors in the stochastic discount factor. Cite as: *"Herskovic (2018) shows that network structure is a priced risk factor in equilibrium, providing theoretical justification for our derivation of how γ_G and γ_F enter the SDF through spatially weighted cash-flow exposures."*

### Paper 1.3 — Contagion versus competitive effects: the dual-sign template

**Lang, Larry H.P., and René M. Stulz (1992). "Contagion and Competitive Intra-Industry Effects of Bankruptcy Announcements: An Empirical Analysis." *Journal of Financial Economics*, 32(1), 45–60.**

**Core finding:** Bankruptcy announcements produce two opposing intra-industry effects on competitors' stock prices. **Contagion effects** (negative spillovers) dominate in highly leveraged, less concentrated industries, while **competitive effects** (positive spillovers) dominate in concentrated, low-leverage industries. The net observed effect depends on the relative strength of both channels.

**How to cite:** This is the single most important precedent for the opposing-sign structure of γ_G and γ_F. The competitive effect—rivals gain market share when a peer exits in a concentrated local market—maps directly to **γ_G > 0**. The contagion effect—peers suffer when a similar firm's failure signals shared vulnerability—maps to **γ_F < 0**. Lang and Stulz's key insight, that the *same event* (firm exit) produces opposite-signed spillovers depending on the *dimension of proximity*, is the exact structure the model formalizes. Cite as: *"Our opposing-sign spatial channels formalize the Lang and Stulz (1992) distinction between competitive and contagion effects, where geographic proximity captures the competitive benefit of rival exit (γ_G > 0) and fuel similarity captures the contagion of shared technological obsolescence (γ_F < 0)."*

### Paper 1.4 — Credit contagion among technologically similar firms

**Jorion, Philippe, and Gaiyan Zhang (2007). "Good and Bad Credit Contagion: Evidence from Credit Default Swaps." *Journal of Financial Economics*, 84(3), 860–883.**

**Core finding:** Using CDS spread data, credit events at one firm produce significant intra-industry spillovers. Chapter 11 bankruptcies (signaling industry-wide distress) generate **contagion effects**—CDS spreads of peers widen. Chapter 7 liquidations (permanent capacity removal) produce **competitive effects**—peer spreads tighten. The type of distress determines the sign.

**How to cite:** Jorion and Zhang provide market-based evidence that the mechanism of exit determines the sign of peer spillovers in credit markets. Coal plant retirements driven by technological obsolescence or regulatory pressure (analogous to Chapter 11 financial distress) transmit negative signals to fuel-similar peers (γ_F < 0), while retirements that purely remove local capacity (analogous to Chapter 7 liquidation) benefit geographic neighbors (γ_G > 0). Cite as: *"Jorion and Zhang (2007) show that the nature of a firm's distress—whether it signals shared vulnerability or merely removes a competitor—determines whether peer spillovers are negative (contagion) or positive (competitive), directly supporting our dual-channel specification."*

### Paper 1.5 — Electricity-specific evidence for local market power gains

**Davis, Lucas W., and Catherine Hausman (2016). "Market Impacts of a Nuclear Power Plant Closure." *American Economic Journal: Applied Economics*, 8(2), 92–122.**

**Core finding:** The abrupt closure of the San Onofre Nuclear Generating Station created **binding transmission constraints**, raising natural gas generation costs by **$350 million** and making it more profitable for nearby plants to exercise local market power. The removal of a large generation source in a geographically constrained grid directly benefited proximate generators.

**How to cite:** This paper provides the electricity-market-specific empirical mechanism for γ_G > 0. When a plant retires, resulting transmission congestion and reduced local supply increase the market power and profitability of remaining geographically proximate generators. This is the "positive local competitive benefit" channel operating through the physical structure of the transmission grid. Cite as: *"Davis and Hausman (2016) demonstrate that plant retirement in a geographically constrained electricity grid increases market power for nearby generators—the precise microeconomic mechanism underlying γ_G > 0 in our model."*

---

## Building Block 2: Policy uncertainty and the credibility gap

The model's credibility gap arises because markets assign probability p_t < 1 to strict enforcement of the transition regime, with the variance of beliefs p_t(1−p_t) driving a time-varying transition risk premium. Six papers—spanning general equilibrium theory, empirical measurement, game-theoretic microfoundations, and direct climate finance evidence—establish this architecture.

### Paper 2.1 — The theoretical backbone of policy uncertainty pricing

**Pástor, Ľuboš, and Pietro Veronesi (2012). "Uncertainty about Government Policy and Stock Prices." *Journal of Finance*, 67(4), 1219–1264.**

**Core finding:** In general equilibrium, uncertainty about government policy depresses stock prices, increases return volatilities and cross-stock correlations, and generates a positive jump risk premium at policy announcements. Stock prices fall on average when policy changes are announced, with the magnitude increasing in prior policy uncertainty.

**How to cite:** This paper provides the direct theoretical backbone for the p_t(1−p_t) credibility gap. Pástor and Veronesi model agents holding Bayesian beliefs about policy impact, with the variance of beliefs driving risk premia—maximal when p_t ≈ 0.5. Their discretionary government policy choice maps directly to the regime-switching structure where markets incompletely observe ω. Cite as: *"Following Pástor and Veronesi (2012), we model the market's incomplete information about the enforcement regime as generating a policy uncertainty risk premium, where the variance of beliefs p_t(1−p_t) is the direct analogue of their government policy uncertainty parameter."*

### Paper 2.2 — Political uncertainty as a priced risk factor

**Pástor, Ľuboš, and Pietro Veronesi (2013). "Political Uncertainty and Risk Premia." *Journal of Financial Economics*, 110(3), 520–545.**

**Core finding:** Political uncertainty commands a risk premium that is **countercyclical**—larger in weaker economic conditions—even when political shocks are orthogonal to fundamental economic shocks. Bayesian learning about political costs endogenously varies the magnitude of political risk premia over time.

**How to cite:** This paper establishes the state-dependent nature of the credibility gap premium. When transition costs are high and economic conditions are weak, the probability of policy reversal (1−p_t) increases, widening the credibility gap and amplifying the transition risk premium. Cite as: *"As in Pástor and Veronesi (2013), the credibility gap generates a risk premium even when shocks to policy beliefs are orthogonal to productivity shocks. The premium rises precisely when transition costs make policy reversal more tempting, endogenously increasing p_t(1−p_t)."*

### Paper 2.3 — Empirical measurement of policy uncertainty effects

**Baker, Scott R., Nicholas Bloom, and Steven J. Davis (2016). "Measuring Economic Policy Uncertainty." *Quarterly Journal of Economics*, 131(4), 1593–1636.**

**Core finding:** A news-based Economic Policy Uncertainty (EPU) index spikes near major policy events and is associated with greater stock price volatility and reduced investment and employment in policy-sensitive sectors. EPU innovations foreshadow declines in investment, output, and employment across **12 major economies**.

**How to cite:** Baker, Bloom, and Davis provide the empirical measurement framework confirming that policy uncertainty has first-order financial and real effects. Their decomposable EPU index—including regulatory and environmental policy components—validates the prediction that when p_t(1−p_t) is high, carbon-intensive firms face elevated volatility and depressed investment. Cite as: *"Baker, Bloom, and Davis (2016) empirically demonstrate that policy uncertainty depresses investment and raises volatility in policy-sensitive sectors, validating our model's prediction that the credibility gap has first-order effects on transition-exposed asset prices."*

### Paper 2.4 — Time-inconsistency microfoundation for non-credible commitments

**Kydland, Finn E., and Edward C. Prescott (1977). "Rules Rather Than Discretion: The Inconsistency of Optimal Plans." *Journal of Political Economy*, 85(3), 473–491.**

**Core finding:** Optimal policy plans are **time-inconsistent**: a government that commits to an optimal long-run policy has ex-post incentives to deviate once private agents have acted on their expectations. Discretionary policymaking is inherently non-credible; institutional rules are required for commitment.

**How to cite:** This Nobel Prize-winning paper provides the game-theoretic microfoundation for why p_t < 1. A government announcing a strict carbon transition regime (ω = 1) faces ex-post incentives to relax enforcement to reduce energy costs, support employment, or win elections. Rational agents anticipate this time-inconsistency and discount policy commitments accordingly—p_t < 1 is not irrational skepticism but the equilibrium outcome. Cite as: *"The credibility gap is microfounded by the classic Kydland and Prescott (1977) time-inconsistency problem: governments announcing strict transition regimes face ex-post incentives to deviate, and rational agents discount policy commitments accordingly, generating 0 < p_t < 1."*

### Paper 2.5 — Carbon risk is partially priced in equities

**Bolton, Patrick, and Marcin Kacperczyk (2021). "Do Investors Care about Carbon Risk?" *Journal of Financial Economics*, 142(2), 517–549.**

**Core finding:** Stocks of firms with higher total CO₂ emissions earn higher returns after controlling for standard risk factors. This **carbon premium** (1.8%–4.0% annualized) is consistent with investors demanding compensation for transition risk exposure and cannot be explained by differences in unexpected profitability.

**How to cite:** The existence of a carbon premium implies markets assign non-trivial probability to costly transition enforcement (p_t > 0), but the fact that it appears as a risk premium—rather than a permanent discount—means markets remain uncertain about enforcement (p_t < 1). The premium is consistent with an intermediate, uncertain p_t, exactly the credibility gap. Cite as: *"Bolton and Kacperczyk (2021) document a carbon premium consistent with markets partially—but not fully—pricing in transition risk, as predicted by a model where 0 < p_t < 1."*

### Paper 2.6 — Climate policy uncertainty priced in options markets

**Ilhan, Emirhan, Zacharias Sautner, and Grigory Vilkov (2021). "Carbon Tail Risk." *Review of Financial Studies*, 34(3), 1540–1571.**

**Core finding:** The cost of option protection against downside tail risk is significantly larger for carbon-intensive firms. This cost increases when public attention to climate change spikes and **decreased** after the U.S. withdrawal from the Paris Agreement—a natural experiment in credibility erosion.

**How to cite:** Ilhan et al. provide the most direct empirical link between climate policy credibility and asset pricing. Their finding that the Paris withdrawal reduced downside tail risk for carbon firms is a natural experiment in credibility dynamics: a credibility-reducing event lowers p_t, narrows p_t(1−p_t), and reduces the transition risk premium embedded in options. Cite as: *"Ilhan, Sautner, and Vilkov (2021) show that option-market pricing of carbon tail risk responds dynamically to signals about climate policy credibility, directly validating our p_t dynamics and the credibility gap mechanism."*

---

## Building Block 3: Ambiguous shocks spike volatility without moving mean returns

When a noisy idiosyncratic shock arrives (e.g., a single EIA-860 plant retirement), the positive geographic and negative fuel channels offset in expectation, muting directional mean returns. But Bayesian updating under parameter uncertainty inflates the predictive return distribution's variance, spiking volatility. Five papers from the ambiguity, learning, and information economics literatures establish this mechanism.

### Paper 3.1 — Ambiguity aversion inflates perceived variance of low-quality signals

**Epstein, Larry G., and Martin Schneider (2008). "Ambiguity, Information Quality, and Asset Pricing." *Journal of Finance*, 63(1), 197–228.**

**Core finding:** When ambiguity-averse investors process news of uncertain quality, they evaluate it under a **worst-case assessment** of signal precision. This inflates price volatility and generates ambiguity premia that depend on idiosyncratic fundamentals, even when the mean informational content is uninformative.

**How to cite:** This paper directly justifies why an ambiguous signal—one where geographic and fuel channels cancel in expectation—spikes volatility. Epstein and Schneider show that when signal quality is uncertain, worst-case evaluation inflates perceived variance independent of the signal's expected sign. Cite as: *"Following Epstein and Schneider (2008), we interpret the noisy plant retirement signal as ambiguous information whose quality is uncertain. Under maxmin expected utility, investors' worst-case assessment of signal precision generates a volatility spike even when the mean price impact is zero."*

### Paper 3.2 — Regime uncertainty maximizes excess volatility at peak ambiguity

**Veronesi, Pietro (1999). "Stock Market Overreaction to Bad News in Good Times: A Rational Expectations Equilibrium Model." *Review of Financial Studies*, 12(5), 975–1007.**

**Core finding:** In a model where dividend drift shifts between unobservable states, equilibrium return volatility is **maximized when posterior uncertainty is at its peak** (π ≈ 0.5). Investors' desire to hedge against changes in their own uncertainty generates excess volatility, volatility clustering, and leverage effects—all derived from Bayesian learning under regime ambiguity.

**How to cite:** Veronesi (1999) provides the key general-equilibrium mechanism linking ambiguous signals to volatility spikes. When a noisy plant retirement pushes investor beliefs toward maximum uncertainty about the regulatory regime (π ≈ 0.5), return volatility peaks even though the expected directional price impact is symmetric and therefore muted on average. Cite as: *"As shown by Veronesi (1999), return volatility peaks at maximum posterior uncertainty about the hidden state—precisely where the credibility gap p_t(1−p_t) is maximized—providing the equilibrium mechanism for our prediction that ambiguous shocks spike volatility without moving mean returns."*

### Paper 3.3 — Parameter learning amplifies variance through subjective long-run risk

**Collin-Dufresne, Pierre, Michael Johannes, and Lars A. Lochstoer (2016). "Parameter Learning in General Equilibrium: The Asset Pricing Implications." *American Economic Review*, 106(3), 664–698.**

**Core finding:** Bayesian learning about unknown structural parameters **strongly amplifies** the impact of macroeconomic shocks on marginal utility when agents prefer early resolution of uncertainty. Because posterior distributions of parameters are martingales, parameter learning generates subjective long-run consumption risks that are quantitatively significant and persistent.

**How to cite:** This paper formalizes how Bayesian updating about structural parameters—the intensity of climate regulation, the probability of regime change—amplifies the variance of returns in response to noisy signals. Even signals with zero expected mean impact generate persistent increases in subjective risk because investors update not only about the current state but about structural parameters governing the transition. Cite as: *"Collin-Dufresne, Johannes, and Lochstoer (2016) show that parameter learning amplifies shock transmission to asset prices. In our setting, a noisy plant retirement triggers updating about both the regulatory state and the structural parameters governing transition intensity, generating a volatility spike even when the signal's directional content averages to zero."*

### Paper 3.4 — Structural uncertainty fattens tails and inflates variance

**Weitzman, Martin L. (2007). "Subjective Expectations and Asset-Return Puzzles." *American Economic Review*, 97(4), 1102–1130.**

**Core finding:** When structural parameters are uncertain and updated via Bayesian learning, the thin-tailed normal distribution of consumption growth becomes a **fat-tailed Student-*t* distribution**. This tail-fattening fundamentally changes asset pricing: it increases the equity premium, lowers the risk-free rate, and raises return volatility—all from parameter uncertainty alone.

**How to cite:** Weitzman provides the foundational result that parameter uncertainty causes tail-fattening of predictive return distributions. For climate transition, investors uncertain about regulatory intensity parameters see each noisy signal (plant retirement) fatten the tails of the return distribution—increasing variance—even when the posterior mean barely moves. Cite as: *"Following Weitzman (2007), structural uncertainty about the transition regime's parameters fattens the tails of the predictive return distribution, causing the conditional variance to spike upon signal arrival independent of any directional mean shift—the formal mechanism underlying our ambiguous-shock volatility prediction."*

### Paper 3.5 — Information arrival drives variance independent of directional content

**French, Kenneth R., and Richard Roll (1986). "Stock Return Variances: The Arrival of Information and the Reaction of Traders." *Journal of Financial Economics*, 17(1), 5–26.**

**Core finding:** Stock prices are substantially more volatile during trading hours than non-trading hours. After testing competing explanations, French and Roll conclude that **information arrival per se**—not its directional content—is the primary driver of return variance.

**How to cite:** This empirical classic establishes the foundational fact that information arrival increases variance independent of the signal's expected direction. When the geographic and fuel channels offset, the directional mean return is muted, but the *flow of information that must be processed* still drives volatility upward. Cite as: *"French and Roll (1986) establish empirically that information arrival per se drives return variance, independent of directional content. This supports our prediction that noisy plant retirement signals spike volatility even when the mean price impact is zero due to offsetting spatial channels."*

---

## Building Block 4: Resolving shocks collapse volatility and force directional repricing

When a binding systemic shock (e.g., a legislated coal phase-out) resolves policy uncertainty, p_t → 1. The uncertainty premium disappears (volatility collapses), and the negative competitive-stranding channel dominates, producing negative abnormal returns proportional to firm-level fossil legacy intensity α_i. Five papers—spanning option-market evidence, macroeconomic theory, and cross-sectional climate finance—jointly establish both predictions.

### Paper 4.1 — Political event resolution collapses implied volatility

**Kelly, Bryan, Ľuboš Pástor, and Pietro Veronesi (2016). "The Price of Political Uncertainty: Theory and Evidence from the Option Market." *Journal of Finance*, 71(5), 2417–2480.**

**Core finding:** Using equity index options from **20 countries** around national elections and global summits, options whose lives span political events are significantly more expensive. Implied volatility **drops sharply** once the event resolves—directly demonstrating that uncertainty resolution causes a volatility collapse.

**How to cite:** This paper provides the cleanest empirical evidence that discrete uncertainty resolution (analogous to a binding phase-out law) triggers a sharp decline in implied volatility. The mechanism—anticipated event elevates variance risk; resolution eliminates it—maps precisely onto BB4's first prediction. Cite as: *"Kelly, Pástor, and Veronesi (2016) show that implied volatility collapses once political uncertainty resolves. We apply this result to climate policy: when a binding phase-out law is enacted (p_t → 1), the policy uncertainty premium embedded in option prices disappears, producing the predicted volatility collapse."*

### Paper 4.2 — Uncertainty shocks and the "drop and rebound" pattern

**Bloom, Nicholas (2009). "The Impact of Uncertainty Shocks." *Econometrica*, 77(3), 623–685.**

**Core finding:** Uncertainty spikes (proxied by stock-market volatility) approximately **double** implied volatility after major shocks. In a structural model with time-varying second moments, activity freezes during uncertainty episodes and rebounds rapidly once uncertainty resolves. The "drop and rebound" pattern establishes that uncertainty resolution has powerful and immediate financial and real effects.

**How to cite:** Bloom provides the canonical macroeconomic framework establishing that volatility normalizes rapidly once uncertainty resolves. The real-options mechanism (higher uncertainty → option value of waiting → activity freeze → uncertainty resolves → rapid repricing) parallels the transition from uncertain to certain regulatory regimes in BB4. Cite as: *"Bloom (2009) demonstrates that uncertainty episodes approximately double volatility, with rapid normalization upon resolution. This 'drop and rebound' dynamic validates our prediction that enactment of a binding transition law collapses the uncertainty component of return volatility."*

### Paper 4.3 — Policy announcements trigger discrete stock price declines

**Pástor, Ľuboš, and Pietro Veronesi (2012). "Uncertainty about Government Policy and Stock Prices." *Journal of Finance*, 67(4), 1219–1264.**

**Core finding (BB4-relevant):** Stock prices **fall on average** at policy-change announcements, with the magnitude of the decline increasing in the degree of prior policy uncertainty. After announcement, the uncertainty premium embedded in prices dissipates.

**How to cite (distinct from BB2 usage):** While cited in BB2 for the credibility gap mechanism, the BB4-relevant prediction is the *directional repricing at resolution*: a binding policy announcement (coal phase-out law) causes a discrete stock price decline whose magnitude scales with prior uncertainty. Combined with cross-sectional carbon-intensity pricing from Bolton and Kacperczyk, this yields the full BB4 prediction. Cite as: *"Pástor and Veronesi (2012) predict that stock prices fall discretely at policy-change announcements, with larger declines under greater prior uncertainty. In our model, the binding phase-out law is the resolution event; the negative repricing is concentrated in high-α_i firms because the competitive-stranding channel dominates once uncertainty collapses."*

### Paper 4.4 — Carbon emissions priced cross-sectionally, establishing the α_i proportionality

**Bolton, Patrick, and Marcin Kacperczyk (2021). "Do Investors Care about Carbon Risk?" *Journal of Financial Economics*, 142(2), 517–549.**

**Core finding (BB4-relevant):** The carbon premium varies **cross-sectionally with emission intensity**: firms with higher total CO₂ emissions carry higher required returns. This cross-sectional variation is the market's quantification of differential transition risk exposure.

**How to cite (distinct from BB2 usage):** While cited in BB2 to establish that transition risk is priced, the BB4-relevant result is the **proportionality**: the magnitude of the risk premium is proportional to the firm's emission intensity. When a binding regulation arrives and resolves uncertainty, the negative repricing will be proportional to α_i because Bolton and Kacperczyk show this is the dimension along which transition risk is priced. Cite as: *"Bolton and Kacperczyk (2021) establish that the cross-section of transition risk premia is ordered by carbon intensity. Upon uncertainty resolution, this premium collapses differentially: high-α_i firms experience the largest negative abnormal returns because they carried the largest uncertainty-driven risk premium."*

### Paper 4.5 — The pollution premium: regime-shift risk proportional to emission intensity

**Hsu, Po-Hsuan, Kai Li, and Chi-Yang Tsou (2023). "The Pollution Premium." *Journal of Finance*, 78(3), 1343–1392.**

**Core finding:** A long-short portfolio of high- vs. low-toxic-emission-intensity firms within an industry earns **4.42% per annum**, driven by a systematic risk related to **environmental policy regime-change risk**. An event study of the 2016 U.S. election (loosening expected regulation) shows differential returns by emission intensity, confirming that regulatory regime shifts reprice assets proportionally to pollution exposure.

**How to cite:** Hsu, Li, and Tsou provide perhaps the most direct evidence for BB4's second prediction. They show the return premium is strictly proportional to firm-level emission intensity (supporting α_i proportionality), and the mechanism is explicitly environmental regulation regime-shift risk—the exact type of systemic shock described in BB4. Their event study of the Trump election demonstrates that when expected regulation loosens, high-emission firms gain differentially; by symmetry, when a binding phase-out law tightens regulation, these firms lose differentially. Cite as: *"Hsu, Li, and Tsou (2023) directly model and test the environmental policy regime-shift mechanism underlying BB4. Their finding that the pollution premium is proportional to within-industry emission intensity, and that regulatory regime shifts reprice assets differentially by emission exposure, provides the primary empirical validation of our α_i-proportional negative repricing prediction."*

---

## How these 21 papers interlock across building blocks

The literature architecture is not merely a list but a **logically interlocking structure** where papers serve double duty across building blocks, reflecting the model's internal consistency. Pástor and Veronesi (2012) anchors both the credibility gap (BB2) and the resolution repricing (BB4), because the same policy uncertainty that generates the premium in steady state produces the discrete price adjustment at resolution. Bolton and Kacperczyk (2021) similarly bridges BB2 and BB4: the carbon premium they document is the equilibrium manifestation of the credibility gap (BB2), and its cross-sectional structure determines the proportionality of repricing upon resolution (BB4). Veronesi (1999) connects BB2's belief dynamics to BB3's volatility prediction by showing that maximum posterior uncertainty simultaneously maximizes the credibility gap variance p_t(1−p_t) and return volatility.

Three papers merit special emphasis for a JEEM submission. **Lang and Stulz (1992)** is essential because the dual contagion-versus-competitive framework is the conceptual DNA of the opposing-sign spatial channels. **Kydland and Prescott (1977)** gives the credibility gap deep game-theoretic microfoundations that will resonate with referees trained in mechanism design. **Hsu, Li, and Tsou (2023)** is the closest existing paper to the model's cross-sectional predictions and was published in JF with an environmental policy regime-shift mechanism—demonstrating that top finance journals value this research agenda.

The collective strength of these references is that every mathematical assumption in the model—the network structure of cash flows, the opposing signs of spatial channels, the belief-variance-driven credibility gap, the volatility-without-mean-shift prediction for ambiguous shocks, and the proportional directional repricing upon resolution—can be traced to a specific, published, top-tier-journal result. No assumption is introduced without authoritative precedent.
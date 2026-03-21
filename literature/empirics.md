# Empirical precedents for spatial climate transition risk econometrics

This literature review identifies **30+ authoritative papers** across four empirical blocks that collectively provide robust methodological justification for every major econometric choice in the paper—from spatial weight matrices and peer-effect identification to volatility event studies and single-sector research design. The strongest precedents come from Hoberg and Phillips (2016, *JPE*) for cosine-similarity networks, Leary and Roberts (2014, *JF*) for causal peer-effect identification, Bolton and Kacperczyk (2021, *JFE*; 2023, *JF*) for carbon-exposure pricing, Beaver (1968, *JAR*) and Hann, Kim, and Zheng (2019, *RAS*) for volatility information transfer, and Cicala (2022, *AER*) for the electricity-sector laboratory defense. Below, each block presents 5–7 papers with exact citations, core findings, and precise citation guidance.

---

## Block 1: Identifying network and peer effects in asset returns

The spatial regression of firm CARs on weighted neighbor shocks faces Manski's reflection problem. The following papers collectively justify the identification strategy, the construction of multiple non-nested weight matrices, and the use of cosine similarity for the fuel-mix network.

### 1. Manski, Charles F. (1993). "Identification of Endogenous Social Effects: The Reflection Problem." *Review of Economic Studies*, 60(3), 531–542.

**Core finding.** Demonstrates that in the linear-in-means model with uniform group membership, it is generically impossible to separately identify endogenous peer effects (behavior responds to group behavior), exogenous contextual effects (behavior responds to group characteristics), and correlated effects (common shocks).

**How to cite.** This paper frames the central identification threat. The paper should cite Manski (1993) to establish that naive regression of firm *i*'s CARs on spatially-weighted neighbor CARs could conflate true spillovers with common climate-policy exposure or correlated fundamentals. The paper's use of **three distinct, non-nested spatial weight matrices** (geographic, fuel-mix, regulatory) is then positioned as a strategy to break the reflection problem, since the same firm has different neighbors on each dimension, creating the intransitivity needed for identification.

### 2. Bramoullé, Yann, Habiba Djebbari, and Bernard Fortin (2009). "Identification of Peer Effects through Social Networks." *Journal of Econometrics*, 150(1), 41–55.

**Core finding.** Shows that when interactions are structured through a network (rather than uniform groups), peer effects are generically identified even in the linear-in-means model. Network intransitivity—friends-of-friends who are not direct friends—generates excluded instruments (G²x, G³x). The necessary and sufficient condition is linear independence of I, G, and G².

**How to cite.** This is the formal identification result that underwrites the entire spatial regression strategy. Unlike SIC-code peer groups where all firms share identical neighbors, the paper's continuous spatial weight matrices (inverse-distance, cosine-similarity) create **firm-specific, asymmetric, non-nested peer groups**. Firm A's geographic neighbor B has a fuel-mix neighbor C who is not A's neighbor on either dimension—exactly the intransitivity Bramoullé et al. require. Cite this to argue that network structure itself provides identification.

### 3. Leary, Mark T. and Michael R. Roberts (2014). "Do Peer Firms Affect Corporate Financial Policy?" *Journal of Finance*, 69(1), 139–178.

**Core finding.** Peer firms' financing decisions significantly influence own-firm capital structure, with a one-standard-deviation change in peer leverage associated with ~11% change in own-firm leverage. They instrument using **idiosyncratic equity return shocks** to peer firms—residuals from a factor model that strips out common industry/market variation—explicitly addressing Manski (1993).

**How to cite.** Provides a direct methodological template. The strategy of using idiosyncratic (residualized) shocks to neighbors as instruments applies directly: the paper can instrument spatially-weighted neighbor CARs using non-climate, non-market return shocks to neighbors. Leary and Roberts demonstrate that peer effects in equity-return-linked outcomes survive rigorous controls for common factors, establishing that this class of models recovers true spillovers in financial settings.

### 4. Hoberg, Gerard and Gordon Phillips (2016). "Text-Based Network Industries and Endogenous Product Differentiation." *Journal of Political Economy*, 124(5), 1423–1465.

**Core finding.** Develops the Text-Based Network Industry Classification (TNIC) system using **cosine similarity of 10-K product descriptions**. Each firm has a unique, time-varying set of competitors defined by proximity in "product space." TNIC classifications significantly outperform SIC/NAICS codes in explaining firm characteristics, profitability, and competitive dynamics.

**How to cite.** This is the conceptual foundation for the w_fuel weight matrix. The paper's construction of fuel-mix similarity using cosine similarity of generation portfolio vectors is **methodologically identical** to Hoberg and Phillips's approach—substituting fuel-type capacity shares for word frequencies. Cite to justify that (a) cosine similarity is the established metric for continuous firm-pair similarity, (b) continuous pairwise similarity creates richer spatial structure than discrete group membership, and (c) firm-specific, time-varying peer groups break the symmetry that causes the reflection problem in traditional industry-group models. See also Hoberg and Phillips (2010, *RFS*) for the earlier application to M&A.

### 5. Pirinsky, Christo and Qinghai Wang (2006). "Does Corporate Headquarters Location Matter for Stock Returns?" *Journal of Finance*, 61(4), 1991–2015.

**Core finding.** Documents strong comovement in stock returns of firms headquartered in the same MSA. Firms that change headquarters experience decreased comovement with the old location and increased comovement at the new one. This local comovement is **not explained by economic fundamentals** and is stronger for smaller firms with more retail investors.

**How to cite.** Establishes the foundational empirical fact that geographic proximity generates return comovement, justifying the inverse-distance w_geo matrix. For power utilities, geographic proximity captures shared grid infrastructure, weather exposure, local regulatory regimes, and overlapping investor bases. Pirinsky and Wang's finding that geographic comovement persists after controlling for industry is especially relevant—w_geo captures information beyond industry groupings.

### 6. Parsons, Christopher A., Riccardo Sabbatucci, and Sheridan Titman (2020). "Geographic Lead-Lag Effects." *Review of Financial Studies*, 33(10), 4721–4770.

**Core finding.** Documents lead-lag effects in stock returns between firms co-headquartered in the same city but operating in different sectors. Geographic lead-lags yield risk-adjusted returns of **5–6% per year** and are linked to the structure of the analyst business (organized by sector, not geography), not to standard investor-inattention proxies.

**How to cite.** Directly justifies including w_geo alongside w_fuel as **non-redundant spatial dimensions**. Parsons et al. show that geographic proximity generates return predictability independent of industry membership, confirming that location-based and product-based (fuel-mix-based) channels capture distinct spillover mechanisms. Cite to defend the simultaneous estimation of effects from both geographic and fuel-mix weight matrices.

### 7. Ahern, Kenneth R. and Jarrad Harford (2014). "The Importance of Industry Links in Merger Waves." *Journal of Finance*, 69(2), 527–576.

**Core finding.** Represents the economy as a network of industries linked through customer-supplier flows. Mergers propagate in waves along these links—transmitting to close industries quickly and distant ones with delay. **Network centrality** determines exposure and transmission speed.

**How to cite.** Establishes that economic shocks propagate along inter-industry networks with centrality amplifying transmission speed. This justifies the hypothesis that a regulatory shock (coal phase-out mandate) propagates through the spatial network of utilities—first to geographic neighbors and fuel-mix-similar peers, then further along the network. Ahern and Harford's centrality finding supports any result that well-connected utilities serve as transmission hubs for transition risk.

---

## Block 2: Event studies on climate policy and asset stranding

The event-study design around binding phase-out laws, interacting event-time dummies with coal-share intensity (α_i), is supported by the following papers that collectively establish climate policy events as valid quasi-experiments, carbon exposure as a priced characteristic, and the DID-within-event-study specification.

### 1. Bolton, Patrick and Marcin Kacperczyk (2021). "Do Investors Care about Carbon Risk?" *Journal of Financial Economics*, 142(2), 517–549.

**Core finding.** Stocks of firms with higher total CO₂ emissions earn **higher returns**, controlling for standard predictors—a "carbon premium." This premium is not explained by unexpected profitability or known risk factors. Institutional investors implement exclusionary screening based on emission intensity, particularly in salient industries including utilities.

**How to cite.** Establishes that carbon exposure is a priced firm characteristic in the cross-section of returns—the foundational result justifying coal-share intensity (α_i) as a continuous exposure variable. If carbon-transition risk is priced unconditionally, then binding phase-out laws should produce directional negative repricing proportional to exposure, exactly as the paper hypothesizes.

### 2. Bolton, Patrick and Marcin Kacperczyk (2023). "Global Pricing of Carbon-Transition Risk." *Journal of Finance*, 78(6), 3677–3754.

**Core finding.** Extends the carbon premium globally across **14,400 firms in 77 countries**. Carbon premia related to emission levels are higher in countries with **stricter domestic climate policies**, and premia increase with investor awareness of climate risk.

**How to cite.** Critical for the paper's global sample design. The finding that carbon premia intensify in stricter policy environments directly supports the hypothesis that binding national phase-out legislation should produce measurable repricing. Cite to justify that the event-study design around Tier-1 laws captures incremental tightening of the already-priced carbon premium.

### 3. Engle, Robert F., Stefano Giglio, Bryan Kelly, Heebum Lee, and Johannes Stroebel (2020). "Hedging Climate Change News." *Review of Financial Studies*, 33(3), 1184–1216.

**Core finding.** Proposes a method to dynamically hedge climate change risk using a mimicking-portfolio approach. Environmental scores serve as exposure weights; portfolios constructed on these weights generate positive returns when climate news is negative, performing well both in-sample and out-of-sample.

**How to cite.** Establishes the principle that firm-level environmental characteristics (analogous to coal-share α_i) serve as valid proxies for climate risk exposure and can be used to model heterogeneous responses to climate news innovations. The paper's use of environmental scores as weights in a portfolio/regression context directly parallels the interaction of α_i with event-time dummies.

### 4. Ilhan, Emirhan, Zacharias Sautner, and Grigory Vilkov (2021). "Carbon Tail Risk." *Review of Financial Studies*, 34(3), 1540–1571.

**Core finding.** Climate policy uncertainty is priced in the options market. The cost of downside tail-risk protection (OTM put prices) is **larger for firms with more carbon-intense business models**. This cost is magnified when public attention to climate change spikes and decreased after Trump's election.

**How to cite.** Provides direct evidence that carbon intensity drives heterogeneous exposure to regulatory downside risk—precisely the mechanism the event study captures. Cite to argue that negative CARs for high-coal-share utilities around phase-out legislation reflect the market pricing increased stranded-asset tail risk, consistent with Ilhan et al.'s finding that carbon-intense firms carry larger downside risk premia.

### 5. Koch, Nicolas, Godefroy Grosjean, Sabine Fuss, and Ottmar Edenhofer (2016). "Politics Matters: Regulatory Events as Catalysts for Price Formation under Cap-and-Trade." *Journal of Environmental Economics and Management*, 78, 121–139.

**Core finding.** Using event studies on **29 hand-collected political announcements** affecting the EU ETS supply schedule, documents high market responsiveness to regulatory events—price drops during backloading debates, positive reactions to 2020/2030 policy packages. Political uncertainty fundamentally shapes carbon market pricing.

**How to cite.** Published in JEEM, this is the key precedent for event-study methodology applied to climate regulatory announcements. Validates treating discrete legislative events (coal phase-out laws) as identifiable shocks whose price impact can be isolated via standard event-study techniques. Their hand-collection of political events parallels the paper's careful identification of Tier-1 binding phase-out laws.

### 6. Seltzer, Lee, Laura Starks, and Qifei Zhu (2022). "Climate Regulatory Risk and Corporate Bonds." *NBER Working Paper No. 29994 / FRB New York Staff Report No. 1014*.

**Core finding.** Uses the Paris Agreement as a shock to expected climate regulation in a **DID specification: Spread_it = β₁(EnvProfile_j × AfterParis_t) + γ_i + κ_t + ε_it**. Bonds of firms with poor environmental profiles experienced credit rating downgrades and wider yield spreads. A structural credit model attributes wider spreads to changes in asset volatilities, not asset values.

**How to cite.** The single closest methodological template for the interaction specification. Their DID approach—interacting a continuous exposure variable with a post-event dummy—is **isomorphic** to the paper's coal-share × post-legislation specification. Substitute α_i for environmental profile and the binding phase-out date for the Paris Agreement date, and the estimating equation is identical. Their finding that spreads reflect volatility changes (not value changes) also connects to Block 3.

### 7. Monasterolo, Irene and Luca de Angelis (2020). "Blind to Carbon Risk? An Analysis of Stock Market Reaction to the Paris Agreement." *Ecological Economics*, 170, 106571.

**Core finding.** Low-carbon indices' systematic risk decreased post-Paris Agreement, while carbon-intensive indices showed milder, partial repricing—suggesting markets only partially priced the transition risk signal from a voluntary international agreement.

**How to cite.** Motivates the paper's focus on **binding national legislation** rather than voluntary international agreements. The finding of weak/partial repricing after the Paris Agreement establishes the baseline against which the paper's results should be read: if even voluntary agreements produce some repricing, legally binding coal phase-out laws with enforcement mechanisms should produce stronger, more detectable CARs.

---

## Block 3: Measuring volatility changes around information shocks

The EIA-860 volatility result—increased realized volatility for exposed neighbors when mean CARs are zero—requires justifying variance-based event studies as an independent information channel. These papers collectively establish that second-moment changes are informative, that they transmit to peers, and that ambiguous signals generate volatility without moving means.

### 1. Beaver, William H. (1968). "The Information Content of Annual Earnings Announcements." *Journal of Accounting Research*, 6 (Supplement), 67–92.

**Core finding.** Both return variance and trading volume spike dramatically during earnings announcement weeks relative to non-announcement weeks. The **variance ratio was approximately 1.67**, even in cases where mean price changes were zero.

**How to cite.** The foundational paper establishing that changes in return variance around information events are themselves informative, independent of the sign or magnitude of mean returns. When neighboring firm CARs average to zero but realized volatility increases, this follows directly from Beaver's framework: the announcement is informationally relevant even without directional price movement. Cite as the conceptual origin of using second-moment changes as evidence of information content.

### 2. Boehmer, Ekkehart, Jim Musumeci, and Annette B. Poulsen (1991). "Event-Study Methodology under Conditions of Event-Induced Variance." *Journal of Financial Economics*, 30(2), 253–272.

**Core finding.** Standard event-study test statistics (e.g., Patell test) severely over-reject the null when events induce variance changes. They propose the **BMP cross-sectional standardized test** that adjusts for event-induced variance.

**How to cite.** The methodological cornerstone for any paper measuring variance changes around events. BMP (1991) establishes that (a) events routinely induce variance changes, (b) ignoring this leads to incorrect inference, and (c) variance changes are an essential feature of the event, not a nuisance parameter. Cite to justify treating the change in realized return volatility (post vs. pre) as a formal test statistic and to defend the statistical framework for testing whether event-induced variance is significantly different from zero across spatially-exposed neighbors.

### 3. Epstein, Larry G. and Martin Schneider (2008). "Ambiguity, Information Quality, and Asset Pricing." *Journal of Finance*, 63(1), 197–228.

**Core finding.** When ambiguity-averse investors process signals of uncertain quality, they adopt worst-case assessments. The model predicts that **poor information quality (ambiguity) increases price volatility** and generates ambiguity premia dependent on idiosyncratic risk. Ambiguous signals contribute to "excess volatility" without necessarily moving mean prices.

**How to cite.** The primary theoretical reference for interpreting why EIA-860 plant retirements increase neighbor volatility without moving mean CARs. Epstein and Schneider predict precisely this pattern: when signal quality is poor—as when a utility cannot easily assess what a peer plant's retirement implies for its own operations and asset values—disagreement increases, elevating return volatility (second moment) without a clear directional effect (first moment). Cite as the theoretical foundation for the "ambiguity-driven volatility" interpretation.

### 4. Hann, Rebecca N., Heedong Kim, and Yue Zheng (2019). "Intra-Industry Information Transfers: Evidence from Changes in Implied Volatility Around Earnings Announcements." *Review of Accounting Studies*, 24(3), 927–971.

**Core finding.** Documents a significantly positive association between changes in **implied volatility of an industry's first earnings announcer and its peer firms**, even after controlling for first-moment (return) information transfer. This "second-moment information transfer" is stronger for bellwether firms and during periods of greater macroeconomic uncertainty.

**How to cite.** The most directly analogous empirical precedent. Hann et al. demonstrate that corporate announcements transmit volatility—not just returns—to industry peers, and that this **second-moment information transfer operates as a distinct channel** from first-moment CAR spillovers. This validates the interpretation that increased realized volatility among spatially-exposed neighbors reflects genuine information transmission about uncertainty, not just correlated return shocks. The spatial weight matrices simply replace the "same industry" peer definition with geographic and fuel-mix proximity.

### 5. Patell, James M. and Mark A. Wolfson (1984). "The Intraday Speed of Adjustment of Stock Prices to Earnings and Dividend Announcements." *Journal of Financial Economics*, 13(2), 223–252.

**Core finding.** Using intraday data, variance spikes within minutes of announcements and **disturbances in variance persist for several hours** into the following trading day, even after mean returns have reverted. Variance is a more sensitive and persistent detector of information arrival than mean returns.

**How to cite.** Establishes that variance disturbances from information shocks persist beyond mean-return effects and that variance is a more sensitive detector of information arrival. Justifies using multi-day post-announcement realized volatility windows around EIA-860 retirement dates: if variance persists longer than mean CAR effects, the volatility test captures information transmission that a standard CAR test misses.

### 6. Dubinsky, Andrew, Michael Johannes, Andreas Kaeck, and Norman J. Seeger (2019). "Option Pricing of Earnings Announcement Risks." *Review of Financial Studies*, 32(2), 646–687.

**Core finding.** Develops models to separate price uncertainty generated by scheduled announcements from normal day-to-day volatility. The anticipated price uncertainty around earnings dates is quantitatively large, varies substantially, and is informative about future return volatility. Option-implied volatility exhibits **deterministic jumps on announcement dates**.

**How to cite.** Provides a direct precedent for isolating the variance component attributable to a specific scheduled information event. Their methodology of measuring excess variance associated with a known event date parallels computing post- vs. pre-announcement realized volatility around EIA-860 retirement dates. Demonstrates that event-associated variance spikes are a well-identified, priced quantity in financial markets—not statistical noise.

---

## Block 4: The power utility sector as an empirical laboratory

The single-sector restriction is a strength, not a weakness. A distinguished tradition in applied economics uses electricity as a laboratory precisely because assets are observable, output is homogeneous, and costs are directly measurable—eliminating the measurement noise that plagues cross-sector climate finance studies.

### 1. Cicala, Steve (2022). "Imperfect Markets versus Imperfect Regulation in US Electricity Generation." *American Economic Review*, 112(2), 409–441.

**Core finding.** Exploiting staggered transition from regulated to market dispatch (1999–2012), constructs a virtually complete **hourly characterization of US electric grid supply and demand**. Markets reduce production costs by 5%, with aggregate savings of ~$3 billion/year.

**How to cite.** Provides the strongest defense of single-sector design. Cicala explicitly argues (p. 410) that *"this study focuses on a single industry that has undergone a profound reorganization"* as an advantage, not a limitation. His use of EIA plant-level data (Forms 860 and 923) combined with EPA CEMS emissions data at hourly resolution demonstrates that no other industry offers comparable granularity. Cite to establish that restricting analysis to power utilities enables measurement precision impossible in cross-sector studies, and to justify EIA-860 as a primary data source.

### 2. Fabrizio, Kira R., Nancy L. Rose, and Catherine D. Wolfram (2007). "Do Markets Reduce Costs? Assessing the Impact of Regulatory Restructuring on US Electric Generation Efficiency." *American Economic Review*, 97(4), 1250–1277.

**Core finding.** Using plant-level data, shows investor-owned plants in restructured states achieved the largest efficiency gains from deregulation, while publicly owned plants had the smallest gains.

**How to cite.** The canonical paper demonstrating that electricity generation offers a uniquely clean empirical laboratory because **plant-level inputs and outputs are directly observable**, EIA data provides comprehensive coverage, output is homogeneous (a MWh is a MWh), and regulatory variation is quasi-experimental. The authors explicitly argue (p. 1253) that electricity generation is particularly well-suited for studying institutional change because plant-level production data are publicly available and cross-plant regulatory variation is exogenous.

### 3. Davis, Lucas W. and Catherine Hausman (2016). "Market Impacts of a Nuclear Power Plant Closure." *American Economic Journal: Applied Economics*, 8(2), 92–122.

**Core finding.** The abrupt closure of San Onofre Nuclear Generating Station was met by increased in-state natural gas generation at ~$68,000/hour, creating binding transmission constraints and short-run inefficiencies.

**How to cite.** Demonstrates that **plant closures are discrete, datable events** providing sharp identification, and that the physical grid imposes observable geographic constraints (transmission congestion) creating exploitable spatial variation. A single plant closure has traceable, system-wide consequences—directly analogous to studying how one plant's retirement transmits risk to neighboring utilities through the grid. The methodology of tracing spatial propagation of a plant-level shock through grid infrastructure parallels the paper's use of geographic weight matrices.

### 4. Borenstein, Severin, James B. Bushnell, and Frank A. Wolak (2002). "Measuring Market Inefficiencies in California's Restructured Wholesale Electricity Market." *American Economic Review*, 92(5), 1376–1405.

**Core finding.** Decomposes wholesale electricity payments into production costs, competitive rents, and market power. In summer 2000, 59% of surging expenditures were attributable to market power.

**How to cite.** The gold standard for demonstrating that the electricity sector permits **direct observation of marginal costs**—something impossible in most industries. Because generation technology is well understood (heat rate × fuel cost = marginal cost), researchers can construct competitive benchmarks from engineering data. This transparency is precisely why the power sector avoids ESG-score noise: researchers verify physical asset characteristics, actual fuel consumption, and output independently using the same EIA data infrastructure.

### 5. Fowlie, Meredith (2010). "Emissions Trading, Electricity Restructuring, and Investment in Pollution Abatement." *American Economic Review*, 100(3), 837–869.

**Core finding.** Under NOx emissions trading, deregulated plants were less likely to adopt capital-intensive environmental compliance options. Pollution concentrates in states with worse air quality due to regulatory heterogeneity.

**How to cite.** Shows how **regulatory heterogeneity across electricity markets creates quasi-experimental variation** for studying investment responses to environmental policy—precisely the variation the climate transition risk paper exploits. Fowlie uses unit-level variation in compliance costs and exogenous state-level restructuring to identify causal effects, demonstrating that the electricity sector's combination of detailed plant-level data, heterogeneous regulatory environments, and observable investment decisions makes it an ideal laboratory for studying firm responses to environmental policy.

### 6. Wolfram, Catherine D. (1999). "Measuring Duopoly Power in the British Electricity Spot Market." *American Economic Review*, 89(4), 805–826.

**Core finding.** Using direct measures of marginal cost, prices in the British electricity market were above marginal cost but not nearly as high as oligopoly theory predicts. Regulatory constraints, entry threats, and contracts moderate pricing.

**How to cite.** Demonstrates that electricity markets provide the rare setting where **direct measures of marginal cost are available to researchers**, enabling precise measurement without indirect estimation approaches. The same cost transparency makes the sector ideal for climate transition risk research—stranding risk can be directly observed from fuel type, technology, and cost position rather than inferred from opaque ESG scores or self-reported emissions data.

### 7. Kellogg, Ryan (2014). "The Effect of Uncertainty on Investment: Evidence from Texas Oil Drilling." *American Economic Review*, 104(6), 1698–1734.

**Core finding.** Firms' responses to changes in oil price volatility match the magnitude prescribed by real options theory—firms reduce drilling when expected price volatility is high, consistent with optimal investment delay under uncertainty.

**How to cite.** While focused on oil rather than electricity, Kellogg demonstrates the broader principle that **energy sectors provide uniquely clean empirical settings for testing investment theories** because investment decisions are discrete and observable, output prices are publicly traded, and production functions are well understood. The same logic applies to power plant investment and retirement decisions under climate transition uncertainty: EIA-860 records observable entry/exit decisions, and electricity price uncertainty can be measured from futures markets. Supports using energy sectors as natural laboratories for investment under policy uncertainty.

---

## How these precedents map to each methodological choice

The table below synthesizes the strongest citation for each contested methodological element, providing a quick-reference guide for drafting the paper's methodology section.

| Methodological element | Primary precedent | Supporting precedents |
|---|---|---|
| Inverse-distance geographic weight matrix (w_geo) | Pirinsky and Wang (2006, *JF*) | Parsons et al. (2020, *RFS*); Davis and Hausman (2016, *AEJ*) |
| Cosine-similarity fuel-mix weight matrix (w_fuel) | Hoberg and Phillips (2016, *JPE*) | Hoberg and Phillips (2010, *RFS*) |
| Multiple non-nested network layers | Bramoullé et al. (2009, *J. Econometrics*) | Fracassi (2017, *Mgmt Science*); Parsons et al. (2020, *RFS*) |
| Reflection-problem identification | Manski (1993, *REStud*) | Leary and Roberts (2014, *JF*); Bramoullé et al. (2009) |
| Climate policy event studies | Koch et al. (2016, *JEEM*) | Ramelli et al. (2021, *RCFS*); Monasterolo and de Angelis (2020) |
| Coal-share × event-time DID interaction | Seltzer et al. (2022, *NBER WP*) | Bolton and Kacperczyk (2021, *JFE*; 2023, *JF*) |
| Stranded-asset negative repricing | Ilhan et al. (2021, *RFS*) | Bolton and Kacperczyk (2023, *JF*) |
| Volatility (not mean) event studies | Beaver (1968, *JAR*) | Boehmer et al. (1991, *JFE*); Patell and Wolfson (1984, *JFE*) |
| Second-moment peer information transfer | Hann et al. (2019, *RAS*) | Barth and So (2014, *Accounting Review*) |
| Ambiguity → volatility without mean shift | Epstein and Schneider (2008, *JF*) | Dubinsky et al. (2019, *RFS*) |
| Single-sector (electricity) defense | Cicala (2022, *AER*) | Fabrizio et al. (2007, *AER*); Borenstein et al. (2002, *AER*) |
| EIA plant-level data in top journals | Fabrizio et al. (2007, *AER*) | Cicala (2015, 2022, *AER*); Fowlie (2010, *AER*) |

## Two gaps worth noting

First, no paper in a top-5 finance journal was found using **EIA-860 retirement announcements as event dates** for a financial markets event study. This appears to be a genuinely novel contribution and should be framed as such. Second, the Ramelli et al. (2021) paper on COVID-era climate repricing was published in the *Review of Corporate Finance Studies*, not the *Review of Financial Studies*—an important distinction when claiming top-tier journal precedents. The strongest top-5 finance citations for the climate event-study block remain Bolton and Kacperczyk (2021, *JFE*; 2023, *JF*), Ilhan et al. (2021, *RFS*), and Engle et al. (2020, *RFS*).
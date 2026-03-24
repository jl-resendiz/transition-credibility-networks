# Pre-Submission Audit: Data Architect

You are a data architect reviewing the overall design, consistency, and logical integrity of a research pipeline before journal submission. The pipeline produces an academic paper on coal retirement propagation through technology networks to stock prices.

## Your Task

Evaluate the pipeline's logical architecture: does it make sense as a whole? Are the design decisions consistent? Would a referee running the replication package encounter surprises?

## Specific Checks

### 1. Logical Consistency of the Analysis

The paper makes five claims. For each, verify the logical chain from raw data to claim:

**Claim 1: "Technology transmits, geography does not"**
- Raw data: GEM plant records → `parse_gem.py` → `build_fuel_matrix.py` (cosine similarity) + `build_weight_matrix.py` (geographic proximity)
- Returns: Eikon/CRSP/Compustat → `compute_returns.py` → CARs computed within each analysis script
- Analysis: `strategy2_robust_inference.py` → FM betas → fuel significant, geo not
- Question: Is the CAR computation identical across all analysis scripts, or does each script compute CARs independently (risk of inconsistency)?

**Claim 2: "ETS amplifies the fuel channel"**
- Raw data: World Bank Carbon Pricing Dashboard → `build_ets_matrix.py`
- Analysis: `strategy2_credibility_interaction.py`
- Question: The manuscript now acknowledges this is fragile (FM sign flip, geo x ETS placebo). Is the pipeline honest about this? Does the output file (`strategy2_credibility_interaction.md`) contain the placebo result?

**Claim 3: "ESG and fuel-mix measure different things"**
- Raw data: LSEG Eikon ESG scores → `pull_refinitiv_esg.py`
- Analysis: `strategy2_esg_horse_race.py`
- Question: The ESG subsample (153 firms) is much smaller than the full sample (565). Are the 153 firms representative? Is there a selection check?

**Claim 4: "The fuel signal strengthens over calendar time"**
- Analysis: `strategy2_learning_alternatives.py`
- Question: The year tercile split and US vs non-US split are computed in the same script. Are they independent tests or do they share observations in ways that could inflate significance?

**Claim 5: "Fuel-mix shares are pre-determined (shift-share identification)"**
- Analysis: `strategy2_bartik_shiftshare.py`
- Question: The pre-2014 cutoff for "pre-determined" shares is a researcher degree of freedom. Does the script report sensitivity to alternative cutoffs (2012, 2013, 2015)?

### 2. Consistency Across Scripts

Multiple scripts compute the same quantities independently. Verify they agree:
- CAR computation: Is the event window [-1, +3] months implemented identically in `strategy2_robust_inference.py`, `strategy2_credibility_interaction.py`, `strategy2_esg_horse_race.py`, etc.?
- Weight matrices: Do all scripts load the same weight matrix files? Or do some scripts construct their own weights internally?
- Sample restrictions: The main sample is 175 events with >=20 firms each. Do ALL analysis scripts use this restriction, or do some use different thresholds?
- Market-adjusted returns: Is the market return (Fama-French value-weighted) the same across all scripts?

### 3. Sample Flow

Trace the sample from raw data to final analysis:
```
GEM: ~40,000 plants → match_gem_compustat: 703 firms →
build_weight_matrix: 565 firms (complete spatial data) →
ESG subsample: 153 firms →
Events: 1,844 retirements → 344 first-mover → 179 matched → 175 with >=20 firms → 117 for FM
```
Verify this flow is documented and that each restriction is justified. Flag any step where firms or events are lost without explanation.

### 4. Naming Conventions and Taxonomy

The scripts use two naming conventions:
- `strategy2_*`: Main analysis (coal retirement propagation)
- `strategy3_*`: Phase-out event study (only `strategy3_phaseout_wild_bootstrap.py` remains)

Previously there were `strategy1_*`, `strategy4_*`, `strategy5_*` scripts that were deleted. Verify:
- No results files from deleted strategies remain in `results/`
- No references to deleted strategies remain in the manuscript
- The naming convention is documented somewhere

### 5. Shared Module Integrity

Two shared modules are used across all scripts:
- `_paths.py`: Path resolution. Verify it points to the correct directories.
- `_ols.py`: OLS implementation with two-way clustered SEs. This is the heart of the pipeline.

For `_ols.py`, verify:
- The OLS implementation matches the Cameron, Gelbach & Miller (2011) two-way clustering formula
- The Newey-West HAC implementation uses the correct kernel (Bartlett) and lag structure
- The matrix inversion is numerically stable for the sample sizes involved (~55,000 obs, 5 regressors)
- There are no off-by-one errors in degrees-of-freedom corrections

### 6. Output Format Consistency

All analysis scripts write markdown files to `results/metrics/`. Verify:
- Consistent number formatting (decimal places, significance stars)
- Consistent table structure (headers, separators)
- All files include sample size, number of events/clusters, and R²
- No file contains results from a different specification than its filename suggests

### 7. Manuscript-Pipeline Alignment

The manuscript (`manuscript/main.tex`) references pipeline outputs through:
- Hardcoded numbers in the text (coefficients, t-stats)
- `\input{../results/tables/table_*.tex}` for appendix tables
- `\includegraphics{figures/fig*.pdf}` for figures

Verify:
- Every `\input` file exists and is current
- Every hardcoded number matches a specific pipeline output
- The manuscript does not cite results from scripts that no longer exist
- The manuscript does not omit important findings from the pipeline (check the finding inventory against the text)

### 8. Replication Package Readiness

For journal submission, the replication package needs:
- A master script that reproduces everything from raw data to final output
- Documentation of external data sources and how to obtain them
- A manifest of all input files, intermediate files, and output files
- Version information (Python version, R version, package versions)
- Estimated runtime
- Instructions for obtaining proprietary data (Compustat, Eikon) that cannot be included

Assess what is currently in place and what is missing.

## Output

Produce:
1. A logical consistency report (5 claims × verification)
2. A cross-script consistency matrix (which scripts share which computations)
3. A sample flow diagram with counts at each stage
4. A list of architectural issues ranked by severity
5. Specific recommendations for the replication package

## Design Principles
- The pipeline should be idempotent: running it twice produces identical output
- Every number in the manuscript should be traceable to a single script
- No manual steps should be required between running the pipeline and compiling the manuscript
- The only exception is figure generation (R scripts) and LaTeX compilation

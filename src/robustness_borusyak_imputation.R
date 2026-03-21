# ──────────────────────────────────────────────────────────────────────
# Borusyak, Jaravel & Spiess (2024) Imputation Estimator
#
# The paper pools 179 first-mover retirement events across different
# dates and geographies. The standard pooled event-study estimate
# could suffer from negative weights if treatment effects are
# heterogeneous across cohorts (early vs late retirements).
#
# This script implements three estimators for comparison:
#   1. Pooled two-way FE (baseline)
#   2. Sun & Abraham (2021) interaction-weighted via sunab()
#   3. Imputation estimator (Borusyak, Jaravel & Spiess 2024)
#
# References:
#   Borusyak, Jaravel & Spiess (2024) Review of Economic Studies
#   Sun & Abraham (2021) Journal of Econometrics
# ──────────────────────────────────────────────────────────────────────
library(fixest)
library(data.table)

# Resolve root: works with Rscript and source()
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("--file=", .args, value = TRUE)
if (length(.file_arg) > 0) {
  .script_dir <- dirname(normalizePath(sub("--file=", "", .file_arg)))
} else {
  .script_dir <- getwd()
}
root <- normalizePath(file.path(.script_dir, ".."))
derived <- file.path(root, "data", "derived")
raw     <- file.path(root, "data", "raw")
results <- file.path(root, "results", "metrics")

cat("=== Borusyak et al. (2024) Imputation Estimator ===\n\n")
cat("Loading data...\n")

# ── Monthly returns ──────────────────────────────────────────────────
ret <- fread(file.path(derived, "returns", "monthly_returns.csv"))
ret[, ym := substr(datadate, 1, 7)]

# Fama-French monthly factors -> vwretd
ff_lines <- readLines(file.path(raw, "factors", "F-F_Research_Data_Factors.csv"))
ff_rows <- list()
for (line in ff_lines) {
  parts <- trimws(strsplit(line, ",")[[1]])
  if (length(parts) >= 5 && nchar(parts[1]) == 6 && grepl("^[0-9]+$", parts[1])) {
    ym <- paste0(substr(parts[1], 1, 4), "-", substr(parts[1], 5, 6))
    mktrf <- as.numeric(parts[2])
    rf    <- as.numeric(parts[5])
    if (!is.na(mktrf) && !is.na(rf)) {
      ff_rows[[length(ff_rows) + 1]] <- data.table(ym = ym, vwretd = (mktrf + rf) / 100)
    }
  }
}
ff <- rbindlist(ff_rows)

# Merge and compute abnormal returns
ret <- merge(ret, ff, by = "ym", all.x = TRUE)
ret[, ar := ret_monthly - vwretd]

cat(sprintf("  Monthly returns: %d obs, %d firms\n", nrow(ret), uniqueN(ret$gvkey)))

# ── Retirement events (first-movers with matched gvkeys) ────────────
events <- fread(file.path(derived, "events", "coal_retirement_events.csv"))

# Filter: first movers only
fm_col <- intersect(c("first_mover", "is_first_mover"), names(events))
if (length(fm_col) > 0) {
  events <- events[tolower(as.character(get(fm_col[1]))) %in% c("1", "true")]
}

# Require matched gvkeys
events <- events[!is.na(matched_gvkeys) & nchar(as.character(matched_gvkeys)) > 0]

# Convert date columns: fread may read dates as integer (days since epoch)
for (col in c("event_date", "announcement_date")) {
  if (col %in% names(events)) {
    if (is.integer(events[[col]]) || is.numeric(events[[col]])) {
      events[, (col) := as.character(as.Date(get(col), origin = "1970-01-01"))]
    } else {
      events[, (col) := as.character(get(col))]
    }
  }
}

# Build event_ym: prefer announcement_date > event_date > ret_year-06
events[, event_ym := NA_character_]
if ("announcement_date" %in% names(events)) {
  events[nchar(announcement_date) >= 7 & !is.na(announcement_date),
         event_ym := substr(announcement_date, 1, 7)]
}
if ("event_date" %in% names(events)) {
  events[is.na(event_ym) & !is.na(event_date) & nchar(event_date) >= 7,
         event_ym := substr(event_date, 1, 7)]
}
if ("ret_year" %in% names(events)) {
  events[is.na(event_ym) & !is.na(ret_year),
         event_ym := paste0(ret_year, "-06")]
}
events <- events[!is.na(event_ym)]

# Create unique event_id
events[, event_id := paste0("ret_", seq_len(.N))]

cat(sprintf("  Events: %d first-mover retirements with matched firms and dates\n", nrow(events)))

# ── Weight matrices ──────────────────────────────────────────────────
load_W <- function(path) {
  w <- fread(path)
  setnames(w, c("gvkey_i", "gvkey_j", "w_ij"))
  w[, w_ij := as.numeric(w_ij)]
  w
}
W_fuel <- load_W(file.path(derived, "networks", "weight_matrix_W_fuel.csv"))
W_geo  <- load_W(file.path(derived, "networks", "weight_matrix_W_geo.csv"))

cat(sprintf("  W_fuel: %d edges, W_geo: %d edges\n", nrow(W_fuel), nrow(W_geo)))

# ── Utility: add months to YYYY-MM string ────────────────────────────
add_months <- function(ym, delta) {
  y <- as.integer(substr(ym, 1, 4))
  m <- as.integer(substr(ym, 6, 7)) + delta
  while (m > 12) { y <- y + 1; m <- m - 12 }
  while (m < 1)  { y <- y - 1; m <- m + 12 }
  sprintf("%04d-%02d", y, m)
}

# ── Build event-firm-month panel ─────────────────────────────────────
cat("Building event-firm-month panel (tau = -6 to +12)...\n")
all_gvkeys <- unique(ret$gvkey)
obs_list <- list()

for (i in seq_len(nrow(events))) {
  ev <- events[i]
  ev_ym     <- ev$event_ym
  ev_id     <- ev$event_id
  ev_gvkeys <- trimws(unlist(strsplit(as.character(ev$matched_gvkeys), ";")))
  ev_gvkeys <- ev_gvkeys[nchar(ev_gvkeys) > 0]
  if (length(ev_gvkeys) == 0) next
  ev_gvkey <- ev_gvkeys[1]  # primary event firm for W lookup

  # Event window: tau = -6 to +12
  tau_range <- -6:12
  months_window <- sapply(tau_range, function(d) add_months(ev_ym, d))

  for (gk in all_gvkeys) {
    # Skip the event firm itself
    if (gk %in% ev_gvkeys) next

    # Fuel-similarity weight (the key exposure channel)
    w_f <- W_fuel[gvkey_i == ev_gvkey & gvkey_j == gk, w_ij]
    if (length(w_f) == 0) w_f <- 0
    w_f <- w_f[1]

    # Get returns in the event window
    firm_ret <- ret[gvkey == gk & ym %in% months_window]
    if (nrow(firm_ret) < 6) next

    for (j in seq_len(nrow(firm_ret))) {
      tau_idx <- which(months_window == firm_ret$ym[j])
      if (length(tau_idx) == 0) next
      tau_val <- tau_range[tau_idx[1]]
      post_val <- as.integer(tau_val >= 0)

      obs_list[[length(obs_list) + 1]] <- data.table(
        event_id  = ev_id,
        gvkey     = gk,
        ym        = firm_ret$ym[j],
        tau       = tau_val,
        ar        = firm_ret$ar[j],
        w_fuel    = w_f,
        post      = post_val,
        w_fuel_post = w_f * post_val,
        cohort    = ev_ym  # cohort = event timing for Sun & Abraham
      )
    }
  }
  if (i %% 20 == 0 || i == nrow(events)) {
    cat(sprintf("  Processed %d/%d events\n", i, nrow(events)))
  }
}

panel <- rbindlist(obs_list)
cat(sprintf("Panel: %d obs, %d events, %d firms, %d months\n",
            nrow(panel), uniqueN(panel$event_id),
            uniqueN(panel$gvkey), uniqueN(panel$ym)))

# Ensure factor types for FE
panel[, gvkey  := as.factor(gvkey)]
panel[, ym     := as.factor(ym)]
panel[, cohort := as.factor(cohort)]
panel[, event_id := as.factor(event_id)]

# ── Estimator 1: Pooled TWFE (baseline) ─────────────────────────────
cat("\n--- Estimator 1: Pooled TWFE ---\n")
est_pooled <- feols(ar ~ w_fuel_post | gvkey + ym,
                    data = panel, cluster = ~event_id)
cat("  Done.\n")
ct_pooled <- coeftable(est_pooled)
cat(sprintf("  w_fuel_post: coeff = %.5f, SE = %.5f, t = %.2f\n",
            ct_pooled["w_fuel_post", 1],
            ct_pooled["w_fuel_post", 2],
            ct_pooled["w_fuel_post", 3]))

# ── Estimator 2: Sun & Abraham (2021) via sunab() ───────────────────
cat("\n--- Estimator 2: Sun & Abraham (2021) ---\n")

# sunab() requires: a cohort variable (timing of treatment) and
# a relative-time variable. Since all firms in an event-cohort
# get "treated" (exposed) at the event date, we construct:
#   - cohort: the event_ym (when the retirement was announced)
#   - tau: relative time to event
# Never-treated units need cohort = Inf or a very large number.
# In our setting every observation is attached to an event,
# so we create a combined panel with untreated cohort markers.

# For sunab, we need to weight by w_fuel to capture exposure.
# The standard approach: interact sunab indicators with w_fuel.
# However, sunab() is designed for binary treatment.
# We use a pragmatic approach: restrict to exposed firms (w_fuel > 0)
# and estimate the dynamic effects.

panel_exposed <- panel[w_fuel > 0]
cat(sprintf("  Exposed panel (w_fuel > 0): %d obs\n", nrow(panel_exposed)))

# For sunab, cohort must indicate *when* treatment starts.
# Never-treated observations need a special cohort marker.
# In our design, all obs are "treated" at tau=0 (exposed to
# an event), so pre-period obs are "not yet treated".
# We convert cohort to numeric year-month for sunab.
panel_exposed[, cohort_num := as.integer(as.factor(cohort))]

est_sa <- tryCatch({
  feols(ar ~ sunab(cohort, tau) | gvkey + ym,
        data = panel_exposed, cluster = ~event_id)
}, error = function(e) {
  cat(sprintf("  sunab() error: %s\n", e$message))
  NULL
})

if (!is.null(est_sa)) {
  cat("  Done.\n")
  # Aggregate the ATT from sunab: average across positive tau
  sa_agg <- summary(est_sa, agg = "ATT")
  ct_sa <- coeftable(sa_agg)
  cat(sprintf("  ATT: coeff = %.5f, SE = %.5f, t = %.2f\n",
              ct_sa[1, 1], ct_sa[1, 2], ct_sa[1, 3]))
} else {
  cat("  sunab() did not converge; attempting manual interaction-weighted estimator...\n")
  # Fallback: manual cohort-specific estimation
  cohorts <- unique(panel_exposed$cohort)
  n_cohorts <- length(cohorts)
  betas <- numeric(n_cohorts)
  weights <- numeric(n_cohorts)

  for (ci in seq_along(cohorts)) {
    c_data <- panel_exposed[cohort == cohorts[ci]]
    if (nrow(c_data) < 10) { betas[ci] <- NA; weights[ci] <- 0; next }
    est_c <- tryCatch(
      feols(ar ~ w_fuel_post | gvkey + ym, data = c_data),
      error = function(e) NULL
    )
    if (is.null(est_c)) { betas[ci] <- NA; weights[ci] <- 0; next }
    betas[ci] <- coef(est_c)["w_fuel_post"]
    weights[ci] <- sum(c_data$post)
  }
  valid <- !is.na(betas) & weights > 0
  sa_coeff <- sum(betas[valid] * weights[valid]) / sum(weights[valid])
  sa_se <- sd(betas[valid]) / sqrt(sum(valid))
  sa_t <- sa_coeff / sa_se
  cat(sprintf("  Manual IW: coeff = %.5f, SE = %.5f, t = %.2f (from %d cohorts)\n",
              sa_coeff, sa_se, sa_t, sum(valid)))
}

# ── Estimator 3: Imputation (Borusyak, Jaravel & Spiess 2024) ───────
cat("\n--- Estimator 3: Imputation (BJS 2024) ---\n")

# Step 1: Estimate counterfactual using only untreated (pre-event) obs
# The counterfactual model includes firm + month FE
cat("  Step 1: Estimating counterfactual from pre-event observations...\n")
panel_pre  <- panel[post == 0]
panel_post <- panel[post == 1]
cat(sprintf("  Pre-event obs: %d, Post-event obs: %d\n",
            nrow(panel_pre), nrow(panel_post)))

est_ctrl <- feols(ar ~ 1 | gvkey + ym, data = panel_pre)

# Step 2: Impute counterfactual for treated (post-event) observations
cat("  Step 2: Imputing counterfactual for post-event observations...\n")
panel_post[, ar_imputed := predict(est_ctrl, newdata = panel_post)]

# Step 3: Treatment effect = actual - imputed counterfactual
panel_post[, tau_hat := ar - ar_imputed]

# Drop any obs where imputation failed (missing FE levels)
n_before <- nrow(panel_post)
panel_post <- panel_post[!is.na(tau_hat)]
n_after <- nrow(panel_post)
if (n_before > n_after) {
  cat(sprintf("  Dropped %d obs with missing imputation (new FE levels)\n",
              n_before - n_after))
}

# Step 4: Regress imputed treatment effect on exposure
cat("  Step 3: Regressing imputed tau_hat on fuel exposure...\n")
est_imp <- feols(tau_hat ~ w_fuel, data = panel_post, cluster = ~event_id)
ct_imp <- coeftable(est_imp)
cat(sprintf("  w_fuel: coeff = %.5f, SE = %.5f, t = %.2f\n",
            ct_imp["w_fuel", 1], ct_imp["w_fuel", 2], ct_imp["w_fuel", 3]))

# ── Collect results ──────────────────────────────────────────────────
cat("\n--- Formatting output ---\n")

# Helper: format a coefficient row
fmt_row <- function(label, coeff, se, tstat, n) {
  sprintf("| %s | %.5f | %.5f | %.2f | %d |", label, coeff, se, tstat, n)
}

# Pooled
row_pooled <- fmt_row(
  "Pooled (baseline)",
  ct_pooled["w_fuel_post", 1],
  ct_pooled["w_fuel_post", 2],
  ct_pooled["w_fuel_post", 3],
  nobs(est_pooled)
)

# Sun & Abraham
if (!is.null(est_sa)) {
  row_sa <- fmt_row(
    "Sun & Abraham (2021)",
    ct_sa[1, 1],
    ct_sa[1, 2],
    ct_sa[1, 3],
    nobs(est_sa)
  )
} else {
  row_sa <- fmt_row(
    "Sun & Abraham (2021)*",
    sa_coeff, sa_se, sa_t,
    nrow(panel_exposed)
  )
}

# Imputation
row_imp <- fmt_row(
  "Imputation (BJS 2024)",
  ct_imp["w_fuel", 1],
  ct_imp["w_fuel", 2],
  ct_imp["w_fuel", 3],
  nrow(panel_post)
)

# Build output
lines <- c(
  "# Staggered Adoption Robustness (Borusyak et al. 2024)",
  "",
  sprintf("Events: %d first-mover retirements, Panel: %d obs, Firms: %d, Cohorts: %d",
          uniqueN(panel$event_id), nrow(panel),
          uniqueN(panel$gvkey), uniqueN(panel$cohort)),
  "",
  "| Estimator | Coeff (w_fuel) | SE | t-stat | N |",
  "|---|---|---|---|---|",
  row_pooled,
  row_sa,
  row_imp,
  ""
)

# Add notes
if (is.null(est_sa)) {
  lines <- c(lines,
    "*Sun & Abraham estimated via manual interaction-weighted average across cohorts.",
    "")
}

lines <- c(lines,
  "Notes: All estimates include firm + month FE. SEs clustered by event.",
  "Pooled estimate uses w_fuel x post as treatment intensity.",
  "Sun & Abraham (2021) reports the aggregated ATT across cohorts and relative-time periods.",
  "Imputation (BJS 2024) estimates counterfactual from pre-event obs, then regresses",
  "the imputed treatment effect (actual minus counterfactual AR) on fuel-similarity exposure.",
  "Exposure variable: w_fuel (fuel-similarity weight from the spatial network)."
)

# ── Write output ─────────────────────────────────────────────────────
dir.create(results, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(results, "robustness_borusyak_imputation.md")
writeLines(lines, outfile)

cat("\n=== Results ===\n")
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nWritten to: %s\n", outfile))

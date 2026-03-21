# ──────────────────────────────────────────────────────────────────────
# Alternative Volatility Measures (Beaver 1968 / Squared Returns)
#
# The paper uses SD-based volatility change. This script adds:
#   1. |AR|-based (Beaver U-statistic): robust to GARCH clustering
#   2. Squared-return ratio: captures variance without normality
#
# Both are estimated via fixest with event-clustered SEs.
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

cat("Loading data...\n")

# Daily returns
daily <- fread(file.path(derived, "returns", "daily_returns.csv"))
daily[, date := as.character(datadate)]

# Fama-French daily → vwretd
ff_lines <- readLines(file.path(raw, "factors", "F-F_Research_Data_Factors_daily.csv"))
ff_rows <- list()
for (line in ff_lines) {
  parts <- trimws(strsplit(line, ",")[[1]])
  if (length(parts) >= 5 && nchar(parts[1]) == 8 && grepl("^[0-9]+$", parts[1])) {
    dt <- paste0(substr(parts[1],1,4), "-", substr(parts[1],5,6), "-", substr(parts[1],7,8))
    mktrf <- as.numeric(parts[2]); rf <- as.numeric(parts[5])
    if (!is.na(mktrf) && !is.na(rf))
      ff_rows[[length(ff_rows)+1]] <- data.table(date = dt, vwretd = (mktrf + rf)/100)
  }
}
ff_d <- rbindlist(ff_rows)
daily <- merge(daily, ff_d, by = "date", all.x = TRUE)
daily[, ar := ret_daily - vwretd]
daily <- daily[!is.na(ar)]

# Build sorted date list per firm
setkey(daily, gvkey, date)

# Load events (EIA-860 + phase-out)
eia <- fread(file.path(derived, "events", "eia860_announcement_events.csv"))
eia[, event_type := "eia860"]
eia <- eia[!is.na(event_date) & nchar(event_date) >= 10]
phaseout <- fread(file.path(derived, "events", "coal_phaseout_shocks_events.csv"))
tier_col <- intersect(c("tier", "exogeneity_tier"), names(phaseout))
if (length(tier_col) > 0) phaseout <- phaseout[get(tier_col[1]) == 1 | get(tier_col[1]) == "1"]
if ("binding" %in% names(phaseout)) phaseout <- phaseout[tolower(binding) %in% c("1", "true", "yes")]
phaseout <- phaseout[!is.na(event_date) & nchar(event_date) >= 10]
phaseout[, event_type := "phaseout"]

# Weight matrix (geographic)
W <- fread(file.path(derived, "networks", "weight_matrix_W_geo.csv"))
setnames(W, c("gvkey_i", "gvkey_j", "w_ij"))

# Alpha panel
alpha <- fread(file.path(derived, "fundamentals", "firm_alpha_panel.csv"))

cat("Computing volatility measures for each event-firm pair...\n")

compute_vol <- function(events_dt, exposure_type = "w_geo") {
  all_gvkeys <- unique(daily$gvkey)
  obs <- list()

  for (i in seq_len(nrow(events_dt))) {
    ev <- events_dt[i]
    ev_date <- ev$event_date
    ev_gvkey <- ev$gvkey
    ev_id <- if ("event_id" %in% names(ev)) ev$event_id else paste0("e_", i)

    for (gk in all_gvkeys) {
      if (gk == ev_gvkey) next
      firm_d <- daily[gvkey == gk]
      if (nrow(firm_d) < 50) next

      # Find event index
      dates_vec <- firm_d$date
      ar_vec    <- firm_d$ar
      ev_idx <- which(dates_vec >= ev_date)[1]
      if (is.na(ev_idx) || ev_idx < 22 || ev_idx + 20 > length(dates_vec)) next

      # Pre window: [-21, -1], Post window: [+1, +20]
      pre_ar  <- ar_vec[(ev_idx - 21):(ev_idx - 1)]
      post_ar <- ar_vec[(ev_idx + 1):(ev_idx + 20)]
      if (length(pre_ar) < 8 || length(post_ar) < 8) next

      # CAR [-1, +20]
      car_ar <- ar_vec[(ev_idx - 1):(ev_idx + 20)]
      car <- sum(car_ar)

      # 1. SD-based (baseline)
      sd_pre  <- sd(pre_ar)
      sd_post <- sd(post_ar)
      vol_sd  <- sd_post - sd_pre

      # 2. |AR|-based (Beaver)
      abs_pre  <- mean(abs(pre_ar))
      abs_post <- mean(abs(post_ar))
      vol_abs  <- if (abs_pre > 1e-10) abs_post / abs_pre - 1 else NA

      # 3. Squared return ratio
      sq_pre  <- mean(pre_ar^2)
      sq_post <- mean(post_ar^2)
      vol_sq  <- if (sq_pre > 1e-10) sq_post / sq_pre - 1 else NA

      # Exposure
      if (exposure_type == "w_geo") {
        exp_val <- W[gvkey_i == ev_gvkey & gvkey_j == gk, w_ij]
        if (length(exp_val) == 0) exp_val <- 0 else exp_val <- exp_val[1]
      } else {
        # coal_share for phase-out
        cs <- alpha[gvkey == gk, alpha]
        if (length(cs) == 0) next
        exp_val <- cs[1]
        if (is.na(exp_val)) next
      }

      obs[[length(obs) + 1]] <- data.table(
        event_id = ev_id, gvkey = gk,
        exposure = exp_val, car = car,
        vol_sd = vol_sd, vol_abs = vol_abs, vol_sq = vol_sq
      )
    }
    if (i %% 10 == 0) cat(sprintf("  %s: %d/%d events\n", exposure_type, i, nrow(events_dt)))
  }
  rbindlist(obs)
}

# EIA-860 events
cat("\n--- EIA-860 Retirement Announcements ---\n")
panel_eia <- compute_vol(eia, "w_geo")
cat(sprintf("  Pairs: %d\n", nrow(panel_eia)))

# Phase-out events
cat("\n--- Binding Phase-Out Laws (Tier-1) ---\n")
panel_po <- compute_vol(phaseout, "coal_share")
cat(sprintf("  Pairs: %d\n", nrow(panel_po)))

# ── Regressions ──────────────────────────────────────────────────────

run_reg <- function(dt, depvar) {
  if (nrow(dt) < 10) return(data.table(coeff=NA, se=NA, tstat=NA, n=0))
  dt_clean <- dt[!is.na(get(depvar)) & is.finite(get(depvar))]
  if (nrow(dt_clean) < 10) return(data.table(coeff=NA, se=NA, tstat=NA, n=0))
  fml <- as.formula(paste0(depvar, " ~ exposure"))
  est <- feols(fml, data = dt_clean, cluster = ~event_id)
  ct <- coeftable(est)
  data.table(coeff = ct["exposure", 1], se = ct["exposure", 2],
             tstat = ct["exposure", 3], n = nrow(dt_clean))
}

format_table <- function(dt, measures, labels) {
  rows <- character(0)
  for (j in seq_along(measures)) {
    r <- dt[[j]]
    rows <- c(rows, sprintf("| %s | %.4f | %.4f | %.2f | %d |",
      labels[j], r$coeff, r$se, r$tstat, r$n))
  }
  rows
}

measures <- c("vol_sd", "vol_abs", "vol_sq", "car")
labels   <- c("SD-based (baseline)", "|AR|-based (Beaver)", "Squared-return ratio", "CAR (level)")

eia_res <- lapply(measures, function(m) run_reg(panel_eia, m))
po_res  <- lapply(measures, function(m) run_reg(panel_po, m))

# ── Output ───────────────────────────────────────────────────────────

lines <- c(
  "# Volatility Robustness: Alternative Measures",
  "",
  "## EIA-860 Retirement Announcements",
  "| Measure | Coeff | SE | t-stat | N |",
  "|---|---|---|---|---|",
  format_table(eia_res, measures, labels),
  "",
  "## Binding Phase-Out Laws (Tier-1)",
  "| Measure | Coeff | SE | t-stat | N |",
  "|---|---|---|---|---|",
  format_table(po_res, measures, labels),
  "",
  "Notes: The |AR|-based measure follows Beaver (1968) and is robust to",
  "GARCH-type volatility clustering. The squared-return measure captures",
  "variance changes without assuming normality. All SEs are event-clustered.",
  "Pre window: [-21, -1] trading days. Post window: [+1, +20] trading days.",
  "CAR window: [-1, +20] trading days."
)

dir.create(results, recursive = TRUE, showWarnings = FALSE)
writeLines(lines, file.path(results, "strategy2_volatility_robustness.md"))

cat("\n=== Results ===\n")
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nWritten to: %s\n", file.path(results, "strategy2_volatility_robustness.md")))
